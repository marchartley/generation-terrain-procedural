# Imports
import argparse
import datetime
import time
from collections.abc import Sequence

import numpy as np
import matplotlib.pyplot as plt
import random, os, math
from functools import lru_cache
from matplotlib.widgets import Button
import curves
from Vectors import Vector2D, Vector3D, line_intersection
import Vectors
from typing import Tuple, List, Callable, Union, Optional, Any
from sketch_maker import LineBuilder, LineBuilder1D, LineBuilder2D, LineBuilderRadial, SketchManagement
from noise import perlin
import noise
import coralize_my_island
from scipy.ndimage import map_coordinates
from matplotlib.widgets import RangeSlider, Slider, CheckButtons
import cProfile
import pstats

# Global Variables
outputImageDims = [256, 256]
# distortionMaps: List[List[List[Vector2D]]] = []
singleDistortionMap = np.zeros((outputImageDims[0], outputImageDims[1], 2)) #: List[List[Vector2D]] = []
fig, fig2 = None, None
axesResults = []
dataset_path = "/media/marc/Data/synthetic_terrains_dataset/"


def numpyIndicesToCoords(_x: int, _y: int, sizeX: int, sizeY: int) -> Tuple[float, float] :
    x, y = intsToCoords(_y, (sizeX - _x), sizeY, sizeX)
    return x, y
def coordsToNumpyIndicesFloat(x: float, y: float, sizeX: int, sizeY: int) -> Tuple[float, float]:
    _x, _y = coordsToFloats(-y, x, sizeX, sizeY)
    return _x, _y
def coordsToNumpyIndices(x: float, y: float, sizeX: int, sizeY: int) -> Tuple[float, float]:
    _x, _y = coordsToInts(-y, x, sizeX, sizeY)
    return _x, _y


def intsToCoords(x: int, y: int, sizeX: int, sizeY: int) -> Tuple[float, float]:
    newX, newY = ((x / sizeX) - .5) * 2.0, ((y / sizeY) - .5) * 2.0
    return newX, newY

def coordsToInts(x: float, y: float, sizeX: int, sizeY: int) -> Tuple[int, int]:
    return math.floor((x + 1)*0.5*(sizeX)), math.floor((y + 1)*0.5*(sizeY))

def coordsToFloats(x: float, y: float, sizeX: int, sizeY: int) -> Tuple[float, float]:
    return (x + 1)*0.5*(sizeX), (y + 1)*0.5*(sizeY)

def wyvill(x: float):
    return ((1.0 - x)**2)**3

def initialDistoMap(sizeX: int, sizeY: int) -> np.ndarray: #List[List[Vector2D]]:
    return np.zeros((sizeX, sizeY, 2)) # [[Vector2D(0, 0) for _ in range(sizeY)] for _ in range(sizeX)]


import numpy as np
from typing import Union, Tuple, Sequence
from scipy.ndimage import map_coordinates


def deform_image(
        image: np.ndarray,
        vector_field: np.ndarray,
        scale: Union[float, Tuple[float, float]] = 1.0,
        interpolation: str = 'linear'
) -> np.ndarray:
    """
    Deform an image using a displacement vector field.

    Parameters:
        image (np.ndarray): Input image. Can be 2D (scalar), 2D with channels (e.g., RGB), or label map.
        vector_field (np.ndarray): Displacement field of shape (H, W, 2), containing (dx, dy) per pixel.
        scale (float or Tuple[float, float]): Scaling for displacement vectors.
        interpolation (str): 'linear' or 'nearest'

    Returns:
        np.ndarray: Deformed image.
    """
    assert interpolation in ('linear', 'nearest'), "Interpolation must be 'linear' or 'nearest'"
    assert vector_field.shape[-1] == 2, "Vector field must have shape (H, W, 2)"

    is_multichannel = image.ndim == 3 and image.shape[-1] == 3
    is_integer_type = np.issubdtype(image.dtype, np.integer)

    if interpolation == 'linear' and is_integer_type:
        raise ValueError("Linear interpolation is not valid for integer-type data")

    if isinstance(scale, Sequence):
        scaleX, scaleY = scale
    else:
        scaleX = scaleY = scale

    repeats = 4
    scaleX /= repeats
    scaleY /= repeats
    deformed_image = np.array(image)
    for _ in range(repeats):
        image = np.array(deformed_image)
        H, W = image.shape[:2]
        coords = np.indices((H, W), dtype=np.float32)

        # Invert and scale vector field
        dx = vector_field[..., 0] * scaleX
        dy = vector_field[..., 1] * scaleY
        displaced_coords = np.array([
            coords[0] - dy,  # y
            coords[1] - dx  # x
        ])

        order = 1 if interpolation == 'linear' else 0
        mode = 'reflect'

        if is_multichannel:
            channels = []
            for c in range(3):
                channel = image[..., c]
                warped = map_coordinates(channel, displaced_coords, order=order, mode=mode)
                channels.append(warped)
            deformed_image = np.stack(channels, axis=-1)
        else:
            deformed_image = map_coordinates(image, displaced_coords, order=order, mode=mode)

        if is_integer_type:
            # Make sure to cast back to original integer type
            deformed_image = np.rint(deformed_image).astype(image.dtype)
    return deformed_image
#
#
# def deform_image(image: np.ndarray, vector_field: np.ndarray, scale: Union[float, Tuple[float, float]] = 1.0) -> np.ndarray:
#     # Check if image and vector_field have the same dimensions
#     assert image.shape == vector_field.shape[:2], "Image and vector field must have the same dimensions"
#     if isinstance(scale, Sequence):
#         scaleX, scaleY = scale[0], scale[1]
#     else:
#         scaleX, scaleY = scale, scale
#     # Create a grid of coordinates in the original image
#     coords = np.indices(image.shape)
#
#     # Adjust the coordinates based on the vector field
#     # Note that we need to invert the vector direction for proper mapping
#     displaced_coords = np.array([
#         coords[0] - (vector_field[..., 1] * scaleY),  # y-coordinates adjust
#         coords[1] - (vector_field[..., 0] * scaleX)   # x-coordinates adjust
#     ])
#
#     # Perform the interpolation
#     deformed_image = map_coordinates(image, displaced_coords, order=1, mode='reflect')
#
#     return deformed_image

def resize_nearest_numpy(image: np.ndarray, new_shape: Tuple[int, int]) -> np.ndarray:
    """
    Resize a 2D or 3D image using nearest-neighbor interpolation in pure NumPy.

    Parameters:
        image (np.ndarray): Input image (H, W) or (H, W, C)
        new_shape (Tuple[int, int]): Desired (new_height, new_width)

    Returns:
        np.ndarray: Resized image
    """
    orig_height, orig_width = image.shape[:2]
    new_height, new_width = new_shape

    # Calculate the ratio of old/new dimensions
    row_scale = orig_height / new_height
    col_scale = orig_width / new_width

    # Compute source indices (nearest neighbor)
    row_indices = (np.arange(new_height) * row_scale).astype(int)
    col_indices = (np.arange(new_width) * col_scale).astype(int)

    # Clip indices to stay within bounds
    row_indices = np.clip(row_indices, 0, orig_height - 1)
    col_indices = np.clip(col_indices, 0, orig_width - 1)

    if image.ndim == 3:
        return image[row_indices[:, None], col_indices[None, :], :]
    else:
        return image[row_indices[:, None], col_indices[None, :]]


def getDisto(x: float, y: float, factor: float = 1.0) -> np.ndarray: # Vector2D:
    # distortion = Vector2D(x, y)
    disto1 = numpyBilinearInterpolation(singleDistortionMap, x, y) * factor
    disto2 = numpyBilinearInterpolation(singleDistortionMap, x - disto1[0], y - disto1[1]) * factor
    return (disto1 + disto2) / 2 #Vector2D(disto.x, disto.y)
    # for disto in reversed(distortionMaps):
    #     d = bilinearInterpolation(disto, distortion.x, distortion.y)
    #     d.y *= -1
    #     distortion += d
    # return distortion - Vector2D(x, y)

def getDistoFromIndices(_x: int, _y: int) -> Vector2D:
    return Vector2D(singleDistortionMap[_x, _y, 0], singleDistortionMap[_x, _y, 1]) # Possibly in the wrong order
    # sizeX, sizeY = len(distortionMaps[0][0]), len(distortionMaps[0])
    # vec: Vector2D = Vector2D(0, 0)
    # if _x < 0 or _x >= sizeX or _y < 0 or _y >= sizeY:
    #     return vec
    # for distoMap in distortionMaps:
    #     vec += distoMap[_x][_y]
    # return vec

def evaluatePosAfterDistortion(prevX: float, prevY: float, factor: float = 1.0) -> np.ndarray: # Tuple[float, float]:
    newPos = np.array([prevX, prevY]) # Vector2D(prevX, prevY)
    nbSteps = 10
    _factor = (factor / nbSteps)
    for _ in range(nbSteps):
        distor = getDisto(newPos[0], newPos[1], _factor)
        newPos -= distor
    return newPos #.x, newPos.y
    # pos = Vector2D(prevX, prevY) - getDisto(prevX, prevY) * factor
    # return pos.x, pos.y

def numpyBilinearInterpolation(arr: np.ndarray, x: float, y: float) -> float:
    sizeX, sizeY = arr.shape[:2]
    nbChannels = arr.shape[2] if arr.ndim == 3 else 0

    # Convert to integer indices and fractional parts
    _x, _y = coordsToNumpyIndices(x, y, sizeX, sizeY)
    __x, __y = coordsToNumpyIndicesFloat(x, y, sizeX, sizeY)
    dx, dy = __x - _x, __y - _y

    # Clip indices to be within array bounds
    _x0 = np.clip(_x, 0, sizeX - 1)
    _y0 = np.clip(_y, 0, sizeY - 1)
    _x1 = np.clip(_x + 1, 0, sizeX - 1)
    _y1 = np.clip(_y + 1, 0, sizeY - 1)

    # Extract corner values efficiently
    d00 = arr[_x0, _y0]
    d10 = arr[_x1, _y0]
    d01 = arr[_x0, _y1]
    d11 = arr[_x1, _y1]

    # Perform bilinear interpolation
    interpolated_value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
    return interpolated_value if nbChannels == 0 else np.clip(interpolated_value, 0, 1)


def numpyNearestNeighbor(arr: np.ndarray, x: float, y: float) -> float:
    sizeX, sizeY = arr.shape[0], arr.shape[1]
    nbChannels = 0 if len(arr.shape) == 2 else arr.shape[2]
    nullArray = 0 if nbChannels == 0 else np.zeros(nbChannels)
    __x, __y = coordsToNumpyIndicesFloat(x, y, sizeX, sizeY)
    _x, _y = min(sizeX - 1, max(0, round(__x))), min(sizeY - 1, max(0, round(__y)))

    value = arr[_x, _y]
    return value

def bilinearInterpolation(arr: np.ndarray, x: float, y: float) -> Vector2D:
    sizeX, sizeY = arr.shape[0], arr.shape[1] #len(arr[0]), len(arr)
    _x, _y = coordsToInts(x, y, sizeX, sizeY)
    __x, __y = coordsToFloats(x, y, sizeX, sizeY)
    dx, dy = __x - _x, __y - _y
    d00 = arr[_x  , _y  ] if _x >= 0 and _x < sizeX and _y >= 0 and _y < sizeY else np.array([0, 0])
    d01 = arr[_x  , _y+1] if _x >= 0 and _x < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else np.array([0, 0])
    d10 = arr[_x+1, _y  ] if (_x+1) >= 0 and (_x+1) < sizeX and _y >= 0 and _y < sizeY else np.array([0, 0])
    d11 = arr[_x+1, _y+1] if (_x+1) >= 0 and (_x+1) < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else np.array([0, 0])
    value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
    return Vector2D(value[0], value[1])

def clamp(x:float, mini: float, maxi: float) -> float:
    return max(min(maxi, x), mini)


def valueAsHSV(value: float, mini: float, maxi: float, L: float = 0.5, S: float = 1.0) -> Tuple[float, float, float]:
    H = 360.0 * (value - mini) / (maxi - mini)
    H = H % 360.0
    H /= 60.0

    C = (1 - abs(2 * L - 1)) * S
    X = C * (1 - abs(H % 2 - 1))

    if 0 <= H < 1:
        R1, G1, B1 = C, X, 0
    elif 1 <= H < 2:
        R1, G1, B1 = X, C, 0
    elif 2 <= H < 3:
        R1, G1, B1 = 0, C, X
    elif 3 <= H < 4:
        R1, G1, B1 = 0, X, C
    elif 4 <= H < 5:
        R1, G1, B1 = X, 0, C
    elif 5 <= H < 6:
        R1, G1, B1 = C, 0, X
    else:
        R1, G1, B1 = 0, 0, 0  # fallback

    m = L - C * 0.5
    R, G, B = R1 + m, G1 + m, B1 + m
    return clamp(R, 0, 1), clamp(G, 0, 1), clamp(B, 0, 1)


def getSequences(profileSketch: LineBuilder):
    centerMarker, islandCenterPos = -1, 0.0
    islandMarker, islandBorderPos = 1, 0.5
    coralMarker, coralsBorderPos = 0, 0.75
    abyssMarker, abyssesPos = 3, 1.0

    def subcurve(curve, t0, t1):
        return curve[int(t0 * len(curve)) : int(t1 * len(curve))]

    curve = [p.y for p in profileSketch.getCurve()]
    sequences = [
        (centerMarker, centerMarker, [curve[0]]),
        (centerMarker, islandMarker, subcurve(curve, islandCenterPos, islandBorderPos)),
        (islandMarker, coralMarker, subcurve(curve, islandBorderPos, coralsBorderPos)),
        (coralMarker, abyssMarker, subcurve(curve, coralsBorderPos, abyssesPos)),
        (abyssMarker, abyssMarker, [curve[-1]]),
    ]
    # for s in list(reversed(sequences)): # Use a copy of it, just to be able to modify it
    #     sequences.append((s[1], s[0], list(reversed(s[2]))))
    return sequences



def splitProfileOnMarkers(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> List[Tuple[int, int, List[float]]]:
    """Extract the curves made by the profile depending on the distance between each island sketch"""
    return getSequences(profileSketch)


def interpolateOnCurve(curve: List[Any], t: float) -> Any:
    if t >= 1.0:
        return curve[-1]
    if t <= 0:
        return curve[0]
    n = len(curve)
    segment = t * (n - 1)
    i = int(segment)
    t = segment - i

    def get_point(index: int) -> Any:
        # Ensure we handle boundary cases: extrapolate linearly at the boundaries
        if index < 0:
            return 2 * curve[0] - curve[1]  # Linear extrapolation
        elif index >= n:
            return 2 * curve[-1] - curve[-2]  # Linear extrapolation
        return curve[index]

    p0 = get_point(i - 1)
    p1 = get_point(i)
    p2 = get_point(i + 1)
    p3 = get_point(i + 2)

    # Catmull-Rom Spline Formula:
    # P(t) = 0.5 * [(2 * P1) + (-P0 + P2) * t + (2*P0 - 5*P1 + 4*P2 - P3) * t^2 + (-P0 + 3*P1 - 3*P2 + P3) * t^3]
    return 0.5 * ((2 * p1) +
                  (-p0 + p2) * t +
                  (2*p0 - 5*p1 + 4*p2 - p3) * t**2 +
                  (-p0 + 3*p1 - 3*p2 + p3) * t**3)

def getDistancesToCurves(x: float, y: float, islandSketches: SketchManagement, profileSlice: Vector2D) -> List[float]:
    distancesToCurves = [(i, curve.intersection(Vector2D(), profileSlice)) for i, curve in enumerate(islandSketches.lineBuilders)]
    distancesToCurves = [vecs[0].norm() if len(vecs) > 0 else 1000.0 for i, vecs in distancesToCurves]
    # distancesToCurves.append(1.5)  # approximation of sqrt(2)
    if distancesToCurves[1] > distancesToCurves[0]:
        tmp = distancesToCurves[0]
        distancesToCurves[0] = distancesToCurves[1]
        distancesToCurves[1] = tmp
    distanceToAbysses = max(1.5, max(distancesToCurves[:2]) + 0.001) # Force the abyss to be just behind the corals
    distancesToCurves.append(distanceToAbysses)
    return distancesToCurves

def getSequenceID(sequence: Tuple[int, int, Any], sequences: List[Tuple[int, int, List[float]]]) -> int:
    for i, seq in enumerate(sequences):
        if seq[0] == sequence[0] and seq[1] == sequence[1]:
            return i
    # We are exactly at the center of the graph
    return max([i for i, seq in enumerate(sequences) if seq[0] == -1 or seq[1] == -1])

def getSequence(sequences: List[Tuple[int, int, List[float]]], distancesToCurves: List[float], distToCenter: float) -> Tuple[int, int, float, List[float]]:
    currSeq = [-1, -1, 0]
    t = 0
    indices = [i for i in range(len(distancesToCurves))]
    indices = [x for _, x in sorted(zip(distancesToCurves, indices))]
    if distToCenter <= distancesToCurves[indices[0]]:
        t = (distToCenter / distancesToCurves[indices[0]] if distancesToCurves[indices[0]] > 0 else 0)
        currSeq = [-1, indices[0], t]
    else:
        found = False
        for i in range(len(distancesToCurves) - 1):
            if distancesToCurves[indices[i]] <= distToCenter < distancesToCurves[indices[i + 1]]:
                startProfile = distancesToCurves[indices[i]]
                endProfile = distancesToCurves[indices[i + 1]]
                t = (distToCenter - startProfile) / (endProfile - startProfile)
                currSeq = [indices[i], indices[i + 1], t]
                found = True
                break
        if not found:
            t = (distToCenter - distancesToCurves[indices[-1]]) / (2.0 - distancesToCurves[indices[-1]])
            currSeq = [indices[-1], len(indices)+1, t]

    curve = [0]
    for seq in sequences:
        if seq[0] == currSeq[0] and seq[1] == currSeq[1]:
            curve = seq[2]
            break

    marker1, marker2, t = currSeq
    return marker1, marker2, interpolateOnCurve(curve, t), curve


def smoothCurve(curve: List[Vector2D]) -> List[Vector2D]:
    res = []
    res.append(curve[0])
    for i in range(len(curve) - 1):
        res.append((curve[i] + curve[i + 1]) * .5)
    res.append(curve[-1])
    return res


def centeredCircle(radius: float, randomness: float = 1.2, seed=None, returnRadial: float = False) -> List[Vector2D]:
    if seed is None:
        seed = getRandom(0, 1000)
    perlin_noise = perlin.SimplexNoise(10)
    points: List[Vector2D] = []
    nbPoints = 30
    for i in range(nbPoints):
        angle = i * math.tau / (nbPoints - 1)
        vertexUnitPos = Vector2D(math.cos(angle), math.sin(angle))
        noiseValue = 1.0 + perlin_noise.noise2(vertexUnitPos.x + seed, vertexUnitPos.y + seed) * (randomness - 1)
        p = vertexUnitPos * (radius * noiseValue)
        if returnRadial: p = p.to_polar()
        points.append(p)
    return points

def randomDistortionCurve() -> List[Vector2D]:
    n = noise.perlin.SimplexNoise()
    p = Vector2D(random.random() * 2 - 1, random.random() * 2 - 1)
    positions: List[Vector2D] = []
    freq = 1.0/3.0
    for _ in range(random.randint(20, 40)):
        x = n.noise2(p.x / freq, p.y / freq)
        y = n.noise2(p.x / freq, p.y / freq + 1000.0)
        # p += Vector2D(x, y) * .5
        p.x, p.y = p.x + x * .5, p.y + y * .5
        toCenter = (p * -1)
        # p += 0.1 * toCenter  # /(1.5 - toCenter.norm())
        p.x, p.y = p.x + toCenter.x * .5, p.y + toCenter.y * .5
        positions.append(p.copy())
    for _ in range(5):
        positions = smoothCurve(positions)
    return curves.catmull_rom_chain(positions)

def getRandom(mini: float, maxi: float) -> float:
    return mini + random.random() * (maxi - mini)


class IslandSketch:
    def __init__(self, topViewAx: plt.Axes, profileAx: plt.Axes, distortionAx: plt.Axes, resistanceAx: plt.Axes, waterLevelSlider: Slider, coralMinMaxSlider: RangeSlider, subsidenceSlider: Slider, names: List[str], colors: List[str]):
        self.topView = SketchManagement(topViewAx)
        # self.sketches = [
        #     # "islandBorders":
        #     LineBuilderRadial(topViewAx, color="green"),
        #     # "beachBorders":
        #     LineBuilderRadial(topViewAx, color="yellow"),
        #     # "lagoonBorders":
        #     LineBuilderRadial(topViewAx, color="blue"),
        #     # "reefBorders":
        #     LineBuilderRadial(topViewAx, color="orange")
        # ]
        self.sketches = [LineBuilderRadial(topViewAx, color=col) for name, col in zip(names, colors)]
        self.topView.addSketchs(self.sketches)

        self.profileView = SketchManagement(profileAx)
        self.profileSketch = LineBuilder1D(profileAx, color="green")
        self.profileView.addSketch(line=self.profileSketch)

        self.resistanceView = SketchManagement(resistanceAx)
        self.resistanceSketch = LineBuilder1D(resistanceAx, color="red")
        self.resistanceView.addSketch(line=self.resistanceSketch)

        self.distortionView = SketchManagement(distortionAx)
        self.distortionSketch = LineBuilder2D(distortionAx, color="blue")
        self.distortionView.addSketch(line=self.distortionSketch)

        def updateWaterLevel(val: float):
            self.waterLevel = val
            self.update()
        def updateSubsidence(val: float):
            self.subsidenceFactor = val
        def updateCoralMinMax(valMinMax: Tuple[float, float]):
            self.coralMin, self.coralMax = valMinMax

        waterLevelSlider.on_changed(updateWaterLevel)
        subsidenceSlider.on_changed(updateSubsidence)
        coralMinMaxSlider.on_changed(updateCoralMinMax)
        self.waterLevel = waterLevelSlider.val
        self.subsidenceFactor = subsidenceSlider.val
        self.coralMin, self.coralMax = coralMinMaxSlider.val

        self.sequences = [
            # (-1, -1, 0.0, 0.0),
            (-1, 0, 0.0, 0.2),
            (0, 1, 0.2, 0.4),
            (1, 2, 0.4, 0.6),
            (2, 3, 0.6, 0.8),
            (3, 1000, 0.8, 1.0)
        ]

        self.topView.onChangeEnded(self.update)
        self.profileView.onChangeEnded(self.update)
        self.distortionView.onChangeEnded(self.update)
        self.resistanceView.onChangeEnded(self.update)

        self.distortionStrength = 0.1
        self.distortionWidth = 0.5

        self.includeCoralSimulation = True

    @lru_cache
    def getSequenceID(self, pos: Vector2D):
        polar = pos.to_polar()
        minDist, minID = math.inf, -1
        maxDist, maxID = math.inf, -1
        distances = self.getDistancesOfBorders(pos)
        for i, sketch in enumerate(self.sketches):
            dist = polar.y - distances[i]
            if dist < 0:
                if abs(dist) < minDist:
                    minDist = abs(dist)
                    minID = i
            if 0 <= dist < maxDist:
                maxDist = dist
                maxID = i

        for i, s in enumerate(self.sequences):
            if s[0] == maxID and s[1] == minID:
                return i
        return -1

    def setActiveTopView(self, sketchID):
        self.topView.activate(sketchID)
        self.update()

    def update(self):
        self.updateTopViewAx()
        self.updateSideViewMarkers()
        self.updateDistortionsAx()
        self.updateResistanceViewMarkers()
        self.profileSketch.line.figure.canvas.draw()

    def updateTopViewAx(self):
        self.topView.draw(False)

    def updateDistortionsAx(self):
        # sizeX = len(distortionMaps[0][0])
        # sizeY = len(distortionMaps[0])
        sizeX = 20 # len(singleDistortionMap[0])
        sizeY = 20 # len(singleDistortionMap)
        addDistortionFromSketch(self.distortionView, self.distortionStrength, self.distortionWidth)
        self.distortionSketch.reset()
        wholeMap = resize_nearest_numpy(singleDistortionMap, (sizeX, sizeY)) # [[getDistoFromIndices(x, y) for y in range(sizeY)] for x in range(sizeX)]

        ax: plt.Axes = self.distortionView.ax
        for l in ax.lines:
            if l is not self.distortionSketch.line:
                l.remove()
        for _y in range(sizeY):
            for _x in range(sizeX):
                y, x = intsToCoords(_x, _y, sizeX, sizeY)
                # x = - x
                y = -y
                ax.plot([x, x + wholeMap[_x, _y, 0]], [y, y - wholeMap[_x, _y, 1]], c="blue")
        self.distortionView.draw(False)

    def updateSideViewMarkers(self):
        ax: plt.Axes = self.profileView.ax
        for line in ax.lines:
            line.set_data([], [])

        ax.axhline(self.waterLevel, color="blue", linestyle="--")
        ax.vlines([0.2, 0.4, 0.6, 0.8], colors=[l.color for l in self.sketches], ymin=0.0, ymax=1.0,
                  linestyles="--")
        self.profileView.draw(False)


    def updateResistanceViewMarkers(self):
        ax: plt.Axes = self.resistanceView.ax
        for line in ax.lines:
            line.set_data([], [])

        ax.vlines([0.2, 0.4, 0.6, 0.8], colors=[l.color for l in self.sketches], ymin=0.0, ymax=1.0,
                  linestyles="--")
        self.resistanceView.draw(False)

    @lru_cache()
    def getDistancesOfBorders(self, pos: Vector2D):
        polar = pos.to_polar()
        distances: List[float] = []
        for s in self.sketches:
            dist = s.getValue(polar.x)
            distances.append(dist)

        for i in range(len(distances) - 1):
            # distances[i] = min(distances[i], distances[i + 1] - 0.05)
            distances[i + 1] = max(distances[i + 1], distances[i] + 0.05)
        return distances


    def evaluateResistance(self, pos: Vector2D) -> float:
        allDistances = self.getDistancesOfBorders(pos)
        polar = pos.to_polar()
        distFromCenter = polar.y
        sequenceID = self.getSequenceID(pos)
        fullResistance = self.resistanceSketch.getCurve()
        marker1, marker2, t1, t2 = self.sequences[sequenceID]
        resistance = curves.subcurve(fullResistance, t1, t2)
        distMin, distMax = (allDistances[marker1] if marker1 != -1 else 0), (
            allDistances[marker2] if marker2 != 1000 else 1.5)
        t = (distFromCenter - distMin) / (distMax - distMin)
        return interpolateOnCurve(resistance, t).y

    def evaluateHeightAndResistance(self, pos: Vector2D) -> Tuple[float, float]:
        allDistances = self.getDistancesOfBorders(pos)
        polar = pos.to_polar()
        distFromCenter = polar.y
        sequenceID = self.getSequenceID(pos)
        fullProfile = self.profileSketch.getCurve()
        fullResistance = self.resistanceSketch.getCurve()
        marker1, marker2, t1, t2 = self.sequences[sequenceID]
        profile = curves.subcurve(fullProfile, t1, t2)
        resistance = curves.subcurve(fullResistance, t1, t2)
        distMin, distMax = (allDistances[marker1] if marker1 != -1 else 0), (
            allDistances[marker2] if marker2 != 1000 else 1.5)
        t = (distFromCenter - distMin) / (distMax - distMin)
        return interpolateOnCurve(profile, t).y, interpolateOnCurve(resistance, t).y

    def createMapsFromSketch(self, path: str = "./", filePrefix: str = "result", randomize_parameters: bool = True):
        featuresFolder = "features/"
        distoFolder = "distortions/"
        heightFolder = "heightmaps/"

        pathHeightmap = path + heightFolder
        pathFeatures = path + featuresFolder
        pathDisto = path + distoFolder
        os.makedirs(pathHeightmap, exist_ok=True)
        os.makedirs(pathFeatures, exist_ok=True)
        os.makedirs(pathDisto, exist_ok=True)

        resistanceMap = self.computeResistanceMap()
        waterLevel = self.waterLevel
        # coralMinHeight = self.coralMin
        # coralMaxHeight = self.coralMax
        heightmap, features, distortions = self.heightFeatsAndDistoFromSketches()

        def getReefDistances(features):
            _distMap = np.zeros((features.shape[0], features.shape[1]))
            # fullProfile = self.profileSketch.getCurve()
            # fullResistance = self.resistanceSketch.getCurve()
            for x in range(features.shape[0]):
                for y in range(features.shape[1]):
                    # sequenceID = features[x, y, 0] # self.getSequenceID(pos)
                    # if sequenceID != 3.0: continue #We are just interested in the reef biome (=3)
                    X = (2*y / float(features.shape[0]) - 1.0)
                    Y = -(2*x / float(features.shape[0]) - 1.0)
                    pos = Vector2D( X, Y)
                    polar = pos.to_polar()
                    distFromCenter = polar.y

                    distMin = self.sketches[2].getValue(polar.x)
                    distMax = self.sketches[3].getValue(polar.x)
                    if distMin > distMax:
                        distMin, distMax = distMax, distMin
                    t = ((distFromCenter - distMin) / (distMax - distMin) if distMax - distMin > 0.01 else -10000 if distFromCenter < distMin else 10000)
                    _distMap[x, y] = t
            return _distMap

        def distortionsToRGB(disto) :
            d = np.zeros((disto.shape[0], disto.shape[1], 3))
            d[:,:,:2] = 0.5 + np.clip(disto, -1, 1) * 0.5
            d[:,:,2] = 0.5
            return d

        # heightmap = coralize_my_island.bw2rgb(np.clip(heightmap, 0.0, 1.0))
        # plt.imsave(f"{pathHeightmap}{filePrefix}.png", heightmap)
        # features = np.array(features)
        # for x in range(features.shape[0]):
        #     for y in range(features.shape[1]):
        #         features[x, y] = valueAsHSV(features[x, y, 0], features[x, y, 1], features[x, y, 2], L=1.0 * 0.5)
        # plt.imsave(f"{pathFeatures}{filePrefix}.png", features)
        # plt.imsave(f"{pathDisto}{filePrefix}.png", distortions)
        # print(f"Heightmap saved at {pathHeightmap}{filePrefix}.png")
        # return heightmap, features, distortions

        distances = getReefDistances(features)
        # _singleDistortionMap = np.zeros((heightmap.shape[0], heightmap.shape[1], 2))
        # for _x in range(_singleDistortionMap.shape[0]):
        #     for _y in range(_singleDistortionMap.shape[1]):
        #         x, y = numpyIndicesToCoords(_x, _y, _singleDistortionMap.shape[0], _singleDistortionMap.shape[1])
        #         d = getDisto(x, y)
        #         _singleDistortionMap[_x, _y] = [d.x, d.y]
        _singleDistortionMap = singleDistortionMap * (1.0 - np.reshape(resistanceMap, (resistanceMap.shape[0], resistanceMap.shape[1], 1)))


        n = noise.perlin.SimplexNoise(1000)
        nbSamples = 5 if randomize_parameters else 1
        for i in range(nbSamples):
            if randomize_parameters:
                subsidence = getRandom(0.1, 1.0)
            else:
                subsidence = self.subsidenceFactor
            _resistanceMap = np.clip(resistanceMap + np.array([[n.noise2(x / 30, y / 30 + 100 * i) * 0.1 * n.noise2(x / 100 + 100 * i, y / 100 + 500 * i) for x in range(heightmap.shape[0])] for y in range(heightmap.shape[1])]), 0.0, 1.0)
            _heightmap = np.array(heightmap)
            _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap*0+1, iterations=100, talus_angle=1.0, erosion_factor=0.2)
            _heightmap = coralize_my_island.apply_hydraulic_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=50, water=1.0, solubility=0.99, evaporation=0.05, capacity=1.0)
            _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=100, talus_angle=0.5, erosion_factor=0.1)
            _features = np.array(features)
            for x in range(features.shape[0]):
                for y in range(features.shape[1]):
                    _features[x, y] = valueAsHSV(_features[x, y, 0], _features[x, y, 1], _features[x, y, 2], L=subsidence*.5)

            _distances = deform_image(distances, _singleDistortionMap * (1.0 - _resistanceMap).reshape((256, 256, 1)), scale = 128.0)
            terrain = _heightmap
            if self.includeCoralSimulation:
                _heightmap, terrain, corals = coralize_my_island.method2Create(_heightmap, _distances, subsidence, waterLevel)
            _heightmap = self.fractal_erosion(_heightmap, _resistanceMap, radialFreq=getRandom(3.0, 6.0), strength=getRandom(5.0, 15.0), randomnessStrength=getRandom(1.0, 3.0))
            _heightmap = coralize_my_island.apply_hydraulic_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=50, water=0.5, solubility=0.99, evaporation=0.05, capacity=1.0)
            _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=100, talus_angle=0.5, erosion_factor=0.01)
            _heightmap = coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0))
            # return _heightmap, _features, distortionsToRGB(distortions)

            # plt.imsave(f"/media/marc/Data/NN Datasets/1/result_height.png", _heightmap)
            plt.imsave(f"{pathHeightmap}{filePrefix}-{i}.png", _heightmap)
            plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
            plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortionsToRGB(distortions))
            # plt.imsave(f"{pathHeightmap}{filePrefix}-{i}.png", coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0)))
            # plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
            # plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortionsToRGB(distortions))
            if not randomize_parameters:
                plt.imsave(f"{pathHeightmap}{filePrefix}-demo-{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.png", _heightmap)
                # plt.imsave(f"/media/marc/Data/NN Datasets/1/result_height.png", _heightmap)
            print(f"Heightmap saved at {pathHeightmap}{filePrefix}-{i}.png")
        return _heightmap, _features, distortionsToRGB(distortions) #, coralize_my_island.bw2rgb(np.clip(terrain, 0.0, 1.0)), coralize_my_island.bw2rgb(np.clip(corals, 0.0, 1.0))

    def fractal_erosion(self, heightmap: np.ndarray, resistance: np.ndarray, radialFreq: float = 5.0, strength: float = 10.0, randomnessStrength: float = 2.0) -> np.ndarray:
        sX, sY = heightmap.shape[0], heightmap.shape[1] #20, 20
        disto = np.zeros((heightmap.shape[0], heightmap.shape[1], 2))
        center = Vector2D(sX / 2, sY / 2)
        n = noise.perlin.SimplexNoise(400)
        for x in range(sX):
            for y in range(sY):
                p = Vector2D(x, y) - center
                polar = p.to_polar()

                v = n.noise3(polar.x * radialFreq, x / 100.0, y / 100.0) * (strength - randomnessStrength)
                d = (p.normalized() * v + Vectors.randomVec2() * randomnessStrength) * (1 - resistance[x, y])
                disto[x, y] = [d.x, d.y]
        return deform_image(heightmap, disto)




    # def createMapsFromSketch_timing(self, path: str = "./", filePrefix: str = "result", coralMinHeight: float = 0.5,
    #                          coralMaxHeight: float = 0.6):
    #     resistanceMap = self.computeResistanceMap()
    #     waterLevel = coralMaxHeight
    #     nTrials = 10
    #     t0 = time.time()
    #     for _ in range(nTrials):
    #         heightmap, features, distortions = self.heightFeatsAndDistoFromSketches()
    #     t1 = time.time()
    #     print(f"{((t1 - t0) * 1000) / nTrials}ms")
    #     return
    #     n = noise.perlin.SimplexNoise(1000)
    #     for i in range(5):
    #         subsidence = getRandom(0.1, 1.0)
    #         _resistanceMap = np.clip(resistanceMap + np.array([[n.noise2(x / 30, y / 30 + 100 * i) * 0.1 * n.noise2(
    #             x / 100 + 100 * i, y / 100 + 500 * i) for x in range(heightmap.shape[0])] for y in
    #                                                            range(heightmap.shape[1])]), 0.0, 1.0)
    #         _heightmap = coralize_my_island.method1Create(heightmap, subsidence, waterLevel, coralMaxHeight,
    #                                                       coralMinHeight)
    #         _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap,
    #                                                               iterations=100, talus_angle=1.0,
    #                                                               erosion_factor=0.2)
    #         _heightmap = coralize_my_island.apply_hydraulic_erosion(_heightmap, resistanceMap=_resistanceMap,
    #                                                                 iterations=50, water=1.0, solubility=0.99,
    #                                                                 evaporation=0.05, capacity=1.0)
    #         _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap,
    #                                                               iterations=100, talus_angle=0.5,
    #                                                               erosion_factor=0.1)
    #
    #         _features = np.array(features)
    #         for x in range(features.shape[0]):
    #             for y in range(features.shape[1]):
    #                 _features[x, y] = valueAsHSV(_features[x, y, 0], _features[x, y, 1], _features[x, y, 2],
    #                                              L=subsidence * 0.5)
    #         plt.imsave(f"{pathHeightmap}{filePrefix}-{i}.png",
    #                    coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0)))
    #         plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
    #         plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortions)
    #         print(f"Heightmap saved at {pathHeightmap}{filePrefix}-{i}.png")

    def computeResistanceMap(self) -> np.ndarray:
        dims = outputImageDims
        resistances: np.ndarray = np.zeros((dims[0], dims[1]))
        disto_resistances: np.ndarray = np.zeros((dims[0], dims[1]))

        for _y in range(dims[1]):
            for _x in range(dims[0]):
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                pos = Vector2D(x, y)
                resistance = self.evaluateResistance(pos)
                resistances[_x, _y] = min(max(resistance, 0.0), 1.0)

        return deform_image(resistances, singleDistortionMap)
        # distortion part :
        for _y in range(dims[1]):
            for _x in range(dims[0]):
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                newX, newY = evaluatePosAfterDistortion(x, y, float(1.0 - resistances[_x, _y]))
                disto_resistances[_x, _y] = numpyBilinearInterpolation(resistances, newX, newY)
        return disto_resistances


    def heightFeatsAndDistoFromSketches(self, featureColorVivid: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        dims = outputImageDims
        heights: np.ndarray = np.zeros((dims[0], dims[1]))
        features: np.ndarray = np.zeros((dims[0], dims[1], 3))
        resistances: np.ndarray = np.zeros((dims[0], dims[1], 1))

        disto_heights: np.ndarray = np.zeros((dims[0], dims[1]))
        disto_features: np.ndarray = np.zeros((dims[0], dims[1], 3))
        disto_forces = np.zeros((dims[0], dims[1], 3))

        n = noise.perlin.SimplexNoise()
        for _x in range(dims[0]):
            for _y in range(dims[1]):
                noiseVal = 1.0 + (n.noise2(_x / 50, _y / 50) / 10.0)
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                pos = Vector2D(x, y)
                h, resistance = self.evaluateHeightAndResistance(pos)
                height = h * noiseVal
                heights[_x, _y] = min(max(height, 0.0), 1.0)
                features[_x, _y] = np.array([self.getSequenceID(pos), float(-1), float(len(self.sequences))])
                resistances[_x, _y] = min(max(resistance, 0.0), 1.0)

        deformation = singleDistortionMap * (1.0 - resistances)
        return deform_image(heights, deformation, 128), deform_image(features, deformation, 128, interpolation="nearest"), singleDistortionMap
        # distortion part :
        for _x in range(dims[0]):
            for _y in range(dims[1]):
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                newX, newY = evaluatePosAfterDistortion(x, y, float(1.0 - resistances[_x, _y]))
                diffX = x - newX
                diffY = y - newY
                disto_forces[_x, _y] = [diffX * .25 + .5, diffY * .25 + .5, 0.5]
                disto_heights[_x, _y] = numpyBilinearInterpolation(heights, newX, newY)
                disto_features[_x, _y] = (float(numpyNearestNeighbor(features, newX, newY)), float(-1), float(len(self.sequences))) # Just add the info, the "toHSV" will be used later
        return disto_heights, disto_features, np.clip(disto_forces, 0.0, 1.0)


def addDistortionFromCurve(curve: List[Vector2D], distortionStrength: float, lineWidth: float = 0.5) -> None:
    global singleDistortionMap
    distortions = singleDistortionMap
    sizeY, sizeX = distortions.shape[:2]

    # Generate full coordinate grid, normalized to [-1, 1] (or whatever `intsToCoords` does)
    xs = np.linspace(1, -1, sizeX)
    ys = np.linspace(1, -1, sizeY)
    grid_y, grid_x = np.meshgrid(ys, xs, indexing='ij')
    positions = np.stack((grid_x, grid_y), axis=-1)  # shape: (H, W, 2)
    pos_flat = positions.reshape(-1, 2)  # shape: (N, 2)

    # Prepare data
    num_points = pos_flat.shape[0]
    best_distances = np.full(num_points, lineWidth**2)
    best_directions = np.zeros((num_points, 2))

    # Evaluate distances from all points to each curve segment
    for i in range(len(curve) - 1):
        p1 = np.array([-curve[i].x, curve[i].y])
        p2 = np.array([-curve[i + 1].x, curve[i + 1].y])
        if np.allclose(p1, p2):
            continue

        # Compute vectorized distances
        v = p2 - p1
        v_len2 = np.sum(v**2)
        w = pos_flat - p1
        t = np.clip(np.sum(w * v, axis=1) / v_len2, 0.0, 1.0)
        proj = p1 + np.outer(t, v)
        dist2 = np.sum((pos_flat - proj)**2, axis=1)

        # Use masking to update only closer distances
        mask = dist2 < best_distances
        best_distances[mask] = dist2[mask]

        direction = v / np.linalg.norm(v)
        best_directions[mask] = direction * np.array([-1, -1])  # same for all, so no need to normalize again

    # Apply distortion where distances are within threshold
    within_range = best_distances < lineWidth**2
    influence = wyvill(np.sqrt(best_distances[within_range]) / lineWidth)
    motions = best_directions[within_range] * distortionStrength * influence[:, None]

    # Reshape and accumulate
    flat_distortions = distortions.reshape(-1, 2)
    flat_distortions[within_range] += motions
    singleDistortionMap = flat_distortions.reshape(distortions.shape)

# def addDistortionFromCurve(curve: List[Vector2D], distortionStrength: float, lineWidth: float = 0.5) -> None:
#     global singleDistortionMap
#     distortions = singleDistortionMap
#     sizeX, sizeY = len(distortions[0]), len(distortions)
#     for _x in range(sizeX):
#         for _y in range(sizeY):
#             x, y = intsToCoords(_x, _y, sizeX, sizeY)
#             closestLineIndex: int = -1
#             closestDistance: float = lineWidth
#             pos = Vector2D(x, y)
#             for i in range(len(curve) - 1):
#                 if curve[i] == curve[i + 1]: continue
#                 # Check original and wrapped distances
#                 candidates = [
#                     pos,  # Original position
#                     Vector2D(pos.x - 2.0, pos.y),  # Wrapped left
#                     Vector2D(pos.x + 2.0, pos.y),  # Wrapped right
#                     Vector2D(pos.x, pos.y - 2.0),  # Wrapped up
#                     Vector2D(pos.x, pos.y + 2.0),  # Wrapped down
#                     Vector2D(pos.x - 2.0, pos.y - 2.0),  # Diagonal wrap: top-left
#                     Vector2D(pos.x + 2.0, pos.y - 2.0),  # Diagonal wrap: top-right
#                     Vector2D(pos.x - 2.0, pos.y + 2.0),  # Diagonal wrap: bottom-left
#                     Vector2D(pos.x + 2.0, pos.y + 2.0),  # Diagonal wrap: bottom-right
#                 ]
#
#                 for candidate in candidates:
#                     distToLine = Vectors.distance2ToLine(candidate, curve[i], curve[i + 1])
#                     if distToLine < closestDistance:
#                         closestDistance = distToLine
#                         closestLineIndex = i
#                 #
#                 # distToLine = Vectors.distanceToLine(pos, curve[i], curve[i + 1])
#                 # if distToLine < closestDistance:
#                 #     closestDistance = distToLine
#                 #     closestLineIndex = i
#             if closestLineIndex > -1:
#                 closestDistance = math.sqrt(closestDistance)
#                 distToLine = (closestDistance - 0) / (lineWidth - 0)
#                 mouseMotion = (curve[closestLineIndex + 1] - curve[closestLineIndex]).normalize() * distortionStrength
#                 distortions[_x, _y] += np.array([mouseMotion.x, mouseMotion.y]) * wyvill(distToLine)


def addDistortionFromSketch(distortionSketcher: SketchManagement, strength: float, width: float) -> None:
    mousePath = distortionSketcher.lineBuilders[0].getCurve()
    addDistortionFromCurve(mousePath, strength, width)
    distortionSketcher.lineBuilders[0].reset()


def updateResultsFigure(images: List[np.ndarray], _axes: List[plt.Axes]):
    global fig2
    axes = fig2.axes
    for i, img in enumerate(images):
        axes[i].imshow(img)
    fig2.canvas.draw()

def randomizeUserCircle(previousOutlines, randomness, scaling: float = 1.0):
    n = noise.perlin.SimplexNoise()
    seed = getRandom(0, 1000)
    newOutlines = []
    for i in range(len(previousOutlines)):
        t = math.tau * i / float(len(previousOutlines) + 1)
        l = previousOutlines[i].x
        p = Vector2D(l * (1 + n.noise2(math.cos(t), math.sin(t) + seed) * (randomness - 1)) * scaling, previousOutlines[i].y)
        newOutlines.append(p)
    return newOutlines


def createDatasetOfRandomIslands(islandSketch: IslandSketch, nbSamples: int = 1000, use_user_input=False):
    global singleDistortionMap
    radiusRandomMin, radiusRandomMax = 0.8, 1.1
    nbCurvesMin, nbCurvesMax = 2, 6
    distoStrengthMin, distoStrengthMax = 0.05, 0.3
    distoWidthMin, distoWidthMax = 0.1, 0.3
    profileRandomMin, profileRandomMax = 0.0, 0.2
    resistanceRandomMin, resistanceRandomMax = 0.0, 0.2

    original_islandBorders = islandSketch.sketches[0].getCurve()
    original_beachBorders  = islandSketch.sketches[1].getCurve()
    original_lagoonBorders = islandSketch.sketches[2].getCurve()
    original_reefBorders   = islandSketch.sketches[3].getCurve()

    original_profile = islandSketch.profileSketch.getCurve()
    original_resistance = islandSketch.resistanceSketch.getCurve()

    for iSample in range(nbSamples):
        if os.path.exists(f"{dataset_path}heightmaps/{iSample}.png"):
            continue
        profileRandomness = getRandom(profileRandomMin, profileRandomMax)
        resistanceRandomness = getRandom(resistanceRandomMin, resistanceRandomMax)
        n = noise.perlin.SimplexNoise(1000)
        singleDistortionMap = initialDistoMap(outputImageDims[0], outputImageDims[1]) #initialDistoMap(30, 30)
        for _ in range(int(getRandom(nbCurvesMin, nbCurvesMax))):
            addDistortionFromCurve(randomDistortionCurve(), getRandom(distoStrengthMin, distoStrengthMax), getRandom(distoWidthMin, distoWidthMax))
        radiusFactor = getRandom(radiusRandomMin, radiusRandomMax)

        if use_user_input:
            islandBorders = randomizeUserCircle(original_islandBorders, 1.2, radiusFactor)
            beachBorders = randomizeUserCircle(original_beachBorders, 1.2, radiusFactor)
            lagoonBorders  = randomizeUserCircle(original_lagoonBorders, 1.2, radiusFactor)
            reefBorders = randomizeUserCircle(original_reefBorders, 1.2, radiusFactor)
        else:
            islandBorders = centeredCircle(0.4 * radiusFactor, 1.2)
            beachBorders = centeredCircle(0.45 * radiusFactor, 1.2)
            lagoonBorders = centeredCircle(0.7 * radiusFactor, 1.2)
            reefBorders = centeredCircle(0.8 * radiusFactor, 1.2)

        # islandTranslate = Vectors.randomVec2() * 0.2
        # islandBorders = [p + islandTranslate for p in islandBorders]
        # beachTranslate = Vectors.randomVec2() * 0.2
        # beachBorders = [p + beachTranslate for p in beachBorders]

        for i in range(len(islandBorders)):
            if beachBorders[i].norm2() < islandBorders[i].norm2():
                beachBorders[i] = islandBorders[i] + islandBorders[i].normalized() * 0.05
            if lagoonBorders[i].norm2() < beachBorders[i].norm2():
                lagoonBorders[i] = beachBorders[i] + beachBorders[i].normalized() * 0.05
            if reefBorders[i].norm2() < lagoonBorders[i].norm2():
                reefBorders[i] = lagoonBorders[i] + lagoonBorders[i].normalized() * 0.05

        islandSketch.sketches[0].setCurve(islandBorders)
        islandSketch.sketches[1].setCurve(beachBorders)
        islandSketch.sketches[2].setCurve(lagoonBorders)
        islandSketch.sketches[3].setCurve(reefBorders)

        _randomProfileCurve = [0.9627757352941178, 0.897671568627451, 0.8425245098039216, 0.775888480392157, 0.6893382352941178, 0.6479779411764706, 0.5974264705882353, 0.5606617647058822, 0.5392156862745099, 0.508578431372549, 0.4986213235294116, 0.48406862745098034, 0.4810049019607843, 0.4779411764705883, 0.4756433823529411, 0.47411151960784303, 0.4595588235294117, 0.4564950980392156, 0.44117647058823517, 0.431985294117647, 0.42585784313725483, 0.315563725490196, 0.2481617647058823, 0.015318627450980338, 0.015318627450980338, 0.015, 0.01, 0.001, 0.0, 0.0] #resizeArray([1.0, 0.0], 10)
        _randomResistanceCurve = [0.7429534313725492, 0.7153799019607845, 0.6755514705882353, 0.626531862745098, 0.5821078431372548, 0.5238970588235294, 0.42585784313725483, 0.27037377450980393, 0.1470588235294117, 0.0850183823529411, 0.036764705882352866, 0.03063725490196073, 0.027573529411764608, 0.027573529411764608, 0.027573529411764608, 0.039828431372548934, 0.06127450980392152, 0.08272058823529405, 0.10569852941176466, 0.15624999999999994, 0.22212009803921565, 0.3117340686274509, 0.4840686274509804, 0.6441482843137254, 0.7467830882352942, 0.8620689655172414, 0.04901960784313719, 0.0337009803921568, 0.021446078431372473, 0.003063725490196012]
        if use_user_input:
            _randomProfileCurve = [p.y for p in original_profile]
            _randomResistanceCurve = [p.y for p in original_resistance]
        randomProfileCurve = []
        randomResistanceCurve = []
        for i in range(len(_randomProfileCurve)):
            randomProfileCurve.append(_randomProfileCurve[i] * (1.0 + n.noise2(i / 15, iSample * 10000) * profileRandomness))
            randomResistanceCurve.append(_randomResistanceCurve[i] * (1.0 + n.noise2(i / 15, (iSample + 100) * 10000) * resistanceRandomness))
            # print(f"Diff {i}: {randomProfileCurve[-1] - _randomProfileCurve[i]} (rand: {profileRandomness}) => from {_randomProfileCurve[i]} to {randomProfileCurve[-1] }")
        islandSketch.profileSketch.setCurve(randomProfileCurve)
        islandSketch.resistanceSketch.setCurve(randomResistanceCurve)


        h, f, d = islandSketch.createMapsFromSketch(path=f"{dataset_path}", filePrefix=str(iSample)) # , coralMinHeight=coralMinHeight, coralMaxHeight=coralMaxHeight)

    if use_user_input:
        islandSketch.sketches[0].setCurve(original_islandBorders)
        islandSketch.sketches[1].setCurve(original_beachBorders)
        islandSketch.sketches[2].setCurve(original_lagoonBorders)
        islandSketch.sketches[3].setCurve(original_reefBorders)
        islandSketch.profileSketch.setCurve(original_profile)
        islandSketch.resistanceSketch.setCurve(original_resistance)
        islandSketch.update()




def main():
    # random.seed(42)
    global fig2
    global axesResults
    # global distortionMaps
    global singleDistortionMap
    # fig2 = plt.figure()
    # axHeight, axFeatures, axDisto = fig2.subplots(1, 3)
    fig2, axesResults = plt.subplots(1, 3, squeeze=True)
    axHeight, axFeatures, axDisto = axesResults
    axHeight.set_title("Height field")
    axFeatures.set_title("Label map")
    axDisto.set_title("Distortion map")
    for ax in axesResults:
        ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)

    singleDistortionMap = initialDistoMap(outputImageDims[0], outputImageDims[1]) # initialDistoMap(20, 20)

    sketch_names = ["Island", "Beach", "Lagoon", "Reef"]
    sketch_colors = ["green", "yellow", "blue", "orange"]

    # Creates a top view for island sketching and a side view for profile editing
    fig, axes = plt.subplots(1, 4, squeeze=True)
    topViewAx: plt.Axes = axes[0]
    sideViewAx: plt.Axes = axes[1]
    distortionsAx: plt.Axes = axes[2]
    resistanceAx: plt.Axes = axes[3]
    fig.subplots_adjust(bottom=0.2)

    topViewAx.set_title('Top view')
    topViewAx.set_xlim(-1, 1)
    topViewAx.set_ylim(-1, 1)
    topViewAx.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    sideViewAx.set_title('Side view')
    sideViewAx.set_xlim(0, 1)
    sideViewAx.set_ylim(0, 1)
    distortionsAx.set_title('Distortions')
    distortionsAx.set_xlim(-1, 1)
    distortionsAx.set_ylim(-1, 1)
    distortionsAx.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    resistanceAx.set_title('Resistance')
    resistanceAx.set_xlim(0, 1)
    resistanceAx.set_ylim(0, 1)

    axWater = fig.add_axes((0.92, 0.2, 0.0225, 0.68))
    water_slider = Slider(
        ax=axWater,
        label="Water\nlevel",
        valmin=0,
        valmax=1,
        valinit=0.4,
        orientation="vertical"
    )

    axCoral = fig.add_axes((95, 0.2, 0.0225, 0.68)) #(0.95, 0.2, 0.0225, 0.68))
    coral_slider = RangeSlider(
        ax=axCoral,
        label="Corals",
        valmin=0,
        valmax=1,
        valinit=(0.2, 0.3),
        orientation="vertical"
    )
    axSubsid = fig.add_axes((0.97, 0.2, 0.0225, 0.68))
    subsid_slider = Slider(
        ax=axSubsid,
        label="Subsidence\nfactor",
        valmin=0,
        valmax=1,
        valinit=0.8,
        orientation="vertical"
    )

    # islandSketch = IslandSketch(topViewAx, sideViewAx, topViewAx, sideViewAx, water_slider, coral_slider, subsid_slider, sketch_names, sketch_colors)
    islandSketch = IslandSketch(topViewAx, sideViewAx, distortionsAx, resistanceAx, water_slider, coral_slider, subsid_slider, sketch_names, sketch_colors)


    buttons: List[Button] = []
    button_size = 0.5 / len(sketch_names)
    button_height = 0.5 / len(sketch_names)
    for i_sketch in range(len(sketch_names)):
        # islandSketches.addSketch(color = sketch_colors[i_sketch], sketch_type = sketch_types[i_sketch])
        sketchID = i_sketch
        # ax_button = fig.add_axes((0.5 + button_size * sketchID, 0.05, button_size - 0.01, 0.075))
        ax_button = fig.add_axes((0.005, 0.25 + (button_height + 0.01) * i_sketch, 0.05, button_height))
        buttons.append(Button(ax_button, sketch_names[sketchID]))
        def activationFunction(id):
            def action(event):
                islandSketch.setActiveTopView(id)
            return action
        buttons[-1].on_clicked(activationFunction(sketchID))

    # Add button for generating a distance map / heightmap
    ax_button = fig.add_axes((0.0, 0.05, 0.15, 0.075))
    distance_button = Button(ax_button, "Gen height map")
    ax_button_dataset = fig.add_axes((0.16, 0.05, 0.15, 0.075))
    dataset_button = Button(ax_button_dataset, "Gen dataset")

    def resetWind():
        global singleDistortionMap
        singleDistortionMap = initialDistoMap(outputImageDims[0], outputImageDims[1])
        islandSketch.update()

    def setWindWidth(newWidth: float):
        islandSketch.distortionWidth = newWidth
    def setWindStrength(newStrength: float):
        islandSketch.distortionStrength = newStrength

    ax_wind_reset = fig.add_axes((0.5, 0.05, 0.08, 0.075))
    wind_reset_button = Button(ax_wind_reset, "Reset")
    wind_reset_button.on_clicked(lambda e: resetWind())

    ax_wind_width = fig.add_axes((0.6, 0.025, 0.08, 0.075))
    wind_width_slider = Slider(ax=ax_wind_width, label="Width", valmin=0.0, valmax=1.0, valinit=islandSketch.distortionWidth, orientation="horizontal")
    wind_width_slider.on_changed(setWindWidth)
    wind_width_slider.valtext.set_visible(False)
    wind_width_slider.label.set_visible(False)
    fig.text(0.64, 0.04, "Width", ha='center', va='top')

    ax_wind_strength = fig.add_axes((0.6, 0.1, 0.08, 0.075))
    wind_strength_slider = Slider(ax=ax_wind_strength, label="Strength", valmin=0.0, valmax=0.5, valinit=islandSketch.distortionStrength, orientation="horizontal")
    wind_strength_slider.on_changed(setWindStrength)
    wind_strength_slider.valtext.set_visible(False)
    wind_strength_slider.label.set_visible(False)
    fig.text(0.64, 0.125, "Strength", ha='center', va='top')

    def setCoralSimulationActive(checkboxes):
        islandSketch.includeCoralSimulation = checkboxes.get_status()[0]

    ax_coral_checkbox = fig.add_axes((0.8, 0.04, 0.15, 0.075))
    coral_checkbox = CheckButtons(ax=ax_coral_checkbox, labels=["Include coral simulation"], actives=[islandSketch.includeCoralSimulation])
    coral_checkbox.on_clicked(lambda l: setCoralSimulationActive(coral_checkbox))
    # coral_checkbox.on_changed(setWindStrength)
    # coral_checkbox.valtext.set_visible(False)
    # coral_checkbox.label.set_visible(False)
    # fig.text(0.64, 0.125, "Strength", ha='center', va='top')

    def genIsland():
        # sequences = getSequences(islandSketch.profileSketch)
        # lagoonSequence = sequences[2][2]
        # reefSequence = sequences[3][2]
        # coralMaxHeight = lagoonSequence[len(lagoonSequence) // 2]
        # coralMinHeight = lagoonSequence[-1] #reefSequence[len(reefSequence) // 2]
        # subsidence = 0.8
        # print(coralMinHeight, coralMaxHeight, subsidence, lagoonSequence + reefSequence)
        # print("Profile:", [p.y for p in islandSketch.profileSketch.getCurve()])
        # print("Resistance", [p.y for p in islandSketch.resistanceSketch.getCurve()])
        h, f, d = islandSketch.createMapsFromSketch(path="test_island_heightmap/", filePrefix="result_height", randomize_parameters=False)
        axesResults[0].imshow(h)
        axesResults[1].imshow(f)
        axesResults[2].imshow(d)
        # for ax in axesResults:
        #     ax.update()
        fig2.canvas.draw()
        fig2.canvas.flush_events()

    def genDataset():
        print(f"Generating dataset at {dataset_path}heightmaps/")
        createDatasetOfRandomIslands(islandSketch, 2, use_user_input=True)


    distance_button.on_clicked(lambda e: genIsland()) # genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut))
    dataset_button.on_clicked(lambda e: genDataset())
    # Add button for splitting sequences
    # ax_button = fig.add_axes((0.1, 0.05, 0.09, 0.075))
    # splitting_button = Button(ax_button, "Split profile")
    # splitting_button.on_clicked(lambda e: splitProfileOnMarkers(profileSketching.lineBuilders[0], islandSketches, sliceCut))

    _randomProfileCurve = [0.9627757352941178, 0.897671568627451, 0.8425245098039216, 0.775888480392157,
                           0.6893382352941178, 0.6479779411764706, 0.5974264705882353, 0.5606617647058822,
                           0.5392156862745099, 0.508578431372549, 0.4986213235294116, 0.48406862745098034,
                           0.4810049019607843, 0.4779411764705883, 0.4756433823529411, 0.47411151960784303,
                           0.4595588235294117, 0.4564950980392156, 0.44117647058823517, 0.431985294117647,
                           0.42585784313725483, 0.315563725490196, 0.2481617647058823, 0.015318627450980338,
                           0.015318627450980338, 0.015, 0.01, 0.001, 0.0, 0.0]  # resizeArray([1.0, 0.0], 10)
    _randomResistanceCurve = [0.7429534313725492, 0.7153799019607845, 0.6755514705882353, 0.626531862745098,
                              0.5821078431372548, 0.5238970588235294, 0.42585784313725483, 0.27037377450980393,
                              0.1470588235294117, 0.0850183823529411, 0.036764705882352866, 0.03063725490196073,
                              0.027573529411764608, 0.027573529411764608, 0.027573529411764608, 0.039828431372548934,
                              0.06127450980392152, 0.08272058823529405, 0.10569852941176466, 0.15624999999999994,
                              0.22212009803921565, 0.3117340686274509, 0.4840686274509804, 0.6441482843137254,
                              0.7467830882352942, 0.8620689655172414, 0.04901960784313719, 0.0337009803921568,
                              0.021446078431372473, 0.003063725490196012]
    islandSketch.profileSketch.setCurve(_randomProfileCurve)
    islandSketch.resistanceSketch.setCurve(_randomResistanceCurve)
    # for i in range(len(sketch_names)):
    #     islandSketch.sketches[i].setCurve(centeredCircle(float(i + 1) / float(len(sketch_names) + 1), 1.1))
    islandSketch.sketches[0].setCurve(centeredCircle(0.25, 1.0))
    islandSketch.sketches[1].setCurve(centeredCircle(0.50, 1.0))
    islandSketch.sketches[2].setCurve(centeredCircle(0.70, 1.0))
    islandSketch.sketches[3].setCurve(centeredCircle(0.80, 1.0))
    # createDatasetOfRandomIslands(islandSketch, 100)
    islandSketch.update()
    # with cProfile.Profile() as profiler:
    #     genIsland()
    # stats = pstats.Stats(profiler)
    # stats.sort_stats(pstats.SortKey.TIME)  # Sorts by time spent in each function
    # stats.print_stats(20)  # Prints top 20 lines of profiling results
    plt.show()

    # import dataAugmentationIslandGeneration
    # dataAugmentationIslandGeneration.main()
    # print("All done!!!")



if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Program for the coral reef island dataset generation, either interactive or automatic")
    # parser.add_argument('path', nargs='?', default='default/path/', help='Path to save the dataset')
    # parser.add_argument('--auto', action='store_true', help='Creates dataset with random layout')
    #
    # args = parser.parse_args()
    #
    # my_path = args.path
    # my_variable = args.test
    main()
