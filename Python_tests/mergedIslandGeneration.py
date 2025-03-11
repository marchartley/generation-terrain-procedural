# Imports
import numpy as np
import matplotlib.pyplot as plt
import random, os, math
from functools import lru_cache
from matplotlib.widgets import Button
import curves
from Python_tests.coralize_my_island import bw2rgb
from Python_tests.heightmap_from_skeleton import readImage
from Python_tests.islandFromSketch import heightmapAndFeatureMapFromSketches
from Python_tests.sketch_maker import resizeArray
from Vectors import Vector2D, Vector3D, line_intersection
import Vectors
from typing import Tuple, List, Callable, Union, Optional, Any
from sketch_maker import LineBuilder, LineBuilder1D, LineBuilder2D, LineBuilderRadial, SketchManagement
from noise import perlin
import noise
from noise.perlin import SimplexNoise
import time
import coralize_my_island
from scipy.ndimage import map_coordinates
from matplotlib.widgets import RangeSlider, Slider
import cProfile
import pstats

# Global Variables
outputImageDims = [256, 256]
# distortionMaps: List[List[List[Vector2D]]] = []
singleDistortionMap: List[List[Vector2D]] = []
fig, fig2 = None, None
dataset_path = "./test_synthetic_terrains_dataset/"


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

def initialDistoMap(sizeX: int, sizeY: int) -> List[List[Vector2D]]:
    return [[Vector2D(0.1, 0) for _ in range(sizeY)] for _ in range(sizeX)]

def deform_image(image: np.ndarray, vector_field: np.ndarray) -> np.ndarray:
    # Check if image and vector_field have the same dimensions
    assert image.shape == vector_field.shape[:2], "Image and vector field must have the same dimensions"

    # Create a grid of coordinates in the original image
    coords = np.indices(image.shape)

    # Adjust the coordinates based on the vector field
    # Note that we need to invert the vector direction for proper mapping
    displaced_coords = np.array([
        coords[0] - vector_field[..., 1],  # y-coordinates adjust
        coords[1] - vector_field[..., 0]   # x-coordinates adjust
    ])

    # Perform the interpolation
    deformed_image = map_coordinates(image, displaced_coords, order=1, mode='reflect')

    return deformed_image

def getDisto(x: float, y: float, factor: float = 1.0) -> Vector2D:
    # distortion = Vector2D(x, y)
    disto1 = bilinearInterpolation(singleDistortionMap, x, y) * factor
    disto2 = bilinearInterpolation(singleDistortionMap, x - disto1.x, y - disto1.y) * factor
    return (disto1 + disto2) / 2 #Vector2D(disto.x, disto.y)
    # for disto in reversed(distortionMaps):
    #     d = bilinearInterpolation(disto, distortion.x, distortion.y)
    #     d.y *= -1
    #     distortion += d
    # return distortion - Vector2D(x, y)

def getDistoFromIndices(_x: int, _y: int) -> Vector2D:
    return singleDistortionMap[_x][_y] # Possibly in the wrong order
    # sizeX, sizeY = len(distortionMaps[0][0]), len(distortionMaps[0])
    # vec: Vector2D = Vector2D(0, 0)
    # if _x < 0 or _x >= sizeX or _y < 0 or _y >= sizeY:
    #     return vec
    # for distoMap in distortionMaps:
    #     vec += distoMap[_x][_y]
    # return vec

def evaluatePosAfterDistortion(prevX: float, prevY: float, factor: float = 1.0) -> Tuple[float, float]:
    newPos = Vector2D(prevX, prevY)
    nbSteps = 10
    _factor = (factor / nbSteps)
    for _ in range(nbSteps):
        distor = getDisto(newPos.x, newPos.y, _factor)
        newPos -= distor
    return newPos.x, newPos.y
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

def bilinearInterpolation(arr: List[List[Vector2D]], x: float, y: float) -> Vector2D:
    sizeX, sizeY = len(arr[0]), len(arr)
    _x, _y = coordsToInts(x, y, sizeX, sizeY)
    __x, __y = coordsToFloats(x, y, sizeX, sizeY)
    dx, dy = __x - _x, __y - _y
    d00 = arr[_x  ][_y  ] if _x >= 0 and _x < sizeX and _y >= 0 and _y < sizeY else Vector2D(0, 0)
    d01 = arr[_x  ][_y+1] if _x >= 0 and _x < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else Vector2D(0, 0)
    d10 = arr[_x+1][_y  ] if (_x+1) >= 0 and (_x+1) < sizeX and _y >= 0 and _y < sizeY else Vector2D(0, 0)
    d11 = arr[_x+1][_y+1] if (_x+1) >= 0 and (_x+1) < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else Vector2D(0, 0)
    value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
    return value

def clamp(x:float, mini: float, maxi: float) -> float:
    return mini if x <= mini else maxi if x >= maxi else x

def valueAsHSV(value: float, mini: float, maxi: float, L: float = 0.5, S: float = 1.0) -> Tuple[float, float, float]:
    H = (360.0 * (value - mini) / (maxi - mini)) / 60
    # L = 0.5
    # S = 1.0

    C = (1 - abs(2 * L - 1)) * S
    try:
        X = C * (1 - abs(math.fmod(H, 2) - 1))
    except Exception as e:
        a = 0
        raise e

    R1, G1, B1 = (C, X, 0) if 0 <= H < 1 else (X, C, 0) if 1 <= H < 2 else (0, C, X) if 2 <= H < 3 else (0, X, C) if 3 <= H < 4 else (X, 0, C) if 4 <= H < 5 else (C, 0, X)
    m = L - C * .5
    R, G, B = R1 + m, G1 + m, B1 + m
    return clamp(R,0,1), clamp(G,0,1), clamp(B,0,1)


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


def centeredCircle(radius: float, randomness: float = 1.2) -> List[Vector2D]:
    perlin_noise = perlin.SimplexNoise(10)
    points: List[Vector2D] = []
    nbPoints = 30
    for i in range(nbPoints):
        angle = i * math.tau / (nbPoints - 1)
        vertexUnitPos = Vector2D(math.cos(angle), math.sin(angle))
        noiseValue = 1.0 + perlin_noise.noise2(vertexUnitPos.x, vertexUnitPos.y) * (randomness - 1)
        points.append(vertexUnitPos * (radius * noiseValue))
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
        sizeX = len(singleDistortionMap[0])
        sizeY = len(singleDistortionMap)
        addDistortionFromSketch(self.distortionView)
        self.distortionSketch.reset()
        wholeMap = [[getDistoFromIndices(x, y) for y in range(sizeY)] for x in range(sizeX)]

        ax: plt.Axes = self.distortionView.ax
        for l in ax.lines:
            if l is not self.distortionSketch.line:
                l.remove()
        for _y in range(sizeY):
            for _x in range(sizeX):
                x, y = intsToCoords(_x, _y, sizeX, sizeY)
                ax.plot([x, x + wholeMap[_x][_y].x], [y, y + wholeMap[_x][_y].y], c="blue")
        self.distortionView.draw(False)

    def updateSideViewMarkers(self):
        ax: plt.Axes = self.profileView.ax
        for line in ax.lines:
            line.set_data([], [])

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

    def createMapsFromSketch(self, path: str = "./", filePrefix: str = "result"):
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
        coralMinHeight = self.coralMin
        coralMaxHeight = self.coralMax
        heightmap, features, distortions = self.heightFeatsAndDistoFromSketches()
        n = noise.perlin.SimplexNoise(1000)
        for i in range(1):
            # subsidence = 0.8 # getRandom(0.1, 1.0)
            subsidence = self.subsidenceFactor
            _resistanceMap = np.clip(resistanceMap + np.array([[n.noise2(x / 30, y / 30 + 100 * i) * 0.1 * n.noise2(x / 100 + 100 * i, y / 100 + 500 * i) for x in range(heightmap.shape[0])] for y in range(heightmap.shape[1])]), 0.0, 1.0)
            _heightmap = coralize_my_island.method1Create(heightmap, subsidence, waterLevel, coralMaxHeight, coralMinHeight)
            # _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=100, talus_angle=1.0, erosion_factor=0.2)
            # _heightmap = coralize_my_island.apply_hydraulic_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=50, water=1.0, solubility=0.99, evaporation=0.05, capacity=1.0)
            # _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap, iterations=100, talus_angle=0.5, erosion_factor=0.1)
            _heightmap = coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0))
            _features = np.array(features)
            for x in range(features.shape[0]):
                for y in range(features.shape[1]):
                    _features[x, y] = valueAsHSV(_features[x, y, 0], _features[x, y, 1], _features[x, y, 2], L=subsidence * 0.5)
            plt.imsave(f"/media/marc/Data/NN Datasets/1/result_height.png", _heightmap)
            plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
            plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortions)
            # plt.imsave(f"{pathHeightmap}{filePrefix}-{i}.png", coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0)))
            # plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
            # plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortions)
            print(f"Heightmap saved at {pathHeightmap}{filePrefix}-{i}.png")
            return _heightmap, _features, distortions

    def createMapsFromSketch_timing(self, path: str = "./", filePrefix: str = "result", coralMinHeight: float = 0.5,
                             coralMaxHeight: float = 0.6):
        resistanceMap = self.computeResistanceMap()
        waterLevel = coralMaxHeight
        nTrials = 10
        t0 = time.time()
        for _ in range(nTrials):
            heightmap, features, distortions = self.heightFeatsAndDistoFromSketches()
        t1 = time.time()
        print(f"{((t1 - t0) * 1000) / nTrials}ms")
        return
        n = noise.perlin.SimplexNoise(1000)
        for i in range(5):
            subsidence = getRandom(0.1, 1.0)
            _resistanceMap = np.clip(resistanceMap + np.array([[n.noise2(x / 30, y / 30 + 100 * i) * 0.1 * n.noise2(
                x / 100 + 100 * i, y / 100 + 500 * i) for x in range(heightmap.shape[0])] for y in
                                                               range(heightmap.shape[1])]), 0.0, 1.0)
            _heightmap = coralize_my_island.method1Create(heightmap, subsidence, waterLevel, coralMaxHeight,
                                                          coralMinHeight)
            _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap,
                                                                  iterations=100, talus_angle=1.0,
                                                                  erosion_factor=0.2)
            _heightmap = coralize_my_island.apply_hydraulic_erosion(_heightmap, resistanceMap=_resistanceMap,
                                                                    iterations=50, water=1.0, solubility=0.99,
                                                                    evaporation=0.05, capacity=1.0)
            _heightmap = coralize_my_island.apply_thermal_erosion(_heightmap, resistanceMap=_resistanceMap,
                                                                  iterations=100, talus_angle=0.5,
                                                                  erosion_factor=0.1)

            _features = np.array(features)
            for x in range(features.shape[0]):
                for y in range(features.shape[1]):
                    _features[x, y] = valueAsHSV(_features[x, y, 0], _features[x, y, 1], _features[x, y, 2],
                                                 L=subsidence * 0.5)
            plt.imsave(f"{pathHeightmap}{filePrefix}-{i}.png",
                       coralize_my_island.bw2rgb(np.clip(_heightmap, 0.0, 1.0)))
            plt.imsave(f"{pathFeatures}{filePrefix}-{i}.png", _features)
            plt.imsave(f"{pathDisto}{filePrefix}-{i}.png", distortions)
            print(f"Heightmap saved at {pathHeightmap}{filePrefix}-{i}.png")

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
        features: np.ndarray = np.zeros((dims[0], dims[1]))
        resistances: np.ndarray = np.zeros((dims[0], dims[1]))

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
                features[_x, _y] = self.getSequenceID(pos)
                resistances[_x, _y] = min(max(resistance, 0.0), 1.0)
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
    sizeX, sizeY = len(distortions[0]), len(distortions)
    for _x in range(sizeX):
        for _y in range(sizeY):
            x, y = intsToCoords(_x, _y, sizeX, sizeY)
            closestLineIndex: int = -1
            closestDistance: float = lineWidth
            pos = Vector2D(x, y)
            for i in range(len(curve) - 1):
                if curve[i] == curve[i + 1]: continue
                # Check original and wrapped distances
                candidates = [
                    pos,  # Original position
                    Vector2D(pos.x - 2.0, pos.y),  # Wrapped left
                    Vector2D(pos.x + 2.0, pos.y),  # Wrapped right
                    Vector2D(pos.x, pos.y - 2.0),  # Wrapped up
                    Vector2D(pos.x, pos.y + 2.0),  # Wrapped down
                    Vector2D(pos.x - 2.0, pos.y - 2.0),  # Diagonal wrap: top-left
                    Vector2D(pos.x + 2.0, pos.y - 2.0),  # Diagonal wrap: top-right
                    Vector2D(pos.x - 2.0, pos.y + 2.0),  # Diagonal wrap: bottom-left
                    Vector2D(pos.x + 2.0, pos.y + 2.0),  # Diagonal wrap: bottom-right
                ]

                for candidate in candidates:
                    distToLine = Vectors.distance2ToLine(candidate, curve[i], curve[i + 1])
                    if distToLine < closestDistance:
                        closestDistance = distToLine
                        closestLineIndex = i
                #
                # distToLine = Vectors.distanceToLine(pos, curve[i], curve[i + 1])
                # if distToLine < closestDistance:
                #     closestDistance = distToLine
                #     closestLineIndex = i
            if closestLineIndex > -1:
                closestDistance = math.sqrt(closestDistance)
                distToLine = (closestDistance - 0) / (lineWidth - 0)
                mouseMotion = (curve[closestLineIndex + 1] - curve[closestLineIndex]).normalize() * distortionStrength
                distortions[_x][_y] += mouseMotion * wyvill(distToLine)


def addDistortionFromSketch(distortionSketcher: SketchManagement) -> None:
    mousePath = distortionSketcher.lineBuilders[0].getCurve()
    addDistortionFromCurve(mousePath, 0.1)


def updateResultsFigure(images: List[np.ndarray], _axes: List[plt.Axes]):
    global fig2
    axes = fig2.axes
    for i, img in enumerate(images):
        axes[i].imshow(img)
    fig2.canvas.draw()



def createDatasetOfRandomIslands(islandSketch: IslandSketch, nbSamples: int = 1000):
    global singleDistortionMap
    radiusRandomMin, radiusRandomMax = 0.8, 1.1
    nbCurvesMin, nbCurvesMax = 2, 6
    distoStrengthMin, distoStrengthMax = 0.05, 0.8
    distoWidthMin, distoWidthMax = 0.1, 0.8
    profileRandomMin, profileRandomMax = 0.0, 0.2
    resistanceRandomMin, resistanceRandomMax = 0.0, 0.2

    for iSample in range(nbSamples):
        if os.path.exists(f"{dataset_path}heightmaps/{iSample}.png"):
            continue
        profileRandomness = getRandom(profileRandomMin, profileRandomMax)
        resistanceRandomness = getRandom(resistanceRandomMin, resistanceRandomMax)
        n = noise.perlin.SimplexNoise(1000)
        for _ in range(int(getRandom(nbCurvesMin, nbCurvesMax))):
            singleDistortionMap = initialDistoMap(30, 30)
            addDistortionFromCurve(randomDistortionCurve(), getRandom(distoStrengthMin, distoStrengthMax),
                                   getRandom(distoWidthMin, distoWidthMax))
        radiusFactor = getRandom(radiusRandomMin, radiusRandomMax)
        islandBorders = centeredCircle(0.4 * radiusFactor, 1.2)
        beachBorders = centeredCircle(0.4 * radiusFactor, 1.2)
        lagoonBorders = centeredCircle(0.7 * radiusFactor, 1.2)
        reefBorders = centeredCircle(0.8 * radiusFactor, 1.2)

        islandTranslate = Vectors.randomVec2() * 0.2
        islandBorders = [p + islandTranslate for p in islandBorders]
        beachTranslate = Vectors.randomVec2() * 0.2
        beachBorders = [p + beachTranslate for p in beachBorders]

        islandSketch.sketches[0].setCurve(islandBorders)
        islandSketch.sketches[1].setCurve(beachBorders)
        islandSketch.sketches[2].setCurve(lagoonBorders)
        islandSketch.sketches[3].setCurve(reefBorders)

        _randomProfileCurve = [0.9627757352941178, 0.897671568627451, 0.8425245098039216, 0.775888480392157, 0.6893382352941178, 0.6479779411764706, 0.5974264705882353, 0.5606617647058822, 0.5392156862745099, 0.508578431372549, 0.4986213235294116, 0.48406862745098034, 0.4810049019607843, 0.4779411764705883, 0.4756433823529411, 0.47411151960784303, 0.4595588235294117, 0.4564950980392156, 0.44117647058823517, 0.431985294117647, 0.42585784313725483, 0.315563725490196, 0.2481617647058823, 0.015318627450980338, 0.015318627450980338, 0.015, 0.01, 0.001, 0.0, 0.0] #resizeArray([1.0, 0.0], 10)
        _randomResistanceCurve = [0.7429534313725492, 0.7153799019607845, 0.6755514705882353, 0.626531862745098, 0.5821078431372548, 0.5238970588235294, 0.42585784313725483, 0.27037377450980393, 0.1470588235294117, 0.0850183823529411, 0.036764705882352866, 0.03063725490196073, 0.027573529411764608, 0.027573529411764608, 0.027573529411764608, 0.039828431372548934, 0.06127450980392152, 0.08272058823529405, 0.10569852941176466, 0.15624999999999994, 0.22212009803921565, 0.3117340686274509, 0.4840686274509804, 0.6441482843137254, 0.7467830882352942, 0.8620689655172414, 0.04901960784313719, 0.0337009803921568, 0.021446078431372473, 0.003063725490196012]
        randomProfileCurve = []
        randomResistanceCurve = []
        for i in range(len(_randomProfileCurve)):
            randomProfileCurve.append(_randomProfileCurve[i] * (1.0 + n.noise2(i / 100, iSample * 10000) * profileRandomness))
            randomResistanceCurve.append(_randomResistanceCurve[i] * (1.0 + n.noise2(i / 100, (iSample + 100) * 10000) * resistanceRandomness))
        islandSketch.profileSketch.setCurve(randomProfileCurve)
        islandSketch.resistanceSketch.setCurve(randomResistanceCurve)

        coralMaxHeight = getSequences(islandSketch.profileSketch)[2][2][-1]
        coralMinHeight = getSequences(islandSketch.profileSketch)[3][2][0]

        islandSketch.createMapsFromSketch(path="new_synthetic_terrains_dataset/", filePrefix=str(iSample)) # , coralMinHeight=coralMinHeight, coralMaxHeight=coralMaxHeight)

        # distortionMaps = distortionMaps[:1]
        # singleDistortionMap = initialDistoMap(20, 20)



def main():
    # random.seed(42)
    global fig2
    # global distortionMaps
    global singleDistortionMap
    # fig2 = plt.figure()
    # axHeight, axFeatures, axDisto = fig2.subplots(1, 3)
    fig2, axesResults = plt.subplots(1, 3, squeeze=True)
    axHeight, axFeatures, axDisto = axesResults

    singleDistortionMap = initialDistoMap(20, 20)

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
    sideViewAx.set_title('Side view')
    sideViewAx.set_xlim(0, 1)
    sideViewAx.set_ylim(0, 1)
    distortionsAx.set_title('Distortions')
    distortionsAx.set_xlim(-1, 1)
    distortionsAx.set_ylim(-1, 1)
    resistanceAx.set_title('Resistance')
    resistanceAx.set_xlim(0, 1)
    resistanceAx.set_ylim(0, 1)

    axWater = fig.add_axes((0.92, 0.2, 0.0225, 0.68))
    water_slider = Slider(
        ax=axWater,
        label="Water",
        valmin=0,
        valmax=1,
        valinit=0.4,
        orientation="vertical"
    )

    axCoral = fig.add_axes((0.95, 0.2, 0.0225, 0.68))
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
        label="Subsidence",
        valmin=0,
        valmax=1,
        valinit=0.8,
        orientation="vertical"
    )

    islandSketch = IslandSketch(topViewAx, sideViewAx, distortionsAx, resistanceAx, water_slider, coral_slider, subsid_slider, sketch_names, sketch_colors)


    buttons: List[Button] = []
    button_size = 0.5 / len(sketch_names)
    for i_sketch in range(len(sketch_names)):
        # islandSketches.addSketch(color = sketch_colors[i_sketch], sketch_type = sketch_types[i_sketch])
        sketchID = i_sketch
        ax_button = fig.add_axes((0.5 + button_size * sketchID, 0.05, button_size - 0.01, 0.075))
        buttons.append(Button(ax_button, sketch_names[sketchID]))
        def activationFunction(id):
            def action(event):
                islandSketch.setActiveTopView(id)
            return action
        buttons[-1].on_clicked(activationFunction(sketchID))

    # Add button for generating a distance map / heightmap
    ax_button = fig.add_axes((0.0, 0.05, 0.09, 0.075))
    distance_button = Button(ax_button, "Gen height map")

    def genIsland():
        sequences = getSequences(islandSketch.profileSketch)
        lagoonSequence = sequences[2][2]
        reefSequence = sequences[3][2]
        coralMaxHeight = lagoonSequence[len(lagoonSequence) // 2]
        coralMinHeight = lagoonSequence[-1] #reefSequence[len(reefSequence) // 2]
        # subsidence = 0.8
        # print(coralMinHeight, coralMaxHeight, subsidence, lagoonSequence + reefSequence)
        # print("Profile:", [p.y for p in islandSketch.profileSketch.getCurve()])
        # print("Resistance", [p.y for p in islandSketch.resistanceSketch.getCurve()])
        h, f, d = islandSketch.createMapsFromSketch(path="test_island_heightmap/", filePrefix="test")
        axesResults[0].imshow(h)
        axesResults[1].imshow(f)
        axesResults[2].imshow(d)
        # for ax in axesResults:
        #     ax.update()
        fig2.canvas.draw()
        fig2.canvas.flush_events()


    distance_button.on_clicked(lambda e: genIsland()) # genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut))
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
    for i in range(len(sketch_names)):
        islandSketch.sketches[i].setCurve(centeredCircle(float(i + 1) / float(len(sketch_names) + 1), 1.1))
    # createDatasetOfRandomIslands(islandSketch, 500)
    islandSketch.update()
    plt.show()




if __name__ == "__main__":
    main()
