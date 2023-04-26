import itertools
import time
from copy import deepcopy
import random

import PIL
import matplotlib.backend_bases
import noise.perlin
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path
from typing import Tuple, List, Callable, Union, Optional, Any

import PIL.Image
from matplotlib.widgets import Button
from Vectors import Vector2D, Vector3D, line_intersection
import Vectors
from noise import perlin
from sketch_maker import LineBuilder1D, LineBuilder2D, LineBuilderRadial, SketchManagement, LineBuilder, resizeArray
from itertools import chain

outputImageDims = [128, 128]
distortionMaps: List[List[List[Vector2D]]] = []
fig2: plt.Figure

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
    return [[Vector2D(0, 0) for _ in range(sizeY)] for _ in range(sizeX)]

def getDisto(x: float, y: float) -> Vector2D:
    distortion = Vector2D(x, y)
    for disto in reversed(distortionMaps):
        d = bilinearInterpolation(disto, distortion.x, distortion.y)
        d.y *= -1
        distortion += d
    return distortion - Vector2D(x, y)

def getDistoFromIndices(_x: int, _y: int) -> Vector2D:
    sizeX, sizeY = len(distortionMaps[0][0]), len(distortionMaps[0])
    vec: Vector2D = Vector2D(0, 0)
    if _x < 0 or _x >= sizeX or _y < 0 or _y >= sizeY:
        return vec
    for distoMap in distortionMaps:
        vec += distoMap[_x][_y]
    return vec

def evaluatePosAfterDistortion(prevX: float, prevY: float) -> Tuple[float, float]:
    pos = Vector2D(prevX, prevY) - getDisto(prevX, prevY)
    return pos.x, pos.y

def numpyBilinearInterpolation(arr: np.ndarray, x: float, y: float) -> float:
    sizeX, sizeY = arr.shape[0], arr.shape[1]
    nbChannels = 0 if len(arr.shape) == 2 else arr.shape[2]
    nullArray = 0 if nbChannels == 0 else np.zeros(nbChannels)
    _x, _y = coordsToNumpyIndices(x, y, sizeX, sizeY)
    __x, __y = coordsToNumpyIndicesFloat(x, y, sizeX, sizeY)

    dx, dy = __x - _x, __y - _y
    d00 = arr[_x  , _y  ] if _x >= 0 and _x < sizeX and _y >= 0 and _y < sizeY else nullArray
    d01 = arr[_x  , _y+1] if _x >= 0 and _x < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else nullArray
    d10 = arr[_x+1, _y  ] if (_x+1) >= 0 and (_x+1) < sizeX and _y >= 0 and _y < sizeY else nullArray
    d11 = arr[_x+1, _y+1] if (_x+1) >= 0 and (_x+1) < sizeX and (_y+1) >= 0 and (_y+1) < sizeY else nullArray
    value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
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

def valueAsHSV(value: float, mini: float, maxi: float) -> Tuple[float, float, float]:
    H = (360.0 * (value - mini) / (maxi - mini)) / 60
    L = 0.5
    S = 1.0

    C = (1 - abs(2 * L - 1)) * S
    X = C * (1 - abs(math.fmod(H, 2) - 1))

    R1, G1, B1 = (C, X, 0) if 0 <= H < 1 else (X, C, 0) if 1 <= H < 2 else (0, C, X) if 2 <= H < 3 else (0, X, C) if 3 <= H < 4 else (X, 0, C) if 4 <= H < 5 else (C, 0, X)
    m = L - C * .5
    R, G, B = R1 + m, G1 + m, B1 + m
    return R, G, B

def distanceMapFromSketches(islandSketch: SketchManagement) -> np.ndarray:
    dims = (100, 100, len(islandSketch.lineBuilders) + 1) # Add the "borders" as last dimension
    heightmap = np.zeros(dims)

    for _y in range(dims[0]):
        for _x in range(dims[1]):
            x = (_x / dims[1]) * 2.0 - 1.0  # Transferring on [-1, 1]^2
            y = (_y / dims[0]) * 2.0 - 1.0
            pos = Vector2D(x, y)
            for iChannel, line in enumerate(islandSketch.lineBuilders):
                minDist = max(dims)
                curve = line.getCurve()
                for i in range(len(curve) - 1):
                    minDist = min(minDist, Vectors.distanceToLine(pos, curve[i], curve[i + 1]))
                heightmap[dims[0] - _y - 1, _x, iChannel] = minDist
            distToBorder = min([abs(x), abs(y), abs(1-x), abs(1-y)])
            heightmap[dims[0] - _y - 1, _x, -1] = distToBorder
    return heightmap

def genAndSaveDistanceMap(islandSketch: SketchManagement):
    image = distanceMapFromSketches(islandSketch)
    image = image[:, :, :3]
    plt.imsave("test.png", image)

def splitProfileOnMarkers(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> List[Tuple[int, int, List[float]]]:
    """Extract the curves made by the profile depending on the distance between each island sketch"""
    markers: List[Tuple[float, int]] = []

    # sliceDirection = (sliceCut[1] - sliceCut[0]).normalize()
    for sketchID, sketch in enumerate(islandSketches.lineBuilders):
        intersections = sketch.intersection(*sliceCut)
        for intersect in intersections:
            # distanceOnSlice = sliceDirection.dot(intersect - sliceCut[0])
            # print()
            distanceOnSlice = intersect.x
            markers.append((distanceOnSlice, sketchID))

    centerMarker = -1
    borderMarker = len(islandSketches.lineBuilders)
    markers.append((0.0, centerMarker))
    markers.sort(key = lambda a: a[0])

    profile = resizeArray(profileSketch.getCurve(), 200)
    curves: List[Tuple[int, int, List[float]]] = []
    begin = (-1.0, borderMarker)  # markers.pop(0)
    ending = markers[0]

    startingIndex = 0
    for i, p in enumerate(profile):
        x, y = p.x, p.y
        if x >= ending[0]:  # End of the curve, beginning of a new one
            curves.append((begin[1], ending[1], [_p.y for _p in profile[startingIndex:i]]))
            if len(markers) > 0:
                startingIndex = i
                begin = markers.pop(0)
                if len(markers) == 0:
                    ending = (1.0, borderMarker)
                else:
                    ending = markers[0]
    return curves


def interpolateOnCurve(curve: List[Any], t: float) -> Any:
    if t >= 1.0:
        return curve[-1]
    if t <= 0:
        return curve[0]
    i = t * (len(curve) - 1)
    a = i - math.floor(i)

    p0 = curve[math.floor(i)]
    p1 = curve[math.ceil(i)]
    return p0 * a + p1 * (1 - a)

def getDistancesToCurves(x: float, y: float, islandSketches: SketchManagement, profileSlice: Vector2D) -> List[float]:
    distancesToCurves = [(i, curve.intersection(Vector2D(), profileSlice)) for i, curve in enumerate(islandSketches.lineBuilders)]
    distancesToCurves = [vecs[0].norm() if len(vecs) > 0 else (Vector2D(x, y) * 1000).norm() for i, vecs in distancesToCurves]
    distancesToCurves.append(1.5)  # approximation of sqrt(2)
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
    if distToCenter < distancesToCurves[indices[0]]:
        t = distToCenter / distancesToCurves[indices[0]]
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


def heightmapFromSketches(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> np.ndarray:
    sequences = splitProfileOnMarkers(profileSketch, islandSketches, sliceCut)
    dims = outputImageDims
    heights: np.ndarray = np.zeros((dims[0], dims[1]))

    for _x in range(dims[0]):
        for _y in range(dims[1]):
            x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
            pos = Vector2D(x, y)
            distToCenter = pos.norm()
            profileSlice = pos.normalized() * 2# if distToCenter > 0.1 else Vector2D(1, 0)

            distancesToCurves = getDistancesToCurves(x, y, islandSketches, profileSlice)
            marker1, marker2, height, curve = getSequence(sequences, distancesToCurves, distToCenter)
            if distToCenter < 0.01:
                height = interpolateOnCurve(profileSketch.getCurve(), 0.5).y
            heights[_x, _y] = height

    # Distortion part :
    distorted = np.zeros((dims[0], dims[1]))
    for _x in range(dims[0]):
        for _y in range(dims[1]):
            x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
            newX, newY = evaluatePosAfterDistortion(x, y)
            val = numpyBilinearInterpolation(heights, newX, newY)
            distorted[_x, _y] = val
    return distorted


def featuresFromSketches(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> np.ndarray:
    sequences = splitProfileOnMarkers(profileSketch, islandSketches, sliceCut)
    dims = outputImageDims
    features: np.ndarray = np.zeros((dims[0], dims[1]))

    for _x in range(dims[0]):
        for _y in range(dims[1]):
            x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
            pos = Vector2D(x, y)
            distToCenter = pos.norm()
            profileSlice = pos.normalized() * 2# if distToCenter > 0.1 else Vector2D(1, 0)

            distancesToCurves = getDistancesToCurves(x, y, islandSketches, profileSlice)
            marker1, marker2, height, curve = getSequence(sequences, distancesToCurves, distToCenter)
            seqID = getSequenceID((marker1, marker2, None), sequences)
            features[_x, _y] = seqID

    # Distortion part :
    distorted = np.zeros((dims[0], dims[1], 3))
    for _x in range(dims[0]):
        for _y in range(dims[1]):
            x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
            newX, newY = evaluatePosAfterDistortion(x, y)
            val = numpyBilinearInterpolation(features, newX, newY)
            distorted[_x, _y] = valueAsHSV(val, 0, len(sequences))
    return distorted


def distortionsFromSketches(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> np.ndarray:
    dims = outputImageDims
    distorted = np.zeros((dims[0], dims[1], 3))
    for _x in range(dims[0]):
        for _y in range(dims[1]):
            x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
            newX, newY = evaluatePosAfterDistortion(x, y)
            diffX = x - newX
            diffY = y - newY
            distorted[_x, _y] = [diffX * .25 + .5, diffY * .25 + .5, 0.5]
    return distorted

def genAndSaveHeightMap(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D], path="./", featuresFolder="features/", distoFolder="distortions/", heightFolder="heightmaps/", filePrefix = "result"):
    pathHeightmap = path + heightFolder
    pathFeatures = path + featuresFolder
    pathDisto = path + distoFolder
    os.makedirs(pathHeightmap, exist_ok=True)
    os.makedirs(pathFeatures, exist_ok=True)
    os.makedirs(pathDisto, exist_ok=True)
    heightmap = heightmapFromSketches(profileSketch, islandSketches, sliceCut)
    features = featuresFromSketches(profileSketch, islandSketches, sliceCut)
    distortions = distortionsFromSketches(profileSketch, islandSketches, sliceCut)
    rgb = np.zeros((heightmap.shape[0], heightmap.shape[1], 3))
    rgb[:, :, 0] = rgb[:, :, 1] = rgb[:, :, 2] = heightmap
    plt.imsave(pathHeightmap + filePrefix + ".png", rgb)
    plt.imsave(pathFeatures + filePrefix + ".png", features)
    plt.imsave(pathDisto + filePrefix + ".png", distortions)

    updateResultsFigure([heightmap, features, distortions], [fig2.axes[0]])

    print("Heightmaps saved at " + os.path.abspath(pathHeightmap + filePrefix + ".png") + ", " +
          os.path.abspath(pathFeatures + filePrefix + ".png") + " and " +
          os.path.abspath(pathDisto + filePrefix + ".png"))


mousePressed = False
def onMousePressed(e: matplotlib.backend_bases.FigureCanvasBase):
    global mousePressed
    mousePressed = True
def onMouseReleased(e: matplotlib.backend_bases.FigureCanvasBase, topViewAx, sideViewAx, distortionsAx, waterLevel, islandSketches, profileSketching, sliceCut, fig, distortionSketcher):
    global mousePressed
    mousePressed = False
    if e.inaxes in [topViewAx, sideViewAx]:
        updateSideViewMarkers(sideViewAx, waterLevel, islandSketches, profileSketching, sliceCut, fig)
    elif e.inaxes == distortionsAx:
        updateDistortionsAx(distortionsAx, distortionSketcher)
def onMouseMove(e: matplotlib.backend_bases.FigureCanvasBase, topViewAx, sideViewAx, distortionsAx, waterLevel, islandSketches, profileSketching, sliceCut, fig, distortionSketcher):
    if mousePressed:
        if e.inaxes in [topViewAx, sideViewAx]:
            updateSideViewMarkers(sideViewAx, waterLevel, islandSketches, profileSketching, sliceCut, fig)
        elif e.inaxes == distortionsAx:
            updateDistortionsSketches(distortionsAx, distortionSketcher)

# Each time the top view sketch is modified, reposition markers on the side view
def updateSideViewMarkers(sideViewAx: plt.Axes, waterLevel: float, islandSketches: SketchManagement, profileSketching: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D], fig: plt.Figure):
    for line in sideViewAx.lines:
        line.set_data([], [])
    # Add water level to profile editing
    sideViewAx.axhline(y = waterLevel, color = "blue")
    for sketch in islandSketches.lineBuilders:
        intersections = sketch.intersection(*sliceCut)
        for intersect in intersections:
            sideViewAx.axvline(x = intersect.x, color = sketch.color)
    islandSketches.draw()
    profileSketching.draw()

def addDistortionFromCurve(curve: List[Vector2D], distortionStrength: float) -> None:
    distortions = deepcopy(distortionMaps[0]) # Get the first distortion map, where everything is null vectors
    sizeX, sizeY = len(distortions[0]), len(distortions)
    lineWidth: float = 0.5
    for _y in range(sizeY):
        for _x in range(sizeX):
            x, y = intsToCoords(_x, _y, sizeX, sizeY)
            closestLineIndex: int = -1
            closestDistance: float = lineWidth
            for i in range(len(curve) - 1):
                distToLine = Vectors.distanceToLine(Vector2D(x, y), curve[i], curve[i + 1])
                if distToLine < closestDistance:
                    closestDistance = distToLine
                    closestLineIndex = i
            if closestLineIndex > -1:
                distToLine = (closestDistance - 0) / (lineWidth - 0)
                mouseMotion = (curve[closestLineIndex + 1] - curve[closestLineIndex]).normalize() * distortionStrength
                distortions[_x][_y] = mouseMotion * wyvill(distToLine)
    distortionMaps.append(distortions)

def addDistortionFromSketch(distortionSketcher: SketchManagement) -> None:
    mousePath = distortionSketcher.lineBuilders[0].getCurve()
    addDistortionFromCurve(mousePath, 0.1)

def updateDistortionsSketches(distortionsAx: plt.Axes, distortionSketcher: SketchManagement):
    distortionSketcher.draw()

def updateDistortionsAx(distortionsAx: plt.Axes, distortionSketcher: SketchManagement):
    sizeX = len(distortionMaps[0][0])
    sizeY = len(distortionMaps[0])
    addDistortionFromSketch(distortionSketcher)
    distortionSketcher.lineBuilders[0].reset()
    wholeMap = [[getDistoFromIndices(x, y) for y in range(sizeY)] for x in range(sizeX)]
    for l in distortionsAx.lines:
        if l is not distortionSketcher.lineBuilders[0].line:
            l.remove()
    for _y in range(sizeY):
        for _x in range(sizeX):
            x, y = intsToCoords(_x, _y, sizeX, sizeY)
            distortionsAx.plot([x, x + wholeMap[_x][_y].x], [y, y + wholeMap[_x][_y].y], c="blue")
    distortionSketcher.draw()

def smoothCurve(curve: List[Vector2D]) -> List[Vector2D]:
    res = []
    res.append(curve[0])
    for i in range(len(curve) - 1):
        res.append((curve[i] + curve[i + 1]) * .5)
    res.append(curve[-1])
    return res


def randomDistortionCurve() -> List[Vector2D]:
    n = noise.perlin.SimplexNoise()
    p = Vector2D(random.random() * 2 - 1, random.random() * 2 - 1)
    positions: List[Vector2D] = []
    freq = 1.0/3.0
    for _ in range(random.randint(20, 40)):
        x = n.noise2(p.x / freq, p.y / freq)
        y = n.noise2(p.x / freq, p.y / freq + 1000.0)
        p += Vector2D(x, y) * .5
        toCenter = (p * -1)
        p += 0.1 * toCenter  # /(1.5 - toCenter.norm())
        positions.append(p.copy())
    for _ in range(5):
        positions = smoothCurve(positions)
    return positions

def createDatasetOfRandomIslands(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D], filePrefix = "result"):
    global distortionMaps
    distortionMaps.append(initialDistoMap(20, 20))
    maxDistortionCurves = 10
    maxDistortionStrength = 0.5
    folder = "correct_synthetic_terrains_dataset/"
    if not os.path.exists(folder):
        os.makedirs(folder)

    for iSample in range(1000000):
        if os.path.exists(folder + "features/" + filePrefix + str(iSample) + ".png"):
            continue
        nbDistortions = random.randint(1, maxDistortionCurves)
        for _ in range(nbDistortions):
            strength = random.random() * maxDistortionStrength
            addDistortionFromCurve(randomDistortionCurve(), strength)
        genAndSaveHeightMap(profileSketch, islandSketches, sliceCut, folder, "features/", "distortions/", "heightmaps/", filePrefix + str(iSample))
        distortionMaps = distortionMaps[0:1]
    print("Finished!")

def updateResultsFigure(images: List[np.ndarray], _axes: List[plt.Axes]):
    global fig2
    axes = fig2.axes
    for i, img in enumerate(images):
        axes[i].imshow(img)
    fig2.canvas.draw()

def main():
    global fig2

    distortionMaps.append(initialDistoMap(20, 20))
    # for _ in range(3):
        # addDistortionFromCurve(randomDistortionCurve(), 0.1)
    # distortionMaps.append([[Vector2D(math.cos(x / 5), math.sin(y / 2)) * 0.1 for y in range(20)] for x in range(20)])
    # distortionMaps.append([[Vector2D(0, 1) if x < 10 else Vector2D(0, 0) for y in range(20)] for x in range(20)])

    sketch_names = ["Reef", "Island", "Passes"]
    sketch_colors = ["orange", "green", "blue"]
    sketch_types = [LineBuilderRadial, LineBuilderRadial, LineBuilder2D]
    waterLevel = 0.5

    # Creates a top view for island sketching and a side view for profile editing
    fig, axes = plt.subplots(1, 3, squeeze=False)
    topViewAx: plt.Axes = axes[0, 0]
    sideViewAx: plt.Axes = axes[0, 1]
    distortionsAx: plt.Axes = axes[0, 2]
    fig.subplots_adjust(bottom=0.2)

    topViewAx.set_title('Top view')
    topViewAx.set_xlim(-1, 1)
    topViewAx.set_ylim(-1, 1)
    sideViewAx.set_title('Side view')
    sideViewAx.set_xlim(-1, 1)
    sideViewAx.set_ylim(0, 1)
    distortionsAx.set_title('Distortions')
    distortionsAx.set_xlim(-1, 1)
    distortionsAx.set_ylim(-1, 1)

    fig2 = plt.figure()
    axHeight, axFeatures, axDisto = fig2.subplots(1, 3)

    # Represent the center of the map with an ellipse
    # topViewAx.add_patch(matplotlib.patches.Circle((0, 0), 0.1))

    # Add sketching for top island sketching
    islandSketches = SketchManagement(topViewAx)

    buttons: List[Button] = []
    button_size = 0.5 / len(sketch_names)
    for i_sketch in range(len(sketch_names)):
        islandSketches.addSketch(color = sketch_colors[i_sketch], sketch_type = sketch_types[i_sketch])
        sketchID = i_sketch
        ax_button = fig.add_axes([0.5 + button_size * sketchID, 0.05, button_size - 0.01, 0.075])
        buttons.append(Button(ax_button, sketch_names[sketchID]))
        def activationFunction(id):
            def action(event):
                islandSketches.activate(id)
            return action
        buttons[-1].on_clicked(activationFunction(sketchID))

    # Add button for generating a distance map / heightmap
    ax_button = fig.add_axes([0.0, 0.05, 0.09, 0.075])
    distance_button = Button(ax_button, "Gen height map")
    distance_button.on_clicked(lambda e: genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut))
    # Add button for splitting sequences
    ax_button = fig.add_axes([0.1, 0.05, 0.09, 0.075])
    splitting_button = Button(ax_button, "Split profile")
    splitting_button.on_clicked(lambda e: splitProfileOnMarkers(profileSketching.lineBuilders[0], islandSketches, sliceCut))

    # Profile will be defined by the horizontal plane passing through center
    sliceCut = (Vector2D(-1, 0), Vector2D(1, 0))
    topViewAx.plot([sliceCut[0].x, sliceCut[1].x], [sliceCut[0].y, sliceCut[1].y], linestyle="--", color="blue")

    # Add sketching for profile editing
    profileSketching: SketchManagement = SketchManagement(sideViewAx)
    sketchPro = LineBuilder1D(sideViewAx, 40, "blue")
    # sketchPro.xMin = 0
    # sketchPro.xMax = 1
    profileSketching.addSketch(line = sketchPro)

    # Testing part
    def centeredCircle(radius:float, randomness: float = 0.2) -> List[Vector2D]:
        perlin_noise = perlin.SimplexNoise(10)
        points: List[Vector2D] = []
        nbPoints = 30
        for i in range(nbPoints):
            angle = i * math.tau / (nbPoints - 1)
            vertexUnitPos = Vector2D(math.cos(angle), math.sin(angle))
            noiseValue = perlin_noise.noise2(vertexUnitPos.x, vertexUnitPos.y) * randomness
            points.append(vertexUnitPos * (radius * (1 - noiseValue)))
        return points
    randomCoralCurve = centeredCircle(0.5)
    islandSketches.lineBuilders[0].setCurve(randomCoralCurve)
    randomIslandCurve = centeredCircle(0.25)
    islandSketches.lineBuilders[1].setCurve(randomIslandCurve)

    # _randomProfileCurve = [.0, .9, .0, .9, .0, .9, .0, .9, .0]
    _randomProfileCurve = [0.02, 0.04, 0.05, 0.06, 0.07, 0.09, 0.13, 0.17, 0.19, 0.7, 0.57, 0.53, 0.61, 0.7, 0.75, 0.86, 0.9, 0.96, 0.99, 0.94, 0.88, 0.82, 0.75, 0.71, 0.52, 0.52, 0.52, 0.65, 0.69, 0.7, 0.3, 0.22, 0.21, 0.2, 0.2, 0.19, 0.19, 0.18, 0.18, 0.05]
    # _randomProfileCurve = [1.0, 0.9, 0.7, 0.3, 0.3, 0.5, 0.1, 0.0]
    randomProfileCurve = []
    for i in range(len(_randomProfileCurve)):
        randomProfileCurve.append(Vector2D((2 * i/(len(_randomProfileCurve) - 1)) - 1, _randomProfileCurve[i]))
    profileSketching.lineBuilders[0].setCurve(randomProfileCurve)

    cidPress = fig.canvas.mpl_connect('button_press_event', onMousePressed)
    cidRelease = fig.canvas.mpl_connect('button_release_event', lambda e: onMouseReleased(e, topViewAx, sideViewAx, distortionsAx, waterLevel, islandSketches, profileSketching, sliceCut, fig, distortionSketcher))
    cidMove = fig.canvas.mpl_connect('motion_notify_event', lambda e: onMouseMove(e, topViewAx, sideViewAx, distortionsAx, waterLevel, islandSketches, profileSketching, sliceCut, fig, distortionSketcher))

    # Todo:
    #  - From profile, get heightmap depending on coral position and island position
    #  - Add passes in the process
    #  - Create heightmap from input and distance map
    #  - Profile differs depending on "r" and "d" (island width, island-coral distance)

    distortionSketcher = SketchManagement(distortionsAx)
    distortionSketcher.addSketch(sketch_type=LineBuilder2D)
    # distortionSketcher.lineBuilders[0].setCurve([Vector2D(-.8, -.8), Vector2D(.8, .3)])
    # distortionSketcher.onChange(lambda: updateDistortionsSketches(distortionsAx, distortionSketcher))
    # distortionSketcher.onChangeEnded(lambda: updateDistortionsAx(distortionsAx, distortionSketcher))
    updateSideViewMarkers(sideViewAx, waterLevel, islandSketches, profileSketching, sliceCut, fig)
    updateDistortionsAx(distortionsAx, distortionSketcher)
    # genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut)

    # createDatasetOfRandomIslands(profileSketching.lineBuilders[0], islandSketches, sliceCut, filePrefix="sample")

    plt.show()


if __name__ == "__main__":
    main()
