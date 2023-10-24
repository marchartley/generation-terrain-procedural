from copy import deepcopy
import random
import math
import os.path
import random
from copy import deepcopy
from functools import lru_cache
from typing import Tuple, List, Any

import matplotlib.pyplot as plt
import noise.perlin
import numpy as np
from noise import perlin

import Vectors
import curves
from Vectors import Vector2D
from sketch_maker import LineBuilder1D, LineBuilder2D, LineBuilderRadial, SketchManagement, resizeArray

outputImageDims = [256, 256]
distortionMaps: List[List[List[Vector2D]]] = []
fig2: plt.Figure
dataset_path = "./correct_synthetic_terrains_dataset/" if os.name == "nt" else "/data/correct_synthetic_terrains_dataset/"

class IslandSketch:
    def __init__(self, topViewAx: plt.Axes, profileAx: plt.Axes, distortionAx: plt.Axes):
        self.topView = SketchManagement(topViewAx)
        self.sketches = [
            # "islandBorders":
            LineBuilderRadial(topViewAx, color="green"),
            # "beachBorders":
            LineBuilderRadial(topViewAx, color="yellow"),
            # "lagoonBorders":
            LineBuilderRadial(topViewAx, color="blue"),
            # "reefBorders":
            LineBuilderRadial(topViewAx, color="orange")
        ]
        self.topView.addSketchs(self.sketches)

        self.profileView = SketchManagement(profileAx)
        self.profileSketch = LineBuilder1D(profileAx, color="green")
        self.profileView.addSketch(line=self.profileSketch)

        self.distortionView = SketchManagement(distortionAx)
        self.distortionSketch = LineBuilder2D(distortionAx, color="blue")
        self.distortionView.addSketch(line=self.distortionSketch)

        self.sequences = [
            # (-1, -1, 0.0, 0.0),
            (-1, 0, 0.0, 0.2),
            (0, 1, 0.2, 0.3),
            (1, 2, 0.3, 0.5),
            (2, 3, 0.5, 0.8),
            (3, 1000, 0.8, 1.0)
        ]

        self.topView.onChangeEnded(self.update)
        self.profileView.onChangeEnded(self.update)
        self.distortionView.onChangeEnded(self.update)

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

    def update(self):
        self.updateTopViewAx()
        self.updateSideViewMarkers()
        self.updateDistortionsAx()

    def updateTopViewAx(self):
        self.topView.draw()

    def updateDistortionsAx(self):
        sizeX = len(distortionMaps[0][0])
        sizeY = len(distortionMaps[0])
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
        self.distortionView.draw()

    def updateSideViewMarkers(self):
        ax: plt.Axes = self.profileView.ax
        for line in ax.lines:
            line.set_data([], [])

        ax.vlines([0.2, 0.3, 0.5, 0.8, 1.0], colors=[l.color for l in self.sketches], ymin=0.0, ymax=1.0,
                  linestyles="--")
        self.profileView.draw()

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

    def evaluateHeight(self, pos: Vector2D):
        allDistances = self.getDistancesOfBorders(pos)
        polar = pos.to_polar()
        distFromCenter = polar.y
        sequenceID = self.getSequenceID(pos)
        fullProfile = self.profileSketch.getCurve()
        marker1, marker2, t1, t2 = self.sequences[sequenceID]
        profile = curves.subcurve(fullProfile, t1, t2)
        distMin, distMax = (allDistances[marker1] if marker1 != -1 else 0), (
            allDistances[marker2] if marker2 != 1000 else 1.5)
        t = (distFromCenter - distMin) / (distMax - distMin)
        return interpolateOnCurve(profile, t).y ** 2.0

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

        heightmap, features, distortions = self.heightFeatsAndDistoFromSketches()
        rgb = np.zeros((heightmap.shape[0], heightmap.shape[1], 3))
        rgb[:, :, 0] = rgb[:, :, 1] = rgb[:, :, 2] = heightmap
        heightmap = rgb
        plt.imsave(pathHeightmap + filePrefix + ".png", heightmap)
        plt.imsave(pathFeatures + filePrefix + ".png", features)
        plt.imsave(pathDisto + filePrefix + ".png", distortions)

        updateResultsFigure([heightmap, features, distortions], [fig2.axes[0]])

        print("Heightmaps saved at " + os.path.abspath(pathHeightmap + filePrefix + ".png") + ", " +
              os.path.abspath(pathFeatures + filePrefix + ".png") + " and " +
              os.path.abspath(pathDisto + filePrefix + ".png"))

    def heightFeatsAndDistoFromSketches(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        dims = outputImageDims
        heights: np.ndarray = np.zeros((dims[0], dims[1]))
        features: np.ndarray = np.zeros((dims[0], dims[1]))
        disto_heights: np.ndarray = np.zeros((dims[0], dims[1]))
        disto_features: np.ndarray = np.zeros((dims[0], dims[1], 3))
        distorted = np.zeros((dims[0], dims[1], 3))

        n = noise.perlin.SimplexNoise()

        for _y in range(dims[1]):
            for _x in range(dims[0]):
                noiseVal = 1.0 + (n.noise2(_x / 20, _y / 20) / 10.0)
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                pos = Vector2D(x, y)
                height = self.evaluateHeight(pos) * noiseVal
                heights[_x, _y] = min(max(height, 0.0), 1.0)
                features[_x, _y] = self.getSequenceID(pos)
        # distortion part :
        for _y in range(dims[1]):
            for _x in range(dims[0]):
                x, y = numpyIndicesToCoords(_x, _y, dims[0], dims[1])
                newX, newY = evaluatePosAfterDistortion(x, y)
                diffX = x - newX
                diffY = y - newY
                distorted[_x, _y] = [diffX * .25 + .5, diffY * .25 + .5, 0.5]
                disto_heights[_x, _y] = numpyBilinearInterpolation(heights, newX, newY)
                disto_features[_x, _y] = valueAsHSV(numpyNearestNeighbor(features, newX, newY), -1, len(self.sequences))
        return disto_heights, disto_features, distorted


def numpyIndicesToCoords(_x: int, _y: int, sizeX: int, sizeY: int) -> Tuple[float, float]:
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
    return math.floor((x + 1) * 0.5 * (sizeX)), math.floor((y + 1) * 0.5 * (sizeY))


def coordsToFloats(x: float, y: float, sizeX: int, sizeY: int) -> Tuple[float, float]:
    return (x + 1) * 0.5 * (sizeX), (y + 1) * 0.5 * (sizeY)


def wyvill(x: float):
    return ((1.0 - x) ** 2) ** 3


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
    d00 = arr[_x, _y] if _x >= 0 and _x < sizeX and _y >= 0 and _y < sizeY else nullArray
    d01 = arr[_x, _y + 1] if _x >= 0 and _x < sizeX and (_y + 1) >= 0 and (_y + 1) < sizeY else nullArray
    d10 = arr[_x + 1, _y] if (_x + 1) >= 0 and (_x + 1) < sizeX and _y >= 0 and _y < sizeY else nullArray
    d11 = arr[_x + 1, _y + 1] if (_x + 1) >= 0 and (_x + 1) < sizeX and (_y + 1) >= 0 and (
            _y + 1) < sizeY else nullArray
    value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
    return value


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
    d00 = arr[_x][_y] if _x >= 0 and _x < sizeX and _y >= 0 and _y < sizeY else Vector2D(0, 0)
    d01 = arr[_x][_y + 1] if _x >= 0 and _x < sizeX and (_y + 1) >= 0 and (_y + 1) < sizeY else Vector2D(0, 0)
    d10 = arr[_x + 1][_y] if (_x + 1) >= 0 and (_x + 1) < sizeX and _y >= 0 and _y < sizeY else Vector2D(0, 0)
    d11 = arr[_x + 1][_y + 1] if (_x + 1) >= 0 and (_x + 1) < sizeX and (_y + 1) >= 0 and (
            _y + 1) < sizeY else Vector2D(0, 0)
    value = (d00 * (1 - dx) + d10 * dx) * (1 - dy) + (d01 * (1 - dx) + d11 * dx) * dy
    return value


def valueAsHSV(value: float, mini: float, maxi: float) -> Tuple[float, float, float]:
    H = (360.0 * (value - mini) / (maxi - mini)) / 60
    L = 0.5
    S = 1.0

    C = (1 - abs(2 * L - 1)) * S
    X = C * (1 - abs(math.fmod(H, 2) - 1))

    R1, G1, B1 = (C, X, 0) if 0 <= H < 1 else (X, C, 0) if 1 <= H < 2 else (0, C, X) if 2 <= H < 3 else (
        0, X, C) if 3 <= H < 4 else (X, 0, C) if 4 <= H < 5 else (C, 0, X)
    m = L - C * .5
    R, G, B = R1 + m, G1 + m, B1 + m
    return R, G, B


def interpolateOnCurve(curve: List[Any], t: float) -> Any:
    if t >= 1.0:
        return curve[-1]
    if t <= 0:
        return curve[0]
    i = t * (len(curve) - 1)
    a = i - math.floor(i)

    p0 = curve[math.floor(i)]
    p1 = curve[math.ceil(i)]
    return p0 * (1 - a) + p1 * a


def addDistortionFromCurve(curve: List[Vector2D], distortionStrength: float, lineWidth: float = 0.5) -> None:
    distortions = deepcopy(distortionMaps[0])  # Get the first distortion map, where everything is null vectors
    sizeX, sizeY = len(distortions[0]), len(distortions)
    for _y in range(sizeY):
        for _x in range(sizeX):
            x, y = intsToCoords(_x, _y, sizeX, sizeY)
            closestLineIndex: int = -1
            closestDistance: float = lineWidth
            pos = Vector2D(x, y)
            for i in range(len(curve) - 1):
                if curve[i] == curve[i + 1]: continue
                distToLine = Vectors.distanceToLine(pos, curve[i], curve[i + 1])
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


def randomDistortionCurve() -> List[Vector2D]:
    n = noise.perlin.SimplexNoise()
    p = Vector2D(random.random() * 2 - 1, random.random() * 2 - 1)
    positions: List[Vector2D] = []
    freq = 1.0 / 3.0
    for _ in range(random.randint(20, 40)):
        x = n.noise2(p.x / freq, p.y / freq)
        y = n.noise2(p.x / freq, p.y / freq + 1000.0)
        p.x, p.y = p.x + x * .5, p.y + y * .5
        positions.append(p.copy())
    # for _ in range(5):
    #     positions = curves.smoothCurve(positions)
    positions = curves.catmull_rom_chain(resizeArray(positions, 10))
    return positions


def getRandom(mini: float, maxi: float) -> float:
    return mini + random.random() * (maxi - mini)


def updateResultsFigure(images: List[np.ndarray], _axes: List[plt.Axes]):
    global fig2
    axes = fig2.axes
    for i, img in enumerate(images):
        axes[i].imshow(img)
    fig2.canvas.draw()


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


def main():
    global fig2
    global distortionMaps
    fig2 = plt.figure()
    axHeight, axFeatures, axDisto = fig2.subplots(1, 3)

    distortionMaps.append(initialDistoMap(20, 20))
    # for _ in range(3):
    #     addDistortionFromCurve(randomDistortionCurve(), 0.1)
    # distortionMaps.append([[Vector2D(math.cos(x / 5), math.sin(y / 2)) * 0.1 for y in range(20)] for x in range(20)])
    # distortionMaps.append([[Vector2D(0, 1) if x < 10 else Vector2D(0, 0) for y in range(20)] for x in range(20)])

    sketch_names = ["Island", "Beach", "Lagoon", "Reef"]
    sketch_colors = ["green", "yellow", "blue", "orange"]
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
    sideViewAx.set_xlim(0, 1)
    sideViewAx.set_ylim(0, 1)
    distortionsAx.set_title('Distortions')
    distortionsAx.set_xlim(-1, 1)
    distortionsAx.set_ylim(-1, 1)

    radiusRandomMin, radiusRandomMax = 0.8, 1.1
    nbCurvesMin, nbCurvesMax = 2, 6
    distoStrengthMin, distoStrengthMax = 0.05, 0.4
    distoWidthMin, distoWidthMax = 0.2, 0.8
    profileRandomMin, profileRandomMax = 0.0, 0.2

    islandSketch = IslandSketch(topViewAx, sideViewAx, distortionsAx)
    # for iSample in range(100000):
    for iSample in range(1):
        if os.path.exists(f"{dataset_path}heightmaps/{iSample}.png"):
            continue
        profileRandomness = getRandom(profileRandomMin, profileRandomMax)
        n = noise.perlin.SimplexNoise(1000)
        for _ in range(int(getRandom(nbCurvesMin, nbCurvesMax))):
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

        _randomProfileCurve = [0.9175971555111971, 0.8957644563657912, 0.8053147027633958, 0.7413760838375647,
                               0.6041419749235855, 0.574511883226249, 0.5557981411016155, 0.5308464849354375,
                               0.4980974362173289, 0.4669078660096063, 0.4372777743122699, 0.41544507516686413,
                               0.40296924708377513, 0.38737446197991376, 0.3651518932069115, 0.3608633273033497,
                               0.3608633273033497, 0.3612531969309462, 0.3671012413448942, 0.3729492857588421,
                               0.3807466783107727, 0.4045287255941612, 0.4357182958018837, 0.48406212962385375,
                               0.5635955336535462, 0.028694404591104672, 0.03493231863264917, 0.024015969059946296,
                               0.025575447570332418, 0.011540140976857238]
        # _randomProfileCurve = [1.0, 0.9, 0.7, 0.3, 0.3, 0.5, 0.1, 0.01]
        randomProfileCurve = []
        for i in range(len(_randomProfileCurve)):
            randomProfileCurve.append(
                _randomProfileCurve[i] * (1.0 + n.noise2(i / 100, iSample * 10000) * profileRandomness))
        islandSketch.profileSketch.setCurve(randomProfileCurve)
        islandSketch.createMapsFromSketch(path="test_island_heightmap", filePrefix=str(iSample))

        distortionMaps = distortionMaps[:1]
    islandSketch.update()
    plt.show()


if __name__ == "__main__":
    main()
