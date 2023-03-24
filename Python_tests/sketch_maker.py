import math
from typing import Tuple, List, Callable, Union, Optional, Any

import matplotlib.patches
import numpy as np
from matplotlib import pyplot as plt

import threading
import time

from matplotlib.widgets import Button

import Vectors
from Vectors import Vector2D, Vector3D, line_intersection

from noise import perlin

def clamp(x, mini, maxi):
    return mini if x < mini else maxi if x > maxi else x

class RepeatedTimer(object):
    def __init__(self, interval, function, *args, **kwargs):
        self._timer = None
        self.interval = interval
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.next_call = time.time()
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self) -> 'RepeatedTimer':
        if not self.is_running:
            self.next_call += self.interval
            self._timer = threading.Timer(self.next_call - time.time(), self._run)
            self._timer.start()
            self.is_running = True

        return self

    def stop(self) -> 'RepeatedTimer':
        self._timer.cancel()
        self.is_running = False
        return self


class LineBuilder:
    def __init__(self, ax: plt.Figure, color: Optional[str] = None, active: bool = True, multiline: bool = False):
        self.line, = ax.plot([0], [0], color = color)
        self.xs = []
        self.ys = []
        self.ax = ax
        self.cid_press = self.line.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_release = self.line.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_move = self.line.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.currentlyPressed = False
        self.active = active
        self.color = color
        self.callbacksOnChange: List[Callable] = []
        self.previousMousePos: Optional[Vector2D] = None
        # self.multiline = multiline

    def addPoint(self, x, y) -> 'LineBuilder':
        if x is None or y is None:
            return self
        self.xs.append(x)
        self.ys.append(y)
        self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def press_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        self.addPoint(event.xdata, event.ydata)
        self.currentlyPressed = True
        self.previousMousePos = Vector2D(event.xdata, event.ydata)

    def release_event(self, _event):
        self.currentlyPressed = False
        self.previousMousePos = None

    def move_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        if self.currentlyPressed:
            self.addPoint(event.xdata, event.ydata)
            self.previousMousePos = Vector2D(event.xdata, event.ydata)

    def reset(self) -> 'LineBuilder':
        self.xs = []
        self.ys = []
        self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def getCurve(self) -> List[Vector2D]:
        return [Vector2D(x, y) for x, y in zip(self.xs, self.ys)]

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilder':
        self.reset()
        for p in points:
            self.addPoint(p.x, p.y)
        return self

    def close(self) -> 'LineBuilder':
        self.xs.append(self.xs[0])
        self.ys.append(self.ys[0])
        return self

    def intersection(self, limitA: Vector2D, limitB: Vector2D) -> List[Vector2D]:
        intersections: List[Vector2D] = []
        curve = self.getCurve()
        for i in range(1, len(curve)):
            p1, p2 = curve[i - 1], curve[i]
            intersect = line_intersection(limitA, limitB, p1, p2)
            if intersect is None:
                continue
            intersections.append(Vector2D(intersect))
        return intersections

    def addCallbackOnChange(self, function: Callable):
        self.callbacksOnChange.append(function)

    def draw(self):
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()


class SketchManagement:
    def __init__(self, ax):
        self.lineBuilders: List[LineBuilder] = []
        self.ax = ax

    def addSketch(self, color: Optional[str] = None, line: Optional[LineBuilder] = None) -> int:
        # isFirstLineCreated = len(self.lineBuilders) == 0  # Activate if first line
        if line is None:
            line = LineBuilder(self.ax, color = color) # , active = isFirstLineCreated)
        self.lineBuilders.append(line)
        self.activate(0)
        return len(self.lineBuilders) - 1

    def activate(self, builder: int):
        for lb in self.lineBuilders:
            lb.active = False
        self.lineBuilders[builder].active = True

    def onChange(self, function: Callable):
        for line in self.lineBuilders:
            line.addCallbackOnChange(function)

    def draw(self):
        for line in self.lineBuilders:
            line.draw()

class LineBuilder1D(LineBuilder):
    def __init__(self, ax: plt.Axes, precision: float = .1, color: Optional[str] = None, active: bool = True):
        super().__init__(ax, color, active)
        self.previousMousePos: Optional[Vector2D] = None
        self.precision = precision
        self.reset()

    def addPoint(self, x, y) -> 'LineBuilder1D':
        if x is None:
            return self
        index = -1
        if x in self.xs:
            index = self.xs.index(x)
        else:
            index = self.getClosestPoint(x)
        self.ys[index] = y
        self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def move_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        if self.currentlyPressed:
            currentMousePos: Vector2D = Vector2D(event.xdata, event.ydata)
            self.setPartialCurve(self.previousMousePos, currentMousePos)
            self.previousMousePos = currentMousePos

    def getClosestPoint(self, x) -> int:
        a, b = self.ax.get_xlim()
        lin = (x - a) / (b - a)
        index = clamp(round(lin * len(self.xs)), 0, len(self.xs) - 1)
        return index

    def setPartialCurve(self, p0: Vector2D, p1: Vector2D) -> 'LineBuilder1D':
        startX = min(p0.x, p1.x)
        endX = max(p0.x, p1.x)
        minIndex = self.getClosestPoint(startX)
        maxIndex = self.getClosestPoint(endX)
        for index in range(minIndex, maxIndex + 1):
            intersection = line_intersection(Vector2D(self.xs[index], 0), Vector2D(self.xs[index], 1), p0, p1, epsilon=self.precision*2)
            if intersection is not None:
                self.ys[index] = intersection[1]
            else:
                self.ys[index] = p0.y  # Last chance: if the  segment is too small, just give it a "random" y-value
        self.draw()
        return self

    def reset(self) -> 'LineBuilder1D':
        nbPoints = math.ceil((self.ax.get_xlim()[1] - self.ax.get_xlim()[0]) / self.precision) + 1
        self.xs = [self.ax.get_xlim()[0] + self.precision * i for i in range(nbPoints)]
        self.ys = [0.0 for _ in range(nbPoints)]
        self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilder1D':
        self.reset()
        intersections: List[Vector2D] = []
        for x in self.xs:
            for iPoint in range(len(points) - 1):
                p0, p1 = points[iPoint], points[iPoint + 1]
                inter = line_intersection(Vector2D(x, 0), Vector2D(x, 1), p0, p1)
                if inter is not None:
                    intersections.append(Vector2D(inter))

        for p in intersections:
            self.addPoint(p.x, p.y)
        return self

    def getPoint(self, x: float) -> float:
        if x in self.xs:
            return self.ys[self.xs.index(x)]
        return self.ys[self.getClosestPoint(x)]

    def close(self) -> 'LineBuilder1D':
        raise NotImplementedError("Cannot close a 1D line")


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
    # image = 1 - np.clip(image * 5.0, 0, 1)
    plt.imsave("test.png", image)

def splitProfileOnMarkers(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]) -> List[Tuple[int, int, List[float]]]:
    """Extract the curves made by the profile depending on the distance between each island sketch"""
    markers: List[Tuple[float, int]] = []

    sliceDirection = (sliceCut[1] - sliceCut[0]).normalize()
    for sketchID, sketch in enumerate(islandSketches.lineBuilders):
        intersections = sketch.intersection(*sliceCut)
        for intersect in intersections:
            # distanceOnSlice = sliceDirection.dot(intersect - sliceCut[0])
            # print()
            distanceOnSlice = intersect.x
            markers.append((distanceOnSlice, sketchID))

    markers.sort(key = lambda a: a[0])

    borderMarker = len(islandSketches.lineBuilders)

    profile = profileSketch.getCurve()
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

def interpolateOnCurve(curve: List[Any], t: float) -> Vector2D:
    if t >= 1.0:
        return curve[-1]
    i = t * len(curve)
    a = i - math.floor(i)
    if math.floor(i) < 0 or math.ceil(i) >= len(curve):
        b = 0
    return curve[math.floor(i)] * a + curve[math.ceil(i)] * (1 - a)

def heightmapFromSketches(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]):
    sequences = splitProfileOnMarkers(profileSketch, islandSketches, sliceCut)
    distances = distanceMapFromSketches(islandSketches)
    heights = np.zeros((distances.shape[0], distances.shape[1]))

    for _y in range(heights.shape[0]):
        for _x in range(heights.shape[1]):
            closestFeatures = [i[0] for i in sorted(enumerate(distances[_y, _x, :]), key=lambda x:x[1])]  # Get indices of sorted array
            top2 = closestFeatures[:2]
            valid_sequences = [s[2] for s in sequences if list(s[:2]) == list(top2)]
            if valid_sequences:
                sequence = valid_sequences[0]
                dist0 = distances[_y, _x, top2[0]]
                dist1 = distances[_y, _x, top2[1]]
                if dist0 + dist1 != 0.0:
                    heights[_y, _x] = interpolateOnCurve(sequence, dist0 / (dist0 + dist1))
                else:
                    heights[_y, _x] = sequence[0]
    print("Heightmap generated")
    return heights

def genAndSaveHeightMap(profileSketch: LineBuilder, islandSketches: SketchManagement, sliceCut: Tuple[Vector2D, Vector2D]):
    image = heightmapFromSketches(profileSketch, islandSketches, sliceCut)
    rgb = np.zeros((image.shape[0], image.shape[1], 3))
    rgb[:, :, 0] = rgb[:, :, 1] = rgb[:, :, 2] = image
    plt.imsave("test.png", rgb)

def main():
    sketch_names = ["Reef", "Island", "Passes"]
    sketch_colors = ["orange", "green", "blue"]
    waterLevel = 0.7

    # Creates a top view for island sketching and a side view for profile editing
    fig, axes = plt.subplots(1, 2, squeeze=False)
    topViewAx: plt.Axes = axes[0, 0]
    sideViewAx: plt.Axes = axes[0, 1]
    fig.subplots_adjust(bottom=0.2)

    topViewAx.set_title('Top view')
    topViewAx.set_xlim(-1, 1)
    topViewAx.set_ylim(-1, 1)
    sideViewAx.set_title('Side view')
    sideViewAx.set_xlim(-1, 1)
    sideViewAx.set_ylim(0, 1)

    # Represent the center of the map with an ellipse
    # topViewAx.add_patch(matplotlib.patches.Circle((0, 0), 0.1))

    # Add sketching for top island sketching
    islandSketches = SketchManagement(topViewAx)

    buttons: List[Button] = []
    button_size = 0.5 / len(sketch_names)
    for i_sketch in range(len(sketch_names)):
        islandSketches.addSketch(color = sketch_colors[i_sketch])
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

    # Each time the top view sketch is modified, reposition markers on the side view
    def updateSideViewMarkers():
        for line in sideViewAx.lines:
            line.set_data([], [])
        # Add water level to profile editing
        sideViewAx.axhline(y = waterLevel, color = "blue")
        for sketch in islandSketches.lineBuilders:
            intersections = sketch.intersection(*sliceCut)
            for intersect in intersections:
                sideViewAx.axvline(x = intersect.x, color = sketch.color)
        profileSketching.draw()
        fig.canvas.draw()
        fig.canvas.flush_events()

    # Add sketching for profile editing
    profileSketching: SketchManagement = SketchManagement(sideViewAx)
    profileSketching.addSketch(line = LineBuilder1D(sideViewAx, .1, "blue"))

    # Testing part
    def centeredCircle(radius:float, randomness: float = 0.2) -> List[Vector2D]:
        noise = perlin.SimplexNoise(10)
        points: List[Vector2D] = []
        nbPoints = 20
        for i in range(nbPoints + 1):
            angle = i * 2 * 3.141592 / nbPoints
            vertexUnitPos = Vector2D(math.cos(angle), math.sin(angle))
            noiseValue = noise.noise2(vertexUnitPos.x, vertexUnitPos.y) * randomness
            points.append(vertexUnitPos * (radius * (1 - noiseValue)))
        return points
    randomCoralCurve = centeredCircle(0.5)
    islandSketches.lineBuilders[0].setCurve(randomCoralCurve)
    randomIslandCurve = centeredCircle(0.25)
    islandSketches.lineBuilders[1].setCurve(randomIslandCurve)

    _randomProfileCurve = [.0, .5, .6, .7, .7, .6, .6, .9, .9, .6, .6, .7, .7, .6, .5, .0]
    randomProfileCurve = []
    for i in range(len(_randomProfileCurve)):
        randomProfileCurve.append(Vector2D((2 * i/(len(_randomProfileCurve) - 1)) - 1, _randomProfileCurve[i]))
    profileSketching.lineBuilders[0].setCurve(randomProfileCurve)

    islandSketches.onChange(updateSideViewMarkers)
    # Todo:
    #  - From profile, get heightmap depending on coral position and island position
    #  - Add passes in the process
    #  - Create heightmap from input and distance map
    #  - Profile differs depending on "r" and "d" (island width, island-coral distance)

    updateSideViewMarkers()
    # genAndSaveHeightMap(profileSketching.lineBuilders[0], islandSketches, sliceCut)
    plt.show()


if __name__ == "__main__":
    main()
