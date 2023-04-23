import math
import timeit
from typing import Tuple, List, Callable, Union, Optional, Any
from matplotlib import pyplot as plt
import threading
import time
from Vectors import Vector2D, Vector3D, line_intersection


def clamp(x, mini, maxi):
    return mini if x < mini else maxi if x > maxi else x

def mapTo(x, prevMin, prevMax, newMin, newMax):
    t = (x - prevMin) / (prevMax - prevMin)
    return newMin + t * (newMax - newMin)

def frange(start: float, end: float, step: float = 1.0):
    while start < end:
        yield start
        start += step

def curveInterpolation(initialCurve: List[Any], desiredLength: int):
    initialLength = len(initialCurve)
    if initialLength == desiredLength:
        return initialCurve

    finalCurve: List[Vector2D] = []
    ratio = (initialLength - 1) / (desiredLength - 1)

    for i in range(desiredLength):
        ii = i * ratio
        t = ii - int(ii)
        finalCurve.append(initialCurve[math.floor(ii)] * (1 - t) + initialCurve[min(math.ceil(ii), initialLength - 1)] * t)

    return finalCurve

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
    def __init__(self, ax: plt.Axes, color: Optional[str] = None, active: bool = True, multiline: bool = False):
        self.line, = ax.plot([0], [0], color = color)
        self.points: List[Vector2D] = []
        self.ax = ax
        self.cid_press = self.line.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_release = self.line.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_move = self.line.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.currentlyPressed = False
        self.active = active
        self.color = color
        self.callbacksOnChange: List[Callable] = []
        self.previousMousePos: Optional[Vector2D] = None
        self.xMin, self.xMax = ax.get_xlim()
        self.yMin, self.yMax = ax.get_ylim()

    def addPoint(self, p: Vector2D) -> 'LineBuilder':
        self.points.append(p.copy())
        # self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def press_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        self.addPoint(Vector2D(event.xdata, event.ydata))
        self.currentlyPressed = True
        self.previousMousePos = Vector2D(event.xdata, event.ydata)
        # print(self, f"Pressed (active = {self.active}, pressed = {self.currentlyPressed})")
        # self.draw()

    def release_event(self, event):
        if event.inaxes != self.line.axes:
            return
        self.currentlyPressed = False
        self.previousMousePos = None
        # print(self, f"Release (active = {self.active}, pressed = {self.currentlyPressed})")

    def move_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        if self.currentlyPressed:
            self.addPoint(Vector2D(event.xdata, event.ydata))
            self.previousMousePos = Vector2D(event.xdata, event.ydata)
            # self.draw()

    def reset(self) -> 'LineBuilder':
        self.points = []
        # self.draw()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def getCurve(self) -> List[Vector2D]:
        return self.points

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilder':
        self.reset()
        for p in points:
            self.addPoint(p)
        return self

    def close(self) -> 'LineBuilder':
        self.addPoint(self.points[0])
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
        curve = self.getCurve()
        self.line.set_data([v.x for v in curve], [v.y for v in curve])
        self.line.figure.canvas.draw()


class LineBuilder1D(LineBuilder):
    def __init__(self, ax: plt.Axes, nb_points: int = 30, color: Optional[str] = None, active: bool = True):
        super().__init__(ax, color, active)
        self.previousMousePos: Optional[Vector2D] = None
        self.nb_points = nb_points
        self.reset()

    def addPoint(self, p: Vector2D) -> 'LineBuilder1D':
        index = self.getClosestPointIndex(p)
        self.points[index] = p.copy()
        for callback in self.callbacksOnChange:
            callback()
        return self

    def move_event(self, event):
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        if self.currentlyPressed:
            mousePos = Vector2D(event.xdata, event.ydata)
            nbPointsToSet = 4*math.ceil(self.nb_points * abs(mousePos.x - self.previousMousePos.x) / (self.xMax - self.xMin))
            for i in range(nbPointsToSet):
                t = i / nbPointsToSet
                p = self.previousMousePos + (mousePos - self.previousMousePos) * t
                self.addPoint(p)
            self.previousMousePos = Vector2D(event.xdata, event.ydata)

    def getClosestPointIndex(self, v: Vector2D) -> int:
        index = -1
        minDist = math.inf
        for i, p in enumerate(self.points):
            if abs(p.x - v.x) < minDist:
                index = i
                minDist = abs(p.x - v.x)
        return index

    def reset(self) -> 'LineBuilder1D':
        nbPoints = self.nb_points
        self.points = [Vector2D(self.xMin + (i / (nbPoints-1)) * (self.xMax - self.xMin), 0) for i in range(nbPoints)]
        for callback in self.callbacksOnChange:
            callback()
        return self

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilder1D':
        self.reset()
        newPoints = curveInterpolation(points, self.nb_points)
        for i, p in enumerate(newPoints):
            self.points[i].y = p.y
        return self

    def close(self) -> 'LineBuilder1D':
        raise NotImplementedError("Cannot close a 1D line")


class LineBuilder2D(LineBuilder):
    def __init__(self, ax: plt.Axes, color: Optional[str] = None, active: bool = True):
        super().__init__(ax, color, active)


def angleDistance(v0: Vector2D, v1: Vector2D) -> float:
    d0 = abs(v0.x - v1.x)
    d1 = abs(v0.x - (v1.x + math.tau))
    return min(d0, d1)


class LineBuilderRadial(LineBuilder1D):
    def __init__(self, ax: plt.Axes, nb_points = 30, color: Optional[str] = None, active: bool = True, center = Vector2D(0, 0)):
        self.center = center
        self.xMin, self.xMax = 0, math.tau
        super().__init__(ax, nb_points = nb_points, color = color, active = active)
        self.center = center
        self.xMin, self.xMax = 0, math.tau
        self.reset()

    def close(self) -> 'LineBuilderRadial':
        return self

    def cartesianToPolar(self, v: Vector2D) -> Vector2D:
        return (v - self.center).to_polar()

    def polarToCartesian(self, v: Vector2D) -> Vector2D:
        return v.to_cartesian() + self.center

    def getClosestPointIndex(self, v: Vector2D) -> int:
        index = -1
        minDist = math.inf
        for i, p in enumerate(self.points):
            if angleDistance(p, v) < minDist:
                index = i
                minDist = angleDistance(p, v)
        return index

    def addPoint(self, p: Vector2D) -> 'LineBuilderRadial':
        index = self.getClosestPointIndex(self.cartesianToPolar(p))
        self.points[index].y = self.cartesianToPolar(p).y
        for callback in self.callbacksOnChange:
            callback()
        return self

    def move_event(self, event):
        pPos = None if self.previousMousePos is None else self.previousMousePos.copy()
        if not self.active:
            return
        if event.inaxes != self.line.axes:
            return
        if self.currentlyPressed and pPos is not None:
            mousePos = Vector2D(event.xdata, event.ydata)
            nbPointsToSet = math.ceil(self.nb_points * abs(self.cartesianToPolar(mousePos).x - self.cartesianToPolar(pPos).x) / math.tau)
            for i in range(nbPointsToSet):
                t = i / nbPointsToSet
                p = pPos + (mousePos - pPos) * t
                self.addPoint(p)
            self.previousMousePos = Vector2D(event.xdata, event.ydata)

    def reset(self) -> 'LineBuilderRadial':
        nbPoints = self.nb_points
        self.points = [Vector2D(mapTo(i / nbPoints, 0, 1, self.xMin, self.xMax), 1) for i in range(nbPoints)]
        for callback in self.callbacksOnChange:
            callback()
        return self

    def getCurve(self) -> List[Vector2D]:
        curve = [self.polarToCartesian(p) for p in self.points] + [self.polarToCartesian(self.points[0])]
        return curve

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilderRadial':
        self.reset()
        newPoints = [p.to_polar() for p in curveInterpolation(points, self.nb_points)]
        for i, p in enumerate(newPoints):
            self.points[i].y = p.y
        return self

class SketchManagement:
    def __init__(self, ax):
        self.lineBuilders: List[LineBuilder] = []
        self.ax = ax

    def addSketch(self, color: Optional[str] = None, sketch_type = LineBuilder, line: Optional[LineBuilder] = None) -> int:
        if line is None:
            line = sketch_type(self.ax, color = color)
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

def main():
    fig, ax = plt.subplots(1, 1)
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    # rad = LineBuilderRadial(ax, 0.01, center=Vector2D(0.5, 0.0))
    rad = LineBuilder1D(ax, nb_points=30)

    rad.draw()
    plt.show()


if __name__ == "__main__":
    main()