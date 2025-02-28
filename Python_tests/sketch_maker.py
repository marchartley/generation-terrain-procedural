import math
import numpy as np
import timeit
from typing import Tuple, List, Callable, Union, Optional, Any
from matplotlib import pyplot as plt
import threading
import time

import curves
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

def resizeArray(initialCurve: List[Any], desiredLength: int):
    initialLength = len(initialCurve)
    if initialLength == desiredLength or initialLength == 0:
        return initialCurve

    # Convert to NumPy array for vectorized operations
    initialCurve = np.array(initialCurve)

    # Compute the indices for the resampled array
    indices = np.linspace(0, initialLength - 1, desiredLength)

    # Get the integer parts (floor) and fractional parts of the indices
    lower_indices = np.floor(indices).astype(int)
    upper_indices = np.ceil(indices).astype(int)
    upper_indices = np.clip(upper_indices, 0, initialLength - 1)  # Clip to avoid out-of-bound indexing
    fractional_parts = indices - lower_indices

    # Perform linear interpolation
    finalCurve = (
        initialCurve[lower_indices] * (1 - fractional_parts) +
        initialCurve[upper_indices] * fractional_parts
    )

    return finalCurve.tolist()

# def resizeArray(initialCurve: List[Any], desiredLength: int):
#     initialLength = len(initialCurve)
#     if initialLength == desiredLength or initialLength == 0:
#         return initialCurve
#
#     finalCurve: List[Vector2D] = []
#     ratio = (initialLength - 1) / (desiredLength - 1)
#
#     for i in range(desiredLength):
#         ii = i * ratio
#         t = ii - int(ii)
#         finalCurve.append(initialCurve[math.floor(ii)] * (1 - t) + initialCurve[min(math.ceil(ii), initialLength - 1)] * t)
#
#     return finalCurve


class LineBuilder:
    def __init__(self, ax: plt.Axes, color: Optional[str] = None, active: bool = True, multiline: bool = False):
        _data = ax.plot([0], [0], color = color)
        self.line: plt.Line2D = _data[0]
        self.points: List[Vector2D] = []
        self.ax = ax
        self.cid_press = self.line.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_release = self.line.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_move = self.line.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.currentlyPressed = False
        self.active = active
        self.color = color
        self.callbacksOnChange: List[Callable] = []
        self.callbacksOnChangeEnded: List[Callable] = []
        self.previousMousePos: Optional[Vector2D] = None
        self.xMin, self.xMax = ax.get_xlim()
        self.yMin, self.yMax = ax.get_ylim()
        self.drawableCurve: List[Vector2D] = []

    def addPoint(self, p: Vector2D) -> 'LineBuilder':
        self.points.append(p.copy())
        self.computeCachedCurve()
        # self.draw()
        # for callback in self.callbacksOnChange:
        #     callback()
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
        if self.currentlyPressed:
            for callback in self.callbacksOnChangeEnded:
                callback()
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
        return self.drawableCurve

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilder':
        self.reset()
        for p in points:
            self.addPoint(p)
        return self

    def close(self) -> 'LineBuilder':
        self.addPoint(self.points[0])
        return self

    def intersection(self, limitA: Vector2D, limitB: Vector2D) -> List[Vector2D]:
        # Retrieve and resample the curve to reduce complexity
        curve = self.getCurve()
        if len(curve) < 2:
            return []
        if len(curve) > 30:  # Avoid redundant resizing if curve is already small
            curve = resizeArray(curve, 30)

        # Preallocate list size conservatively to reduce dynamic resizing
        intersections = []

        # Iterate through curve segments
        previous_point = curve[0]
        for current_point in curve[1:]:
            # Check for intersection between the line segment and the curve segment
            intersect = line_intersection(limitA, limitB, previous_point, current_point)
            if intersect and intersect not in intersections:  # Avoid duplicates
                intersections.append(intersect)
            previous_point = current_point  # Update for next segment

        return intersections
    # def intersection(self, limitA: Vector2D, limitB: Vector2D) -> List[Vector2D]:
    #     intersections: List[Vector2D] = []
    #     curve = self.getCurve()
    #     curve = resizeArray(curve, 30)  # Just to gain in time... should be removed later
    #     for i in range(1, len(curve)):
    #         p1, p2 = curve[i - 1], curve[i]
    #         intersect = line_intersection(limitA, limitB, p1, p2)
    #         if intersect is None or intersect in intersections:
    #             continue
    #         intersections.append(intersect)
    #     return intersections

    def addCallbackOnChange(self, function: Callable):
        self.callbacksOnChange.append(function)

    def addCallbackOnChangeEnded(self, function: Callable):
        self.callbacksOnChangeEnded.append(function)

    def draw(self, forceRedraw = True):
        curve = self.getCurve()
        c1 = curve
        c2 = curve  # self.points
        self.line.set_data([v.x for v in c1], [v.y for v in c1])
        self.line.set_linewidth(3.0 if self.active else 2.0)
        # ax: plt.Axes = self.line.axes
        # ax.scatter([p.x for p in c2], [p.y for p in c2], marker="x")
        if forceRedraw:
            self.line.figure.canvas.draw()

    def computeCachedCurve(self):
        self.drawableCurve = resizeArray(curves.catmull_rom_chain(self.points, num_points=5), 30)
        for callback in self.callbacksOnChange:
            callback()


class LineBuilder1D(LineBuilder):
    def __init__(self, ax: plt.Axes, nb_points: int = 30, color: Optional[str] = None, active: bool = True):
        super().__init__(ax, color, active)
        self.previousMousePos: Optional[Vector2D] = None
        self.nb_points = nb_points
        self.reset()

    def addPoint(self, p: Vector2D) -> 'LineBuilder1D':
        index = self.getClosestPointIndex(p)
        self.points[index] = p.copy()
        self.computeCachedCurve()
        # for callback in self.callbacksOnChange:
        #     callback()
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
        return self

    def setCurve(self, points: Union[List[float], List[Vector2D]]) -> 'LineBuilder1D':
        self.reset()
        points = resizeArray(points, self.nb_points)
        # self.points = []
        # for i, p in enumerate(points):
        #     if isinstance(p, Vector2D):
        #         self.points.append(p.copy())
        #     else:
        #         self.points.append(Vector2D(self.xMin + (i * (self.xMax - self.xMin) / (len(points) + 1)), p))
        self.points = [p if isinstance(p, Vector2D) else Vector2D(self.xMin + (i * (self.xMax - self.xMin) / (len(points) + 1)), p) for i, p in enumerate(points)]
        self.computeCachedCurve()
        return self

    def close(self) -> 'LineBuilder1D':
        raise NotImplementedError("Cannot close a 1D line")

    def draw(self, forceRedraw = True):
        return super().draw(forceRedraw)


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

    def computeCachedCurve(self):
        self.drawableCurve = resizeArray(curves.catmull_rom_chain([self.polarToCartesian(p) for p in self.points], closed=True), 20)
        for callback in self.callbacksOnChange:
            callback()

    def addPoint(self, p: Vector2D) -> 'LineBuilderRadial':
        index = self.getClosestPointIndex(self.cartesianToPolar(p))
        self.points[index].y = self.cartesianToPolar(p).y
        self.computeCachedCurve()
        # for callback in self.callbacksOnChange:
        #     callback()
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

            for callback in self.callbacksOnChangeEnded:
                callback()

    def reset(self) -> 'LineBuilderRadial':
        nbPoints = self.nb_points
        self.points = [Vector2D(mapTo(i / nbPoints, 0, 1, self.xMin, self.xMax), 1) for i in range(nbPoints)]
        self.computeCachedCurve()
        # for callback in self.callbacksOnChange:
        #     callback()
        return self

    def getCurve(self) -> List[Vector2D]:
        return self.drawableCurve
        # curve = curves.catmull_rom_chain([self.polarToCartesian(p) for p in self.points], closed=True)  # + [self.polarToCartesian(self.points[0])], closed=True)
        # return curve

    def setCurve(self, points: List[Vector2D]) -> 'LineBuilderRadial':
        self.reset()
        newPoints = [p.to_polar() for p in resizeArray(points, self.nb_points)]
        newPoints.sort(key = lambda v: v.x)
        for i, p in enumerate(newPoints):
            self.points[i].y = p.y
        self.computeCachedCurve()
        return self

    def getValue(self, angle: float) -> float:
        angle = math.fmod(angle + 2 * math.pi, 2 * math.pi)
        closestBefore, closestAfter = self.points[-1], None
        closestBeforeID = -1
        for i, p in enumerate(self.points):  # Points are defined as (theta, r)
            if p.x <= angle:
                closestBefore = p
                closestBeforeID = i
        closestAfter = self.points[(closestBeforeID + 1) % len(self.points)]
        t = (angle - closestBefore.x) / ((closestAfter.x if closestBeforeID+1 < len(self.points) else 2*math.pi) - closestBefore.x)
        return closestBefore.y + t * (closestAfter.y - closestBefore.y)

    def draw(self, forceRedraw = True):
        return super().draw(forceRedraw)

class SketchManagement:
    def __init__(self, ax):
        self.lineBuilders: List[LineBuilder] = []
        self.ax = ax

    def addSketch(self, color: Optional[str] = None, sketch_type = LineBuilder, line: Optional[LineBuilder] = None) -> LineBuilder:
        if line is None:
            line = sketch_type(self.ax, color = color)
        self.lineBuilders.append(line)
        self.activate(0)
        return line

    def activate(self, builder: int):
        for lb in self.lineBuilders:
            lb.active = False
        self.lineBuilders[builder].active = True

    def onChange(self, function: Callable):
        for line in self.lineBuilders:
            line.addCallbackOnChange(function)

    def onChangeEnded(self, function: Callable):
        for line in self.lineBuilders:
            line.addCallbackOnChangeEnded(function)

    def draw(self, forceRedraw = True):
        for line in self.lineBuilders:
            line.draw(forceRedraw)

    def addSketchs(self, sketches: List[LineBuilder]):
        for sketch in sketches:
            self.addSketch(line = sketch)


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