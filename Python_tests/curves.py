from copy import copy
from typing import List
import numpy as np
from Vectors import Vector2D

def flatten(list_of_lists) -> list:
    # E.g. mapping [[1, 2], [3], [4, 5]] to  [1, 2, 3, 4, 5]
    return [elem for lst in list_of_lists for elem in lst]


def catmull_rom_spline(P0: Vector2D, P1: Vector2D, P2: Vector2D, P3: Vector2D, num_points: int, alpha: float = 1.0,
):
    """
    Compute the points in the spline segment
    :param P0, P1, P2, and P3: The (x,y) point pairs that define the Catmull-Rom spline
    :param num_points: The number of points to include in the resulting curve segment
    :param alpha: 0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.
    :return: The points
    """

    # Calculate t0 to t4. Then only calculate points between P1 and P2.
    # Reshape linspace so that we can multiply by the points P0 to P3
    # and get a point for each value of t.
    def tj(ti: float, pi: Vector2D, pj: Vector2D) -> float:
        delta = pj - pi
        l = (delta.x ** 2 + delta.y ** 2) ** 0.5
        return ti + l ** alpha

    t0: float = 0.0
    t1: float = tj(t0, P0, P1)
    t2: float = tj(t1, P1, P2)
    t3: float = tj(t2, P2, P3)
    t = np.linspace(t1, t2, num_points).reshape(num_points, 1)

    A1 = (t1 - t) / (t1 - t0) * P0 + (t - t0) / (t1 - t0) * P1
    A2 = (t2 - t) / (t2 - t1) * P1 + (t - t1) / (t2 - t1) * P2
    A3 = (t3 - t) / (t3 - t2) * P2 + (t - t2) / (t3 - t2) * P3
    B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
    B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
    points = (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2
    return [p[0] for p in points] # My conversion for using Vector2Ds


def catmull_rom_chain(points: List[Vector2D], num_points: int = -1, closed: bool = False) -> List[Vector2D]:
    """
    Calculate Catmull-Rom for a sequence of initial points and return the combined curve.
    :param points: Base points from which the quadruples for the algorithm are taken
    :param num_points: The number of points to include in each curve segment
    :return: The chain of all points (points of all segments)
    """
    if num_points == -1:
        num_points = 5
    if closed:
        usedList = [points[-1]] + points + [points[0], points[1]]
    else:
        if len(points) < 3:
            return points
        usedList = [points[0] - (points[1] - points[0])*0.01] + points + [points[-1] + (points[-2] - points[-1])*0.01]
    point_quadruples = [  # Prepare function inputs
        [usedList[idx_segment_start + d] for d in range(4)]
        for idx_segment_start in range(len(usedList) - 3)
    ]
    all_splines = (catmull_rom_spline(pq[0], pq[1], pq[2], pq[3], num_points) for pq in point_quadruples)
    return flatten(all_splines)


def smoothCurve(curve: List[Vector2D]) -> List[Vector2D]:
    res = []
    res.append(curve[0])
    for i in range(len(curve) - 1):
        res.append((curve[i] + curve[i + 1]) * .5)
    res.append(curve[-1])
    return res

def subcurve(curve: List[Vector2D], t0: float, t1: float) -> List[Vector2D]:
    return curve[int(t0 * len(curve)): int(t1 * len(curve))]

def main():
    import matplotlib.pyplot as plt
    # points = [Vector2D(x, y) for x, y in [(0, 1.5), (2, 2), (3, 1), (4, 0.5), (5, 1), (6, 2), (7, 3)]]
    points = [Vector2D(x, y) for x, y in [(0, 0), (1, 1), (2, 0)]]
    # numPoints = 100
    chain = catmull_rom_chain(points, closed = True) #, numPoints)
    plt.plot([p.x for p in chain], [p.y for p in chain], c="blue")
    plt.plot([p.x for p in points], [p.y for p in points], c="red", linestyle="none", marker="o")
    plt.show()


if __name__ == "__main__":
    main()
