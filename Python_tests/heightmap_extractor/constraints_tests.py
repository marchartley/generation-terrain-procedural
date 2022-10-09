import copy
import random
import time
from collections.abc import Iterable
from typing import List, Union

import matplotlib.pyplot as plt
import numpy as np
import math


class Vector2D:
    """A two-dimensional vector with Cartesian coordinates."""

    def __init__(self, x: Union[float, Iterable] = None, y: float = None):
        if y is None:
            if x is None:
                x = 0
                y = 0
            elif isinstance(x, Iterable):
                y = x[1]
                x = x[0]
        self.x, self.y = x, y

    def __str__(self):
        """Human-readable string representation of the vector."""
        return '{:g}i + {:g}j'.format(self.x, self.y)

    def __repr__(self):
        """Unambiguous string representation of the vector."""
        return repr((self.x, self.y))

    def dot(self, other):
        """The scalar (dot) product of self and other. Both must be vectors."""

        if not isinstance(other, Vector2D):
            raise TypeError('Can only take dot product of two Vector2D objects')
        return self.x * other.x + self.y * other.y

    # Alias the __matmul__ method to dot so we can use a @ b as well as a.dot(b).
    __matmul__ = dot

    def alignedWith(self, vecA: 'Vector2D', vecB: 'Vector2D', tol: float = 0.01) -> bool:
        AB = (vecB - vecA).normalize()
        AC = (self - vecA).normalize()
        return AB.dot(AC) > 1 - tol

    def __sub__(self, other):
        """Vector subtraction."""
        return Vector2D(self.x - other.x, self.y - other.y)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __add__(self, other):
        """Vector addition."""
        return Vector2D(self.x + other.x, self.y + other.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __mul__(self, scalar):
        """Multiplication of a vector by a scalar."""

        if isinstance(scalar, int) or isinstance(scalar, float):
            return Vector2D(self.x * scalar, self.y * scalar)
        raise NotImplementedError('Can only multiply Vector2D by a scalar')

    def __imul__(self, scalar):
        """Multiplication of a vector by a scalar."""

        if isinstance(scalar, int) or isinstance(scalar, float):
            self.x *= scalar
            self.y *= scalar
            return self
        raise NotImplementedError('Can only multiply Vector2D by a scalar')

    def __rmul__(self, scalar):
        """Reflected multiplication so vector * scalar also works."""
        return self.__mul__(scalar)

    def __neg__(self):
        """Negation of the vector (invert through origin.)"""
        return Vector2D(-self.x, -self.y)

    def __truediv__(self, scalar):
        """True division of the vector by a scalar."""
        return Vector2D(self.x / scalar, self.y / scalar)

    def __mod__(self, scalar):
        """One way to implement modulus operation: for each component."""
        return Vector2D(self.x % scalar, self.y % scalar)

    def norm2(self) -> float:
        return self.x ** 2 + self.y ** 2

    def __abs__(self):
        """Absolute value (magnitude) of the vector."""
        return math.sqrt(self.norm2())

    def norm(self):
        return math.sqrt(self.norm2())

    def normalize(self):
        if self.norm() == 0:
            return self
        return self / self.norm()

    def normalized(self):
        a = Vector2D(self.x, self.y)
        return a.normalize()

    def distance_to(self, other):
        """The distance between vectors self and other."""
        return abs(self - other)

    def to_polar(self):
        """Return the vector's components in polar coordinates."""
        return self.__abs__(), math.atan2(self.y, self.x)

    def copy(self):
        return copy.deepcopy(self)

    def asarray(self) -> List[float]:
        return [self.x, self.y]

    def rotate(self, radians: float):
        newX = math.cos(radians) * self.x - math.sin(radians) * self.y
        newY = math.sin(radians) * self.x + math.cos(radians) * self.y
        self.x = newX
        self.y = newY
        return self


def randomVec2():
    return Vector2D(random.random(), random.random())


def str_distances(nodes, distance_target):
    out = "Distances:\n"
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if distance_target[i][j] > 0:
                out += str(chr(ord("A") + i)) + str(chr(ord("A") + j)) + ": " + str(
                    round((nodes[i] - nodes[j]).norm(), 3)) + "/" + str(round(distance_target[i][j], 3)) + "\n"
    return out


def circle_coords(pos: Vector2D, radius: float):
    xs, ys = [], []
    nb_points = 100
    for i in range(nb_points + 1):
        theta = (i / nb_points) * 3.141592 * 2
        xs.append(math.cos(theta) * radius + pos.x)
        ys.append(math.sin(theta) * radius + pos.y)
    return xs, ys


def id_to_str(id: int) -> str:
    return chr(ord('A') + id)


def check_feasibility(distances: list, verbose: bool = False):
    for i in range(len(distances) - 1):  # No need to check last nodes constraints
        for j in range(i + 2, len(distances)):  # No need to check self-constraints
            if distances[i][j - 1] is None or distances[j - 1][j] is None or distances[i][j] is None or \
                    distances[i][j - 1] < 0 or distances[j - 1][j] < 0 or distances[i][j] < 0:
                continue
            distA = max(distances[i][j - 1], 0)
            distB = max(distances[j - 1][j], 0)
            distC = max(distances[i][j], 0)
            if verbose:
                print(
                    f"{id_to_str(i)}{id_to_str(j - 1)}({distA}) + {id_to_str(j - 1)}{id_to_str(j)}({distB}) >= {id_to_str(i)}{id_to_str(j)}({distC}) ? {distA + distB >= distC}")
            if distA + distB < distC:
                return False
    return True


def main(display: bool = True):
    if display:
        fig, ax = plt.subplots()
        plt.ion()

    nb_nodes = 6
    # samples = poisson_disc.Bridson_sampling() # Approx. 0.2s just for Poisson Disk Sampling
    # indices = random.sample(list(range(len(samples))), nb_nodes)
    # nodes = [Vector2D(samples[i, 0], samples[i, 1]) * 100 for i in indices]
    nodes = [randomVec2() * 100 for _ in range(nb_nodes)]

    distances = [
        [0.000, 2.000, 2.828, 4.472, None, None],
        [None, 0.000, 4.472, 4.000, None, None],
        [None, None, 0.000, 4.472, None, 8.100],
        [None, None, None, 0.000, 1.000, None],
        [None, None, None, None, 0.000, 3.100],
        [None, None, None, None, None, 0.000],
    ]
    for i in range(len(distances)):
        for j in range(len(distances)):
            if distances[i][j] is None:
                distances[i][j] = -1

    for i in range(len(distances)):
        for j in range(len(distances)):
            distances[i][j] = max(distances[i][j], distances[j][i])
        for j in range(i):
            distances[i][j] = -1

    if not check_feasibility(distances, verbose=False):
        print("No config possible to get the constraints done :")
        check_feasibility(distances, verbose=True)
        return
    epsilon = 1e-2
    tries = 0
    deltaMove = .5
    deltaMoveDamping = 1.

    while True:
        dists = np.zeros((len(nodes), len(nodes)))
        nodes_pairs = []
        for i in range(len(nodes)):
            nodes_pairs.append([])
            for j in range(len(nodes)):
                nodes_pairs[i].append([])
        moves = []
        for i in range(len(nodes)):
            for j in range(len(nodes)):
                nodes_pairs[i][j] = nodes[j] - nodes[i]
                dists[i, j] = (nodes[i] - nodes[j]).norm()

        for i in range(len(nodes)):
            move = Vector2D(0, 0)
            divisor = 0
            for j in range(i, len(nodes)):
                if distances[i][j] >= 0:
                    move += nodes_pairs[i][j].normalized() * (nodes_pairs[i][j].norm() - distances[i][j])
                    divisor += 1
            move /= divisor
            if move.norm() > 0.01 and random.random() < 0.001:
                move *= deltaMove
            moves.append(move)
            nodes[i] += move

        if display:
            plt.cla()
            plt.title(f"Step {tries + 1}, Delta = {deltaMove}")
            for i, node in enumerate(nodes):
                circle_distances = distances[i][i:]
                valid_distances = [abs(distances[i][j] - dists[i][j]) < epsilon * 10 for j in
                                   range(i, len(distances[i]))]
                for j, dist in enumerate(circle_distances):
                    if dist <= 0:
                        continue
                    circle_x, circle_y = circle_coords(node, dist)
                    color = 'green' if valid_distances[j] else 'red'
                    ax.plot(circle_x, circle_y, "-", c=color)

            ax.scatter([vec.x for vec in nodes], [vec.y for vec in nodes])
            for i, node in enumerate(nodes):
                ax.text(node.x, node.y, chr(ord("A") + i))
                if moves[i].norm() > 0.00001:
                    ax.arrow(node.x, node.y, moves[i].normalize().x, moves[i].normalize().y)

            ax.text(0, .5, str_distances(nodes, distances), transform=ax.transAxes)
            ax.set_ylabel("Y")
            ax.set_xlabel("X")
            ax.set_aspect('equal')
            plt.pause(0.05)

        tries += 1
        if np.all(np.array([vec.norm() for vec in moves]) < epsilon) or tries >= 1000:
            break

        deltaMove *= deltaMoveDamping

    if display:
        print(str_distances(nodes, distances))
        plt.ioff()
        plt.show()


if __name__ == "__main__":
    main()

    nb_tests = 10
    start = time.time()
    for _ in range(nb_tests):
        main(display=False)
    end = time.time()
    print(f"Process takes {(end - start) / nb_tests}s")
