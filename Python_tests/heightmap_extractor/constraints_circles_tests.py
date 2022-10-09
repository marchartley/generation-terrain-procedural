import random
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math

class Vector2D:
    """A two-dimensional vector with Cartesian coordinates."""

    def __init__(self, x, y):
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
            return Vector2D(self.x*scalar, self.y*scalar)
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

    def __abs__(self):
        """Absolute value (magnitude) of the vector."""
        return math.sqrt(self.getX()**2 + self.getY()**2)

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def norm(self):
        return abs(self)

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
        return self.__abs__(), math.atan2(self.getY(), self.getX())

    @staticmethod
    def from_polar(length, theta):
        return Vector2D(math.cos(theta), math.sin(theta)) * length

class PolarVector2D(Vector2D):
    def __init__(self, x, y, length, theta):
        super().__init__(x, y)
        self.length = length
        self.theta = theta

    def getX(self):
        return self.x + math.cos(self.theta) * self.length

    def getY(self):
        return self.y + math.sin(self.theta) * self.length

    def angular_distance_to(self, other):
        # Return the angle only
        translatedA = Vector2D(self.getX() - self.x, self.getY() - self.y)
        translatedB = Vector2D(other.getX() - self.x, other.getY() - self.y)
        return (translatedB - translatedA).to_polar()[1]

    def distance_to(self, other):
        return math.sqrt((self.getX() - other.getX())**2 + (self.getY() - other.getY())**2)

def randomVec2():
    return PolarVector2D(random.random(), random.random(), 0, 0)

def str_distances(nodes, distance_target):
    out = "Distances:\n"
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if distance_target[i][j] > 0:
                out += str(chr(ord("A") + i)) + str(chr(ord("A") + j)) + ": " + str(round(math.sqrt((nodes[i].getX() - nodes[j].getX())**2 + (nodes[i].getY() - nodes[j].getY())**2), 3)) + "/" + str(distance_target[i][j]) + "\n"
    return out

def circle_coords(pos: Vector2D, radius: float):
    xs, ys = [], []
    nb_points = 100
    for i in range(nb_points + 1):
        theta = (i / nb_points) * 3.141592 * 2
        xs.append(math.cos(theta) * radius + pos.getX())
        ys.append(math.sin(theta) * radius + pos.getY())
    return xs, ys

def main():
    fig, ax = plt.subplots()
    plt.ion()

    distances = [
        [0.000, 2.000, 2.828, 4.472, -.100],
        [2.000, 0.000, 4.472, 4.000, -.100],
        [2.828, 4.472, 0.000, 4.472, -.100],
        [4.472, 4.000, 4.472, 0.000, 1.000],
        [-.100, -.100, -.100, -.100, 0.000],
    ]
    nodes = [randomVec2() * 100 for _ in range(len(distances))]
    for i in range(len(nodes)):
        if i == 0:
            nodes[i] = PolarVector2D(nodes[i].x, nodes[i].y, 0, 0)
        else:
            nodes[i] = PolarVector2D(nodes[i - 1].x, nodes[i - 1].y, max(distances[i][i - 1], distances[i - 1][i]), random.random() * 3.141592 * 2)

    for i in range(len(distances)):
        for j in range(i):
            distances[i][j] = -1

    epsilon = 1e-5
    tries = 0
    deltaMove = .01

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
                if nodes[i].distance_to(nodes[j]) > distances[i][j]:
                    dists[i, j] = nodes[i].angular_distance_to(nodes[j])

        for i in range(len(nodes)):
            move = 0
            for j in range(len(nodes)):
                if distances[i][j] >= 0:
                    move += dists[i][j]
            move *= deltaMove
            moves.append(move)
            nodes[i].theta += move

        plt.cla()
        plt.title(f"Step {tries + 1}, Delta = {deltaMove}")
        for i, node in enumerate(nodes):
            circle_distances = list(filter(lambda d: d > 0, sorted(distances[i])))
            for j, dist in enumerate(circle_distances):
                alpha = (dist / circle_distances[-1]) * .5 + .5
                circle_x, circle_y = circle_coords(node, dist)
                ax.plot(circle_x, circle_y, "-", c = (.0, .0, .0, alpha))

        ax.scatter([vec.getX() for vec in nodes], [vec.getY() for vec in nodes])
        for i, node in enumerate(nodes):
            ax.text(node.getX(), node.getY(), chr(ord("A") + i))
            if moves[i] > 0.00001:
                moveVec = PolarVector2D(0, 0, 1, moves[i])
                ax.arrow(node.x, node.y, moveVec.getX(), moveVec.getY())

        ax.text(0, .5, str_distances(nodes, distances), transform=ax.transAxes)
        ax.set_ylabel("Y")
        ax.set_xlabel("X")
        plt.pause(0.01)

        tries += 1
        if np.all(np.array([abs(vec) for vec in moves]) < epsilon) or tries >= 1000:
            break

        # deltaMove *= .90

    print(str_distances(nodes, distances))
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
