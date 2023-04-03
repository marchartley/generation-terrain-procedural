import random
from typing import Union, Iterable, List, Optional
import math
import copy

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




class Vector3D:
    """A 3-dimensional vector with Cartesian coordinates."""

    def __init__(self, x: Union[float, Iterable] = None, y: float = None, z: float = None):
        if y is None and z is None:
            if x is None:
                x = 0
                y = 0
                z = 0
            elif isinstance(x, Iterable):
                x = x[0]
                y = x[1]
                z = x[2]
        self.x, self.y, self.z = x, y, z

    def __str__(self):
        """Human-readable string representation of the vector."""
        return f"[{self.x}, {self.y}, {self.z}]"

    def __repr__(self):
        """Unambiguous string representation of the vector."""
        return repr((self.x, self.y, self.z))

    def dot(self, other):
        """The scalar (dot) product of self and other. Both must be vectors."""

        if not isinstance(other, Vector2D):
            raise TypeError('Can only take dot product of two Vector2D objects')
        return self.x * other.x + self.y * other.y + self.z * other.z

    # Alias the __matmul__ method to dot so we can use a @ b as well as a.dot(b).
    __matmul__ = dot

    def alignedWith(self, vecA: 'Vector3D', vecB: 'Vector3D', tol: float = 0.01) -> bool:
        AB = (vecB - vecA).normalize()
        AC = (self - vecA).normalize()
        return AB.dot(AC) > 1 - tol

    def __sub__(self, other):
        """Vector subtraction."""
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __add__(self, other):
        """Vector addition."""
        return Vector3D(self.x + other.x, self.y + other.y, self.z + other.z)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __mul__(self, scalar):
        """Multiplication of a vector by a scalar."""

        if isinstance(scalar, int) or isinstance(scalar, float):
            return Vector3D(self.x * scalar, self.y * scalar, self.z * scalar)
        raise NotImplementedError('Can only multiply Vector3D by a scalar')

    def __imul__(self, scalar):
        """Multiplication of a vector by a scalar."""

        if isinstance(scalar, int) or isinstance(scalar, float):
            self.x *= scalar
            self.y *= scalar
            self.z *= scalar
            return self
        raise NotImplementedError('Can only multiply Vector3D by a scalar')

    def __rmul__(self, scalar):
        """Reflected multiplication so vector * scalar also works."""
        return self.__mul__(scalar)

    def __neg__(self):
        """Negation of the vector (invert through origin.)"""
        return Vector3D(-self.x, -self.y, -self.z)

    def __truediv__(self, scalar):
        """True division of the vector by a scalar."""
        return Vector3D(self.x / scalar, self.y / scalar, self.z / scalar)

    def __mod__(self, scalar):
        """One way to implement modulus operation: for each component."""
        return Vector3D(self.x % scalar, self.y % scalar, self.z % scalar)

    def norm2(self) -> float:
        return self.x ** 2 + self.y ** 2 + self.z ** 2

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
        a = Vector3D(self.x, self.y, self.z)
        return a.normalize()

    def distance_to(self, other):
        """The distance between vectors self and other."""
        return abs(self - other)

    def to_polar(self):
        """Return the vector's components in polar coordinates."""
        raise NotImplementedError("Cannot get polar coordinates from a 3D-vector for now")
        # return self.__abs__(), math.atan2(self.y, self.x)

    def copy(self):
        return copy.deepcopy(self)

    def asarray(self) -> List[float]:
        return [self.x, self.y, self.z]

    def rotate(self, radians: float):
        raise NotImplementedError("Cannot apply rotation on a 3D vector for now")
        # newX = math.cos(radians) * self.x - math.sin(radians) * self.y
        # newY = math.sin(radians) * self.x + math.cos(radians) * self.y
        # self.x = newX
        # self.y = newY
        # return self


def randomVec3():
    return Vector3D(random.random(), random.random(), random.random())


def line_intersection(P11, P12, P21, P22, epsilon: float = 0.0) -> Optional[List[float]]:
    if isinstance(P11, Vector2D):
        P11 = P11.asarray()
    if isinstance(P12, Vector2D):
        P12 = P12.asarray()
    if isinstance(P21, Vector2D):
        P21 = P21.asarray()
    if isinstance(P22, Vector2D):
        P22 = P22.asarray()

    divisor = (P11[0] - P12[0]) * (P21[1] - P22[1]) - (P11[1] - P12[1]) * (P21[0] - P22[0])
    if divisor == 0:
        return None
    t = ((P11[0] - P21[0]) * (P21[1] - P22[1]) - (P11[1] - P21[1]) * (P21[0] - P22[0])) / divisor
    u = ((P11[0] - P21[0]) * (P11[1] - P12[1]) - (P11[1] - P21[1]) * (P11[0] - P12[0])) / divisor

    # check if line actually intersect
    if (epsilon != 0.0 and 0 + epsilon <= t <= 1 - epsilon and 0 + epsilon <= u <= 1 - epsilon) or (
            epsilon == 0.0 and 0 <= t <= 1 and 0 <= u <= 1):
        return [P11[0] + t * (P12[0] - P11[0]), P11[1] + t * (P12[1] - P11[1])]
    else:
        return None


def closestPointToLine(point: Vector2D, lineA: Vector2D, lineB: Vector2D, limitedToSegment: bool = True) -> Vector2D:
    a = (point - lineA)
    b = (lineB - lineA)
    t = ((a.dot(b)) / (b.dot(b)))
    return (t if not limitedToSegment else min(1, max(0, t))) * b + lineA

def distanceToLine(point: Vector2D, lineA: Vector2D, lineB: Vector2D, limitedToSegment: bool = True) -> float:
    return (point - closestPointToLine(point, lineA, lineB, limitedToSegment)).norm()
