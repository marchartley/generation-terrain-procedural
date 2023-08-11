"""import numpy as np
from matplotlib import pyplot as plt

import numpy as np


def g(X, p):
    d = np.linalg.norm(p)
    if d < 1:
        return (3 * X / np.pi) * (1 - d)
    else:
        return 0


def monte_carlo_integral(X, num_samples=1000000):
    total_value = 0
    count_inside = 0

    for _ in range(num_samples):
        # Randomly sample a point in the bounding box [-1, 1] x [-1, 1] x [-1, 1]
        p = np.array([2 * np.random.rand() - 1 for _ in range(3)])

        # Check if the point is inside the metaball
        if np.linalg.norm(p) < 1:
            count_inside += 1
            total_value += g(X, p)

    # The volume of the bounding box is 2^3 = 8
    # The volume of the unit sphere (metaball) is (4/3) * pi
    # The ratio of the volumes is used to scale the integral approximation
    volume_ratio = (4 / 3) * np.pi / 8
    average_value = total_value / count_inside

    return average_value * (count_inside / num_samples) * 8


# Test the Monte Carlo approximation
X_test = 10000
approximated_value = monte_carlo_integral(X_test)
print(f"For X = {X_test}, the approximated integral value is: {approximated_value}")

# You can easily test other values of X by calling the monte_carlo_integral function with different X values.

"""


"""def g(X, p):
    d = np.linalg.norm(p)
    if d < 1:
        return (3 * X / np.pi) * (1 - d)
    else:
        return 0


def monte_carlo_integral(X, num_samples=100000):
    total_value = 0

    for _ in range(num_samples):
        # Randomly sample a point in the bounding box [-1, 1] x [-1, 1] x [-1, 1]
        p = np.array([2 * np.random.rand() - 1 for _ in range(3)])

        # Check if the point is inside the metaball
        if np.linalg.norm(p) < 1:
            total_value += g(X, p)

    # The volume of the bounding box is 2^3 = 8
    # The volume of the unit sphere (metaball) is (4/3) * pi
    # The ratio of the volumes is used to scale the integral approximation
    volume_ratio = (4 / 3) * np.pi / 8
    average_value = total_value / num_samples

    return average_value * 8 / volume_ratio  # Scale by the ratio of the volumes


def visualize_g(X, grid_size=100):
    # Create a 2D grid
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    X_grid, Y_grid = np.meshgrid(x, y)

    # Evaluate g on the grid
    Z = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            Z[i, j] = g(X, np.array([X_grid[i, j], Y_grid[i, j], 0]))

    # Display the grid using matplotlib
    plt.imshow(Z, extent=(-1, 1, -1, 1), origin='lower', cmap='viridis')
    plt.colorbar(label='g(X, p)')
    plt.title(f'Visualization of g for X = {X}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

# Test the Monte Carlo approximation
X_test = 100000
visualize_g(X_test)
approximated_value = monte_carlo_integral(X_test)
print(f"For X = {X_test}, the approximated integral value is: {approximated_value}")

# You can easily test other values of X by calling the monte_carlo_integral function with different X values.
"""





"""
import numpy as np


def f(p):
    d = np.linalg.norm(p)
    if d < 1:
        return 1 - d
    return 0


def monte_carlo_integration(samples=1000000):
    count_inside_sphere = 0
    total_value = 0

    for _ in range(samples):
        # Randomly sample a point in the unit cube [-1, 1] x [-1, 1] x [-1, 1]
        p = np.array([2 * np.random.rand() - 1, 2 * np.random.rand() - 1, 2 * np.random.rand() - 1])

        # Check if the point is inside the unit sphere
        if np.linalg.norm(p) < 1:
            count_inside_sphere += 1
            total_value += f(p)

    # The volume of the unit cube is 8 (2x2x2)
    # The estimated volume of the unit sphere is 8 * (count_inside_sphere / samples)
    # Multiply this by the average value of f(p) for points inside the sphere
    integral_estimate = (8 * count_inside_sphere / samples) * (total_value / count_inside_sphere)

    return integral_estimate


# Compute the numerical integral
result = monte_carlo_integration()
print(f"Numerical estimate of the integral: {result}")
"""





import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def g(X, p):
    d = np.linalg.norm(p)
    if d < 1:
        return (3 * X / np.pi) * (1 - d)
    else:
        return 0

def ball(nX) -> np.ndarray:
    center = [nX/2, nX/2, nX/2]
    r = nX / 2
    arr = [[[abs(x - center[0])**2 + abs(y - center[1])**2 + abs(z - center[2]) for z in range(nX)] for x in range(nX)] for y in range(nX)]
    return np.maximum(1 - np.sqrt(np.array(arr)) / r, 0)

def ball2D(nX) -> np.ndarray:
    center = [nX/2, nX/2, nX/2]
    r = nX / 2
    arr = [[abs(x - center[0])**2 + abs(y - center[1])**2 for x in range(nX)] for y in range(nX)]
    return np.maximum(1 - np.sqrt(np.array(arr)) / r, 0)
    # arr = np.array(arr)
    # arr[:, :, :] = (1.0 / np.maximum(arr[:, :, :], 1))**2
    # return arr

def ballMidpoint(nX) -> np.ndarray:
    center = [nX/2, nX/2, nX/2]
    r = nX / 2
    arr = [[[g(1, ((np.array([x, y, z]) - center) / r)) for z in range(nX)] for x in range(nX)] for y in range(nX)]
    return np.array(arr) / (r**3)

def main():
    nX = 30
    b = ball(nX)
    # b /= nX**2
    Q = 100
    # b *= Q
    voxels = (b * Q)
    total = voxels.sum()
    voxels /= total
    print(total)
    midpoint = ballMidpoint(nX) * Q
    diff = voxels - midpoint

    print(voxels.sum(), midpoint.sum(), diff.sum())

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, squeeze=True)
    ax1.imshow(voxels[:, :, nX//3])
    ax2.imshow(midpoint[:, :, nX//3])
    ax3.imshow(diff[:, :, nX//3])
    plt.show()
    return
    x = [i for i in range(2, 50)]
    y = []
    estimation = 3 / math.pi
    for i in x:
        b = ball(i)
        inside = b > 0.0
        y.append(b.sum() / inside.sum())
        if i % 10 == 0:
            print(f"{i}: \t{y[-1]}")
    plt.plot(x, np.array(y) / estimation)
    # plt.yscale("log")
    plt.show()

if __name__ == "__main__":
    main()
