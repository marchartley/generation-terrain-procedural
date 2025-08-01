import math

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


def rotation_matrix(theta, phi, degrees=True):
    """
    Construct rotation matrix that rotates the ellipsoid's z-axis to the (theta, phi) direction.
    theta: azimuthal angle (rotation around z-axis)
    phi: inclination from z-axis (tilt)
    """
    if degrees:
        theta = np.radians(theta)
        phi = np.radians(phi)

    # Rotation around z (azimuth)
    Rz = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])

    # Rotation around y (inclination)
    Ry = np.array([
        [np.cos(phi), 0, np.sin(phi)],
        [0, 1, 0],
        [-np.sin(phi), 0, np.cos(phi)]
    ])

    # Final rotation: R = Rz * Ry
    return Rz @ Ry


def ellipsoid_coefficients(a, b, c, R):
    """
    Returns the coefficients of the quadratic form of the rotated ellipsoid.
    The ellipsoid is defined by x.T @ Q @ x = 1, where Q = R @ A @ R.T
    """
    A = np.diag([1 / a ** 2, 1 / b ** 2, 1 / c ** 2])
    Q = R @ A @ R.T  # The rotated quadric matrix
    return Q

def ellipsoid_roots(x, y, Q, R, z0=0.0, inverse_cut=False):
    """
    Given (x, y), return (z_lower_clipped, z_upper_clipped), where clipping occurs
    at the plane z' = z0 in the local frame, intersecting a global vertical ray.
    """
    Qxx, Qxy, Qxz = Q[0, 0], Q[0, 1], Q[0, 2]
    Qyy, Qyz = Q[1, 1], Q[1, 2]
    Qzz = Q[2, 2]

    # Solve quadratic in z for vertical line at (x, y)
    A = Qzz
    B = 2 * (Qxz * x + Qyz * y)
    C = Qxx * x**2 + 2 * Qxy * x * y + Qyy * y**2 - 1

    D = B**2 - 4 * A * C
    if D < 0:
        return (np.nan, np.nan)

    sqrt_D = np.sqrt(D)
    z1 = (-B + sqrt_D) / (2 * A)
    z2 = (-B - sqrt_D) / (2 * A)

    z_upper = max(z1, z2)
    z_lower = min(z1, z2)


    # Compute intersection of vertical line with the cutting plane z' = z0
    # Find z such that (R^T @ [x, y, z])[2] == z0
    # => row 3 of R^T dotted with [x, y, z] == z0
    # => R[2,0] * x + R[2,1] * y + R[2,2] * z_cut = z0
    # Solve for z_cut
    v = R.T[2]  # third row of R^T == third column of R # Don't know why I have to translate there...
    if abs(v[2]) < 1e-5:
        # Avoid division by near-zero: plane is vertical
        return (z_lower, z_upper)  # no intersection with vertical ray

    z_cut = (z0 - v[0]*x - v[1]*y) / v[2]

    if inverse_cut:
        z_upper, z_lower, z_cut = -z_lower, -z_upper, -z_cut

    # Decide whether to replace upper or lower based on orientation
    if v[2] > 0:
        if z_cut > z_upper:
            z_lower, z_upper = np.nan, np.nan
        else:
            # local z' points upward in global space — keep upper, clip lower
            z_lower = max(z_cut, z_lower)
            # local z' points downward in global space — keep lower, clip upper
            z_upper = max(z_cut, z_upper)
    elif z_cut < z_upper:
        if z_cut < z_lower:
            z_lower, z_upper = np.nan, np.nan
        else:
            # local z' points upward in global space — keep upper, clip lower
            z_lower = min(z_cut, z_lower)
            # local z' points downward in global space — keep lower, clip upper
            z_upper = min(z_cut, z_upper)

    if inverse_cut:
        z_upper, z_lower, z_cut = -z_lower, -z_upper, -z_cut

    return (z_lower, z_upper)

def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm

def rotation_from_target(target_x: np.array, target_y: np.array, target_z: np.array):
    return np.column_stack([target_x, target_y, target_z])


def h_plus(x, y, Q, R, z0, inverse_cut=False):
    return ellipsoid_roots(x, y, Q, R, z0, inverse_cut)[1]
def h_minus(x, y, Q, R, z0, inverse_cut=False):
    return ellipsoid_roots(x, y, Q, R, z0, inverse_cut)[0]
def h(x, y, Q, R, z0, inverse_cut=False):
    plus = h_plus(x, y, Q, R, z0, inverse_cut)
    minus = h_minus(x, y, Q, R, z0, inverse_cut)
    return 0 if plus != plus else plus - minus


def height(x, y):
    octaves = 5
    lacun = 2.0
    gain = 1.5
    return (sum([(math.cos(x * lacun**i * 0.1) + math.sin(15 + y * lacun**i)) / (gain**i) for i in range(octaves)]) / (2.0 * (1 - (1/gain)**(octaves + 1)) * 2.0)) * 0.5 + 0.5

# Ellipsoid and rotation setup
a, b, c = 0.250, 0.250, 0.150  # semi-axes
theta = 0               # fixed azimuthal rotation

# Grid
x_range = np.linspace(-1.1, 1.1, 100)
y_range = np.linspace(-1.1, 1.1, 100)
X, Y = np.meshgrid(x_range, y_range)


# --- Animation Setup ---
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf1 = [None]
surf2 = [None]
surf3 = [None]

cx, cy = 0, 0

speedX = 0*0.05
speedY = 0*0.02

shift1 = 20
shift2 = -25

def update(frame):
    global surf1, surf2, surf3, cx, cy
    ax.clear()
    z0 = -0.0
    phi = frame  # degrees

    cx += speedX
    cy += speedY

    def rotation_and_ellipsoid_at(cx, cy):
        center_height = height(cx, cy)
        dx = (center_height - height(cx + 0.01, cy)) / 0.01
        dy = (center_height - height(cx, cy + 0.01)) / 0.01
        target_X = normalize(np.array([1, 0, -dx])) # - dx / np.linalg.norm(dx)
        target_Y = normalize(np.array([0, 1, -dy])) # - dy / np.linalg.norm(dy)
        target_Z = np.cross(target_X, target_Y)

        R = rotation_from_target(target_X, target_Y, target_Z)
        # R = rotation_matrix(theta, phi)
        Q = ellipsoid_coefficients(a, b, c, R)
        return R, Q

    center_height = height(cx, cy)
    R1, Q1 = rotation_and_ellipsoid_at(cx + shift1, cy)
    R2, Q2 = rotation_and_ellipsoid_at(cx + shift2, cy)
    R, Q = rotation_and_ellipsoid_at(cx, cy)
    # Z1 = np.array([[h_minus(x, y, Q, R, z0, False) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
    # Z2 = np.array([[h_plus(x, y, Q, R, z0, True) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
    Z3 = np.array([[h(x, y, Q1, R1, z0, False) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
    Z4 = np.array([[h(x, y, Q2, R2, z0, True) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])

    # surf1[0] = ax.plot_surface(X, Y, Z1, cmap='viridis', edgecolor='none', alpha=0.9)
    # surf2[0] = ax.plot_surface(X, Y, Z2, cmap='plasma', edgecolor='none', alpha=0.9)

    # ax.plot_surface(X, Y, Z3)


    Surface = np.array([[height(cx + x * 2.0, cy + y * 2.0) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])
    # Z1[Z1 != Z1] = 0
    # Z2[Z2 != Z2] = 0
    Z3[Z3 != Z3] = -np.inf
    Z4[Z4 != Z4] = -np.inf

    # surf3[0] = ax.plot_surface(X, Y, np.maximum(Surface, (Z2 - Z1)) - 1.5)
    surf3[0] = ax.plot_surface(X, Y, (Surface - np.roll(Z3, shift1, axis=0) + np.roll(Z4, shift2, axis=0)) - 1.5)

    # surf3[0] = ax.plot_surface(X, Y, np.maximum(Surface, np.roll(Z3, shift1, axis=0) + height(cx + shift1 * 2.0, cy)) - 1.5)
    # -- Cutting Plane: z' = z0 in local frame, transformed to global
    plane_size = 1.2
    x_plane = np.linspace(-plane_size, plane_size, 20)
    y_plane = np.linspace(-plane_size, plane_size, 20)
    X_plane, Y_plane = np.meshgrid(x_plane, y_plane)
    Z_plane_local = np.full_like(X_plane, z0)

    # Stack local coordinates (x', y', z0)
    plane_points_local = np.stack([X_plane, Y_plane, Z_plane_local], axis=-1)
    plane_points_global = plane_points_local @ R.T  # shape (N, N, 3)

    # Extract transformed (x, y, z)
    X_cut = plane_points_global[:, :, 0]
    Y_cut = plane_points_global[:, :, 1]
    Z_cut = plane_points_global[:, :, 2]

    # Plot the cutting plane
    # ax.plot_surface(X_cut, Y_cut, Z_cut + center_height - 1.5, color='gray', alpha=0.4, edgecolor='none')


    ax.set_title(f"Rotated Ellipsoid — phi = {phi:.1f}°")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)

ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, 360, 60), interval=10)
# update(45)
plt.show()