import numpy as np
import matplotlib.pyplot as plt

def hat_tile_vertices(base_length, height):
    """
    Function to compute the positions of vertices of a hat tile.

    Args:
    - base_length: Length of the base of the hat tile.
    - height: Height of the hat tile (distance from the top vertex to the base).

    Returns:
    - Array of vertex positions (2D coordinates).
    """
    # Top vertex (origin)
    top_vertex = np.array([0.0, 0.0])

    # Bottom vertices of the triangle
    bottom_left_vertex = np.array([-base_length / 2, 0.0])
    bottom_right_vertex = np.array([base_length / 2, 0.0])

    # Apex vertex (top of the hat)
    apex_vertex = np.array([0.0, height])

    return np.array([top_vertex, bottom_left_vertex, bottom_right_vertex, apex_vertex])

# Define parameters for the hat tile
base_length = 4.0
height = 3.0

# Compute vertex positions
vertices = hat_tile_vertices(base_length, height)

# Plot the hat tile
plt.figure(figsize=(6, 4))
plt.plot(vertices[:, 0], vertices[:, 1], 'b-')
plt.plot(vertices[:, 0], vertices[:, 1], 'ro')
plt.title('Vertices of Hat Tile')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
plt.grid(True)
plt.show()
