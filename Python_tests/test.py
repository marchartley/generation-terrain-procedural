import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
grid_size = 50
num_steps = 100
G = 9.81e10 # 6.67430e-11  # gravitational constant

# Initialize grid
grid = np.zeros((grid_size, grid_size))  # grid representing mass distribution
velocity = np.zeros((grid_size, grid_size, 2))  # velocity of each cell (x, y)
acceleration = np.zeros((grid_size, grid_size, 2))  # acceleration of each cell (x, y)

# Add a massive object at the center
massive_object_mass = 1e10
grid[grid_size // 2, grid_size // 2] = massive_object_mass

# Function to calculate gravitational force
def calculate_gravity(m1, m2, r):
    theta = 0
    force_magnitude = G * m1 * m2 / (r ** 2)
    force_x = force_magnitude * np.cos(theta)
    force_y = force_magnitude * np.sin(theta)
    return force_x, force_y

# Function to update the grid based on gravitational forces
def update_grid():
    for i in range(grid_size):
        for j in range(grid_size):
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    ni, nj = i + dx, j + dy
                    if 0 <= ni < grid_size and 0 <= nj < grid_size:
                        # Calculate distance and angle between cells
                        r = np.sqrt(dx ** 2 + dy ** 2)
                        theta = np.arctan2(dy, dx)
                        # Calculate gravitational force
                        force_x, force_y = calculate_gravity(grid[i, j], grid[ni, nj], r)
                        # Update acceleration
                        acceleration[i, j, 0] += force_x #/ grid[i, j]
                        acceleration[i, j, 1] += force_y #/ grid[i, j]

# Function to update cell velocities and positions
def update_simulation():
    global velocity, acceleration, grid
    update_grid()
    # Update velocity
    velocity += acceleration
    # Update position
    grid += (velocity[:, :, 0]**2 + velocity[:, :, 1]**2)**0.5
    # Reset acceleration
    acceleration = np.zeros((grid_size, grid_size, 2))

# Function to animate simulation
def animate(frame):
    print("Goooo")
    update_simulation()
    im.set_array(grid)
    return im,

# Create animation
fig, ax = plt.subplots()
im = ax.imshow(grid, cmap='viridis', animated=True)
ani = animation.FuncAnimation(fig, animate, frames=num_steps, repeat=False, blit=True)
plt.show()
