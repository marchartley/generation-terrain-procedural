# Complete comparison script with numeric checks included
import numpy as np
import random
import matplotlib.pyplot as plt

def bridson_sampling(width, height, radius, k=30):
    cell_size = radius / np.sqrt(2)
    grid_width, grid_height = int(width / cell_size) + 1, int(height / cell_size) + 1
    grid = [[None for _ in range(grid_height)] for _ in range(grid_width)]
    samples, active_list = [], []
    checks = 0

    def get_grid_coords(p):
        return int(p[0] / cell_size), int(p[1] / cell_size)

    def is_valid(pt):
        nonlocal checks
        gx, gy = get_grid_coords(pt)
        for i in range(max(0, gx - 2), min(grid_width, gx + 3)):
            for j in range(max(0, gy - 2), min(grid_height, gy + 3)):
                neighbor = grid[i][j]
                if neighbor:
                    checks += 1
                    if np.linalg.norm(np.array(pt) - np.array(neighbor)) < radius:
                        return False
        return True

    p0 = (random.uniform(0, width), random.uniform(0, height))
    samples.append(p0)
    gx, gy = get_grid_coords(p0)
    grid[gx][gy] = p0
    active_list.append(p0)

    while active_list:
        idx = random.randint(0, len(active_list)-1)
        center = active_list[idx]
        placed = False

        for _ in range(k):
            angle = random.uniform(0, 2 * np.pi)
            distance = random.uniform(radius, 2 * radius)
            pt = (center[0] + distance * np.cos(angle),
                  center[1] + distance * np.sin(angle))

            if not (0 <= pt[0] < width and 0 <= pt[1] < height):
                continue
            if is_valid(pt):
                samples.append(pt)
                active_list.append(pt)
                gx, gy = get_grid_coords(pt)
                grid[gx][gy] = pt
                placed = True

        if not placed:
            active_list.pop(idx)

    return samples, checks

def robust_variable_radii_sampling(width, height, radii, k=30):
    cell_size = max(radii) / np.sqrt(2)
    grid_width, grid_height = int(width / cell_size) + 1, int(height / cell_size) + 1
    grid = [[[] for _ in range(grid_height)] for _ in range(grid_width)]
    samples, active_list = [], []
    checks = 0

    def get_grid_coords(p):
        return int(p[0] / cell_size), int(p[1] / cell_size)

    def is_valid(pt, r):
        nonlocal checks
        gx, gy = get_grid_coords(pt)
        search_range = int(np.ceil((max(radii)*2) / cell_size))
        for i in range(max(0, gx - search_range), min(grid_width, gx + search_range + 1)):
            for j in range(max(0, gy - search_range), min(grid_height, gy + search_range + 1)):
                for neighbor, nr in grid[i][j]:
                    checks += 1
                    if np.linalg.norm(np.array(pt) - np.array(neighbor)) < (r + nr):
                        return False
        return True

    r0 = random.choice(radii)
    p0 = (random.uniform(0, width), random.uniform(0, height))
    samples.append((p0, r0))
    gx, gy = get_grid_coords(p0)
    grid[gx][gy].append((p0, r0))
    active_list.append((p0, r0))

    while active_list:
        idx = random.randint(0, len(active_list)-1)
        center, center_radius = active_list[idx]
        placed = False

        for _ in range(k):
            new_radius = random.choice(radii)
            angle = random.uniform(0, 2 * np.pi)
            min_dist = center_radius + new_radius
            distance = random.uniform(min_dist, 2 * min_dist)
            pt = (center[0] + distance * np.cos(angle),
                  center[1] + distance * np.sin(angle))

            if not (0 <= pt[0] < width and 0 <= pt[1] < height):
                continue
            if is_valid(pt, new_radius):
                samples.append((pt, new_radius))
                active_list.append((pt, new_radius))
                gx, gy = get_grid_coords(pt)
                grid[gx][gy].append((pt, new_radius))
                placed = True

        if not placed:
            active_list.pop(idx)

    return samples, checks

# Run both algorithms
radius = 2
width, height = 100, 100

bridson_samples, bridson_checks = bridson_sampling(width, height, 2 * radius, k=30)
variable_samples, variable_checks = robust_variable_radii_sampling(width, height, [radius, 2*radius, 3*radius], k=30)

# Results summary
print(f"Bridson's Algorithm: {len(bridson_samples)} samples, {bridson_checks} checks")
print(f"Variable Radii Algorithm: {len(variable_samples)} samples, {variable_checks} checks")

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Plot Bridson's Algorithm
ax1 = axes[0]
for x, y in bridson_samples:
    circle = plt.Circle((x, y), radius, edgecolor='blue', facecolor='none', alpha=0.5)
    ax1.add_patch(circle)
ax1.set_xlim([0, width])
ax1.set_ylim([0, height])
ax1.set_aspect('equal')
ax1.set_title(f"Bridson's Algorithm ({len(bridson_samples)} samples)")

# Plot Variable Radii Algorithm
ax2 = axes[1]
for (x, y), r in variable_samples:
    circle = plt.Circle((x, y), r, edgecolor='green', facecolor='none', alpha=0.5)
    ax2.add_patch(circle)
ax2.set_xlim([0, width])
ax2.set_ylim([0, height])
ax2.set_aspect('equal')
ax2.set_title(f'Variable Radii Algorithm ({len(variable_samples)} samples)')

plt.suptitle('Comparison of Bridson vs Variable Radii Algorithm')
plt.show()
