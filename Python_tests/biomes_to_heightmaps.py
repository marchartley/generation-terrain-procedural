import numpy as np
import matplotlib.pyplot as plt
from noise import pnoise2
from scipy.ndimage import gaussian_filter

# Map size
width, height = 256, 256
scale = 200.0

# Perlin noise parameters
octaves = 4
persistence = 0.5
lacunarity = 2.0
seed = np.random.randint(0, 100)

# Generate normalized noise map
def generate_noise_map(width, height, scale, seed_offset=0):
    raw = np.array([
        [pnoise2((x + seed_offset) / scale, (y + seed_offset) / scale,
                 octaves=octaves, persistence=persistence,
                 lacunarity=lacunarity, repeatx=1024, repeaty=1024, base=seed)
         for x in range(width)] for y in range(height)
    ])
    return (raw - raw.min()) / (raw.max() - raw.min())

# Generate environmental maps
elevation = generate_noise_map(width, height, scale)
temperature = generate_noise_map(width, height, scale, seed_offset=100)
moisture = generate_noise_map(width, height, scale, seed_offset=200)

# Classify biome based on thresholds
def classify_biome(e, t, m):
    if e < 0.3:
        return 0  # Ocean
    elif e < 0.35:
        return 1  # Beach
    elif e > 0.8 and t < 0.4:
        return 9  # Snow
    elif e > 0.8:
        return 8  # Mountain
    elif t > 0.8 and m > 0.7:
        return 5  # Tropical Rainforest
    elif t > 0.7 and m < 0.4:
        return 2  # Desert
    elif t < 0.3 and m < 0.6:
        return 7  # Tundra
    elif t < 0.4 and m > 0.4:
        return 6  # Taiga
    elif m > 0.6:
        return 4  # Forest
    else:
        return 3  # Grassland


import numpy as np
from scipy.spatial import KDTree


def generate_biome_map(elevation, temperature, moisture, width, height, num_seeds=100):
    """
    Generates a biome map using Voronoi diagram logic based on seed sampling.

    Args:
        elevation, temperature, moisture: 2D arrays of shape (height, width)
        width, height: dimensions of the map
        num_seeds: number of Voronoi seed points to sample

    Returns:
        biome_map: 2D array of biome IDs, shape (height, width)
    """
    # Step 1: Randomly sample seed points
    # np.random.seed(42)
    seed_coords = np.column_stack((
        np.random.randint(0, width, num_seeds),
        np.random.randint(0, height, num_seeds)
    ))

    # Step 2: Get biome classification for each seed
    seed_biomes = []
    for x, y in seed_coords:
        e = elevation[y, x]
        t = temperature[y, x]
        m = moisture[y, x]
        biome_id = classify_biome(e, t, m)
        seed_biomes.append(biome_id)

    seed_biomes = np.array(seed_biomes)

    # Step 3: Build a KDTree for fast nearest neighbor lookup
    kdtree = KDTree(seed_coords)

    # Step 4: Create biome map by assigning each pixel to the nearest seed
    biome_map = np.zeros((height, width), dtype=np.uint8)

    # Generate grid of coordinates
    grid_x, grid_y = np.meshgrid(np.arange(width), np.arange(height))
    pixel_coords = np.column_stack((grid_x.ravel(), grid_y.ravel()))

    # Query KDTree for nearest seed point for each pixel
    _, nearest_seed_indices = kdtree.query(pixel_coords)

    # Assign biome based on the seed's biome classification
    biome_map = seed_biomes[nearest_seed_indices].reshape((height, width))

    return biome_map

biome_map = generate_biome_map(elevation, temperature, moisture, width, height, num_seeds=50)

# Generate biome map
# biome_map = np.zeros((height, width), dtype=int)
# for y in range(height):
#     for x in range(width):
#         biome_map[y, x] = classify_biome(elevation[y, x], temperature[y, x], moisture[y, x])

# Define biome-specific noise parameters
biome_noise_params = {
    0: {"base": 0.1, "amplitude": 0.1, "frequency": 0.5, "lacunarity": 2.0, "gain": 0.4},  # Ocean
    1: {"base": 0.15, "amplitude": 0.05, "frequency": 1.0, "lacunarity": 2.0, "gain": 0.4},  # Beach
    2: {"base": 0.3, "amplitude": 0.2, "frequency": 1.5, "lacunarity": 2.5, "gain": 0.3},  # Desert
    3: {"base": 0.4, "amplitude": 0.3, "frequency": 1.2, "lacunarity": 2.0, "gain": 0.5},  # Grassland
    4: {"base": 0.5, "amplitude": 0.4, "frequency": 1.0, "lacunarity": 2.1, "gain": 0.5},  # Forest
    5: {"base": 0.55, "amplitude": 0.45, "frequency": 1.1, "lacunarity": 2.0, "gain": 0.6},  # Tropical Rainforest
    6: {"base": 0.45, "amplitude": 0.35, "frequency": 1.3, "lacunarity": 2.0, "gain": 0.5},  # Taiga
    7: {"base": 0.35, "amplitude": 0.2, "frequency": 1.4, "lacunarity": 2.1, "gain": 0.4},  # Tundra
    8: {"base": 0.7, "amplitude": 0.8, "frequency": 2.0, "lacunarity": 2.8, "gain": 0.6},  # Mountain
    9: {"base": 0.6, "amplitude": 0.5, "frequency": 1.5, "lacunarity": 2.5, "gain": 0.4},  # Snow
}


# Smoothstep helper
def smoothstep(edge0, edge1, x):
    x = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return x * x * (3 - 2 * x)

# Compute edge distances for blending
def compute_edge_distance(biome_map):
    edge_distance = np.ones((height, width), dtype=np.float32)
    for y in range(1, height - 1):
        for x in range(1, width - 1):
            center = biome_map[y, x]
            neighbors = biome_map[y-1:y+2, x-1:x+2]
            if np.any(neighbors != center):
                edge_distance[y, x] = 0
    return gaussian_filter(edge_distance, sigma=2)

edge_distances = compute_edge_distance(biome_map)

# Generate biome-aware heightmap
def generate_biome_height_map():
    height_map = np.zeros((height, width), dtype=np.float32)
    for y in range(height):
        for x in range(width):
            biome_id = biome_map[y, x]
            params = biome_noise_params[biome_id]

            blend = smoothstep(0, 1, edge_distances[y, x])
            freq = params["frequency"]
            amp = params["amplitude"]
            lac = params["lacunarity"]
            gain = params["gain"]

            # Fractal noise generation
            nx = x / width
            ny = y / height
            value = 0
            amplitude = amp
            frequency = freq
            for _ in range(4):
                value += pnoise2(nx * frequency, ny * frequency, base=seed) * amplitude
                frequency *= lac
                amplitude *= gain

            height_map[y, x] = value
    # Normalize
    return (height_map - height_map.min()) / (height_map.max() - height_map.min())

# Generate and display the height map
biome_height_map = generate_biome_height_map()

# Helper: Linearly interpolate two parameters
def lerp(a, b, t):
    return a * (1 - t) + b * t

# Get neighboring biome different from current
def find_neighbor_biome(x, y, biome_map, current_biome):
    for dy in range(-1, 2):
        for dx in range(-1, 2):
            nx, ny = x + dx, y + dy
            if 0 <= nx < width and 0 <= ny < height:
                neighbor_biome = biome_map[ny, nx]
                if neighbor_biome != current_biome:
                    return neighbor_biome
    return current_biome  # fallback to self if isolated

# Blend parameters between two biome configs
def blend_params(params_a, params_b, t):
    return {
        "amplitude": lerp(params_a["amplitude"], params_b["amplitude"], t),
        "frequency": lerp(params_a["frequency"], params_b["frequency"], t),
        "lacunarity": lerp(params_a["lacunarity"], params_b["lacunarity"], t),
        "gain": lerp(params_a["gain"], params_b["gain"], t),
        "base": lerp(params_a["base"], params_b["base"], t)
    }

# Updated biome-aware heightmap with parameter interpolation
def generate_blended_biome_height_map():
    height_map = np.zeros((height, width), dtype=np.float32)
    for _ in range(3):
        for y in range(height):
            for x in range(width):
                biome_a = biome_map[y, x]
                biome_b = find_neighbor_biome(x, y, biome_map, biome_a)

                params_a = biome_noise_params[biome_a]
                params_b = biome_noise_params[biome_b]

                blend = smoothstep(0, 1, 1 - edge_distances[y, x])  # near edge = more blend
                params = blend_params(params_a, params_b, blend)

                # Fractal noise generation using blended parameters
                nx = x / width
                ny = y / height
                value = params["base"] if _ > 0 else 0
                amplitude = params["amplitude"]
                frequency = params["frequency"]
                for _ in range(4):
                    value += pnoise2(nx * frequency, ny * frequency, base=seed) * amplitude
                    frequency *= params["lacunarity"]
                    amplitude *= params["gain"]
                # value = blend

                height_map[y, x] = value
        height_map = gaussian_filter(height_map, sigma=2)
    return (height_map - height_map.min()) / (height_map.max() - height_map.min())

# Generate and show the updated height map
blended_biome_height_map = generate_blended_biome_height_map()

# plt.figure(figsize=(8, 8))
fig, axes = plt.subplots(1, 3, squeeze=True)
axes[0].imshow(blended_biome_height_map, cmap='gray')
axes[0].axis("off")

axes[1].imshow(biome_height_map, cmap='gray')
axes[1].axis("off")

axes[2].imshow(biome_map)
axes[2].axis("off")
plt.show()

