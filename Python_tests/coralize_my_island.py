import matplotlib
import skimage.morphology
from jupyterlab.semver import outside
from scipy.ndimage import convolve
from scipy.signal import convolve2d

from heightmap_from_skeleton import *
import numpy as np
import matplotlib.pyplot as plt
import noise.perlin
from skimage.morphology import skeletonize
from Vectors import *

def plot3D(grid: np.ndarray, ax = None, cmap='viridis'):
    grid = grid.copy()
    callShow = ax is None
    x, y = np.meshgrid(range(grid.shape[1]), range(grid.shape[0]))
    grid[0, 0] = 0  # just to block the normalization process
    grid *= min(grid.shape)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(grid)))
    ax.plot_surface(x, y, grid, rstride=4, cstride=4, linewidth=0, cmap=cmap, vmin=0.0, vmax=min(grid.shape), antialiased=False, shade=True)
    if callShow:
        plt.show()

def bw2rgb(bw: np.ndarray) -> np.ndarray:
    rgb = np.zeros((bw.shape[0], bw.shape[1], 3))
    rgb[:, :, 0] = bw
    rgb[:, :, 1] = bw
    rgb[:, :, 2] = bw
    return np.clip(rgb, 0.0, 1.0)

def normalizeImage(img: np.ndarray) -> np.ndarray:
    mini, maxi = scipy.ndimage.minimum(img), scipy.ndimage.maximum(img)
    return (img - mini) / (maxi - mini)

def distortion(img: np.ndarray, frequency: float = 0.01, noiseImpact: float = 0.1) -> np.ndarray:
    # scipy.ndimage.dis
    distort: np.ndarray = img.copy()
    p = noise.perlin.SimplexNoise(200)
    for y in range(distort.shape[0]):
        for x in range(distort.shape[1]):
            distort[y, x] = distort[y, x] + p.noise2(x * frequency, y * frequency) * noiseImpact
    return distort

def binarizeImage(img: np.ndarray, minVal: float, maxVal: float) -> np.ndarray:
    binary = np.zeros_like(img)
    binary[np.logical_and(minVal < img, img < maxVal)] = 1
    return binary

def irange(start, stop, step = 1):
    while start < stop:
        yield start
        start += step
    return

def clamp(x, _min, _max):
    return min(max(x, _min), _max)

def removeParasitesOnSkeleton(skeleton: np.ndarray, minimalSize: int = 10):
    for y in range(skeleton.shape[0]):
        for x in range(skeleton.shape[1]):
            if not get(skeleton, x, y):
                continue
            nb_neighbors = 0
            for dy in irange(-1, 2):
                for dx in irange(-1, 2):
                    nb_neighbors += get(skeleton, x + dx, y + dy)
            if nb_neighbors > 3:
                for dy in irange(-1, 2):
                    for dx in irange(-1, 2):
                        set(skeleton, x + dx, y + dy, 0)

    return skimage.morphology.remove_small_objects(skeleton, min_size = minimalSize, connectivity=2)

def extractInsideAndOutsideOfReef(heights: np.ndarray, minHeight: float, maxHeight: float, minGradient: float, maxGradient: float, subsidence: float) -> Tuple[np.ndarray, np.ndarray]:
    inside = np.zeros_like(heights)
    outside = np.zeros_like(heights)

    # define inside as the area above minHeight
    inside[heights > minHeight] = 1
    outside[heights < minHeight] = 1

    return inside, outside


def extractCoralArea(heights: np.ndarray, minHeight: float, maxHeight: float, minGradient: float, maxGradient: float, subsidence: float) -> np.ndarray:
    coralArea = np.zeros_like(heights)
    gradient = normalizeImage(scipy.ndimage.morphological_gradient(heights, size=(10, 10)))
    coralArea[np.logical_and(heights > minHeight, np.logical_and(heights < maxHeight, np.logical_and(minGradient < gradient, gradient < maxGradient)))] = 1
    # coralArea[np.logical_and(heights > minHeight, heights < maxHeight)] = 1
    skeleton = removeParasitesOnSkeleton(skeletonize(coralArea), 10)
    finalCoralArea = np.zeros_like(coralArea)

    for y in range(heights.shape[0]):
        for x in range(heights.shape[1]):
            gradientValue = get(gradient, x, y)
            if not get(skeleton, x, y) or gradientValue == 0:
                continue
            sizeFactor = 2  # clamp(1 / gradientValue, 0, 100)
            for dy in irange(-sizeFactor, sizeFactor):
                for dx in irange(-sizeFactor, sizeFactor):
                    set(finalCoralArea, x + dx, y + dy, 1.0)
    return finalCoralArea

def extractIslandArea(heights: np.ndarray, waterLevel: float) -> np.ndarray:
    islandArea = np.zeros_like(heights)
    islandArea[heights > waterLevel] = 1
    return islandArea

def lerp(x, mini, maxi):
    return mini + x * (maxi - mini)

def createRandomIslandFromHeightmap(heightmap: np.ndarray, subsidence: float = 0.9) -> np.ndarray:
    waterLevel: float = .7
    small_distortion = distortion(heightmap, 0.05, 0.0)
    big_distortion = distortion(heightmap, 0.05, 0.0)
    corals = extractCoralArea(small_distortion, waterLevel - .1, waterLevel, 0.0, 1.0, subsidence)

    big_distortion *= subsidence
    island = extractIslandArea(big_distortion, waterLevel)

    distanceToCoral = distanceTransform(1 - corals)
    distanceToIsland = distanceTransform(1 - island)
    difference = distanceToIsland - distanceToCoral

    output = np.zeros_like(heightmap)
    # output = big_distortion * .5
    output[:] = waterLevel - 0.01
    # output = (lerp(1 - np.clip(distanceToCoral, 0, 1), waterLevel, waterLevel - 0.5) + lerp(1 - np.clip(distanceToIsland, 0, 1), waterLevel, waterLevel - 0.5)) * .5
    # output[difference > 0] = waterLevel - 0.5
    # output[difference <= 0] = lerp(1 - np.clip(distanceToCoral[difference <= 0], 0, 1), waterLevel, waterLevel - 0.5)
    output[corals == 1] = max(np.max(big_distortion), waterLevel)
    output[island == 1] = big_distortion[island == 1]
    return output


"""
Might be useful at one point :
    higherBand = np.zeros_like(heights)
    surfaceIsland = np.zeros_like(heights)

    surfaceIsland[heights > waterLevel] = 1
    higherBand[np.logical_and(initialCoral != dilatedInitialSeed, heights >= minCoralHeight)] = 1
    islandSkeleton = skeletonize(surfaceIsland)
    distanceToSkeleton = distanceMap(1 - islandSkeleton)
    gradientFromIsland = np.gradient(distanceToSkeleton)  # grad[0] -> towards top, grad[1] -> towards left
    gradient = (gradientFromIsland[0]**2 + gradientFromIsland[1]**2)**.5
    distanceFromIsland = distanceMap(1 - surfaceIsland)
"""

def smin(a: np.ndarray, b: np.ndarray, k: float) -> np.ndarray:
    k *= 6.0
    h = np.maximum( k-np.abs(a-b), 0.0 )/k
    return np.minimum(a,b) - h*h*h*k*(1.0/6.0)

def smax(a: np.ndarray, b: np.ndarray, k: float) -> np.ndarray:
    return -smin(-a, -b, k)

def method1Create(heights: np.ndarray, subsidence: float, waterLevel: float = 0.5, maxCoralHeight:float = None, minCoralHeight: float = None, outsideSlopeFactor: float = 3.0) -> np.ndarray:
    # waterLevel: float = .7
    if maxCoralHeight is None:
        maxCoralHeight = waterLevel
    if minCoralHeight is None:
        minCoralHeight = maxCoralHeight - 0.1
    # minCoralHeight, maxCoralHeight = waterLevel - .1, waterLevel

    verticalScale = 1.0 / 1000.0
    horizontalScale = 1.0 / 1.0
    vh = 1.0 / horizontalScale
    vv = 1.0 / verticalScale
    alpha = 0.1

    initialCoral = np.zeros_like(heights)
    lowerBand = np.zeros_like(heights)

    initialCoral[np.logical_and(minCoralHeight <= heights, heights <= maxCoralHeight)] = 1
    dilatedInitialSeed = dilation(initialCoral)
    lowerBand[np.logical_and(initialCoral != dilatedInitialSeed, heights <= minCoralHeight)] = 1
    distanceFromLowerCorals = 1 - distanceMap(1 - lowerBand)
    distanceFromLowerCorals[heights > minCoralHeight] *= -1

    distanceWhenLowerCoralTouchesWater = (vh / vv) * (maxCoralHeight - minCoralHeight) / verticalScale
    insideCorals = np.zeros_like(heights)
    insideCorals[np.logical_or(heights > waterLevel, np.logical_and(distanceFromLowerCorals < distanceWhenLowerCoralTouchesWater, heights < maxCoralHeight))] = 1
    insideCorals *= 1 - ((distanceWhenLowerCoralTouchesWater - distanceFromLowerCorals) / distanceWhenLowerCoralTouchesWater) * alpha
    insideCorals *= waterLevel

    outsideCorals = np.zeros_like(heights)
    outsideCorals[distanceFromLowerCorals >= distanceWhenLowerCoralTouchesWater] = 1
    outsideCorals *= waterLevel + (distanceWhenLowerCoralTouchesWater - distanceFromLowerCorals) * (outsideSlopeFactor / 3.0) / distanceWhenLowerCoralTouchesWater

    insideOutsideCorals = smax(insideCorals, outsideCorals, 0.01) * clamp((1-subsidence)**0.5 + .8, 0, 1)
    finalMap = smax(insideOutsideCorals, heights * subsidence, 0.05)
    return finalMap


def method1(heights: np.ndarray, subsidenceBetweenFrames: float = 0.05):
    """Find initial seeds and grow a cone from it"""
    currentSubsidence = .85
    finalMap = method1Create(heights, currentSubsidence)
    rgb = bw2rgb(finalMap)
    # plt.imsave("map.png", rgb)

    results = [2, 2]
    fig, axes = plt.subplots(results[0], results[1], squeeze=False)
    plotNumber = 0
    for i in range(results[0]):
        for j in range(results[1]):
            currentSubsidence = 1 - ((plotNumber + 1) * subsidenceBetweenFrames)
            finalMap = method1Create(heights, currentSubsidence)
            axes[i, j].imshow(finalMap, vmin=0, vmax=1)
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            plotNumber += 1
    plt.show()
    plot3D(finalMap)

def method2(heights: np.ndarray, subsidenceBetweenFrames: float = 0.05):
    """Grow a wide band from the initial seed towards the sea"""
    pass

def method3(heights: np.ndarray, subsidenceBetweenFrames: float = 0.05):
    """Find initial seeds and build vertical wall from them"""
    results = [2, 2]
    fig, axes = plt.subplots(results[0], results[1], squeeze=False)
    plotNumber = 0
    currentSubsidence = 0.85
    finalMap = createRandomIslandFromHeightmap(heights, subsidence=currentSubsidence)
    rgb = bw2rgb(finalMap)
    plt.imsave("map_type3.png", rgb)
    for i in range(results[0]):
        for j in range(results[1]):
            currentSubsidence = 1 - ((plotNumber + 1) * subsidenceBetweenFrames)
            finalMap = createRandomIslandFromHeightmap(heights, subsidence= currentSubsidence)
            axes[i, j].imshow(finalMap, vmin=0, vmax=1)
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            plotNumber += 1
    plt.show()
    plot3D(createRandomIslandFromHeightmap(heights, subsidence = currentSubsidence))

def main():
    random.seed(1)
    heights = readImage("../saved_maps/heightmaps/Mt_Ruapehu_Mt_Ngauruhoe.png")
    # heights = readImage("../saved_maps/heightmaps/new_one_slope.png")
    # heights = distanceTransform(readImage("skeleton_1.png"))
    # heights = distanceMap(createRandomWalk())
    # plt.imshow(heights, cmap="gray")
    # plt.show()
    # heights = normalizeImage(readImage("random_heightmaps/Height Map PNG_cropped.png", (200, 200)))
    # heights = normalizeImage(readImage("random_heightmaps/gebco_2022_n-9.16_s-9.5953_w45.9874_e46.7857.png", (200, 200)))
    # heights = distortion(heights, 0.05, 0.03)

    subsidenceBetweenFrames = 0.05

    method1(heights, subsidenceBetweenFrames)
    # method2(heights, subsidenceBetweenFrames)
    # method3(heights, subsidenceBetweenFrames)


def apply_thermal_erosion(heightmap, resistanceMap, iterations=10, talus_angle=0.5, erosion_factor=0.1):
    """
    Apply thermal erosion to a heightmap.

    Parameters:
        heightmap (numpy.ndarray): 2D array representing the heightmap.
        iterations (int): Number of erosion iterations.
        talus_angle (float): Maximum allowable slope angle before erosion occurs.
        erosion_factor (float): Proportion of height difference to erode.

    Returns:
        numpy.ndarray: Eroded heightmap.
    """
    # Ensure the heightmap is a NumPy array
    heightmap = np.array(heightmap, dtype=float) * 100.0
    rows, cols = heightmap.shape

    for _ in range(iterations):
        # Pad heightmap for neighbor calculations
        padded_heightmap = np.pad(heightmap, 1, mode='edge')
        neighbors = [
            padded_heightmap[:-2, 1:-1],  # Up
            padded_heightmap[2:, 1:-1],  # Down
            padded_heightmap[1:-1, :-2],  # Left
            padded_heightmap[1:-1, 2:]  # Right
        ]

        # Calculate height differences and transfer amounts
        delta_height = np.zeros_like(heightmap)
        for neighbor in neighbors:
            height_diff = heightmap - neighbor
            transfer = np.where(height_diff > talus_angle, erosion_factor * (height_diff - talus_angle) * (1 - resistanceMap), 0)
            delta_height -= transfer
            neighbor += transfer

        # Apply calculated changes
        heightmap += delta_height

    return heightmap / 100.0


# def apply_hydraulic_erosion(heightmap, iterations=100, water=1.0, solubility=0.1, evaporation=0.1, capacity=1.0):
#     """
#     Apply simple hydraulic erosion to a heightmap.
#
#     Parameters:
#         heightmap (numpy.ndarray): 2D array representing the heightmap.
#         iterations (int): Number of erosion iterations.
#         water (float): Initial amount of water at each cell.
#         solubility (float): Proportion of material dissolved into the water.
#         evaporation (float): Proportion of water that evaporates per iteration.
#         capacity (float): Maximum sediment a unit of water can carry.
#
#     Returns:
#         numpy.ndarray: Eroded heightmap.
#     """
#     # Ensure the heightmap is a NumPy array
#     heightmap = np.array(heightmap, dtype=float) * 100.0
#     rows, cols = heightmap.shape
#
#     # Initialize water and sediment arrays
#     water_map = np.full_like(heightmap, water, dtype=float)
#     sediment_map = np.zeros_like(heightmap)
#
#     for _ in range(iterations):
#         # Create temporary arrays for changes
#         delta_height = np.zeros_like(heightmap)
#         delta_sediment = np.zeros_like(heightmap)
#
#         for i in range(1, rows - 1):
#             for j in range(1, cols - 1):
#                 # Get neighbors
#                 neighbors = [
#                     (i - 1, j),  # Up
#                     (i + 1, j),  # Down
#                     (i, j - 1),  # Left
#                     (i, j + 1),  # Right
#                 ]
#
#                 # Distribute water flow and sediment
#                 total_outflow = 0
#                 flows = []
#                 for ni, nj in neighbors:
#                     height_diff = (heightmap[i, j] + water_map[i, j]) - (heightmap[ni, nj] + water_map[ni, nj])
#                     flow = max(0, height_diff) / 4  # Divide by 4 for even distribution
#                     flows.append((ni, nj, flow))
#                     total_outflow += flow
#
#                 if total_outflow > 0:
#                     for ni, nj, flow in flows:
#                         # Proportion of outflow to each neighbor
#                         proportion = flow / total_outflow
#
#                         # Transfer water
#                         delta_height[ni, nj] += proportion * water_map[i, j]
#                         delta_height[i, j] -= proportion * water_map[i, j]
#
#                         # Transfer sediment
#                         delta_sediment[ni, nj] += proportion * sediment_map[i, j]
#                         delta_sediment[i, j] -= proportion * sediment_map[i, j]
#
#                 # Dissolve material
#                 dissolved = solubility * water_map[i, j]
#                 sediment_map[i, j] += dissolved
#                 heightmap[i, j] -= dissolved
#
#         # Update water and sediment maps
#         water_map += delta_height
#         sediment_map += delta_sediment
#
#         # Evaporate water and deposit sediment
#         water_map *= (1 - evaporation)
#         excess_sediment = sediment_map - (water_map * capacity)
#         sediment_map = np.maximum(sediment_map - excess_sediment, 0)
#         heightmap += np.maximum(excess_sediment, 0)
#
#     return heightmap / 100.0

def apply_hydraulic_erosion(heightmap, resistanceMap, iterations=100, water=1.0, solubility=0.8, evaporation=0.1, capacity=1.0):
    """
    Optimized hydraulic erosion function with fixed numerical stability.

    Parameters:
        heightmap (numpy.ndarray): 2D array representing the heightmap.
        iterations (int): Number of erosion iterations.
        water (float): Initial amount of water at each cell.
        solubility (float): Proportion of material dissolved into the water.
        evaporation (float): Proportion of water that evaporates per iteration.
        capacity (float): Maximum sediment a unit of water can carry.

    Returns:
        numpy.ndarray: Eroded heightmap.
    """
    heightmap = np.array(heightmap, dtype=float) * 100.0
    rows, cols = heightmap.shape

    water_map = np.full_like(heightmap, water, dtype=float)
    sediment_map = np.zeros_like(heightmap)

    for _ in range(iterations):
        # Calculate slopes and flow direction
        height_plus_water = heightmap + water_map
        padded = np.pad(height_plus_water, 1, mode='edge')
        neighbors = [
            padded[:-2, 1:-1],  # Up
            padded[2:, 1:-1],   # Down
            padded[1:-1, :-2],  # Left
            padded[1:-1, 2:]    # Right
        ]
        diffs = [height_plus_water - neighbor for neighbor in neighbors]
        flows = [np.maximum(0, diff) for diff in diffs]

        total_outflow = sum(flows)
        total_outflow = np.maximum(total_outflow, 1e-8)  # Prevent division by zero

        # Distribute water and sediment based on flow
        for flow, neighbor in zip(flows, neighbors):
            proportion = flow / total_outflow
            transfer = proportion * water_map
            water_map -= transfer

            # Handle sediment safely with water_map > 0
            valid_mask = water_map > 1e-8
            sediment_map[valid_mask] -= transfer[valid_mask] * (sediment_map[valid_mask] / water_map[valid_mask])
            water_map += np.pad(transfer, 1, mode='constant')[1:-1, 1:-1]

        # Dissolve material
        dissolved = solubility * water_map * (1.0 - resistanceMap)
        sediment_map += dissolved
        heightmap -= dissolved

        # Evaporation and sediment deposition
        water_map *= (1 - evaporation)
        excess_sediment = sediment_map - (water_map * capacity)
        sediment_map = np.maximum(sediment_map - excess_sediment, 0)
        heightmap += np.maximum(excess_sediment, 0)

    return heightmap / 100.0

if __name__ == "__main__":
    main()
