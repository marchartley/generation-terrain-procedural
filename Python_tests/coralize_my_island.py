import matplotlib
import skimage.morphology
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
    corals = extractCoralArea(small_distortion, waterLevel - .1, waterLevel, 0.2, 1.0, subsidence)

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
def method1Create(heights: np.ndarray, subsidence: float) -> np.ndarray:
    waterLevel: float = .7
    minCoralHeight, maxCoralHeight = waterLevel - .1, waterLevel

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
    outsideCorals *= waterLevel + (distanceWhenLowerCoralTouchesWater - distanceFromLowerCorals) / (distanceWhenLowerCoralTouchesWater * 5.)
    finalMap = np.maximum(np.maximum(insideCorals, outsideCorals) * clamp((1-subsidence)**0.5 + .8, 0, 1), heights * subsidence)
    return finalMap


def method1(heights: np.ndarray, subsidenceBetweenFrames: float = 0.05):
    """Find initial seeds and grow a cone from it"""
    currentSubsidence = .85
    finalMap = method1Create(heights, currentSubsidence)
    rgb = np.zeros((finalMap.shape[0], finalMap.shape[1], 3))
    rgb[:, :, 0] = finalMap
    rgb[:, :, 1] = finalMap
    rgb[:, :, 2] = finalMap
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
    rgb = np.zeros((finalMap.shape[0], finalMap.shape[1], 3))
    rgb[:, :, 0] = finalMap
    rgb[:, :, 1] = finalMap
    rgb[:, :, 2] = finalMap
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
    # heights = distanceTransform(readImage("skeleton_1.png"))
    # heights = distanceMap(createRandomWalk())q
    # heights = normalizeImage(readImage("random_heightmaps/Height Map PNG_cropped.png", (200, 200)))
    heights = normalizeImage(readImage("random_heightmaps/gebco_2022_n-9.16_s-9.5953_w45.9874_e46.7857.png", (200, 200)))
    # heights = distortion(heights, 0.05, 0.03)

    subsidenceBetweenFrames = 0.05

    method1(heights, subsidenceBetweenFrames)
    # method2(heights, subsidenceBetweenFrames)
    # method3(heights, subsidenceBetweenFrames)


if __name__ == "__main__":
    main()
