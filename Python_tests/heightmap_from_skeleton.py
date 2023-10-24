import random
from typing import Tuple

import numpy as np
import PIL.Image
import matplotlib.pyplot as plt
import scipy.ndimage
import os

def get(img, x, y, default = 0):
    return default if x < 0 or x >= img.shape[1] or y < 0 or y >= img.shape[0] else img[y, x]
def set(img, x, y, val):
    if x < 0 or x >= img.shape[1] or y < 0 or y >= img.shape[0]:
        return
    img[int(y), int(x)] = val

def erosion(img: np.ndarray, useCross: bool = True) -> np.ndarray:
    eroded = scipy.ndimage.grey_erosion(img, size=(3, 3))
    return eroded

def dilation(img: np.ndarray, useCross: bool = True) -> np.ndarray:
    dilated = scipy.ndimage.grey_dilation(img, size=(3, 3))
    return dilated

def createRandomWalk() -> np.ndarray:
    x = 200
    y = 200
    minPos = .2
    maxPos = 1 - minPos
    img_file = np.ones((y, x))
    pos = [y * .5, x * .5]
    for i in range(40):
        pos[0] = (pos[0] + random.randint(-2, 4))
        pos[1] = (pos[1] + random.randint(-2, 4))
        if i % 5 == 0:
            pos[0] = (pos[0] + random.randint(-40, 40))
            pos[1] = (pos[1] + random.randint(-40, 40))
        pos[0] = y * minPos if pos[0] > y * maxPos else y * maxPos if pos[0] < y * minPos else pos[0]
        pos[1] = x * minPos if pos[1] > x * maxPos else x * maxPos if pos[1] < x * minPos else pos[1]
        set(img_file, int(pos[0]), int(pos[1]), 0)
    return img_file

def readImage(filename: str, size: Tuple[int, int] = (200, 100)) -> np.ndarray:
    img_file = PIL.Image.open(filename).resize(size, resample=0)
    img = np.array(img_file) / 255.0
    if len(img.shape) == 3:
        img = img[:, :, 0]
    return img

def distanceTransform(img: np.ndarray) -> np.ndarray:
    img = img.copy()
    # img[:, 0] = .999
    # img[:, -1] = .999
    # img[0, :] = .999
    # img[-1, :] = .999
    forcedPoints = img.copy()
    forcedPoints[forcedPoints < 1.0] = 0

    """forcedCoord = np.where(img != 1)

    def dist(x0, y0, x1, y1):
        # return abs(x1 - x0) + abs(y1 - y0)
        return ((x1 - x0)**2 + (y1 - y0)**2)**.5

    distances: np.ndarray = np.zeros_like(img)
    for x in range(img.shape[1]):
        for y in range(img.shape[0]):
            heightAndWeight = []
            heights = []
            totalDistances = 0
            for _y, _x, val in zip(forcedCoord[0], forcedCoord[1], img[forcedCoord[0], forcedCoord[1]]):
                distanceToForced = dist(x, y, _x, _y)
                forcedValue = val # get(img, _x, _y)
                heightAndWeight.append((forcedValue, distanceToForced))
                heights.append(forcedValue / max(1, distanceToForced))
                totalDistances += distanceToForced
            set(distances, x, y, sum(heights)/len(heights))
            # set(distances, x, y, sum([height * (1 - (distance / totalDistances)**2) for height, distance in heightAndWeight]))

    distances /= max(distances.flatten())

    plt.imshow(distances)
    plt.show()
    return"""

    indices: np.ndarray
    distances, indices = scipy.ndimage.distance_transform_edt(forcedPoints, return_indices=True)
    rgb = np.zeros((img.shape[0], img.shape[1], 3))
    rgb[:, :, 0] = indices[0]
    rgb[:, :, 1] = indices[1]
    rgb[:, :, 2] = max(rgb.flatten())

    rgb /= max(rgb.flatten())

    distances /= max(distances.flatten())
    distances = (1 - distances) * (1 - img[indices[0], indices[1]])
    return distances


def distanceMap(img: np.ndarray) -> np.ndarray:
    dist = scipy.ndimage.distance_transform_edt(img)
    dist = 1 - (dist / scipy.ndimage.maximum(dist))
    return dist

def saveImage(img: np.ndarray, filename: str):
    rgb = np.zeros((img.shape[0], img.shape[1], 3))
    rgb[:, :, 0] = img
    rgb[:, :, 1] = img
    rgb[:, :, 2] = img
    plt.imsave(filename, rgb)

def main():
    skeleton = createRandomWalk()
    distance = distanceMap(skeleton)
    saveImage(distance, "random_result.png")
    print("Done")


if __name__ == "__main__":
    main()
