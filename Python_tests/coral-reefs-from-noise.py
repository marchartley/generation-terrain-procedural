import math
import random
from math import floor
from typing import Tuple

import noise.perlin
import numpy as np
import matplotlib.pyplot as plt

from Python_tests.mergedIslandGeneration import valueAsHSV


def sign(x):
    return 1 if x > 0 else -1

def clamp(x:float, mini: float, maxi: float) -> float:
    return mini if x <= mini else maxi if x >= maxi else x

def valueAsHSV(value: float, mini: float, maxi: float, L: float = 0.5, S: float = 1.0) -> Tuple[float, float, float]:
    H = (360.0 * (math.fmod(value, maxi) - mini) / (maxi - mini)) / 60
    # L = 0.5
    # S = 1.0

    C = (1 - abs(2 * L - 1)) * S
    try:
        X = C * (1 - abs(math.fmod(H, 2) - 1))
    except Exception as e:
        a = 0
        raise e

    R1, G1, B1 = (C, X, 0) if 0 <= H < 1 else (X, C, 0) if 1 <= H < 2 else (0, C, X) if 2 <= H < 3 else (0, X, C) if 3 <= H < 4 else (X, 0, C) if 4 <= H < 5 else (C, 0, X)
    m = L - C * .5
    R, G, B = R1 + m, G1 + m, B1 + m
    return clamp(R,0,1), clamp(G,0,1), clamp(B,0,1)

def falloff(x, y, max_x, max_y, width):
    return math.sin(clamp(min([x, y, max_x - x, max_y - y]) / width, 0, 1) * math.pi * 0.5)

arr = np.zeros((256, 256, 3))

n = noise.perlin.SimplexNoise(period=200)

X = 6
val_array = np.zeros((256, 256))
for y in range(arr.shape[0]):
    for x in range(arr.shape[1]):
        noise_val = n.noise2(x / 150, y / 150) * 0.8 + n.noise2(x / 100, y / 100) * 0.3 + n.noise2(x / 50, y / 50) * 0.2
        val = clamp(sign(noise_val) * abs(noise_val)**1.5 + 0.6, -1, 0.99)
        # val = clamp(n.noise2(x / 100, y / 100) + n.noise2((x + 20) / 50, (y + 20) / 50) * .33 + 0.5, -1, 0.99)
        # val = (X) - min(max(sign(val) * abs(val ** 2.0) * X, 0), X - 1) #max(floor(val * 5), 0)
        val = val * falloff(x, y, arr.shape[1], arr.shape[0], 100)
        if val < 0.4: val = 6   # Abyss
        elif val < 0.5: val = 4 # Reef
        elif val < 0.7: val = 3 # Lagoon
        elif val < 0.8: val = 2 # Beach
        else: val = 1           # Island
        # Abyss  = 0 or 6
        # Island = 1
        # Beach  = 2
        # Lagoon = 3
        # Reef   = 4
        # val = (X) - floor(max(val * X * falloff(x, y, arr.shape[1], arr.shape[0], 100), 0))
        # if val == X - 1:
        #     val = X
        val_array[y, x] = val

nbLevels = 30
for iteration in range(nbLevels):
    L = 0.5 * (1 - (iteration / (nbLevels - 1)))
    for y in range(arr.shape[0]):
        for x in range(arr.shape[1]):
            arr[x, y] = valueAsHSV(val_array[x, y], 0, X, L=L) # np.array(valueAsHSV(val, 0, (X), L=0.05))

    if iteration == 0:
        plt.imshow(arr)
        plt.show()
    plt.imsave(f"/media/marc/Data/NN Datasets/1/noise_labels/input_label-{iteration}.png", arr)
