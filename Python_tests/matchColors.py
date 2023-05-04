import glob
import os.path
import shutil
from scipy.ndimage import median_filter

import PIL
import PIL.Image
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    # file = "correct_synthetic_terrains_dataset_even_other_colors"
    prefix = "___"
    goodFolder = "correct_synthetic_terrains_dataset/features/"
    badFolder = "correct_synthetic_terrains_dataset_even_other_colors/features/"
    # hueValues = {
    #     76: 153,  # abyss
    #     51: 101,  # coral
    #     25: 51  # island
    # }
    hueValues = {
        0: 153,  # abyss
        223: 153,  # abyss (2)
        190: 153,
        220: 153,
        159: 101,  # coral
        140: 101,
        95: 51  # island
    }

    files = glob.glob(badFolder + "*.png")
    for filename in files:
        img = PIL.Image.open(filename)
        img = img.convert('HSV')
        arr = np.array(img)
        arr[:, :, 0] = median_filter(arr[:, :, 0], size=3)
        # plt.imshow(arr)
        # plt.show()
        # exit(0)
        converted = arr.copy()
        for toReplace, replaced in hueValues.items():
            converted[:, :, 0][abs(arr[:, :, 0].astype(int) - int(toReplace)) < 30] = replaced
        # print(goodFolder + prefix + os.path.basename(filename), filename)
        # converted = median_filter(median_filter(median_filter(converted, size=(3, 3, 1)), size=(3, 3, 1)), size=(3, 3, 1)) #, cval=0, mode='constant')
        res = np.array(PIL.Image.fromarray(converted, mode="HSV").convert('RGB'))
        # plt.imshow(res)
        # plt.show()
        # exit(0)
        plt.imsave(goodFolder + prefix + os.path.basename(filename), res)
        # exit(0)

    goodFolder = "correct_synthetic_terrains_dataset/heightmaps/"
    badFolder = "correct_synthetic_terrains_dataset_even_other_colors/heightmaps/"
    files = glob.glob(badFolder + "*.png")
    for filename in files:
        shutil.copyfile(filename, goodFolder + prefix + os.path.basename(filename))
    goodFolder = "correct_synthetic_terrains_dataset/distortions/"
    badFolder = "correct_synthetic_terrains_dataset_even_other_colors/distortions/"
    files = glob.glob(badFolder + "*.png")
    for filename in files:
        shutil.copyfile(filename, goodFolder + prefix + os.path.basename(filename))


if __name__ == "__main__":
    main()
