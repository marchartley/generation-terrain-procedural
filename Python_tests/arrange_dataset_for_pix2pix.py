import glob
import math
import os
import re
from pathlib import Path
import shutil
from random import shuffle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_order(file, file_pattern):
    match = file_pattern.match(Path(file).name)
    if not match:
        return math.inf
    return int(match.groups()[0])


def moveToFolders(path):
    allFiles = glob.glob(path + "*.png")
    allFiles += glob.glob(path + "*/*.png")

    file_pattern = re.compile(r'.*?(\d+).*?')
    sortedFiles = sorted(allFiles, key=lambda f: get_order(f, file_pattern))
    sortedFiles = [(get_order(f, file_pattern), f) for f in sortedFiles]

    # All the following just to reassign the duplicate ID values...
    allIDs = [ID for ID, file in sortedFiles]
    uniques = np.unique(allIDs, return_counts=True)
    duplicates = {uniques[0][i]: uniques[1][i] for i in range(len(uniques[0]))}
    availableIDs = [i for i in range(len(sortedFiles)) if i not in allIDs]

    if len(availableIDs) != 0:
        currentSpot = 0
        for i, (ID, file) in enumerate(sortedFiles):
            if duplicates[ID] > 1:
                duplicates[ID] -= 1
                sortedFiles[i] = (availableIDs[currentSpot], file)
                currentSpot += 1

    os.makedirs(path + "train", exist_ok=True)
    os.makedirs(path + "val", exist_ok=True)
    os.makedirs(path + "test", exist_ok=True)

    trainLimit = int(len(sortedFiles) * 3 / 4)
    valLimit = trainLimit + int(len(sortedFiles) * 1 / 8)

    shuffle(sortedFiles)
    print(f"Arranging {path}")
    for i in range(len(sortedFiles)):
        if i % 50 == 0:
            print(f"{i + 1} / {len(sortedFiles)}")
        previous_filename = sortedFiles[i][1]
        fileID = sortedFiles[i][0]
        newFilename = path + ("train" if i < trainLimit else "val" if i < valLimit else "test") + "/" + str(
            fileID) + ".png"
        shutil.move(previous_filename, newFilename)


def main():
    mainPath = ""
    if os.name == 'nt':
        mainPath = "./correct_synthetic_terrains_dataset/"
    else:
        mainPath = "/data/correct_synthetic_terrains_dataset/"
    pathHeightmap = mainPath + "heightmaps/"
    pathFeatures = mainPath + "features/"
    pathDistortions = mainPath + "distortions/"
    os.makedirs(mainPath + "AB", exist_ok=True)

    moveToFolders(pathHeightmap)
    moveToFolders(pathFeatures)
    moveToFolders(pathDistortions)


if __name__ == "__main__":
    main()
