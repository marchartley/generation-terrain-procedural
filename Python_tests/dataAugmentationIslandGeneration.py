import glob
import os.path
import random
from collections.abc import Callable

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PIL.Image
import PIL.ImageChops

def bb_intersection_over_union(boxA, boxB):
    # determine the (x, y)-coordinates of the intersection rectangle
    xA = max(boxA[0], boxB[0])
    yA = max(boxA[1], boxB[1])
    xB = min(boxA[2], boxB[2])
    yB = min(boxA[3], boxB[3])

    # compute the area of intersection rectangle
    interArea = abs(max((xB - xA, 0)) * max((yB - yA), 0))
    if interArea == 0:
        return 0
    # compute the area of both the prediction and ground-truth
    # rectangles
    boxAArea = abs((boxA[2] - boxA[0]) * (boxA[3] - boxA[1]))
    boxBArea = abs((boxB[2] - boxB[0]) * (boxB[3] - boxB[1]))

    # compute the intersection over union by taking the intersection
    # area and dividing it by the sum of prediction + ground-truth
    # areas - the interesection area
    iou = interArea / float(boxAArea + boxBArea - interArea)

    # return the intersection over union value
    return iou

class Image2:
    def __init__(self, other: PIL.Image.Image, fillColor = 0):
        self.img = other
        self.fillColor = fillColor

    def __getattr__(self, item):
        if item == "fillColor":
            return self.fillColor
        return self.img.__getattribute__(item)

def nbToAlpha(x: int, nbChars: int = 4):
    ret = ""
    if x == 0:
        ret = "A"
    while x > 0:
        x, r = divmod(x, 27)
        ret = chr(ord("A") + r) + ret
    if len(ret) < nbChars:
        ret = ret + ("_" * (nbChars - len(ret)))
    return ret

def main():
    dataset_folder = ("./correct_synthetic_terrains_dataset/" if os.name == "nt" else "/data/correct_synthetic_terrains_dataset/")
    heightmaps_folder = dataset_folder + "heightmaps/"
    features_folder = dataset_folder + "features/"
    distortions_folder = dataset_folder + "distortions/"

    allHeightmaps = glob.glob(heightmaps_folder + "[0-9]*.png")
    for iFile, fullpath in enumerate(allHeightmaps):
        print(f"{iFile + 1}/{len(allHeightmaps)}")
        original_filename = os.path.basename(fullpath)

        folders = [heightmaps_folder, features_folder, distortions_folder]
        filenames = [f + original_filename for f in folders]
        originals = [Image2(PIL.Image.open(filenames[0]).convert("L")), Image2(PIL.Image.open(filenames[1]).convert("RGB")), Image2(PIL.Image.open(filenames[2]).convert("RGB"))]
        defaultColors = [0, (255, 0, 0), (127, 127, 127)]
        for i in range(len(originals)):
            originals[i].fillColor = defaultColors[i]

        def offset(img: Image2, dx: int, dy: int):
            return Image2(PIL.ImageChops.offset(img.img, dx, dy), img.fillColor)

        def rotate(img: Image2, angle: float):
            return Image2(img.rotate(angle, fillcolor=img.fillColor), img.fillColor)

        def resize(img: Image2, newSizeX: int, newSizeY: int):
            return Image2(img.resize((newSizeX, newSizeY), resample=PIL.Image.NEAREST), img.fillColor)

        def getRandom(start, end):
            return start + random.random() * (end - start)

        def copiesFunc(nb):
            usedBoxes = []
            allFunctions = []
            for _ in range(nb):
                newSize = (int(getRandom(30, 100)), int(getRandom(30, 100)))
                functions = [(rotate, getRandom(-180, 180)), (resize, newSize[0], newSize[1])]
                for tries in range(10):
                    ok = True
                    box = [int(getRandom(0, 128)), int(getRandom(0, 128))]
                    box += [box[0] + newSize[0], box[1] + newSize[1]]
                    for otherBox in usedBoxes:
                        if bb_intersection_over_union(box, otherBox) > 0:
                            ok = False
                            break
                    if ok:
                        functions.append((box[0], box[1]))
                        usedBoxes.append(box)
                        allFunctions.append(functions)
                        break
            return allFunctions

        transfos = [] + \
                   [[[(offset, int(getRandom(-100, 100)), int(getRandom(-100, 100)))]] for _ in range(4)] + \
                   [copiesFunc(5) for _ in range(4)]

        for iTransfoCombi, combination in enumerate(transfos):
            results = [PIL.Image.new("L", originals[0].size, defaultColors[0]), PIL.Image.new("RGB", originals[1].size, defaultColors[1]),
                       PIL.Image.new("RGB", originals[2].size, defaultColors[2])]

            for iSubImg in range(len(combination)):
                tmpOrig = [Image2(img.copy(), img.fillColor) for img in originals]
                position = (0, 0)
                for transfo in combination[iSubImg]:
                    allArgs = transfo
                    if isinstance(allArgs[0], Callable):
                        func, *args = allArgs
                        for i in range(len(tmpOrig)):
                            tmpOrig[i] = func(tmpOrig[i], *args)
                    else:
                        position = allArgs

                for i, (orig, res) in enumerate(zip(tmpOrig, results)):
                    res.paste(orig.img, position)

            for i, res in enumerate(results):
                res.save(folders[i] + nbToAlpha(iTransfoCombi) + "a_" + original_filename)
                res.rotate(90).save(folders[i] + nbToAlpha(iTransfoCombi) + "b_" + original_filename)
                res.rotate(180).save(folders[i] + nbToAlpha(iTransfoCombi) + "c_" + original_filename)
                res.rotate(270).save(folders[i] + nbToAlpha(iTransfoCombi) + "d_" + original_filename)
            plt.show()


if __name__ == "__main__":
    main()
