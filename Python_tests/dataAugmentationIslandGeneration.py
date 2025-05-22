import cProfile
import glob
import os.path
import pstats
import random
import time
from collections.abc import Callable
from typing import Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PIL.Image
import PIL.ImageChops

from Python_tests.smoothmax_tests import smoothmax

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
#
# def addHeightmaps(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
#     curr = np.array(currentImg)
#     paste = np.array(addedImg)
#     w = paste.shape[1]
#     h = paste.shape[0]
#     for x in range(curr.shape[1]):
#         for y in range(curr.shape[0]):
#             _x = x - position[0]
#             _y = y - position[1]
#             if _x < 0 or _x >= w or _y < 0 or _y >= h: continue
#
#             prev_val = curr[y, x]
#             paste_val = paste[_y, _x]
#             new_val = (prev_val**5 + paste_val**5)**(1/5) # smoothmax(prev_val, paste_val, k=1.0)
#             curr[y, x] = new_val
#     currentImg.paste(PIL.Image.fromarray(curr.astype(np.uint8)))
#     return currentImg
#
#
# def addFeatures(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
#     curr = np.array(currentImg.copy().convert('HSV'))
#     paste = np.array(addedImg.copy().convert('HSV'))
#     w = paste.shape[1]
#     h = paste.shape[0]
#     for x in range(curr.shape[1]):
#         for y in range(curr.shape[0]):
#             _x = x - position[0]
#             _y = y - position[1]
#             if _x < 0 or _x >= w or _y < 0 or _y >= h: continue
#
#
#             prev_val = curr[y, x]
#             paste_val = paste[_y, _x]
#             new_val = (min(float(prev_val[0]), float(paste_val[0])) if (prev_val[0] != 0 and paste_val[0] != 0) else max(float(paste_val[0]), float(prev_val[0])),
#                        min(float(prev_val[1]), float(paste_val[1])), min(float(prev_val[2]), float(paste_val[2])))
#             curr[y, x] = new_val
#     currentImg.paste(PIL.Image.fromarray(curr.astype(np.uint8), 'HSV'))
#     return currentImg
#
# def addDistortions(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
#     curr = np.array(currentImg, dtype=np.uint8)
#     paste = np.array(addedImg, dtype=np.uint8)
#     w = paste.shape[1]
#     h = paste.shape[0]
#     for x in range(curr.shape[1]):
#         for y in range(curr.shape[0]):
#             _x = x - position[0]
#             _y = y - position[1]
#             if _x < 0 or _x >= w or _y < 0 or _y >= h: continue
#
#             prev_val = curr[y, x]
#             paste_val = paste[_y, _x]
#             new_val = (((prev_val[0] - 127) + (paste_val[0] - 127)) + 127, ((prev_val[1] - 127) + (paste_val[1] - 127)) + 127, ((prev_val[2] - 127) + (paste_val[2] - 127)) + 127)
#             curr[y, x] = new_val
#     currentImg.paste(PIL.Image.fromarray(curr.astype(np.uint8)))
#     return currentImg

def addHeightmaps(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
    curr = np.array(currentImg, dtype=np.float32)
    paste = np.array(addedImg, dtype=np.float32)

    x_off, y_off = position
    h, w = paste.shape[:2]

    # Calculate bounds
    y1, y2 = max(0, y_off), min(curr.shape[0], y_off + h)
    x1, x2 = max(0, x_off), min(curr.shape[1], x_off + w)
    py1, py2 = y1 - y_off, y2 - y_off
    px1, px2 = x1 - x_off, x2 - x_off

    # Apply smooth maximum (fifth root of sum of fifth powers)
    region = curr[y1:y2, x1:x2]
    paste_region = paste[py1:py2, px1:px2]
    result = np.power(np.power(region, 5) + np.power(paste_region, 5), 1/5)

    curr[y1:y2, x1:x2] = result
    currentImg.paste(PIL.Image.fromarray(np.clip(curr, 0, 255).astype(np.uint8)))
    return currentImg

def addFeatures(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
    curr = np.array(currentImg.convert('HSV'), dtype=np.uint8)
    paste = np.array(addedImg.convert('HSV'), dtype=np.uint8)

    x_off, y_off = position
    h, w = paste.shape[:2]

    y1, y2 = max(0, y_off), min(curr.shape[0], y_off + h)
    x1, x2 = max(0, x_off), min(curr.shape[1], x_off + w)
    py1, py2 = y1 - y_off, y2 - y_off
    px1, px2 = x1 - x_off, x2 - x_off

    region = curr[y1:y2, x1:x2].astype(np.float32)
    paste_region = paste[py1:py2, px1:px2].astype(np.float32)

    hue_mask = (region[..., 0] != 0) & (paste_region[..., 0] != 0)
    hue = np.where(hue_mask, np.minimum(region[..., 0], paste_region[..., 0]),
                              np.maximum(region[..., 0], paste_region[..., 0]))

    saturation = np.minimum(region[..., 1], paste_region[..., 1])
    value = np.minimum(region[..., 2], paste_region[..., 2])

    result = np.stack([hue, saturation, value], axis=-1)

    curr[y1:y2, x1:x2] = result.astype(np.uint8)
    currentImg.paste(PIL.Image.fromarray(curr, mode='HSV'))
    return currentImg

def addDistortions(currentImg: PIL.Image.Image, addedImg: PIL.Image.Image, position: Tuple[int, int]) -> PIL.Image.Image:
    curr = np.array(currentImg, dtype=np.int16)
    paste = np.array(addedImg, dtype=np.int16)

    x_off, y_off = position
    h, w = paste.shape[:2]

    y1, y2 = max(0, y_off), min(curr.shape[0], y_off + h)
    x1, x2 = max(0, x_off), min(curr.shape[1], x_off + w)
    py1, py2 = y1 - y_off, y2 - y_off
    px1, px2 = x1 - x_off, x2 - x_off

    region = curr[y1:y2, x1:x2]
    paste_region = paste[py1:py2, px1:px2]

    result = (region - 127) + (paste_region - 127) + 127
    result = np.clip(result, 0, 255)

    curr[y1:y2, x1:x2] = result
    currentImg.paste(PIL.Image.fromarray(curr.astype(np.uint8)))
    return currentImg





def main():
    dataset_folder = "_new_synthetic_terrains_dataset/"
    heightmaps_folder = dataset_folder + "heightmaps/"
    features_folder = dataset_folder + "features/"
    distortions_folder = dataset_folder + "distortions/"

    allHeightmaps = glob.glob(heightmaps_folder + "*.png")
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
                newSize = (int(getRandom(30, 128)), int(getRandom(30, 128)))
                functions = [(rotate, getRandom(-180, 180)), (resize, newSize[0], newSize[1])]
                for tries in range(10):
                    ok = True
                    box = [int(getRandom(0, 256)), int(getRandom(0, 256))]
                    box += [box[0] + newSize[0], box[1] + newSize[1]]
                    for otherBox in usedBoxes:
                        if bb_intersection_over_union(box, otherBox) > 0.1:
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
            results = [PIL.Image.new("L", originals[0].size, defaultColors[0]), PIL.Image.new("RGB", originals[1].size, originals[1].img.getpixel((0, 0))),
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
                    # res.paste(orig.img, position)
                    if i == 0:
                        addHeightmaps(res, orig.img, position)
                    elif i == 1:
                        addFeatures(res, orig.img, position)
                    elif i == 2:
                        addDistortions(res, orig.img, position)
                    else:
                        print("Impossible case")
            for i, res in enumerate(results):
                res.rotate(0).save  (folders[i] + nbToAlpha(iTransfoCombi) + "a_" + original_filename)
                res.rotate(90).save (folders[i] + nbToAlpha(iTransfoCombi) + "b_" + original_filename)
                res.rotate(180).save(folders[i] + nbToAlpha(iTransfoCombi) + "c_" + original_filename)
                res.rotate(270).save(folders[i] + nbToAlpha(iTransfoCombi) + "d_" + original_filename)


                # os.remove(folders[i] + nbToAlpha(iTransfoCombi) + "a_" + original_filename)
                # os.remove(folders[i] + nbToAlpha(iTransfoCombi) + "b_" + original_filename)
                # os.remove(folders[i] + nbToAlpha(iTransfoCombi) + "c_" + original_filename)
                # os.remove(folders[i] + nbToAlpha(iTransfoCombi) + "d_" + original_filename)

            # plt.imshow(results[0], cmap="gray", vmin=0, vmax=255)
            # plt.show()
            # plt.imshow(results[1], cmap="gray", vmin=0, vmax=255)
            # plt.show()
            # plt.imshow(results[2], cmap="gray", vmin=0, vmax=255)
            # plt.show()


if __name__ == "__main__":
    # with cProfile.Profile() as profiler:
        main()
    # stats = pstats.Stats(profiler)
    # stats.sort_stats(pstats.SortKey.TIME)  # Sorts by time spent in each function
    # stats.print_stats(20)  # Prints top 20 lines of profiling results

