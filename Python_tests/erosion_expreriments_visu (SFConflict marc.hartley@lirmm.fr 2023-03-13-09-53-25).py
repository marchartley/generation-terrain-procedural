import os.path
from typing import List, Set, Dict, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PIL.Image as Image
import sklearn.preprocessing
from sklearn.preprocessing import MinMaxScaler
import pickle

parameters: Dict[str, float] = {}
visualizedIteration: int = 99
allImages: List[List[Optional[np.ndarray]]] = []
originalDataFolder = "/data/erosionsTests/"
pickleDataFolder = "/data/_erosionsCachedTests/"


def main():
    global parameters
    global visualizedIteration
    global allImages

    def sliderChanged(attribute: str, newVal: float):
        global parameters
        parameters[attribute] = newVal
        updateImg()

    def setIteration(newIteration: int):
        global visualizedIteration
        visualizedIteration = int(newIteration)
        updateImg()

    def updateImg():
        global visualizedIteration
        global parameters
        displayMap(visualizedIteration, normalizedData, parameters, min_max_scaler, img_ax)

    def displayMap(iteration: int, normalizedData: pd.DataFrame, param: Dict[str, float],
                   scaler: MinMaxScaler, ax: plt.Axes):
        distances: List[float] = []
        data: List[List[float]] = []
        files: List[str] = []
        _img = allImages

        for i in range(normalizedData.shape[0]):
            data.append([])
            files.append(
                originalDataFolder + normalizedData["folder_name"].values[i] + "/screen/" + str(iteration) + ".jpg")

        param_data = []
        for col in param.keys():
            param_data.append(param[col])
            for i, item in enumerate(normalizedData[col]):
                data[i].append(item)
        param_array = scaler.transform(np.array([param_data]))

        for i in range(len(files)):
            dist = np.linalg.norm(param_array - np.array(data[i]))
            if dist == 0:
                dist = 0.000001
            distances.append(1 / dist ** 2)

        images: List[Image] = []
        for i in range(len(files)):
            try:
                images.append(allImages[i][iteration])
            except:
                images.append(None)
                print(f"Dead link : {files[i]}")

        finalImage = None
        i = 0
        while finalImage is None:
            if images[i] is not None:
                finalImage = np.zeros(np.array(images[i]).shape)
            i += 1

        distSum = 0
        finalImage = images[np.argmax(distances)]
        # for i in range(len(distances)):
        #     if images[i] is not None:
        #         finalImage += np.array(images[i]) * distances[i]
        #         distSum += distances[i]
        #     else:
        #         print(f"Image[{i}] is null : {files[i]}")
        # print(distSum, distances)
        # finalImage /= distSum
        ax.clear()
        ax.imshow(finalImage / 255.0)

    data = pd.read_csv(originalDataFolder + "allData.csv", header=0, delimiter=";", decimal=",")
    print(data)
    sliders: Dict[str, Optional[plt.Slider]] = {}

    fig: plt.Figure
    fig = plt.figure()
    # ax = fig.add_axes([.1, .1, .6, .8])
    fig.subplots_adjust(right=.25)
    # sliderAx: plt.Axes = fig.add_axes([0.7, 0.1, 0.3, 0.8])
    for col in data.columns:
        if data[col].dtype.kind in 'biufc' and data[col].min() != data[col].max():
            sliders[col] = None

    sliderHeight = .8 / (len(sliders) + 1)

    def sliderChangeFunction(columnName):
        return lambda v: sliderChanged(columnName, v)

    parameters = {}
    for i, col in enumerate(sliders.keys()):
        parameters[col] = data[col].mean()
        _ax = fig.add_axes([0.8, 0.1 + (sliderHeight * i), 0.2, sliderHeight])
        _ax.set_xticks([])
        _ax.set_yticks([])
        s: plt.Slider = plt.Slider(_ax, col, data[col].min(), data[col].max(), data[col].mean())
        func = sliderChangeFunction(col)
        s.on_changed(func)
        sliders[col] = s

    normalizedData = data[sliders.keys()]
    x = normalizedData.values
    min_max_scaler = MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    normalizedData.loc[:, :] = x_scaled
    # normalizedData = pd.DataFrame(x_scaled, columns=normalizedData.columns, index=normalizedData.index)
    # normalizedData = normalizedData.apply(lambda col: (col - col.min()) / (col.max() - col.min()), axis=0)
    normalizedData.loc[:, "folder_name"] = data["folder_name"]

    _ax = fig.add_axes([0.8, 0.1 + (sliderHeight * len(sliders)), 0.2, sliderHeight])
    _ax.set_xticks([])
    _ax.set_yticks([])
    s: plt.Slider = plt.Slider(_ax, "iteration", 0, 99, 99)
    s.on_changed(lambda x: setIteration(x))
    sliders["iteration"] = s

    img_ax = fig.add_axes([.1, .1, .6, .8])
    for col in parameters.keys():
        parameters[col] = data[col].iloc[1]

    print("Loading images :")
    for i, row in enumerate(normalizedData.iterrows()):
        print(f"{i + 1}/{len(normalizedData.values)}...")
        folder = row[1]["folder_name"]
        allImages.append([])
        storedFolder = pickleDataFolder
        pickleFilename = storedFolder + folder + ".pickle"
        if os.path.exists(pickleFilename):
            allImages[-1] = pickle.load(open(pickleFilename, 'rb'))
        else:
            if os.name == 'nt':
                os.makedirs(storedFolder, exist_ok=True)
            else:
                os.makedirs(storedFolder)

            for i in range(100):
                # try:
                allImages[-1].append(np.array(
                    Image.open(originalDataFolder + folder + "/screen/" + str(i) + ".jpg")))  # .resize((150, 150))))
                # except:
                #     allImages[-1].append(None)
            pickle.dump(allImages[-1], open(pickleFilename, 'wb'))

    displayMap(99, normalizedData, parameters, min_max_scaler, img_ax)
    plt.show()


if __name__ == "__main__":
    main()
