import glob
import os.path
import xml.etree.ElementTree

import numpy as np
import matplotlib.pyplot as plt
import h5py
import binascii

def loadMB41InHeightmap(img_file: str) -> np.ndarray:
    filename = os.path.split(img_file)[-1]

    with open(img_file, "rb") as file:
        i = 0
        for line in file:
            print(line)
            i += 1
            if i > 1000: return
    return


def loadBagInHeightmap(img_file: str) -> np.ndarray:
    filename = os.path.split(img_file)[-1]
    array = np.zeros((0, 0))
    print(f"Working on {img_file}")
    # File is encoded using HDF5 format
    f = h5py.File(img_file, "r")
    datasets = f.keys()
    if len(datasets) == 0:
        print("No data to read...")
        return None
    dset = f[list(datasets)[0]]
    if "elevation" not in list(dset.keys()):
        print("Did not find a 'elevation' dataset in file")
        print("Available datasets are : ", list(dset.keys()))
        return None

    array = np.array(dset["elevation"])

    if "metadata" in list(dset.keys()):
        meta = "".join([byte.decode() for byte in np.array(dset["metadata"]).flatten().tolist()])
        xmlData = xml.etree.ElementTree.fromstring(meta)
        queue = list((e, 0) for e in xmlData.findall("*"))
        while len(queue) > 0:
            elem, level = queue.pop()
            queue += list((e, level + 1) for e in elem.findall("*"))
            if elem.text:
                print(f"{'>' * level}{elem} : {elem.text}")
        return
    resolution = 0.001
    print(f"File {filename} is defined on {nCols * resolution * 100}x{nRows * resolution * 100}m with depths going from {np.max(array)}m to {np.min(array)}m")
    return array

def loadAscInHeightmap(img_file: str) -> np.ndarray:
    filename = os.path.split(img_file)[-1]
    array = np.zeros((0, 0))
    print(f"Working on {img_file}")
    with open(img_file, "r") as file:
        line = str(file.readline())
        nCols = int(line.split(" ")[-1])
        line = str(file.readline())
        nRows = int(line.split(" ")[-1])
        line = str(file.readline())
        xMini = float(line.split(" ")[-1])
        line = str(file.readline())
        yMini = float(line.split(" ")[-1])
        line = str(file.readline())
        resolution = float(line.split(" ")[-1])
        line = str(file.readline())
        nanValue = float(line.split(" ")[-1])

        array = np.zeros((nRows, nCols))
        i = 0
        for line in file:
            values = list(map(lambda x: 0 if x == nanValue else x, map(lambda s: float(s), filter(lambda s: s, line.split(" ")))))
            array[i] = values
            i += 1

        print(f"File {filename} is defined on {nCols * resolution * 1000}x{nRows * resolution * 1000}m with depths going from {np.max(array)}m to {np.min(array)}m")
    return array

def loadGlzInHeightmap(img_file: str) -> np.ndarray:
    filename = os.path.split(img_file)[-1]
    x = []
    y = []
    z = []
    array = np.zeros((0, 0))
    resolution = 1000.0
    print(f"Working on {img_file}")
    with open(img_file) as file:
        for line in file:
            line = line.replace("\t", " ")
            s_line = line.split(" ")
            orig_x, orig_y, orig_z = s_line
            _x, _y, _z = round(float(orig_x) * resolution), round(float(orig_y) * resolution), float(orig_z)
            x.append(_x)
            y.append(_y)
            z.append(_z)

        minX, maxX = min(x), max(x)
        minY, maxY = min(y), max(y)
        minZ, maxZ = min(z), max(z)

        array = np.zeros((maxY - minY + 1, maxX - minX + 1))

        file.seek(0)  # go back to start of file
        for line in file:
            line = line.replace("\t", " ")
            s_line = line.split(" ")
            orig_x, orig_y, orig_z = s_line
            _x, _y, _z = round(float(orig_x) * resolution), round(float(orig_y) * resolution), float(orig_z)
            array[(maxY - minY) - (_y - minY), _x - minX] = _z

        print(f"File {filename} is defined on {maxX - minX}x{maxY - minY}m with depths going from {minZ}m to {maxZ}m")
    return array


def main():
    folder = "C:/codes/Qt/Shared_folder_IRIT_Seafile/generation-terrain-procedural/Other_Data/"
    # folder = "C:/codes/Qt/Shared_folder_IRIT_Seafile/generation-terrain-procedural/DataShomData/"

    # for img_file in glob.glob(folder + "*.mb41"):
    #     heightmap = loadMB41InHeightmap(img_file)

    # for img_file in glob.glob(folder + "*.bag"):
    #     heightmap = loadBagInHeightmap(img_file)
    #     plt.imsave(img_file.replace(".bag", ".png"), heightmap, cmap="gray")
    #
    for img_file in glob.glob(folder + "*.asc"):
        heightmap = loadAscInHeightmap(img_file)
        plt.imsave(img_file.replace(".asc", ".png"), heightmap, cmap="gray")

    # for img_file in glob.glob(folder + "*.glz"):
    #     heightmap = loadGlzInHeightmap(img_file)
    #     plt.imsave(img_file.replace(".glz", ".png"), heightmap, cmap="gray")


if __name__ == "__main__":
    main()