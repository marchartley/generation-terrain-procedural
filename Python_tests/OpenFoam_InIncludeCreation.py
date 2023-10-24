import glob
import os.path
import shutil

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    rootFolder = "C:/codes/Qt/Shared_folder_IRIT_Seafile/generation-terrain-procedural/src/third-party/openFoam/"
    copyFolder = "C:/codes/Qt/Shared_folder_IRIT_Seafile/generation-terrain-procedural/src/third-party/openFoam/InInclude/"

    os.makedirs(copyFolder, exist_ok=True)
    files = glob.glob(rootFolder + "**", recursive=True)
    for f in files:
        if f.endswith(".H") or f.endswith(".C"):
            filename = os.path.basename(f)
            try:
                shutil.copyfile(f, copyFolder + filename)
            except shutil.SameFileError:
                pass
            print(filename)
            # break


if __name__ == "__main__":
    main()
