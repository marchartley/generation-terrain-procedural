import time
from multiprocessing import Process
import os
import glob
import shutil

allFiles = glob.glob("/mnt/nvme0n1p5/LaTex/*/main.tex")

for file in allFiles:
    path = os.path.split(file)
    parent = os.path.split(path[-2])[-1]
    newFile = file.replace("main.tex", parent + ".tex")
    shutil.move(file, newFile)
print("Done")