import glob
import os
import shutil

def main():
    src = "/data/ErosionAnimation/"
    dst = "/data/ErosionAnimation/results10/"
    os.makedirs(dst, exist_ok=True)
    files = glob.glob(src + "*/10.png")

    for filename in files:
        # filename = src + "/" + filename
        filename = filename
        basename = os.path.basename(filename)
        absPath = os.path.abspath(filename)[:-len(basename)]
        parentFolderName = os.path.basename(os.path.abspath(os.path.dirname(filename)))

        newFolder = absPath.replace(src, dst)
        shutil.copy(filename, dst + parentFolderName + ".png")
        # print(newFolder.rstrip("/") + ".png")
        # print(parentFolderName)

    print(files)

if __name__ == "__main__":
    main()