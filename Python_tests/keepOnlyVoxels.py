import os
import shutil
import glob


def main():
    initialFolder = os.path.abspath("/media/simulateurrsm/Marc_perso_data/PreviousComputer/") #"../saved_maps/Geometry/")
    # initialFolder = os.path.abspath("/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/saved_maps/Geometry/") #"../saved_maps/Geometry/")
    dstFolder = os.path.abspath("/media/simulateurrsm/FC3F-140B/DonneesErosion/") #"../saved_maps/onlyVoxels/")

    allFiles = glob.glob("*/*voxels.stl", root_dir=initialFolder, recursive=True) + glob.glob("*.stl", root_dir=initialFolder, recursive=True)
    for i, filename in enumerate(allFiles):
        print(f"{100 * (i + 1) / len(allFiles) }% {(i + 1)} / {len(allFiles)}")
        filename = initialFolder + "/" + filename
        basename = os.path.basename(filename)
        absPath = os.path.abspath(filename)[:-len(basename)]

        newFolder = absPath.replace(initialFolder, dstFolder)
        os.makedirs(newFolder, exist_ok=True)
        shutil.copy(filename, newFolder + basename)
        # print(filename)
        # print(newFolder + basename)
        # print()
    print(initialFolder, dstFolder)

if __name__ == "__main__":
    main()