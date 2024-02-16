import glob
import shutil

pattern = " (SFConflict marc.hartley@lirmm.fr 2024-01-21-15-"
files = glob.glob(f"./*{pattern}*.cpp", recursive=True) +  glob.glob(f"./*/*{pattern}*.cpp", recursive=True) + glob.glob(f"./*/*{pattern}*.h*", recursive=True)

for file in files:
    patternPos = file.find(pattern)
    originalFilename = file.split(pattern)[0]
    originalExtension = file.split(".")[-1]
    originalFilename = originalFilename + "." + originalExtension
    res = shutil.move(file, originalFilename)
    print(f"{originalFilename} -> {res}")

# print(files)