#ifndef VOXEL_H
#define VOXEL_H

enum TerrainTypes {
    AIR,
    CURRENT_MIDDLE,
    CURRENT_TOP,
    CURRENT_BOTTOM,
    WATER,
    CORAL,
    SAND,
    DIRT,
    ROCK,
    BEDROCK,
    LAST
};

enum VOXEL_NEIGHBOR {
    TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5
};


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "DataStructure/Matrix3.h"

class VoxelDataFile {
public:
    VoxelDataFile();
    VoxelDataFile(const GridF& dVec);

    void write(const std::string& filename);
    void load(const std::string& filename);

    GridF data;

protected:

    void loadFromFileBinary(const std::string& filename);
    void loadFromFileOld(const std::string& filename);
};

class VectorFieldDataFile {
public:
    VectorFieldDataFile();
    VectorFieldDataFile(const GridV3& dVec);

    void write(const std::string& filename);
    void load(const std::string& filename);

    GridV3 data;

protected:

    void loadFromFileBinary(const std::string& filename);
    void loadFromFileOld(const std::string& filename);
};

#endif // VOXEL_H
