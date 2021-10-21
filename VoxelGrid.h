#ifndef VOXELGRID_H
#define VOXELGRID_H

class Voxel;
class VoxelGrid;

#include "Grid.h"
#include <vector>

enum TerrainTypes {
    AIR = 0,
    DIRT = 1,
    WATER = 2
};

class Voxel {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize);

    void display();

    int getX() { return this->x; }
    int getY() { return this->y; }
    int getZ() { return this->z; }
protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
};

class VoxelGrid {
public:
    VoxelGrid();
    VoxelGrid(Grid& grid);
    VoxelGrid(int nx, int ny, int nz, float blockSize);

    void from2DGrid(Grid grid);

    void display();

    int getSizeX() { return this->sizeX; }
    int getSizeY() { return this->sizeY; }
    int getSizeZ() { return this->sizeZ; }
    float getBlockSize() { return this->blockSize; }

    std::vector<Voxel> getVoxels() { return this->voxels; }
    int getHeight(int x, int y);

protected:
    int sizeX, sizeY, sizeZ;
    std::vector<Voxel> voxels;
    float blockSize;
};

#endif // VOXELGRID_H
