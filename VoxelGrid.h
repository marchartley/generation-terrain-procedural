#ifndef VOXELGRID_H
#define VOXELGRID_H

class VoxelGrid;

#include "Grid.h"
#include "VoxelChunk.h"
#include "Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Mesh.h"

class VoxelGrid {
public:
    VoxelGrid();
    VoxelGrid(Grid& grid);
    VoxelGrid(int nx, int ny, int nz, float blockSize, float noise_shifting = 0.0);

    void from2DGrid(Grid grid);

    void display();

    void createMesh();

    void makeItFall(int groupId = -1);

    int getSizeX() { return this->sizeX; }
    int getSizeY() { return this->sizeY; }
    int getSizeZ() { return this->sizeZ; }
    float getBlockSize() { return this->blockSize; }

    std::vector<Voxel> getVoxels() { return this->voxels; }
    int getHeight(int x, int y);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);
    Voxel* getVoxel(Vector3 pos);
    Voxel* getVoxel(float x, float y, float z);

    void remeshAll();

    std::string toString();
    std::string toShortString();

//protected:
    int sizeX, sizeY, sizeZ;
    std::vector<Voxel> voxels;
    float blockSize;
    std::vector<VoxelChunk*> chunks;
    float noise_shifting;

    int chunkSize = 40;
    bool displayWithMarchingCubes = false;
};

#endif // VOXELGRID_H
