#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>
#include "MarchingCubes.h"

class VoxelChunk
{
public:
    VoxelChunk();
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector<TerrainTypes>>> data);

    void display(bool apply_marching_cubes = false, bool display_vertices = false, float isolevel = 0.0);
    void createMesh();

//protected:
    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    int height;
    int x, y;
    int sizeX, sizeY;
    int surrounding_heights[4];
    std::vector<std::vector<std::vector<Voxel*>>> voxels;
    std::map<VOXEL_NEIGHBOR, VoxelChunk*> neighboring_chunks;
    bool lastChunkOnX = false, lastChunkOnY = false;
};

#endif // VOXELCHUNK_H
