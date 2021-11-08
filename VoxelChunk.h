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
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< TerrainTypes > > > data);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< float > > > iso_data);

    void display(bool apply_marching_cubes = false, bool display_vertices = false, float isolevel = 0.0);
    void createMesh();

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

//protected:
    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    std::vector<std::vector<std::vector<float>>> iso_data;
    int height;
    int x, y;
    int sizeX, sizeY;
    int surrounding_heights[4];
    std::vector<std::vector<std::vector<Voxel*>>> voxels;
    std::map<VOXEL_NEIGHBOR, VoxelChunk*> neighboring_chunks;
    bool needRemeshing = true;
    bool lastChunkOnX = false, lastChunkOnY = false;
};

#endif // VOXELCHUNK_H
