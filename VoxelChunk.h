#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>

class VoxelChunk
{
public:
    VoxelChunk();
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector<TerrainTypes>>> data);

    void display(bool apply_marching_cubes = false, bool display_vertices = false);
    void createMesh();

//protected:
    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    int height;
    int x, y;
    int sizeX, sizeY;
    int surrounding_heights[4];
    std::vector<std::vector<std::vector<Voxel*>>> voxels;
    std::map<VOXEL_NEIGHBOR, VoxelChunk*> neighboring_chunks;


    void draw_part(int start_height, int end_height, TerrainTypes type, bool draw_bottom);
};

#endif // VOXELCHUNK_H
