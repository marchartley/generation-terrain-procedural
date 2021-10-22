#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>

class VoxelChunk
{
public:
    VoxelChunk();
    VoxelChunk(int x, int y, int height, std::vector<TerrainTypes> data, int surrounding_heights[4]);

    void display();

//protected:
    std::vector<TerrainTypes> data;
    int height;
    int x, y;
    int surrounding_heights[4];

    void draw_part(int start_height, int end_height, TerrainTypes type, bool draw_bottom);
};

#endif // VOXELCHUNK_H
