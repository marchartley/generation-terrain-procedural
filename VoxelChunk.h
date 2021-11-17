#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>
#include "MarchingCubes.h"
#include "Mesh.h"

class VoxelChunk
{
public:
    VoxelChunk();
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< TerrainTypes > > > data, VoxelGrid* parent);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< float > > > iso_data, VoxelGrid* parent);

    void display(bool apply_marching_cubes = false, bool display_vertices = false, float isolevel = 0.0);
    void createMesh(bool updateMesh = true);

    void makeItFall(int groupId = -1);
    void computeGroups();

    std::vector<Vector3> applyMarchingCubes();

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

    VoxelGrid* parent;

    Mesh mesh;
};

#endif // VOXELCHUNK_H
