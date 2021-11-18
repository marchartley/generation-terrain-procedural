#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>
#include "MarchingCubes.h"
#include "Mesh.h"
#include <functional>

class VoxelChunk
{
public:
    VoxelChunk();
//    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< TerrainTypes > > > data, VoxelGrid* parent);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< float > > > iso_data, VoxelGrid* parent);

    void display();
    void createMesh(bool applyMarchingCubes = true, bool updateMesh = true);

    void makeItFall(int groupId = -1);
    void computeGroups();

    std::vector<Vector3> applyMarchingCubes();

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

//protected:
//    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    std::vector<std::vector<std::vector<float>>> iso_data;
    int x, y;
    int sizeX, sizeY, height;
    std::vector<std::vector<std::vector<Voxel*>>> voxels;
    std::map<VOXEL_NEIGHBOR, VoxelChunk*> neighboring_chunks;
    bool lastChunkOnX = false, lastChunkOnY = false;

    template<class R>
    void applyToVoxels(std::function<R(Voxel*)> func) {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    func(this->voxels[v_x][v_y][h]);
                }
            }
        }
    }
    template <class F>
    void applyToVoxels(F func) {
        this->applyToVoxels(std::function<void(Voxel*)>(func));
    }

    VoxelGrid* parent;

    Mesh mesh;
    bool needRemeshing = true;
};

#endif // VOXELCHUNK_H
