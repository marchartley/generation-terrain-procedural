#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>
#include "MarchingCubes.h"
#include "Mesh.h"
#include "Vertex.h"
#include <functional>

class VoxelChunk
{
public:
    VoxelChunk();
//    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< TerrainTypes > > > data, VoxelGrid* parent);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< float > > > iso_data, VoxelGrid* parent);
    ~VoxelChunk();

    void display();
    void createMesh(bool applyMarchingCubes = true, bool updateMesh = true);

    void makeItFall();
    void letGravityMakeSandFall();
    void computeGroups();

    void resetVoxelsNeighbors();

    std::vector<Vector3> applyMarchingCubes(bool useGlobalCoords = false, std::vector<Vector3> *outColors = nullptr);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    std::vector<Vector3> computeMarchingCube(Vertex vertices[8], float isolevel, bool useGlobalCoords = false, std::vector<Vector3> *outColors = nullptr);

//protected:
//    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    std::vector<std::vector<std::vector<float>>> iso_data;
    int x, y;
    int sizeX, sizeY, height;
    std::vector<std::vector<std::vector<Voxel*>>> voxels;
    std::vector<std::vector<std::vector<float>>> voxelValues;
    std::vector<std::vector<std::vector<float>>> originalVoxelValues;
    std::vector<std::vector<std::vector<int>>> voxelGroups;
    int LoDIndex = 1;
    std::vector<int> LoDs;

    std::map<VOXEL_NEIGHBOR, VoxelChunk*> neighboring_chunks;
    bool lastChunkOnX = false, lastChunkOnY = false;

    template<class F>
    void applyToVoxels(F func) {
        /*for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    func(this->voxels[v_x][v_y][h]);
                }
            }
        }*/
        this->applyTo(this->voxels, func);
    }
    template<class F, class T>
    void applyTo(std::vector<std::vector<std::vector<T>>>& array, F func) {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    func(array[v_x][v_y][h]);
                }
            }
        }
    }

    std::vector<std::vector<std::vector<float>>>& toFloat();
    std::vector<std::vector<std::vector<Voxel*>>>& toVoxels();

    VoxelGrid* parent;

    Mesh mesh;
    bool needRemeshing = true;
};

#endif // VOXELCHUNK_H
