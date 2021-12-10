#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "Voxel.h"
#include <vector>
#include "MarchingCubes.h"
#include "Mesh.h"
#include "Vertex.h"
#include <functional>
#include <memory>
#include <optional>

class VoxelChunk : public std::enable_shared_from_this<VoxelChunk>
{
public:
    VoxelChunk();
//    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< TerrainTypes > > > data, std::shared_ptr<VoxelGrid> parent);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, std::vector<std::vector<std::vector< float > > > iso_data, std::shared_ptr<VoxelGrid> parent);
    ~VoxelChunk();

    void display();
    void createVoxels();
    void createMesh(bool applyMarchingCubes = true, bool updateMesh = true);

    void updateLoDsAvailable();

    void makeItFall();
    void letGravityMakeSandFall();
    void computeGroups();

    void resetVoxelsNeighbors();
    void computeFlowfield(Vector3 sea_current = Vector3()); //int blur_iterations = 5);

    std::vector<Vector3> applyMarchingCubes(bool useGlobalCoords = false, std::vector<Vector3>* outColors = nullptr);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    std::vector<Vector3> computeMarchingCube(Vertex vertices[8], float isolevel, bool useGlobalCoords = false, std::vector<Vector3>* outColors = nullptr);

    Vector3 computeNormal(Vector3 pos);
    Vector3 computeNormal(int x, int y, int z);

    void computeDistanceField();


//protected:
//    std::vector<std::vector<std::vector<TerrainTypes>>> data;
    std::vector<std::vector<std::vector<float>>> iso_data;
    int x, y;
    int sizeX, sizeY, height;
    std::vector<std::vector<std::vector<std::shared_ptr<Voxel>>>> voxels;
    std::vector<std::vector<std::vector<float>>> voxelValues;
    std::vector<std::vector<std::vector<float>>> originalVoxelValues;
    std::vector<std::vector<std::vector<int>>> voxelGroups;
    std::vector<std::vector<std::vector<Vector3>>> flowField;
    std::vector<std::vector<std::vector<int>>> distanceField;
    std::vector<std::vector<std::vector<float>>> pressureField;
    int LoDIndex = 1;
    std::vector<int> LoDs;

    std::map<VOXEL_NEIGHBOR, std::shared_ptr<VoxelChunk>> neighboring_chunks;
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
    std::vector<std::vector<std::vector<std::shared_ptr<Voxel>>>>& toVoxels();

    std::shared_ptr<VoxelGrid> parent;

    Mesh mesh;
    bool needRemeshing = true;
};

#endif // VOXELCHUNK_H
