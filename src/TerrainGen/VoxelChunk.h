#ifndef VOXELCHUNK_H
#define VOXELCHUNK_H

class VoxelChunk;

#include "DataStructure/Voxel.h"
#include <vector>
#include "Graphics/MarchingCubes.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Vertex.h"
#include "DataStructure/Matrix3.h"
#include <functional>
#include <memory>
#include <optional>

class VoxelChunk : public std::enable_shared_from_this<VoxelChunk>
{
public:
    VoxelChunk();
//    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, Matrix3<TerrainTypes> data, std::shared_ptr<VoxelGrid> parent);
    VoxelChunk(int x, int y, int sizeX, int sizeY,  int height, Matrix3<float> iso_data, std::shared_ptr<VoxelGrid> parent);
    ~VoxelChunk();

    void display();
    void createVoxels();
    void createMesh(bool applyMarchingCubes = true, bool updateMesh = true);

    void updateLoDsAvailable();

    void makeItFall();
    void letGravityMakeSandFall();
    void computeGroups();

    void applyModification(Matrix3<float> modifications);

    void resetVoxelsNeighbors();
    void computeFlowfield(Vector3 sea_current = Vector3()); //int blur_iterations = 5);

    std::vector<Vector3> applyMarchingCubes(bool useGlobalCoords = false, std::vector<Vector3>* outColors = nullptr);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    float getVoxelValue(int x, int y, int z);
    float getVoxelValue(Vector3 pos);
    Matrix3<float> getVoxelValues();

    std::vector<Vector3> computeMarchingCube(Vertex vertices[8], float isolevel, bool useGlobalCoords = false, std::vector<Vector3>* outColors = nullptr);

    Vector3 computeNormal(Vector3 pos);
    Vector3 computeNormal(int x, int y, int z);

    void computeDistanceField();


//protected:
//    Matrix3<TerrainTypes> data;
    Matrix3<float> iso_data;
    int x, y;
    int sizeX, sizeY, height;
    Matrix3<std::shared_ptr<Voxel>> voxels;
    Matrix3<float> voxelValues;
    Matrix3<float> originalVoxelValues;
    Matrix3<int> voxelGroups;
    Matrix3<Vector3> flowField;
    Matrix3<int> distanceField;
    Matrix3<float> pressureField;
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
    void applyTo(Matrix3<T>& array, F func) {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < height; h++) {
                    func(array(v_x, v_y, h));
                }
            }
        }
    }

    Matrix3<float>& toFloat();
    Matrix3<std::shared_ptr<Voxel>>& toVoxels();

    std::shared_ptr<VoxelGrid> parent;

    Mesh mesh;
    bool needRemeshing = true;

    std::vector<Matrix3<float>> voxelsValuesStack;
};

#endif // VOXELCHUNK_H
