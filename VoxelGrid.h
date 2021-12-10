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
#include "FastNoiseLit.h"

struct NoiseMinMax {
    float min;
    float max;

    void update(float val) {
        min = std::min(min, val);
        max = std::max(max, val + 0.00001f);
    }
    float remap(float x, float newMin, float newMax) {
        /*float a = (x - min);
        float b = a / (max - min);
        float c = b - .0;
        float d = c * (newMax - newMin);
        float e = d + newMin;*/
        return (((x - min) / (max - min)) * (newMax - newMin) + newMin);
    }
};

class VoxelGrid : public std::enable_shared_from_this<VoxelGrid> {
public:
    VoxelGrid();
    VoxelGrid(Grid& grid);
    VoxelGrid(int nx, int ny, int nz, float blockSize, float noise_shifting = 0.0);
    ~VoxelGrid();
    void from2DGrid(Grid grid);
    std::shared_ptr<VoxelGrid> fromIsoData();

    void initMap();

    void display();

    void createMesh();

    void makeItFall(float erosionStrength = 0.0);
    void letGravityMakeSandFall(bool remesh = true);

    int getSizeX() { return this->sizeX; }
    int getSizeY() { return this->sizeY; }
    int getSizeZ() { return this->sizeZ; }
    float getBlockSize() { return this->blockSize; }

    std::vector<Voxel> getVoxels() { return this->voxels; }
    std::vector<std::vector<std::vector<std::vector<float>>>> tempData;
    int getHeight(int x, int y);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);
    std::shared_ptr<Voxel> getVoxel(Vector3 pos);
    std::shared_ptr<Voxel> getVoxel(float x, float y, float z);

    void remeshAll();

    int numberOfChunksX() { return this->sizeX / this->chunkSize; }
    int numberOfChunksY() { return this->sizeY / this->chunkSize; }

    std::string toString();
    std::string toShortString();

    std::vector<std::vector<std::vector<float>>> toFloat();
    void toVoxels(std::vector<std::vector<std::vector<float>>> arr);
    void toVoxels();

    void computeVoxelGroups();

    std::tuple<int, int, int, int> getChunksAndVoxelIndices(Vector3 pos);
    std::tuple<int, int, int, int> getChunksAndVoxelIndices(float x, float y, float z);
    float getVoxelValue(Vector3 pos);
    float getVoxelValue(float x, float y, float z);
    void setVoxelValue(Vector3 pos, float newVal);
    void setVoxelValue(float x, float y, float z, float newVal);
    float getOriginalVoxelValue(Vector3 pos);
    float getOriginalVoxelValue(float x, float y, float z);
    Vector3 getFlowfield(Vector3 pos);
    Vector3 getFlowfield(float x, float y, float z);
    void setFlowfield(Vector3 pos, Vector3 newVal);
    void setFlowfield(float x, float y, float z, Vector3 newVal);
    int getVoxelGroup(Vector3 pos);
    int getVoxelGroup(float x, float y, float z);
    void setVoxelGroup(Vector3 pos, int newVal);
    void setVoxelGroup(float x, float y, float z, int newVal);
    bool getVoxelIsOnGround(Vector3 pos);
    bool getVoxelIsOnGround(float x, float y, float z);
    void setVoxelIsOnGround(Vector3 pos, bool newVal);
    void setVoxelIsOnGround(float x, float y, float z, bool newVal);

    int getMaxLoD();
//protected:
    int sizeX, sizeY, sizeZ;
    std::vector<Voxel> voxels;
    float blockSize;
    std::vector<std::shared_ptr<VoxelChunk>> chunks;
    float noise_shifting;

    int chunkSize = 30;
    bool displayWithMarchingCubes = false;
    FastNoiseLite noise;
    NoiseMinMax noiseMinMax;

    Mesh mesh;
};

#endif // VOXELGRID_H
