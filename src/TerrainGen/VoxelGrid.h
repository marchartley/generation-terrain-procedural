#ifndef VOXELGRID_H
#define VOXELGRID_H

class VoxelGrid;

#include "TerrainGen/Grid.h"
#include "TerrainGen/VoxelChunk.h"
#include "DataStructure/Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include <limits>
#include "Graphics/Mesh.h"
#include "Utils/FastNoiseLit.h"
#include "DataStructure/Matrix3.h"

#include "FluidSimulation/FluidSimulation.h"
#include "src/sim-fluid-ethanjli/fluidsystem.h"

struct NoiseMinMax {
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::min();

    void update(float val) {
        min = std::min(min, val);
        max = std::max(max, val + 0.00001f);
    }
    float remap(float x, float newMin, float newMax) {
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
    void letGravityMakeSandFallWithFlow(bool remesh = true);
    Matrix3<float> shareSandWithNeighbors(); // Doesn't affect the grid directly, but changes are returned to be applied after
    void applyModification(Matrix3<float> modifications);
    void add2DHeightModification(Matrix3<float> heightmapModifier, float factor = 1.f);
    void undo();
    void redo();

    int getSizeX() { return this->sizeX; }
    int getSizeY() { return this->sizeY; }
    int getSizeZ() { return this->sizeZ; }
    Vector3 getDimensions() { return Vector3(getSizeX(), getSizeY(), getSizeZ()); }
    float getBlockSize() { return this->blockSize; }

    std::vector<Matrix3<float>> tempData;
    int getHeight(int x, int y);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    void remeshAll();

    int numberOfChunksX() { return this->sizeX / this->chunkSize; }
    int numberOfChunksY() { return this->sizeY / this->chunkSize; }

    std::string toString();
    std::string toShortString();

    void smoothVoxels();


    void computeVoxelGroups();
    Matrix3<float> getVoxelValues();

    std::tuple<int, int, int, int> getChunksAndVoxelIndices(Vector3 pos);
    std::tuple<int, int, int, int> getChunksAndVoxelIndices(float x, float y, float z);
    float getVoxelValue(Vector3 pos);
    float getVoxelValue(float x, float y, float z);
    void setVoxelValue(Vector3 pos, float newVal);
    void setVoxelValue(float x, float y, float z, float newVal);
    float getOriginalVoxelValue(Vector3 pos);
    float getOriginalVoxelValue(float x, float y, float z);
    Matrix3<Vector3> getFlowfield();
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

    void computeFlowfield();

    void affectFlowfieldAround(Vector3 pos, Vector3 newVal, int kernelSize = 3);
    void affectFlowfieldAround(float x, float y, float z, Vector3 newVal, int kernelSize = 3);
    void affectFlowfieldAround(Vector3 pos, float alphaEffect, int kernelSize = 3);
    void affectFlowfieldAround(float x, float y, float z, float alphaEffect, int kernelSize = 3);

    int getMaxLoD();

    void saveState();

    void saveMap(std::string filename);
    void retrieveMap(std::string filename);

//protected:
    int sizeX, sizeY, sizeZ;
    float blockSize;
    std::vector<std::shared_ptr<VoxelChunk>> chunks;
    float noise_shifting;

    int chunkSize = 30;
    bool displayWithMarchingCubes = false;
    FastNoiseLite noise;
    NoiseMinMax noiseMinMax;

    Mesh mesh;
    FluidSimulation fluidSimulation;
    std::unique_ptr<FluidSystem> fluidSystem;
    std::unique_ptr<VelocityField> constantFlowSource;

    int fluidSimRescale = 4;


    bool _smoothingNeeded = false; // Just used when we come from a 2D grid.

    Matrix3<Vector3> flowField;
    Matrix3<int> distanceField;
    Matrix3<float> pressureField;

    Vector3 sea_current = Vector3(1.f, 0.0, 0.0);

protected:
    float getNoiseValue(int x, int y, int z, float noise_shift = 0.f);
};

#endif // VOXELGRID_H
