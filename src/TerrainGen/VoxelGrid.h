#ifndef VOXELGRID_H
#define VOXELGRID_H

class VoxelGrid;

#include "TerrainGen/Heightmap.h"
#include "TerrainGen/VoxelChunk.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "DataStructure/Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include <limits>
#include "Graphics/Mesh.h"
#include "Utils/FastNoiseLit.h"
#include "DataStructure/Matrix3.h"
#include "TerrainGen/ImplicitPatch.h"

#include "FluidSimulation/StableFluidsFluidSimulation.h"
// #include "src/sim-fluid-ethanjli/fluidsystem.h"
#include "TerrainGen/TerrainModel.h"

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

class VoxelGrid : public TerrainModel {
public:
    VoxelGrid();
    VoxelGrid(Heightmap& grid);
    VoxelGrid(int nx, int ny, int nz, float noise_shifting = 0.0);
    ~VoxelGrid();
    void from2DGrid(Heightmap grid, Vector3 subsectionStart = Vector3(), Vector3 subsectionEnd = Vector3(), float scaleFactor = 1.f);
    void fromLayerBased(LayerBasedGrid layerBased, int fixedHeight = -1);
    void fromImplicit(ImplicitPatch* implicitTerrain, int fixedHeight = -1);
    VoxelGrid *fromCachedData();
    void setVoxelValues(const Matrix3<float>& values);

    void initMap();

    void makeItFall(float erosionStrength = 0.0);
    void letGravityMakeSandFall(bool remesh = true);
    void letGravityMakeSandFallWithFlow(bool remesh = true);
    Matrix3<float> shareSandWithNeighbors(); // Doesn't affect the grid directly, but changes are returned to be applied after
    void applyModification(Matrix3<float> modifications, Vector3 anchor = Vector3());
    void add2DHeightModification(Matrix3<float> heightmapModifier, float factor = 1.f, Vector3 anchor = Vector3());
    bool undo();
    bool redo();
    size_t getCurrentHistoryIndex() const;

//    int numberOfChunksX() { return std::ceil(this->getSizeX() / (float)this->chunkSize); }
//    int numberOfChunksY() { return std::ceil(this->getSizeY() / (float)this->chunkSize); }

//    std::vector<Matrix3<float>> tempData;
    float getHeight(float x, float y);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    virtual bool checkIsInGround(Vector3 position);

//    void remeshAll();
    void limitVoxelValues(float limitedValue);

    std::string toString();
    std::string toShortString();

    void smoothVoxels();


    void computeVoxelGroups();
    Matrix3<float> getVoxelValues();
    Matrix3<float> getVoxelized(Vector3 dimensions = Vector3(false), Vector3 scale = Vector3(1.f, 1.f, 1.f));

    virtual Mesh getGeometry(Vector3 reducedResolution = Vector3(false));

//    std::tuple<int, int, int, int> getChunksAndVoxelIndices(Vector3 pos);
//    std::tuple<int, int, int, int> getChunksAndVoxelIndices(float x, float y, float z);
    float getVoxelValue(Vector3 pos);
    float getVoxelValue(float x, float y, float z);
    void setVoxelValue(Vector3 pos, float newVal);
    void setVoxelValue(float x, float y, float z, float newVal);
//    float getOriginalVoxelValue(Vector3 pos);
//    float getOriginalVoxelValue(float x, float y, float z);
//    Matrix3<Vector3> getFlowfield();
//    Matrix3<Vector3> getFlowfield(size_t flowIndex);
//    Vector3 getFlowfield(Vector3 pos);
//    Vector3 getFlowfield(float x, float y, float z);
//    void setFlowfield(Vector3 pos, Vector3 newVal);
//    void setFlowfield(float x, float y, float z, Vector3 newVal);
//    int getVoxelGroup(Vector3 pos);
//    int getVoxelGroup(float x, float y, float z);
//    void setVoxelGroup(Vector3 pos, int newVal);
//    void setVoxelGroup(float x, float y, float z, int newVal);
//    bool getVoxelIsOnGround(Vector3 pos);
//    bool getVoxelIsOnGround(float x, float y, float z);
//    void setVoxelIsOnGround(Vector3 pos, bool newVal);
//    void setVoxelIsOnGround(float x, float y, float z, bool newVal);

//    void computeFlowfield();
//    void computeMultipleFlowfields(int steps = 30, ImplicitNaryOperator* primitives = nullptr);

//    void affectFlowfieldAround(Vector3 pos, Vector3 newVal, int kernelSize = 3);
//    void affectFlowfieldAround(float x, float y, float z, Vector3 newVal, int kernelSize = 3);
//    void affectFlowfieldAround(Vector3 pos, float alphaEffect, int kernelSize = 3);
//    void affectFlowfieldAround(float x, float y, float z, float alphaEffect, int kernelSize = 3);

//    int getMaxLoD();

    void saveState();

    void saveMap(std::string filename);
    void retrieveMap(std::string filename);

    Vector3 getFirstIntersectingVoxel(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));
    Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));

    float getSizeX() { return _cachedVoxelValues.sizeX; }
    float getSizeY() { return _cachedVoxelValues.sizeY; }
    float getSizeZ() { return _cachedVoxelValues.sizeZ; }

    Vector3 fluidSimRescale;

//    int getChunkSize() const { return this->chunkSize; }

    void saveHeightmap(std::string heightmap_filename);

//protected:
//    std::vector<std::shared_ptr<VoxelChunk>> chunks;
//    float noise_shifting;
//    int chunkSize = 20;
    FastNoiseLite noise;
    NoiseMinMax noiseMinMax;

    bool _smoothingNeeded = false; // Just used when we come from a 2D grid.


    float getNoiseValue(int x, int y, int z, float noise_shift = 0.f);
//    int _cachedHistoryIndex = -1;
    Matrix3<float> _cachedVoxelValues;

    std::vector<Matrix3<float>> voxelsValuesStack;
    std::vector<Vector3> voxelsValuesAnchorStack;
//    void computeFlowfield(FluidSimType type);
//    void computeMultipleFlowfields(FluidSimType type, int steps = 30, ImplicitNaryOperator *primitives = nullptr);
};

#endif // VOXELGRID_H
