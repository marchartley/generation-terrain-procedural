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
    float max = std::numeric_limits<float>::lowest();

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
    void setVoxelValues(const GridF& values);

    void initMap();

    void makeItFall(float erosionStrength = 0.0);
    void letGravityMakeSandFall(bool remesh = true);
    void letGravityMakeSandFallWithFlow(bool remesh = true);
    GridF shareSandWithNeighbors(); // Doesn't affect the grid directly, but changes are returned to be applied after
    void applyModification(GridF modifications, const Vector3& anchor = Vector3());
    void add2DHeightModification(GridF heightmapModifier, float factor = 1.f, const Vector3& anchor = Vector3());
    bool undo();
    bool redo();
    size_t getCurrentHistoryIndex() const;

    GridF getHeights() const;
    float getHeight(float x, float y);

    bool contains(const Vector3& v);
    bool contains(float x, float y, float z);

    virtual bool checkIsInGround(const Vector3& position);

    void limitVoxelValues(float limitedValue, bool increaseHeightIfNeeded = true);

    std::string toString();
    std::string toShortString();

    void smoothVoxels();

    virtual GridV3 getNormals();

    void computeVoxelGroups();
    GridF getVoxelValues();
    GridF getVoxelized(const Vector3& dimensions = Vector3(false), const Vector3& scale = Vector3(1.f, 1.f, 1.f));

    Mesh getGeometry(const Vector3& dimensions = Vector3(false));

    float getVoxelValue(const Vector3& pos);
    float getVoxelValue(float x, float y, float z);
    void setVoxelValue(const Vector3& pos, float newVal);
    void setVoxelValue(float x, float y, float z, float newVal);

    void saveState();

    void saveMap(std::string filename);
    void retrieveMap(std::string filename);

    Vector3 getFirstIntersectingVoxel(const Vector3& origin, const Vector3& dir, const Vector3& minPos = Vector3(false), const Vector3& maxPos = Vector3(false));
    Vector3 getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos = Vector3(false), const Vector3& maxPos = Vector3(false));

    float getSizeX() const { return _cachedVoxelValues.sizeX; }
    float getSizeY() const { return _cachedVoxelValues.sizeY; }
    float getSizeZ() const { return _cachedVoxelValues.sizeZ; }

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
    GridF _cachedVoxelValues;
    GridF _cachedMaxHeights;

    std::vector<GridF> voxelsValuesStack;
    std::vector<Vector3> voxelsValuesAnchorStack;
};

#endif // VOXELGRID_H
