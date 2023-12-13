#ifndef LAYERBASEDGRID_H
#define LAYERBASEDGRID_H

class LayerBasedGrid;
#include <vector>
#include <tuple>

#include "Graphics/Mesh.h"
//#include "TerrainGen/Heightmap.h"
//#include "TerrainGen/VoxelChunk.h"
#include "DataStructure/Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Graphics/Mesh.h"
//#include "Utils/FastNoiseLit.h"
//#include "Utils/ShapeCurve.h"
#include "TerrainGen/ImplicitPatch.h"
#include "TerrainGen/TerrainModel.h"

#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"

class LayerBasedGrid : public TerrainModel
{
public:
    LayerBasedGrid();
    LayerBasedGrid(int nx, int ny, float nz = 1);

    TerrainTypes getValue(const Vector3& pos);
    TerrainTypes getValue(float x, float y, float z);

    void addLayer(const Vector3& position, float height, TerrainTypes material);

    void reorderLayers();

    float getHeight(float x, float y);

    void from2DGrid(Heightmap grid);
    void fromVoxelGrid(VoxelGrid& voxelGrid);
    void fromImplicit(ImplicitPatch* implicitTerrain);

    VoxelGrid toVoxelGrid();
    GridF voxelize(int fixedHeight = -1, float kernelSize = 1.f);
    GridF getVoxelized(const Vector3& dimensions = Vector3(false), const Vector3& scale = Vector3(1.f, 1.f, 1.f));
    std::map<TerrainTypes, float> getKernel(const Vector3& pos, float kernelSize);
    std::pair<TerrainTypes, float> getMaterialAndHeight(const Vector3& pos);

    Vector3 getFirstIntersectingStack(const Vector3& origin, const Vector3& dir, const Vector3& minPos = Vector3(false), const Vector3& maxPos = Vector3(false));
    Vector3 getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos = Vector3(false), const Vector3& maxPos = Vector3(false));

    std::pair<GridI, GridF> getMaterialAndHeightsGrid();

    virtual GridV3 getNormals();

    void thermalErosion();

    void cleanLayer(int x, int y, float minLayerHeight = 0.f);
    void cleanLayers(float minLayerHeight = 0.f);

    LayerBasedGrid *transformLayer(int x, int y, float startZ, float endZ, TerrainTypes material);

    void add(ImplicitPatch* patch);

    Mesh getGeometry(const Vector3& dimensions = Vector3(false));

    static std::map<TerrainTypes, std::pair<float, float>> materialLimits;
    static TerrainTypes materialFromDensity(float density);
    static float densityFromMaterial(TerrainTypes material);
    static float minDensityFromMaterial(TerrainTypes material);
    static float maxDensityFromMaterial(TerrainTypes material);

    static std::set<TerrainTypes> invisibleLayers;
    static std::set<TerrainTypes> instanciableLayers;

    std::vector<std::pair<std::map<TerrainTypes, float>, std::map<TerrainTypes, float>>> transformationRules;

    virtual bool checkIsInGround(const Vector3& position);

    float getSizeX() const { return this->layers.sizeX; }
    float getSizeY() const { return this->layers.sizeY; }
    float getSizeZ() const;

    void initMap() {};

    bool undo() { return false; };
    bool redo() { return false; };

    void saveMap(std::string filename) {};
    void retrieveMap(std::string filename) {};

    std::string toString() { return this->toShortString(); };
    std::string toShortString() { return "Layered terrain : " + std::to_string(this->getDimensions().x) + "x" + std::to_string(this->getDimensions().y) + "x" + std::to_string(this->getDimensions().z); };

    void reset() { this->layers = previousState; }
    Matrix3<std::vector<std::pair<TerrainTypes, float>>>& getLayers() { return this->layers; }
protected:
    std::pair<GridI, GridF > _cachedMaterialAndHeights;

    Matrix3<std::vector<std::pair<TerrainTypes, float>>> layers;
    Matrix3<std::vector<std::pair<TerrainTypes, float>>> previousState;
};

#endif // LAYERBASEDGRID_H
