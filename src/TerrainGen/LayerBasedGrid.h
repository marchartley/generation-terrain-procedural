#ifndef LAYERBASEDGRID_H
#define LAYERBASEDGRID_H

class LayerBasedGrid;
#include <vector>
#include <tuple>

#include "Graphics/Mesh.h"
//#include "TerrainGen/Grid.h"
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

class LayerBasedGrid : public TerrainModel
{
public:
    LayerBasedGrid();
    LayerBasedGrid(int nx, int ny, float nz = 1);

    TerrainTypes getValue(Vector3 pos);
    TerrainTypes getValue(float x, float y, float z);

    void addLayer(Vector3 position, float height, TerrainTypes material);

    void reorderLayers();

    float getHeight(float x, float y);

    void from2DGrid(Heightmap grid);
    void fromVoxelGrid(VoxelGrid& voxelGrid);

    VoxelGrid toVoxelGrid();
    Matrix3<float> voxelize(int fixedHeight = -1, float kernelSize = 1.f);
    std::map<TerrainTypes, float> getKernel(Vector3 pos, float kernelSize);
    std::pair<TerrainTypes, float> getMaterialAndHeight(Vector3 pos);

    Vector3 getFirstIntersectingStack(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));
    Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));

    std::pair<Matrix3<int>, Matrix3<float>> getMaterialAndHeightsGrid();

    void thermalErosion();

    void cleanLayers(float minLayerHeight = 0.1f);

    void add(ImplicitPatch* patch);

    Mesh getGeometry();

    static std::map<TerrainTypes, std::pair<float, float>> materialLimits;
    static TerrainTypes materialFromDensity(float density);
    static float densityFromMaterial(TerrainTypes material);
    static float minDensityFromMaterial(TerrainTypes material);
    static float maxDensityFromMaterial(TerrainTypes material);

    static std::vector<TerrainTypes> invisibleLayers;
    static std::vector<TerrainTypes> instanciableLayers;

    std::vector<std::pair<std::map<TerrainTypes, float>, std::map<TerrainTypes, float>>> transformationRules;


    float getSizeX() { return this->layers.sizeX; }
    float getSizeY() { return this->layers.sizeY; }
    float getSizeZ();

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
    std::pair<Matrix3<int>, Matrix3<float> > _cachedMaterialAndHeights;

    Matrix3<std::vector<std::pair<TerrainTypes, float>>> layers;
    Matrix3<std::vector<std::pair<TerrainTypes, float>>> previousState;
};

#endif // LAYERBASEDGRID_H
