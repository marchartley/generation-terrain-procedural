#ifndef LAYERBASEDGRID_H
#define LAYERBASEDGRID_H

class LayerBasedGrid;
#include <vector>
#include <tuple>

#include "Graphics/Mesh.h"
#include "TerrainGen/Grid.h"
#include "TerrainGen/VoxelChunk.h"
#include "DataStructure/Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Graphics/Mesh.h"
#include "Utils/FastNoiseLit.h"
#include "Utils/ShapeCurve.h"
#include "TerrainGen/ImplicitPatch.h"

class ImplicitPatch;
class Patch2D;
class Patch3D;

class LayerBasedGrid
{
public:
    LayerBasedGrid();
    LayerBasedGrid(int nx, int ny, float nz = 1);

    void createMesh();
    void display();
    TerrainTypes getValue(Vector3 pos);
    TerrainTypes getValue(float x, float y, float z);

    void addLayer(Vector3 position, float height, TerrainTypes material);

    void reorderLayers();

    Matrix3<std::vector<std::pair<TerrainTypes, float>>> layers;
    Matrix3<std::vector<std::pair<TerrainTypes, float>>> previousState;

    int getSizeX() { return layers.sizeX; }
    int getSizeY() { return layers.sizeY; }
    float getSizeZ();
    Vector3 getDimensions() { return Vector3(getSizeX(), getSizeY(), getSizeZ()); }
//    int sizeX, sizeY;
//    float sizeZ;

    Mesh mesh;

    float getHeight(float x, float y);
    float getHeight(Vector3 pos);

    void from2DGrid(Grid grid);
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
/*
    void add(Patch2D patch, TerrainTypes material, bool applyDistanceFalloff = true, float distancePower = 1.f);
    void add(Patch3D patch, TerrainTypes material, bool applyDistanceFalloff = true, float distancePower = 1.f);
    */
    void add(ImplicitPatch* patch, TerrainTypes material, bool applyDistanceFalloff = true, float distancePower = 1.f);

    int currentHistoryIndex = 0;
    int _historyIndex = 0;

    std::pair<Matrix3<int>, Matrix3<float> > _cachedMaterialAndHeights;

    static std::map<TerrainTypes, std::pair<float, float>> materialLimits;
    static TerrainTypes materialFromDensity(float density);
    static float densityFromMaterial(TerrainTypes material);
};

#endif // LAYERBASEDGRID_H
