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

class Patch;
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

    void from2DGrid(Grid grid);
    void fromVoxelGrid(VoxelGrid& voxelGrid);

    VoxelGrid toVoxelGrid();
    Matrix3<float> voxelize(int fixedHeight = -1);

    std::pair<Matrix3<int>, Matrix3<float> > getMaterialAndHeightsGrid();

    void thermalErosion();

    void cleanLayers(float minLayerHeight = 0.1f);

    void add(Patch2D patch, TerrainTypes material, bool applyDistanceFalloff = true, float distancePower = 1.f);
    void add(Patch3D patch, TerrainTypes material, bool applyDistanceFalloff = true, float distancePower = 1.f);

    static std::map<TerrainTypes, std::pair<float, float>> materialLimits;
    static TerrainTypes materialFromDensity(float density);
    static float densityFromMaterial(TerrainTypes material);
};



class Patch
{
public:
    Patch();
    Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, float densityValue = 1.f);
    Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
//    virtual ~Patch();

    virtual float getMaxHeight(Vector3 position);
    virtual float evaluate(Vector3 position);
    float get(Vector3 position) {
        float composition = this->compositionFunction(position);
        return composition * this->densityValue;
    }

    float getAlphaValue(Vector3 pos) {
        float alpha = this->alphaDistanceFunction(pos);
        if (alpha > 1.f) {
            int a = 0;
        }
        return alpha;
    }

    Vector3 position;
    Vector3 boundMin;
    Vector3 boundMax;
    std::function<float(Vector3)> evalFunction;
    float densityValue;
    std::function<float(Vector3)> compositionFunction;
    Matrix3<float> distanceTransform;
    std::function<float(Vector3)> alphaDistanceFunction;

    Matrix3<float> _cachedHeights;
};

class Patch2D: public Patch
{
public:
    Patch2D();
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, float densityValue = 1.f);
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
//    ~Patch2D() {}

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);
};

class Patch3D: public Patch
{
public:
    Patch3D();
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, float densityValue = 1.f);
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, std::function<float(Vector3)> alphaDistanceFunction, float densityValue = 1.f);
    Patch3D(const Patch3D& copy);
//    ~Patch3D() {}

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);

    static Patch3D stack(Patch3D *P1, Patch3D *P2);
    static Patch3D replace(Patch3D* P1, Patch3D* P2);
    static Patch3D blend(Patch3D* P1, Patch3D* P2);
};


#endif // LAYERBASEDGRID_H
