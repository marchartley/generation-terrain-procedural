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

    int getSizeX() { return layers.sizeX; }
    int getSizeY() { return layers.sizeY; }
    float getSizeZ();
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
    virtual ~Patch();

    virtual float getMaxHeight(Vector3 position) = 0;
    virtual float evaluate(Vector3 position) = 0;
    float get(Vector3 position) {
        return compositionFunction(position) * this->densityValue;
    }

    Vector3 position;
    Vector3 boundMin;
    Vector3 boundMax;
    std::function<float(Vector3)> evalFunction;
    float densityValue;
    std::function<float(Vector3)> compositionFunction;
};

class Patch2D: public Patch
{
public:
    Patch2D();
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, float densityValue = 1.f);
    Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> heightFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);

//    ShapeCurve shape;

};

class Patch3D: public Patch
{
public:
    Patch3D();
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, float densityValue = 1.f);
    Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float(Vector3)> evalFunction, std::function<float(Vector3)> compositionFunction, float densityValue = 1.f);

    float getMaxHeight(Vector3 position);
    float evaluate(Vector3 position);

    static Patch3D stack(Patch& P1, Patch& P2);
    static Patch3D replace(Patch& P1, Patch& P2);
    static Patch3D blend(Patch& P1, Patch& P2);
};


/*
class LayerBasedGrid;

#include "TerrainGen/Grid.h"
#include "VoxelChunk.h"
#include "Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Graphics/Mesh.h"
#include "FastNoiseLit.h"

class LayerBasedGrid {
public:
    LayerBasedGrid();
    LayerBasedGrid(Grid& grid);
    LayerBasedGrid(int nx, int ny, int nz, float blockSize = 1.0, float noise_shifting = 0.0);
    ~LayerBasedGrid();
    void from2DGrid(Grid grid);

    void initMap();

    void display();

    void createMesh();

    void makeItFall(int groupId = -1);

    int getSizeX() { return this->sizeX; }
    int getSizeY() { return this->sizeY; }
    int getSizeZ() { return this->sizeZ; }
    float getBlockSize() { return this->blockSize; }

    std::vector<Voxel> getVoxels() { return this->voxels; }
    int getHeight(int x, int y);

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);
    std::shared_ptr<Voxel> getVoxel(Vector3 pos);
    std::shared_ptr<Voxel> getVoxel(float x, float y, float z);

    void remeshAll();

    std::string toString();
    std::string toShortString();

//protected:
    int sizeX, sizeY, sizeZ;
    std::vector<Voxel> voxels;
    float blockSize;
    std::vector<VoxelChunk*> chunks;
    float noise_shifting;

    int chunkSize = 40;
    bool displayWithMarchingCubes = false;
    FastNoiseLite noise;
    NoiseMinMax noiseMinMax;
};
*/
#endif // LAYERBASEDGRID_H
