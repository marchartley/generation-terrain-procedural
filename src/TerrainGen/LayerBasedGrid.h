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

class LayerBasedGrid
{
public:
    LayerBasedGrid();
    LayerBasedGrid(int nx, int ny, float nz);

    void createMesh();
    void display();
    TerrainTypes getValue(Vector3 pos);
    TerrainTypes getValue(float x, float y, float z);

    Matrix3<std::vector<std::pair<TerrainTypes, float>>> layers;

    int getSizeX() { return layers.sizeX; }
    int getSizeY() { return layers.sizeY; }
    float getSizeZ();
//    int sizeX, sizeY;
//    float sizeZ;

    Mesh mesh;

    void from2DGrid(Grid grid);
    void fromVoxelGrid(VoxelGrid voxelGrid);

    static std::map<TerrainTypes, std::pair<float, float>> materialLimits;
    static TerrainTypes materialFromDensity(float density);
    static float densityFromMaterial(TerrainTypes material);
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
