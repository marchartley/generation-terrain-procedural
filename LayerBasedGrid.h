#ifndef LAYERBASEDGRID_H
#define LAYERBASEDGRID_H

class LayerBasedGrid;
#include <vector>
#include <tuple>

#include "Voxel.h"
#include "Mesh.h"
#include "Grid.h"
#include "VoxelChunk.h"
#include "Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Mesh.h"
#include "FastNoiseLit.h"

class LayerBasedGrid
{
public:
    LayerBasedGrid();
    LayerBasedGrid(int nx, int ny, float nz);

    void createMesh();
    void display();

    std::vector<std::vector<std::vector<std::tuple<TerrainTypes, float>>>> layers;
    int sizeX, sizeY;
    int height;
    float sizeZ;

    Mesh mesh;
    std::vector<std::vector<std::vector<Voxel*>>> voxels;


    template<class R>
    void applyToVoxels(std::function<R(Voxel*)> func) {
        for(int v_x = 0; v_x < sizeX; v_x++) {
            for(int v_y = 0; v_y < sizeY; v_y++) {
                for(int h = 0; h < sizeZ; h++) {
                    func(this->voxels[v_x][v_y][h]);
                }
            }
        }
    }
    template <class F>
    void applyToVoxels(F func) {
        this->applyToVoxels(std::function<void(Voxel*)>(func));
    }
};
/*
class LayerBasedGrid;

#include "Grid.h"
#include "VoxelChunk.h"
#include "Voxel.h"
#include <vector>
#include <map>
#include <tuple>
#include "Mesh.h"
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
    Voxel* getVoxel(Vector3 pos);
    Voxel* getVoxel(float x, float y, float z);

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
