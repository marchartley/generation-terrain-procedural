#ifndef VOXEL_H
#define VOXEL_H

class Voxel;

enum TerrainTypes {
    AIR = 0,
    DIRT = 1,
    WATER = 2
};

enum VOXEL_NEIGHBOR {
    TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5
};

#include <map>
#include "Grid.h"
#include "Vertex.h"


class Voxel {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize);

    void display(bool apply_marching_cubes = false, bool display_vertices = false);

    int getX() { return this->x; }
    int getY() { return this->y; }
    int getZ() { return this->z; }

    float getIsosurface();

    void addNeighbor(Voxel& neighbor);
    void removeNeighbor(Voxel& neighbor);
    std::map<VOXEL_NEIGHBOR, bool> has_neighbors;

    operator bool() { return this->type != TerrainTypes::AIR; }

    friend std::ostream& operator<<(std::ostream& io, const Voxel& v);
    friend std::ostream& operator<<(std::ostream& io, Voxel* v);

    Vertex vertices[8];
    float* isosurfaces[8];
//protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
    std::map<VOXEL_NEIGHBOR, Voxel&> neighbors;
    float isosurface;
    float manual_isosurface = 0.0;
    VoxelChunk* parent;

};
#endif // VOXEL_H
