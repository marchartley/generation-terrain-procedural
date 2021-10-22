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


class Voxel {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize);

    void display();

    int getX() { return this->x; }
    int getY() { return this->y; }
    int getZ() { return this->z; }

    void addNeighbor(Voxel& neighbor);
    void removeNeighbor(Voxel& neighbor);
    std::map<VOXEL_NEIGHBOR, bool> has_neighbors;
protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
    std::map<VOXEL_NEIGHBOR, Voxel&> neighbors;
};
#endif // VOXEL_H
