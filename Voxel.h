#ifndef VOXEL_H
#define VOXEL_H

#include <map>
#include <set>

class Voxel;

enum TerrainTypes {
    AIR = 0,
    DIRT = 1,
    WATER = 2,
    SAND = 3
};

enum VOXEL_NEIGHBOR {
    TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5
};

#include "Grid.h"
#include "Vertex.h"


class Voxel {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize, float isosurface);

    void display(bool apply_marching_cubes = false, bool display_vertices = false);

    int getX() { return this->x; }
    int getY() { return this->y; }
    int getZ() { return this->z; }

    float getIsosurface();

    float globalX();
    float globalY();
    float globalZ();
    Vector3 globalPos() { return Vector3(globalX(), globalY(), globalZ()); }

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    void addNeighbor(Voxel* neighbor);
    void removeNeighbor(Voxel* neighbor);
    void resetNeighbors();

    std::vector<Vector3> getMeshVertices();

    operator bool() { return isMatter[this->getType()]; }
    TerrainTypes getType() { return (this->getIsosurface() > 0.0 ? DIRT : AIR); }

    friend std::ostream& operator<<(std::ostream& io, const Voxel& v);
    friend std::ostream& operator<<(std::ostream& io, Voxel* v);

    int shareGroup(Voxel* v);

//protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
    std::map<VOXEL_NEIGHBOR, Voxel*> neighbors;
    float isosurface = 0.0;
    float manual_isosurface = 0.0;
    bool isOnGround = false;
    VoxelChunk* parent;

    static std::map<TerrainTypes, bool> isMatter;

    static std::vector<std::set<int>> voxelGroups;
    static int currentLabelIndex;
    int group = -1;

};
#endif // VOXEL_H
