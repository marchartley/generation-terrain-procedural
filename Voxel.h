#ifndef VOXEL_H
#define VOXEL_H

#include <unordered_map>
#include <set>
#include <tuple>
#include <unordered_set>

class Voxel;

enum TerrainTypes {
    AIR = 0,
    DIRT = 1,
    WATER = 2,
    SAND = 3,
    LAST
};

enum VOXEL_NEIGHBOR {
    TOP = 0, BOTTOM = 1, LEFT = 2, RIGHT = 3, FRONT = 4, BACK = 5
};

#include "Grid.h"
#include "Vertex.h"
#include "LayerBasedGrid.h"


class Voxel {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize, float isosurface, VoxelChunk* parent);

    int getX() { return this->x; }
    int getY() { return this->y; }
    int getZ() { return this->z; }

    float getIsosurface();

    float globalX();
    float globalY();
    float globalZ();
    Vector3 localPos() { return Vector3(getX(), getY(), getZ()); }
    Vector3 globalPos() { return Vector3(globalX(), globalY(), globalZ()); }

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    void addNeighbor(Voxel* neighbor);
    void removeNeighbor(Voxel* neighbor);
    void resetNeighbors();

    std::vector<Vector3> getMeshVertices(bool useGlobalCoords = false);

    operator bool() { return isMatter[this->getType()]; }
    TerrainTypes getType() { return (this->getIsosurface() > 0.0 ? (this->getIsosurface() > 0.5 ? DIRT : SAND) : (this->getIsosurface() < -0.5 ? AIR : WATER)); }

    friend std::ostream& operator<<(std::ostream& io, const Voxel& v);
    friend std::ostream& operator<<(std::ostream& io, Voxel* v);

    int shareGroup(Voxel* v);

//protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
    std::unordered_map<VOXEL_NEIGHBOR, Voxel*> neighbors;
    float isosurface = 0.0;
    float manual_isosurface = 0.0;
    bool isOnGround = false;
    VoxelChunk* parent = nullptr;
    LayerBasedGrid* lParent = nullptr;

    static std::unordered_map<TerrainTypes, bool> isMatter;

    static std::vector<std::set<int>> voxelGroups;
    static int currentLabelIndex;
    int group = -1;

    template <typename F>
    inline void applyToNeighbors(F func) {
        for(auto& v : this->neighbors) {
            func(v.second);
        }
    }
    /*template <typename F>
    inline void applyToNeighbors(F func) {
        for(auto& v : this->neighbors) {
            func(v.second);
        }
    }*/
    /*
    Voxel* getNeighbor(VOXEL_NEIGHBOR dir);
    void InsertNeighborsIfNotOnGround(std::unordered_set<Voxel*>& set);
    int CountNeighbors();
    */
    /*
    enum Func {
        InsertIfNotOnGround = 1,
        DecrementIfIsAir = 2
    };
    template <typename T>
    void applyToNeighbors(Func functionToApply, T& param) {
        for (int x = this->globalX()-1; x <= this->globalX()+1; x++) {
            for (int y = this->globalY()-1; y <= this->globalY()+1; y++) {
                for (int z = this->globalZ()-1; z <= this->globalZ()+1; z++) {
                    Voxel* v = this->parent->parent->getVoxel(x, y, z);
                    if(v != this) {
                        switch(functionToApply) {
                        case InsertIfNotOnGround:
                        if(v != nullptr && (bool)*v && !v->isOnGround) {
                            param.insert(v);
                        }
                        break;
                        case DecrementIfIsAir:
                        if(v == nullptr || !(bool)*v) {
                            param --;
                        }
                        break;
                        }
                    }
                }
            }
        }
    }
    */

};
#endif // VOXEL_H
