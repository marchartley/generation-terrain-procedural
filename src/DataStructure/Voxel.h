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

#include "TerrainGen/Grid.h"
#include "Vertex.h"
#include "TerrainGen/LayerBasedGrid.h"
#include <memory>


class Voxel : public std::enable_shared_from_this<Voxel> {
public:
    Voxel();
    Voxel(int x, int y, int z, TerrainTypes type, float blockSize, float isosurface, std::shared_ptr<VoxelChunk> parent);

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

    void addNeighbor(std::shared_ptr<Voxel> neighbor);
    void removeNeighbor(std::shared_ptr<Voxel> neighbor);
    void resetNeighbors();

    std::vector<Vector3> getMeshVertices(bool useGlobalCoords = false);

    operator bool() { return isMatter[this->getType()]; }
    TerrainTypes getType() { return (this->getIsosurface() > 0.0 ? (this->getIsosurface() > 0.5 ? DIRT : SAND) : (this->getIsosurface() < -0.5 ? AIR : WATER)); }

    friend std::ostream& operator<<(std::ostream& io, const Voxel& v);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Voxel> v);

    int shareGroup(std::shared_ptr<Voxel> v);

//protected:
    int x, y, z;
    TerrainTypes type;
    float blockSize;
    std::unordered_map<VOXEL_NEIGHBOR, std::shared_ptr<Voxel>> neighbors;
    float isosurface = 0.0;
    float manual_isosurface = 0.0;
    bool isOnGround = false;
    std::shared_ptr<VoxelChunk> parent = nullptr;
    std::shared_ptr<LayerBasedGrid> lParent = nullptr;

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
    std::shared_ptr<Voxel> getNeighbor(VOXEL_NEIGHBOR dir);
    void InsertNeighborsIfNotOnGround(std::unordered_set<std::shared_ptr<Voxel>>& set);
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
                    std::shared_ptr<Voxel> v = this->parent->parent->getVoxel(x, y, z);
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
