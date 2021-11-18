#include "Voxel.h"

#include "MarchingCubes.h"

std::vector<std::set<int>> Voxel::voxelGroups;
int Voxel::currentLabelIndex = 0;
std::map<TerrainTypes, bool> Voxel::isMatter = {
    {AIR, false},
    {DIRT, true},
    {WATER, false},
    {SAND, true}
};

Voxel::Voxel(int x, int y, int z, TerrainTypes type, float blockSize, float isosurface, VoxelChunk* parent)
    : x(x), y(y), z(z), type(type), blockSize(blockSize), isosurface(isosurface), parent(parent) {

}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0, 0.0, nullptr) {

}
void Voxel::addNeighbor(Voxel* neighbor)
{
    if(neighbor == nullptr || !(bool)*neighbor || !(bool)*this)
        return;
    if (this->globalX() < neighbor->globalX()) {
        this->neighbors[RIGHT] = neighbor;
        neighbor->neighbors[LEFT] = this;
    }
    if (this->globalX() > neighbor->globalX()) {
        this->neighbors[LEFT] = neighbor;
        neighbor->neighbors[RIGHT] = this;
    }
    if (this->globalY() < neighbor->globalY()) {
        this->neighbors[BACK] = neighbor;
        neighbor->neighbors[FRONT] = this;
    }
    if (this->globalY() > neighbor->globalY()) {
        this->neighbors[FRONT] = neighbor;
        neighbor->neighbors[BACK] = this;
    }
    if (this->globalZ() < neighbor->globalZ()) {
        this->neighbors[TOP] = neighbor;
        neighbor->neighbors[BOTTOM] = this;
    }
    if (this->globalZ() > neighbor->globalZ()) {
        this->neighbors[BOTTOM] = neighbor;
        neighbor->neighbors[TOP] = this;
    }
    this->shareGroup(neighbor);
}
int Voxel::shareGroup(Voxel *v)
{
    if (this->group == -1)
        this->group = v->group;
    this->group = std::min(this->group, v->group);
    this->isOnGround = this->isOnGround || v->isOnGround;
}
void Voxel::removeNeighbor(Voxel* neighbor)
{
    if (neighbor == nullptr)
        return;
    if (this->globalX() < neighbor->globalX()) {
        neighbor->neighbors[RIGHT] = nullptr;
        this->neighbors[LEFT] = nullptr;
    }
    if (this->globalX() > neighbor->globalX()) {
        neighbor->neighbors[LEFT] = nullptr;
        this->neighbors[RIGHT] = nullptr;
    }
    if (this->globalY() < neighbor->globalY()) {
        neighbor->neighbors[FRONT] = nullptr;
        this->neighbors[BACK] = nullptr;
    }
    if (this->globalY() > neighbor->globalY()) {
        neighbor->neighbors[BACK] = nullptr;
        this->neighbors[FRONT] = nullptr;
    }
    if (this->globalZ() < neighbor->globalZ()) {
        neighbor->neighbors[BOTTOM] = nullptr;
        this->neighbors[TOP] = nullptr;
    }
    if (this->globalZ() > neighbor->globalZ()) {
        neighbor->neighbors[TOP] = nullptr;
        this->neighbors[BOTTOM] = nullptr;
    }
}
void Voxel::resetNeighbors() {
    this->neighbors[TOP]    = this->parent->parent->getVoxel(this->globalPos() + Vector3( 0,  0,  1) + Vector3(.5, .5, .5));
    this->neighbors[BOTTOM] = this->parent->parent->getVoxel(this->globalPos() + Vector3( 0,  0, -1) + Vector3(.5, .5, .5));
    this->neighbors[RIGHT]  = this->parent->parent->getVoxel(this->globalPos() + Vector3( 1,  0,  0) + Vector3(.5, .5, .5));
    this->neighbors[LEFT]   = this->parent->parent->getVoxel(this->globalPos() + Vector3(-1,  0,  0) + Vector3(.5, .5, .5));
    this->neighbors[FRONT]  = this->parent->parent->getVoxel(this->globalPos() + Vector3( 0, -1,  0) + Vector3(.5, .5, .5));
    this->neighbors[BACK]   = this->parent->parent->getVoxel(this->globalPos() + Vector3( 0,  1,  0) + Vector3(.5, .5, .5));
//    for(std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it = this->neighbors.begin(); it != this->neighbors.end(); it++)
//        removeNeighbor(it->second);
}

float Voxel::getIsosurface() {
    return this->isosurface + this->manual_isosurface;
}
std::ostream& operator<<(std::ostream& io, const Voxel& v)
{
    io << "Voxel (" << v.x << ", " << v.y << ", " << v.z << ")";
    return io;
}
std::ostream& operator<<(std::ostream& io, Voxel* v)
{
    io << "Voxel (" << v->x << ", " << v->y << ", " << v->z << ")";
    return io;
}

bool Voxel::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool Voxel::contains(float x, float y, float z) {
    return (this->globalX() <= x && x < this->globalX() + 1 && this->globalY() <= y && y < this->globalY() + 1 && this->globalZ() <= z && z < this->globalZ() + 1);
}

float Voxel::globalX()  {
    return this->x + this->parent->x;
}
float Voxel::globalY()  {
    return this->y + this->parent->y;
}
float Voxel::globalZ()  {
    return this->z;
}


std::vector<Vector3> Voxel::getMeshVertices()
{
    std::vector<Vector3> vertex;
    if ((bool)*this) {
        // BOTTOM
        if (!this->neighbors[BOTTOM] || !(bool)*this->neighbors[BOTTOM]) {
            vertex.push_back(Vector3(0, 0, 0));
            vertex.push_back(Vector3(1, 0, 0));
            vertex.push_back(Vector3(1, 1, 0));
            vertex.push_back(Vector3(0, 1, 0));
        }

        // RIGHT
        if (!this->neighbors[RIGHT] || !(bool)*this->neighbors[RIGHT]) {
            vertex.push_back(Vector3(1, 0, 0));
            vertex.push_back(Vector3(1, 1, 0));
            vertex.push_back(Vector3(1, 1, 1));
            vertex.push_back(Vector3(1, 0, 1));
        }

        // TOP
        if (!this->neighbors[TOP] || !(bool)*this->neighbors[TOP]) {
            vertex.push_back(Vector3(1, 0, 1));
            vertex.push_back(Vector3(1, 1, 1));
            vertex.push_back(Vector3(0, 1, 1));
            vertex.push_back(Vector3(0, 0, 1));
        }

        // LEFT
        if (!this->neighbors[LEFT] || !(bool)*this->neighbors[LEFT]) {
            vertex.push_back(Vector3(0, 0, 1));
            vertex.push_back(Vector3(0, 1, 1));
            vertex.push_back(Vector3(0, 1, 0));
            vertex.push_back(Vector3(0, 0, 0));
        }

        // FRONT
        if (!this->neighbors[FRONT] || !(bool)*this->neighbors[FRONT]) {
            vertex.push_back(Vector3(0, 0, 0));
            vertex.push_back(Vector3(0, 0, 1));
            vertex.push_back(Vector3(1, 0, 1));
            vertex.push_back(Vector3(1, 0, 0));
        }

        // BACK
        if (!this->neighbors[BACK] || !(bool)*this->neighbors[BACK]) {
            vertex.push_back(Vector3(0, 1, 0));
            vertex.push_back(Vector3(1, 1, 0));
            vertex.push_back(Vector3(1, 1, 1));
            vertex.push_back(Vector3(0, 1, 1));
        }
    }
    std::vector<Vector3> returnArray;
    for (size_t i = 0; i < vertex.size(); i+=4)
    {
        returnArray.push_back(globalPos() + vertex[i + 0]);
        returnArray.push_back(globalPos() + vertex[i + 1]);
        returnArray.push_back(globalPos() + vertex[i + 2]);
        returnArray.push_back(globalPos() + vertex[i + 2]);
        returnArray.push_back(globalPos() + vertex[i + 3]);
        returnArray.push_back(globalPos() + vertex[i + 0]);
    }
    return returnArray;
}
