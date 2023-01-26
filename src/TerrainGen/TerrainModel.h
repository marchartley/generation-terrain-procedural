#ifndef TERRAINMODEL_H
#define TERRAINMODEL_H

#include "DataStructure/Vector3.h"
#include "Graphics/Mesh.h"

class TerrainModel : public std::enable_shared_from_this<TerrainModel>
{
public:
    TerrainModel();


    virtual void saveMap(std::string filename) = 0;
    virtual void retrieveMap(std::string filename) = 0;

    virtual Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false)) = 0;

    virtual Mesh getGeometry() = 0;

    virtual std::string toString() = 0;
    virtual std::string toShortString();

    virtual int getHeight(int x, int y) = 0;

    bool contains(Vector3 v);
    bool contains(float x, float y, float z);

    virtual bool undo() = 0;
    virtual bool redo() = 0;
    size_t getCurrentHistoryIndex() const;

    int getSizeX() { return this->getDimensions().x; }
    int getSizeY() { return this->getDimensions().y; }
    int getSizeZ() { return this->getDimensions().z; }
    Vector3 getDimensions() { return this->dimensions; }

    virtual void initMap() = 0;

//    int sizeX, sizeY, sizeZ;
    int _cachedHistoryIndex = -1;
    Vector3 dimensions;
    Vector3 scaling = Vector3(1.f, 1.f, 1.f);
    Vector3 translation;
};

#endif // TERRAINMODEL_H
