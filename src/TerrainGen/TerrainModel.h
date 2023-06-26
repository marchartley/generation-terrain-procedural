#ifndef TERRAINMODEL_H
#define TERRAINMODEL_H

#include "DataStructure/Vector3.h"
#include "Graphics/Mesh.h"

class TerrainModel// : public std::enable_shared_from_this<TerrainModel>
{
public:
    TerrainModel();
    virtual ~TerrainModel();

    virtual void initMap() = 0;

    virtual bool undo() = 0;
    virtual bool redo() = 0;

    virtual void saveMap(std::string filename) = 0;
    virtual void retrieveMap(std::string filename) = 0;
    virtual Mesh getGeometry() = 0;

    virtual Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false)) = 0;

    virtual std::string toString() = 0;
    virtual std::string toShortString() = 0;

    virtual float getHeight(float x, float y) = 0;
    virtual float getHeight(Vector3 pos);

    virtual bool contains(Vector3 v);
    virtual bool contains(float x, float y, float z);

    virtual size_t getCurrentHistoryIndex() const;

    virtual float getSizeX() = 0;
    virtual float getSizeY() = 0;
    virtual float getSizeZ() = 0;
    virtual Vector3 getDimensions() { return Vector3(getSizeX(), getSizeY(), getSizeZ()); }

    Vector3 getTerrainPos(Vector3 pos) { return (pos * scaling) - translation; }
    Vector3 getWorldPos(Vector3 pos) { return (pos + translation) / scaling; }
    void setScaling(Vector3 newScale) { this->scaling = newScale; }
    void setScaling(float newScale) { this->scaling = Vector3(newScale, newScale, newScale); }
    void setTranslation(Vector3 newTranslation) { this->translation = newTranslation; }

    virtual bool checkIsInGround(Vector3 position) = 0;


//    int sizeX, sizeY, sizeZ;
    int _cachedHistoryIndex = -1;
    int currentHistoryIndex = 0;
    int _historyIndex = -3;

//    Vector3 dimensions;
    Vector3 scaling = Vector3(1.f, 1.f, 1.f);
    Vector3 translation;
};

#endif // TERRAINMODEL_H
