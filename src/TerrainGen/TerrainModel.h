#ifndef TERRAINMODEL_H
#define TERRAINMODEL_H

class TerrainModel;

#include "DataStructure/Vector3.h"
#include "Graphics/Mesh.h"
#include "TerrainGen/GlobalTerrainProperties.h"

class TerrainModel
{
public:
    TerrainModel();
    virtual ~TerrainModel();

    virtual void initMap() = 0;

    virtual bool undo() = 0;
    virtual bool redo() = 0;

    virtual void saveMap(std::string filename) = 0;
    virtual void retrieveMap(std::string filename) = 0;
    virtual Mesh getGeometry(const Vector3& dimensions = Vector3(false)) = 0;

    virtual Vector3 getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos = Vector3(false), const Vector3& maxPos = Vector3(false)) = 0;

    virtual std::string toString() = 0;
    virtual std::string toShortString() = 0;

    virtual float getHeight(float x, float y) = 0;
    virtual float getHeight(const Vector3& pos);

    virtual bool contains(const Vector3& v);
    virtual bool contains(float x, float y, float z);

    virtual size_t getCurrentHistoryIndex() const;

    virtual float getSizeX() const = 0;
    virtual float getSizeY() const = 0;
    virtual float getSizeZ() const = 0;
    virtual Vector3 getDimensions() const { return Vector3(getSizeX(), getSizeY(), getSizeZ()); }

    Vector3 getTerrainPos(const Vector3& pos) { return (pos * scaling) - translation; }
    Vector3 getWorldPos(const Vector3& pos) { return (pos + translation) / scaling; }
    void setScaling(const Vector3& newScale) { this->scaling = newScale; }
    void setScaling(float newScale) { this->scaling = Vector3(newScale, newScale, newScale); }
    void setTranslation(const Vector3& newTranslation) { this->translation = newTranslation; }

    virtual bool checkIsInGround(const Vector3& position) = 0;

    virtual GridF getVoxelized(const Vector3& dimensions = Vector3(false), const Vector3& scale = Vector3(1.f, 1.f, 1.f)) = 0;

    void initFluidSim();
    void initEnvironmentalDensities();

    virtual GridV3 getNormals() = 0;

    virtual GridV3 getFlowfield(FluidSimType simu = LBM);
//    virtual GridV3 getFlowfield(size_t flowIndex);
//    virtual void computeFlowfield(FluidSimType simu = LBM);
    virtual void computeFlowfield(FluidSimType simu = LBM, int steps = 30, TerrainModel* implicit = nullptr);

    GridF& getEnvironmentalDensities();
    void updateEnvironmentalDensities(float waterLevel);


//    int sizeX, sizeY, sizeZ;
    int _cachedHistoryIndex = -1;
    int currentHistoryIndex = 0;
    int _historyIndex = -3;

//    Vector3 dimensions;
    Vector3 scaling = Vector3(1.f, 1.f, 1.f);
    Vector3 translation;

    GlobalTerrainProperties* properties;


//    StableFluidsSimulation fluidSimulation;
    Vector3 fluidSimRescale;

//    std::vector<StableFluidsSimulation> multipleFluidSimulations;
//    std::vector<Vector3> multipleSeaCurrents;
//    std::vector<GridV3> multipleFlowFields;

//    GridV3 flowField;
    GridI distanceField;
    GridF pressureField;

    Vector3 sea_current = Vector3(1.f, 0.0, 0.0);

//    float storedWaterLevel = 0.f;
};

#endif // TERRAINMODEL_H
