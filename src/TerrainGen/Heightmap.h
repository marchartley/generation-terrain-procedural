#ifndef HEIGHTMAP_H
#define HEIGHTMAP_H

#include "DataStructure/Vector3.h"

class Heightmap;

#include "TerrainGen/VoxelGrid.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Matrix3.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/TerrainModel.h"
#include "TerrainGen/ImplicitPatch.h"

class Heightmap : public TerrainModel {
public:
    Heightmap();
    Heightmap(int nx, int ny, float heightFactor);
    Heightmap(std::string heightmap_filename, int nx = 30, int ny = 30, float heightFactor = 10);

    Matrix3<float> getHeights() { return this->heights; }
    float getHeight(float x, float y) { return this->heights.at(x, y); }
    float getHeight(Vector3 pos) { return this->getHeight(pos.x, pos.y); }

    float getMaxHeight();
    float getHeightFactor() { return this->heightFactor; }
    float getSizeX() { return this->heights.sizeX; }
    float getSizeY() { return this->heights.sizeY; }
    float getSizeZ() { return this->getMaxHeight(); }

    virtual bool checkIsInGround(Vector3 position);

    /// Erosion functions (should be in another class, I guess...)
    std::vector<std::vector<Vector3>> hydraulicErosion(int numIterations = 1000,
                                                       int erosionRadius = 10,
                                                       int maxDropletLifetime = 30,
                                                       float erodeSpeed = .3f,
                                                       float depositSpeed = .3f,
                                                       float evaporateSpeed = .01f,
                                                       float gravity = 4,
                                                       float inertia = .05f,
                                                       float sedimentCapacityFactor = 1,
                                                       bool applyDeposit = true);
    void thermalErosion(float erosionCoef = .1f, float minSlope = .01f);
    std::vector<std::vector<Vector3>> windErosion(int numberOfParticles = 100,
                                                  Vector3 windDirection = Vector3(2.f, 0.f, 0.f),
                                                  float bedrocksProportionInGround = .0f,
                                                  float suspension = .002f,
                                                  float abrasion = .01f,
                                                  float roughness = .005f,
                                                  float settling = .05f,
                                                  float scale = 40.f,
                                                  float dt = .1f);

    void raise(Matrix3<float> elevation);

    Heightmap& fromVoxelGrid(VoxelGrid& voxelGrid);
    Heightmap& fromLayerGrid(LayerBasedGrid& layerGrid);
    Heightmap& fromImplicit(ImplicitPatch *implicitTerrain);
    Matrix3<float> getVoxelized(Vector3 dimensions = Vector3(false), Vector3 scale = Vector3(1.f, 1.f, 1.f));

    void randomFaultTerrainGeneration(int numberOfFaults = 50, int maxNumberOfSubpointsInFaults = 2, float faultHeight = 1.f);

    void saveMap(std::string filename) { return this->saveHeightmap(filename); }
    void retrieveMap(std::string filename) { this->loadFromHeightmap(filename, getSizeX(), getSizeY(), heightFactor); }
    Heightmap& loadFromHeightmap(std::string heightmap_filename, int nx = -1, int ny = -1, float heightFactor = -1);
    void saveHeightmap(std::string heightmap_filename);

    Vector3 getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos = Vector3(false), Vector3 maxPos = Vector3(false));
    Vector3 findSurfaceBetween(Vector3 start, Vector3 end);

    Mesh getGeometry(Vector3 dimensions = Vector3(false));

    void initMap() {};

    bool undo() { return false; };
    bool redo() { return false; };

    std::string toString() { return this->toShortString(); };
    std::string toShortString() { return "Grid: " + std::to_string(getDimensions().x) + "x" + std::to_string(getDimensions().y); };

    Matrix3<std::vector<int>>& getBiomeIndices() { return this->biomeIndices; }

    Matrix3<Vector3> getNormals();

//protected:
    Matrix3<float> heights;
//    float maxHeight;
    float heightFactor = 1.f;
    Matrix3<std::vector<int>> biomeIndices;
};

#endif // HEIGHTMAP_H
