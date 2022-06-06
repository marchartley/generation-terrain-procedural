#ifndef GRID_H
#define GRID_H

#include "DataStructure/Vector3.h"

class Grid;

#include "TerrainGen/VoxelGrid.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Matrix3.h"

class Grid {
public:
    Grid();
    Grid(int nx, int ny, float max_height, float tileSize = 0.1);
    Grid(std::string heightmap_filename, int nx = -1, int ny = -1, float max_height = -1, float tileSize = 0.1);

    void display(bool displayNormals = false);

    int getSizeX() {return heights.sizeX;}
    int getSizeY() {return heights.sizeY;}

    Matrix3<float> getHeights() { return this->heights; }
    float getHeight(int x, int y) { return this->heights.at(x, y); }
    float getMaxHeight();
    float getTileSize() { return this->tileSize; }
    Vector3 getNormal(int x, int y) { return this->normals.at(x, y); }

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

    void createMesh();

    void fromVoxelGrid(VoxelGrid& voxelGrid);

    void randomFaultTerrainGeneration(int numberOfFaults = 50, int maxNumberOfSubpointsInFaults = 2, float faultHeight = 1.f);

    void loadFromHeightmap(std::string heightmap_filename, int nx = -1, int ny = -1, float max_height = -1, float tileSize = 0.1);
    void saveHeightmap(std::string heightmap_filename);

//protected:
    void computeNormals();
    Matrix3<float> heights;
    Matrix3<Vector3> normals;
    float maxHeight;
    float tileSize;
    std::vector<float> vertexArrayFloat;
    Mesh mesh;
};

#endif // GRID_H
