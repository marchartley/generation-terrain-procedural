#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Octree.h"

class FluidSimulation
{
public:
    FluidSimulation() {}
    FluidSimulation(int sizeX, int sizeY, int sizeZ) : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), obstacleGrid(Matrix3<int>(sizeX, sizeY, sizeZ)), obstacleGradient(Matrix3<Vector3>(sizeX, sizeY, sizeZ))
    {}
    virtual ~FluidSimulation() {
    }


    virtual void step() = 0;
    virtual void handleCollisions() = 0;

    virtual Matrix3<Vector3> getVelocities(int newSizeX, int newSizeY, int newSizeZ) = 0;
    virtual void addVelocity(int x, int y, int z, Vector3 amount) = 0;

    void setObstacles(const std::vector<std::vector<Vector3>>& triangles);
    void setObstacles(Matrix3<float> obstacle);
    void addObstacles(const std::vector<std::vector<Vector3>>& triangles);
    void addObstacles(Matrix3<float> obstacle);


    void simulate() {
        step();
    }


    int sizeX, sizeY, sizeZ;
    std::vector<std::vector<Vector3>> triangles;
    Matrix3<int> obstacleGrid;
    Matrix3<Vector3> obstacleGradient;
    Octree* obstacleTrianglesOctree;
};

#endif // FLUIDSIMULATION_H
