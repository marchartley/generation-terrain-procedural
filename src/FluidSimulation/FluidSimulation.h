#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "DataStructure/BVH.h"
#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Octree.h"
#include "DataStructure/Particle.h"
#include "DataStructure/KDTree.h"

class FluidSimulation
{
public:
    FluidSimulation();
    FluidSimulation(int sizeX, int sizeY, int sizeZ);
    virtual ~FluidSimulation() {
    }


    virtual void step() = 0;
    virtual void handleCollisions() = 0;

    Matrix3<Vector3> getVelocities(const Vector3 &dimensions);
    virtual Matrix3<Vector3> getVelocities(int newSizeX, int newSizeY, int newSizeZ) = 0;
    virtual Vector3 getVelocity(const Vector3& pos);
    virtual Vector3 getVelocity(int x, int y, int z) = 0;
    virtual void addVelocity(const Vector3& pos, const Vector3& amount);
    virtual void addVelocity(int x, int y, int z, const Vector3& amount) = 0;
    virtual void setVelocity(const Vector3& pos, const Vector3& amount);
    virtual void setVelocity(int x, int y, int z, const Vector3& amount);

    virtual void setObstacles(const std::vector<std::vector<Vector3>>& triangles);
    virtual void setObstacles(const Matrix3<float>& obstacle);
    virtual void addObstacles(const std::vector<std::vector<Vector3>>& triangles);
    virtual void addObstacles(const Matrix3<float> &obstacle);


    void simulate() {
        step();
    }


    Vector3 dimensions;
//    int sizeX, sizeY, sizeZ;
    std::vector<std::vector<Vector3>> triangles;
    Matrix3<float> obstacleGrid;
    Matrix3<Vector3> obstacleGradient;
//    Octree* obstacleTrianglesOctree;
    BVHTree obstacleTriangleTree;

    int currentStep = 0;
    int _cachedStep = -1;
    Matrix3<Vector3> _cachedVelocity;
};

#endif // FLUIDSIMULATION_H
