#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Octree.h"

class Particle;
class KDTree;

class Particle {
public:
    Vector3 position;
    Vector3 velocity;
    Vector3 force;
    float density;
    float pressure;
    float mass;
    float smoothingRadius;
    float gasConstant;
    float restDensity;
    float viscosity;
    bool isGhost;

    int index;
};

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
    Octree* obstacleTrianglesOctree;

    int currentStep = 0;
    int _cachedStep = -1;
    Matrix3<Vector3> _cachedVelocity;
};


class KDNode {
public:
    size_t particleIndex;
    KDNode* left;
    KDNode* right;
    int axis;

    KDNode(size_t pIndex, int a);
};

class KDTree {
public:
    KDTree();
    ~KDTree();
    KDNode* root;

    KDTree(std::vector<Particle>& particles);

    KDNode* build(std::vector<Particle> particles, int depth);

    void findNeighbors(std::vector<Particle>& particles, KDNode* node, Vector3& position, float maxDistance, std::vector<size_t>& neighbors);
};

#endif // FLUIDSIMULATION_H
