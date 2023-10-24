#ifndef SPHSIMULATION_H
#define SPHSIMULATION_H

#include <vector>
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include "FluidSimulation/FluidSimulation.h"

class SPHSimulation;
//class KDNode;
//class KDTree;

class SPHSimulation : public FluidSimulation
{
public:
    SPHSimulation();
    virtual ~SPHSimulation();

    std::vector<Particle> particles;
//    Vector3 dimensions;
    KDTree* tree = nullptr;

    int nbParticles;
    float dt;
    float damping;
    float t;
    float computeTime;

    std::vector<std::vector<size_t>> precomputedNeighbors;

    void computeNeighbors() ;
    void initialize(std::vector<std::vector<Vector3>> meshBoundaries = {});
    void computeDensityAndPressure();
    void computeForces();
    void integrate();
    void relaxDensity();
    void handleCollisions();
    void step();
    GridV3 getVelocities(int newSizeX, int newSizeY, int newSizeZ);
    Vector3 getVelocity(int x, int y, int z);
    void addVelocity(int x, int y, int z, const Vector3& amount);

    std::vector<size_t> getNeighbors(Vector3& position, float distance);

};

#endif // SPHSIMULATION_H
