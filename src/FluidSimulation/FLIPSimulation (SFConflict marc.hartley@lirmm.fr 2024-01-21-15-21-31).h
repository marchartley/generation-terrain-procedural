#ifndef FLIPSIMULATION_H
#define FLIPSIMULATION_H

#include <iostream>
#include <cmath>
#include <vector>
#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "FluidSimulation/FluidSimulation.h"


class FLIPSimulation : public FluidSimulation
{
public:
    Vector3 fNum;

    float cellSize;
    float fInvSpacing;

    GridF u;
    GridF v;
    GridF w;

    GridF du;
    GridF dv;
    GridF dw;

    GridF prevU;
    GridF prevV;
    GridF prevW;

    GridF p;
    GridF s;
    GridI cellType;

    KDTree particleTree;

    std::vector<Particle> particles;
    std::vector<int> numCellParticles;
    float density;
    int maxParticles;

    float pInvSpacing;
    int pNumX;
    int pNumY;
    int pNumZ;
    int pNumCells;

    std::vector<int> firstCellParticle;
    std::vector<int> cellParticleIds;

    Vector3 dimensions;

    int SOLID_CELL = 0;
    int FLUID_CELL = 1;
    int AIR_CELL = 2;

    float dt = .1f;
    float flipRatio = .9f;
    int averaging = 100;
    int maxAverageSize = 100;

    bool compensateDrift = true; // Whether to compensate for drift
    float gravityValue = 9.81 * .5f; // Acceleration due to gravity
    int numIterations = 10; // Number of iterations for each step
    float overRelaxation = 1.9f; // Over-relaxation parameter

    GridF particleDensity;
    float particleRestDensity;
    float particleRadius = 1.6f;

    FLIPSimulation();
    FLIPSimulation(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt);
    virtual ~FLIPSimulation() {}
    void init(float density, float width, float depth, float height, float spacing, float particleRadius, float maxParticles, float dt);

    void integrateParticles(double dt, const Vector3& gravity);
    void pushParticlesApart(int numIters);
//    void handleParticleCollisions(const Vector3& obstaclePos, float obstacleRadius);
    void updateParticleDensity();
    void transferVelocities(bool toGrid, float flipRatio);
    void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true);
    void step();
    void simulate();
    void respawnLostParticles();

    void storeVelocities();

    void handleCollisions();


    std::vector<GridV3> velocitiesHistory;
    GridV3 getVelocities(int newSizeX, int newSizeY, int newSizeZ);
    Vector3 getVelocity(int x, int y, int z);
    void addVelocity(int x, int y, int z, const Vector3& amount);

    void reset();

    bool useVelocityForCollisionDetection;

    bool averagingVelocities = false;

    std::vector<Particle> savedState;

    int countCells() const { return fNum.x * fNum.y * fNum.z; }

    float getU(const Vector3& p, const GridF& uGrid) const;
    float getV(const Vector3& p, const GridF& vGrid) const;
    float getW(const Vector3& p, const GridF& wGrid) const;

    Vector3 toUcoords(const Vector3& p) const { return p - Vector3(0, .5f, .5f); }
    Vector3 toVcoords(const Vector3& p) const { return p - Vector3(.5f, 0, .5f); }
    Vector3 toWcoords(const Vector3& p) const { return p - Vector3(.5f, .5f, 0); }

    Vector3 fromUcoords(const Vector3& p) const { return p + Vector3(0, .5f, .5f); }
    Vector3 fromVcoords(const Vector3& p) const { return p + Vector3(.5f, 0, .5f); }
    Vector3 fromWcoords(const Vector3& p) const { return p + Vector3(.5f, .5f, 0); }
};

#endif // FLIPSIMULATION_H
