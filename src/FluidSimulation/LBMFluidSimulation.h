#ifndef LBMFLUIDSIMULATION_H
#define LBMFLUIDSIMULATION_H

#include "FluidSimulation.h"
#include <cmath>
#include <limits>

class LBMFluidSimulation : public FluidSimulation
{
public:
    LBMFluidSimulation(bool uses3D = true);

    Vector3 computeMacroscopicVelocity(int x, int y, int z);
    float computeEquilibrium(const Vector3& c_i, const Vector3& u, float rho, float w_i);
    void handleCollisions();
    void stream();
    void step();


    virtual GridV3 getVelocities(int newSizeX, int newSizeY, int newSizeZ);
    Vector3 getVelocity(int x, int y, int z);
    void addVelocity(int x, int y, int z, const Vector3& amount);
    void setVelocity(int x, int y, int z, const Vector3& amount);
    void setObstacles(const GridF &obstacles);
    void addObstacles(const GridF& obstacles);

    float getWi(int i);

    std::vector<Vector3> c2D = {
        // Define the D2Q9 model
        Vector3(0, 0, 0),  // Rest
        Vector3(1, 0, 0),  // East
        Vector3(-1, 0, 0), // West
        Vector3(0, 1, 0),  // North
        Vector3(0, -1, 0), // South
        Vector3(1, 1, 0),  // Northeast
        Vector3(-1, 1, 0), // Northwest
        Vector3(-1, -1, 0),// Southwest
        Vector3(1, -1, 0)  // Southeast
    };
    std::vector<Vector3> c3D = {
        // Define the D3Q19 model
        Vector3(0, 0, 0),  // Rest
        Vector3(1, 0, 0),  // 6 directions along the axes
        Vector3(-1, 0, 0),
        Vector3(0, 1, 0),
        Vector3(0, -1, 0),
        Vector3(0, 0, 1),
        Vector3(0, 0, -1),
        Vector3(1, 1, 0),  // 12 directions along the diagonals
        Vector3(-1, -1, 0),
        Vector3(1, -1, 0),
        Vector3(-1, 1, 0),
        Vector3(1, 0, 1),
        Vector3(-1, 0, -1),
        Vector3(1, 0, -1),
        Vector3(-1, 0, 1),
        Vector3(0, 1, 1),
        Vector3(0, -1, -1),
        Vector3(0, 1, -1),
        Vector3(0, -1, 1)
    };

    std::vector<Vector3> c;

    // Initialize the distribution functions
    Matrix3<GridF> f;
    Matrix3<GridF> f_next;

    bool uses3D = false;
};


#endif // LBMFLUIDSIMULATION_H
