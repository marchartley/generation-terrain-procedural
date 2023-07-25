#ifndef WARPEDFLUIDSIMULATION_H
#define WARPEDFLUIDSIMULATION_H

#include "FluidSimulation.h"

class WarpedFluidSimulation : public FluidSimulation
{
public:
    WarpedFluidSimulation();
    WarpedFluidSimulation(int x, int y, int z);

    virtual void step();
    virtual void handleCollisions();

    virtual Matrix3<Vector3> getVelocities(int newSizeX, int newSizeY, int newSizeZ);
    virtual Vector3 getVelocity(int x, int y, int z);
    virtual void addVelocity(int x, int y, int z, const Vector3& amount);

    Vector3 mainDirection = Vector3(1, 0, 0);

    Matrix3<Vector3> velocities;
};

#endif // WARPEDFLUIDSIMULATION_H
