#include "WarpedFluidSimulation.h"

WarpedFluidSimulation::WarpedFluidSimulation()
{
}

WarpedFluidSimulation::WarpedFluidSimulation(int x, int y, int z)
    : FluidSimulation(x, y, z)
{
    velocities = GridV3(dimensions, mainDirection);
}

void WarpedFluidSimulation::step()
{
    // Nothing to do
    currentStep++;
}

void WarpedFluidSimulation::handleCollisions()
{
    // Nothing to do
}

GridV3 WarpedFluidSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ)
{
    if (_cachedStep != currentStep) {
        _cachedStep = currentStep;
        _cachedVelocity = velocities;
    }
    return _cachedVelocity.resize(newSizeX, newSizeY, newSizeZ);
}

Vector3 WarpedFluidSimulation::getVelocity(int x, int y, int z)
{
    return FluidSimulation::getVelocity(x, y, z);
}

void WarpedFluidSimulation::addVelocity(int x, int y, int z, const Vector3 &amount)
{
    velocities.at(x, y, z) += amount;
}
