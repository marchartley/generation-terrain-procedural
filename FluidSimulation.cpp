#include "FluidSimulation.h"

FluidSimulation::FluidSimulation()
{

}

FluidSimulation::FluidSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations)
    : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), dt(dt), diffusionAmount(diffusionAmount), viscosity(viscosity), iterations(iterations)
{
    obstacles = Matrix3<int>(sizeX, sizeY, sizeZ);
    density_old = Matrix3<float>(sizeX, sizeY, sizeZ);
    density = Matrix3<float>(sizeX, sizeY, sizeZ);

    velocity = Matrix3<Vector3>(sizeX, sizeY, sizeZ);
    velocity_old = Matrix3<Vector3>(sizeX, sizeY, sizeZ);

    divergence = Matrix3<float>(sizeX, sizeY, sizeZ);
    pressure = Matrix3<float>(sizeX, sizeY, sizeZ);
}

void FluidSimulation::setObstacles(Matrix3<int> new_obstacles)
{
    this->obstacles = new_obstacles.resize(sizeX, sizeY, sizeZ);
}

void FluidSimulation::addDensity(int x, int y, int z, float amount)
{
    this->density(x, y, z) += amount;
}

void FluidSimulation::addVelocity(int x, int y, int z, Vector3 amount)
{
    this->velocity(x, y, z) += amount;
}

Matrix3<Vector3> FluidSimulation::getVelocities(int rescaleX, int rescaleY, int rescaleZ)
{
    return this->velocity.resize(rescaleX, rescaleY, rescaleZ);
}

void FluidSimulation::step()
{
    /*
    this->diffuse(this->velocity, this->velocity_old, this->diffusionAmount, true);
    this->project();
    this->advect(this->velocity, this->velocity_old, true);

    this->project();
    this->diffuse();*/
    this->velocityStep();
    this->densityStep();
    this->set_bounds(this->velocity, true, false);
    this->set_bounds(this->density, false, true);
}

void FluidSimulation::velocityStep()
{
    swapArrays(this->velocity, this->velocity_old);
    this->diffuse(this->velocity, this->velocity_old, this->viscosity);
    this->project();
    swapArrays(this->velocity_old, this->velocity);
    this->advect(this->velocity, this->velocity_old, this->velocity_old);
    this->project();
}

void FluidSimulation::densityStep()
{
    swapArrays(this->density, this->density_old);
    this->diffuse(this->density, this->density_old, this->diffusionAmount);
    swapArrays(this->density_old, this->density);
    this->advect(this->density, this->density_old, this->velocity);
}

void FluidSimulation::project()
{
    for (int x = 1; x < this->sizeX - 1; x++) {
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                this->divergence(x, y, z) = -.5 * (
                          velocity(x+1, y  , z  ).x
                        - velocity(x-1, y  , z  ).x
                        + velocity(x  , y+1, z  ).y
                        - velocity(x  , y-1, z  ).y
                        + velocity(x  , y  , z+1).z
                        - velocity(x  , y  , z-1).z) / std::sqrt(this->sizeX * this->sizeY * this->sizeZ);
                this->pressure(x, y, z) = 0.0;
            }
        }
    }
    this->set_bounds(this->divergence);
    this->set_bounds(this->pressure);
    this->solve_linear(this->pressure, this->divergence, 1, 6, false);

    for (int x = 1; x < this->sizeX - 1; x++) {
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                this->velocity(x, y, z).x -= .5 * (pressure(x-1, y  , z  ) - pressure(x+1, y  , z  )) * this->sizeX;
                this->velocity(x, y, z).y -= .5 * (pressure(x  , y-1, z  ) - pressure(x  , y+1, z  )) * this->sizeY;
                this->velocity(x, y, z).z -= .5 * (pressure(x  , y  , z-1) - pressure(x  , y  , z+1)) * this->sizeZ;
            }
        }
    }
    this->set_bounds(this->velocity, true);
}



