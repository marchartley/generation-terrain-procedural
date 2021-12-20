#include "FluidSimulation.h"

FluidSimulation::FluidSimulation()
{

}

FluidSimulation::FluidSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations)
    : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), dt(dt), diffusionAmount(diffusionAmount), viscosity(viscosity), iterations(iterations)
{
    obstacles = Matrix3<float>(sizeX, sizeY, sizeZ);
    density_old = Matrix3<float>(sizeX, sizeY, sizeZ);
    density = Matrix3<float>(sizeX, sizeY, sizeZ);

    velocity = Matrix3<Vector3>(sizeX, sizeY, sizeZ);
    velocity_old = Matrix3<Vector3>(sizeX, sizeY, sizeZ);

    divergence = Matrix3<float>(sizeX, sizeY, sizeZ);
    pressure = Matrix3<float>(sizeX, sizeY, sizeZ);
}

void FluidSimulation::setObstacles(Matrix3<float> new_obstacles)
{
    this->obstacles = new_obstacles.resize(sizeX, sizeY, sizeZ);
/*
    for(int i = 0; i < this->obstacles.data.size(); i++) {
        if (this->obstacles.data[i] < 0.5) {
            int x, y, z;
            std::tie(x, y, z) = this->obstacles.getCoord(i);
            std::cout << x << " " << y << " " << z << " | " << std::endl;
        }
    }*/
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
    /*Matrix3<Vector3> mat(sizeX, sizeY, sizeZ);
    for (size_t i = 0; i < mat.data.size(); i++)
        mat.at(i) = (Vector3(1, 0, 0) * obstacles.at(i) );
    return mat.resize(rescaleX, rescaleY, rescaleZ);*/
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
//    this->densityStep(); // This is actually useless right now. Maybe we could use it to erod the walls depending on the density, but that's not for now.
    this->set_bounds(this->velocity, true, false);
//    this->set_bounds(this->density, false, true);
}

void FluidSimulation::velocityStep()
{
    swapArrays(this->velocity, this->velocity_old); // Originally was on first line
    this->diffuse(this->velocity, this->velocity_old, this->viscosity);
    this->project();
    swapArrays(this->velocity_old, this->velocity);
    this->advect(this->velocity, this->velocity_old, this->velocity_old);
    this->project();
    this->setVelocityBounds();
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
    this->pressure = Matrix3<float>(this->sizeX, this->sizeY, this->sizeZ);

    // What is "h"??? In Josh Stam's code, it's just 1/N
    float h = 1/std::sqrt(this->sizeX * this->sizeY * this->sizeZ);
    this->divergence = velocity.divergence() * (-h/2.0);
    /*
    for (int x = 1; x < this->sizeX - 1; x++) {
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                this->divergence(x, y, z) = -.5 * (
                          velocity(x+1, y  , z  ).x
                        - velocity(x-1, y  , z  ).x
                        + velocity(x  , y+1, z  ).y
                        - velocity(x  , y-1, z  ).y
                        + velocity(x  , y  , z+1).z
                        - velocity(x  , y  , z-1).z) * h;
            }
        }
    }*/
    this->set_bounds(this->divergence);
    this->set_bounds(this->pressure);
    this->solve_linear(this->pressure, this->divergence, 1, false);

    Matrix3<Vector3> pressureGradient = this->pressure.gradient() / h;
    this->velocity -= pressureGradient / 2.f;
    /*
    for (int x = 1; x < this->sizeX - 1; x++) {
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                this->velocity(x, y, z).x -= .5 * (pressure(x-1, y  , z  ) - pressure(x+1, y  , z  )) * this->sizeX;
                this->velocity(x, y, z).y -= .5 * (pressure(x  , y-1, z  ) - pressure(x  , y+1, z  )) * this->sizeY;
                this->velocity(x, y, z).z -= .5 * (pressure(x  , y  , z-1) - pressure(x  , y  , z+1)) * this->sizeZ;
            }
        }
    }*/
    this->set_bounds(this->velocity, true);
}

void FluidSimulation::setVelocityBounds()
{
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3 origin = Vector3(x, y, z);
                Vector3 vel = velocity(x, y, z);
                // This is not the best way, but it's quite fast
                if ((origin + vel).x < 0) vel.x = -origin.x;
                if ((origin + vel).x >= sizeX) vel.x = sizeX-origin.x;
                if ((origin + vel).y < 0) vel.y = -origin.y;
                if ((origin + vel).y >= sizeY) vel.y = sizeY-origin.y;
                if ((origin + vel).z < 0) vel.z = -origin.z;
                if ((origin + vel).z >= sizeZ) vel.z = sizeZ-origin.z;
                velocity(x, y, z) = vel;
            }
        }
    }
}



