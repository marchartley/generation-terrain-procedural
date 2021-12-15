#include "FluidSimulation.h"

FluidSimulation::FluidSimulation()
{

}

FluidSimulation::FluidSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations)
    : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), dt(dt), diffusionAmount(diffusionAmount), viscosity(viscosity), iterations(iterations)
{
    obstacles = std::vector<std::vector<std::vector<bool>>>(sizeX, std::vector<std::vector<bool>>(sizeY, std::vector<bool>(sizeZ)));
    density_old = std::vector<std::vector<std::vector<float>>>(sizeX, std::vector<std::vector<float>>(sizeY, std::vector<float>(sizeZ)));
    density = std::vector<std::vector<std::vector<float>>>(sizeX, std::vector<std::vector<float>>(sizeY, std::vector<float>(sizeZ)));

    velocity = std::vector<std::vector<std::vector<Vector3>>>(sizeX, std::vector<std::vector<Vector3>>(sizeY, std::vector<Vector3>(sizeZ)));
    velocity_old = std::vector<std::vector<std::vector<Vector3>>>(sizeX, std::vector<std::vector<Vector3>>(sizeY, std::vector<Vector3>(sizeZ)));

    divergence = std::vector<std::vector<std::vector<float>>>(sizeX, std::vector<std::vector<float>>(sizeY, std::vector<float>(sizeZ)));
    pressure = std::vector<std::vector<std::vector<float>>>(sizeX, std::vector<std::vector<float>>(sizeY, std::vector<float>(sizeZ)));
}

void FluidSimulation::setObstacles(std::vector<std::vector<std::vector<bool> > > new_obstacles)
{
    if (int(new_obstacles.size()) == this->sizeX && int(new_obstacles[0].size()) == this->sizeY && int(new_obstacles[0][0].size()) == this->sizeZ) {
        this->obstacles = new_obstacles;
        return;
    }
    float rx = int(new_obstacles.size()) / (float)sizeX, ry = int(new_obstacles[0].size()) / (float)sizeY, rz = int(new_obstacles[0][0].size()) / (float)sizeZ;
    std::vector<std::vector<std::vector<float> > > obstacleRate = std::vector<std::vector<std::vector<float>>>(sizeX, std::vector<std::vector<float>>(sizeY, std::vector<float>(sizeZ)));
    std::vector<std::vector<std::vector<int> > > obstacleNumbers = std::vector<std::vector<std::vector<int>>>(sizeX, std::vector<std::vector<int>>(sizeY, std::vector<int>(sizeZ)));
    for (int x = 0; x < int(new_obstacles.size()); x++) {
        for (int y = 0; y < int(new_obstacles[0].size()); y++) {
            for (int z = 0; z < int(new_obstacles[0][0].size()); z++) {
                int _x = std::max(0, std::min(int(new_obstacles.size() - 1), x)), _x_resized = std::max(0, std::min(this->sizeX - 1, int(x / rx)));
                int _y = std::max(0, std::min(int(new_obstacles[0].size() - 1), y)), _y_resized = std::max(0, std::min(this->sizeY - 1, int(y / ry)));
                int _z = std::max(0, std::min(int(new_obstacles[0][0].size() - 1), z)), _z_resized = std::max(0, std::min(this->sizeZ - 1, int(z / rz)));
                obstacleRate[_x_resized][_y_resized][_z_resized] += (new_obstacles[_x][_y][_z] ? 1.0 : 0.0);
                obstacleNumbers[_x_resized][_y_resized][_z_resized] ++;
//                float rate = obstacleRate[_x_resized][_y_resized][_z_resized] / (float)obstacleNumbers[_x_resized][_y_resized][_z_resized];
//                std::cout << rate << std::endl;
            }
        }
    }
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (obstacleRate[x][y][z] != obstacleNumbers[x][y][z])
                    int w = 0;
                this->obstacles[x][y][z] = (obstacleRate[x][y][z] / (float)obstacleNumbers[x][y][z]) > .5f;
            }
        }
    }
}

void FluidSimulation::addDensity(int x, int y, int z, float amount)
{
    this->density[x][y][z] += amount;
}

void FluidSimulation::addVelocity(int x, int y, int z, Vector3 amount)
{
    this->velocity[x][y][z] += amount;
}

std::vector<std::vector<std::vector<Vector3> > > FluidSimulation::getVelocities(int rescaleX, int rescaleY, int rescaleZ)
{
    if (rescaleX == -1) rescaleX = sizeX;
    if (rescaleY == -1) rescaleY = sizeY;
    if (rescaleZ == -1) rescaleZ = sizeZ;
    if (rescaleX == this->sizeX && rescaleY == this->sizeY && rescaleZ == this->sizeZ) return this->velocity;
    float rx = rescaleX / (float)sizeX, ry = rescaleY / (float)sizeY, rz = rescaleZ / (float)sizeZ;
    std::vector<std::vector<std::vector<Vector3> > > out;
    for (int x = 0; x < rescaleX; x++) {
        out.push_back(std::vector<std::vector<Vector3>>());
        for (int y = 0; y < rescaleY; y++) {
            out[x].push_back(std::vector<Vector3>());
            for (int z = 0; z < rescaleZ; z++) {
                int _x = std::max(0, std::min(this->sizeX - 1, int(x/rx)));
                int _y = std::max(0, std::min(this->sizeY - 1, int(y/ry)));
                int _z = std::max(0, std::min(this->sizeZ - 1, int(z/rz)));
                out[x][y].push_back(this->velocity[_x][_y][_z]);
            }
        }
    }
    return out;
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
                this->divergence[x][y][z] = -.5 * (
                          velocity[x+1][y  ][z  ].x
                        - velocity[x-1][y  ][z  ].x
                        + velocity[x  ][y+1][z  ].y
                        - velocity[x  ][y-1][z  ].y
                        + velocity[x  ][y  ][z+1].z
                        - velocity[x  ][y  ][z-1].z) / std::sqrt(this->sizeX * this->sizeY * this->sizeZ);
                this->pressure[x][y][z] = 0.0;
            }
        }
    }
    this->set_bounds(this->divergence);
    this->set_bounds(this->pressure);
    this->solve_linear(this->pressure, this->divergence, 1, 6, false);

    for (int x = 1; x < this->sizeX - 1; x++) {
        for (int y = 1; y < this->sizeY - 1; y++) {
            for (int z = 1; z < this->sizeZ - 1; z++) {
                this->velocity[x][y][z].x -= .5 * (pressure[x-1][y  ][z  ] - pressure[x+1][y  ][z  ]) * this->sizeX;
                this->velocity[x][y][z].y -= .5 * (pressure[x  ][y-1][z  ] - pressure[x  ][y+1][z  ]) * this->sizeY;
                this->velocity[x][y][z].z -= .5 * (pressure[x  ][y  ][z-1] - pressure[x  ][y  ][z+1]) * this->sizeZ;
            }
        }
    }
    this->set_bounds(this->velocity, true);
}



