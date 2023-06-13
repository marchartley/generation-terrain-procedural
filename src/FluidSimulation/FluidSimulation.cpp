#include "FluidSimulation/FluidSimulation.h"
/*
template<class T>
T clamp(T val, T min, T max) {
    val = std::min(std::max(val, min), max);
    return val;
}
*/

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

    maxSpeedSquared = maxSpeed * maxSpeed;
}

void FluidSimulation::setObstacles(Matrix3<float> new_obstacles)
{
    this->obstacles = new_obstacles.resize(sizeX, sizeY, sizeZ).binarize(0.5);
    for(size_t i = 0; i < this->obstacles.data.size(); i++) {
        if (this->obstacles[i] > 0.001) {
            this->velocity[i] = Vector3();
            this->velocity_old[i] = Vector3();
        }
    }
}

void FluidSimulation::addDensity(int x, int y, int z, float amount)
{
    this->density(x, y, z) += amount;
}

void FluidSimulation::addVelocity(int x, int y, int z, Vector3 amount)
{
    this->velocity(x, y, z) += amount;
}

void FluidSimulation::setMaxSpeed(float speed)
{
    this->maxSpeed = speed;
    this->maxSpeedSquared = speed * speed;
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
    this->currentStep ++;
    this->velocityStep();
//    this->densityStep(); // This is actually useless right now. Maybe we could use it to erod the walls depending on the density, but that's not for now.
    this->set_bounds(this->velocity, true, false);
    this->set_bounds(this->density, false, true);
}

void FluidSimulation::velocityStep()
{
    this->velocity += this->velocity_old * this->dt;
//    swapArrays(this->velocity, this->velocity_old); // Originally was on first line
    this->diffuseVelocity(); // Removed for now (considering viscosity = 0)
//    velocity_old = velocity; // Replaced by a simple affectation

    this->projectVelocity();

    swapArrays(this->velocity_old, this->velocity);
    this->advectVelocity();
    this->projectVelocity();
    this->setVelocityBounds();
}

void FluidSimulation::diffuseVelocity()
{
    this->velocity_old.raiseErrorOnBadCoord = false;

    float a = dt * viscosity; // * sizeX * sizeY * sizeZ; // Removing "N^3" following ethanjli code
    for (int i = 0; i < this->iterations; i++) {
#pragma omp parallel for collapse(3)
        for (int x = 0; x < sizeX; x++) {
            for (int y = 0; y < sizeY; y++) {
                for (int z = 0; z < sizeZ; z++) {
                    this->velocity_old(x, y, z) = this->velocity(x, y, z) + ((
                                velocity_old(x    , y    , z - 1) +
                                velocity_old(x    , y    , z + 1) +
                                velocity_old(x    , y - 1, z    ) +
                                velocity_old(x    , y + 1, z    ) +
                                velocity_old(x - 1, y    , z    ) +
                                velocity_old(x + 1, y    , z    )
                                ) * a) / (float)(1 + velocity.getNumberNeighbors(x, y, z) * a);
                }
            }
        }
//        swapArrays(velocity, velocity_old);
//        this->setVelocityBounds();
//        swapArrays(velocity, velocity_old);
    }
    this->velocity_old.raiseErrorOnBadCoord = true;
}

void FluidSimulation::advectVelocity()
{
    Vector3 dt0 = Vector3(1.f, 1.f, 1.f) * dt; // = Vector3(sizeX, sizeY, sizeZ) * dt;

    velocity_old.raiseErrorOnBadCoord = false;
#pragma omp parallel for collapse(3)
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                Vector3 emittingCell = Vector3(x, y, z) - (velocity_old(x, y, z) * dt0);
                float xCell = clamp(emittingCell.x, 0.f, (float)sizeX - 0.0001f);
                float yCell = clamp(emittingCell.y, 0.f, (float)sizeY - 0.0001f);
                float zCell = clamp(emittingCell.z, 0.f, (float)sizeZ - 0.0001f);

//                if (!velocity.checkCoord(emittingCell)) {
//                    velocity.at(emittingCell) -= velocity_old.at(x, y, z);
//                    continue;
//                }
                // Just computing some linear interpolation between the bounds of the emitting cell
                int xFloor = int(xCell); int xCeil = xFloor + 1; float xAlpha = xCeil - x;
                int yFloor = int(yCell); int yCeil = yFloor + 1; float yAlpha = yCeil - y;
                int zFloor = int(zCell); int zCeil = zFloor + 1; float zAlpha = zCeil - z;

                velocity(x, y, z) = ((velocity_old(xFloor, yFloor, zFloor) * zAlpha + velocity_old(xFloor, yFloor, zCeil) * (1 - zAlpha)) * yAlpha +
                                     (velocity_old(xFloor,  yCeil, zFloor) * zAlpha + velocity_old(xFloor,  yCeil, zCeil) * (1 - zAlpha)) * (1 - yAlpha)) * xAlpha +
                                    ((velocity_old( xCeil, yFloor, zFloor) * zAlpha + velocity_old( xCeil, yFloor, zCeil) * (1 - zAlpha)) * yAlpha +
                                     (velocity_old( xCeil,  yCeil, zFloor) * zAlpha + velocity_old( xCeil,  yCeil, zCeil) * (1 - zAlpha)) * (1 - yAlpha)) * (1 - xAlpha);
            }
        }
    }
    velocity_old.raiseErrorOnBadCoord = true;
    this->setVelocityBounds();
}

void FluidSimulation::densityStep()
{
    swapArrays(this->density, this->density_old);
    this->diffuse(this->density, this->density_old, this->diffusionAmount);
    swapArrays(this->density_old, this->density);
    this->advect(this->density, this->density_old, this->velocity);
}

void FluidSimulation::projectVelocity()
{
    this->pressure = Matrix3<float>(this->sizeX, this->sizeY, this->sizeZ);

    // What is "h"??? In Josh Stam's code, it's just 1/N
    // So I think that if "h" is big (.5 < h < 1.), there are alternating attraction/repulsion "poles"
    // These poles are having a distance of N/h between them, so if we set h to 1/N, this effect looks canceleld
//    float h = 1.f / (float)(sizeX * sizeY * sizeZ); // (float)std::sqrt(this->sizeX * this->sizeY * this->sizeZ);
//    this->divergence = velocity.divergence() * (h/2.f);
    // In ethanjli's code, the divergence is *= -1, so we try that here
    float h = 1.f/(float)std::max(this->sizeX, std::max(this->sizeY, this->sizeZ));
    this->divergence = velocity.divergence() * -h;

    this->set_bounds(this->divergence);
    this->set_bounds(this->pressure);

    this->pressure = divergence;
    this->pressure.raiseErrorOnBadCoord = false;
    Matrix3<float> tmp = pressure;
    for (int i = 0; i < this->iterations; i++) {
#pragma omp parallel for collapse(3)
        for (int x = 0; x < sizeX; x++) {
            for (int y = 0; y < sizeY; y++) {
                for (int z = 0; z < sizeZ; z++) {
                    tmp(x, y, z) = (this->divergence(x, y, z) +
                                                pressure(x    , y    , z - 1) +
                                                pressure(x    , y    , z + 1) +
                                                pressure(x    , y - 1, z    ) +
                                                pressure(x    , y + 1, z    ) +
                                                pressure(x - 1, y    , z    ) +
                                                pressure(x + 1, y    , z    )
                                ) / (float)(pressure.getNumberNeighbors(x, y, z));
                }
            }
        }
        pressure = tmp;
        set_bounds(pressure);
    }
    this->pressure.raiseErrorOnBadCoord = true;

//    Matrix3<Vector3> pressureGradient = this->pressure.gradient() / (-h * 2.f);
    // Following ethanjli code :
    Matrix3<Vector3> pressureGradient = this->pressure.gradient() * -h;
    this->velocity += pressureGradient;
    this->setVelocityBounds();
}

void FluidSimulation::setVelocityBounds()
{
    bool nullifyOnBounds = false;
    bool inverseOnBounds = true;
    Matrix3<Vector3> boundariesGradient = this->obstacles.gradient() * (-1.f);

#pragma omp parallel for collapse(3)
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                Vector3 origin(x, y, z);
                if (velocity(x, y, z).norm2() > this->maxSpeedSquared) velocity(x, y, z) = velocity(x, y, z).normalized() * this->maxSpeed;
                bool isGoingThroughObstable = (velocity.checkCoord(origin + velocity.at(origin) * 1.5f) ? obstacles.at(origin + velocity.at(origin) * 1.5f) > .01 : false);
                if (isGoingThroughObstable) {
                    if (inverseOnBounds)
                        velocity.at(x, y, z) = boundariesGradient(x, y, z).normalized() * velocity.at(x, y, z).norm(); //velocity.at(x, y, z) - boundariesGradient(x, y, z).normalized() * (velocity.at(x, y, z).normalized().dot(boundariesGradient(x, y, z).normalized() * -1.f) * .05f) * 2.f;
                    if (nullifyOnBounds)
                        velocity.at(x, y, z) *= 0.f;
                }
            }
        }
    }
}



