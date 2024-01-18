#include "FluidSimulation/StableFluidsFluidSimulation.h"
#include "FluidSimulation/OpenFoamParser.h"

/*
StableFluids::StableFluidsSimulation::StableFluidsSimulation() {}

StableFluids::StableFluidsSimulation::StableFluidsSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusion, float viscosity, int solverIterations)
    : FluidSimulation(sizeX, sizeY, sizeZ), diffusion(diffusion), viscosity(viscosity), dt(dt), solverIterations(solverIterations)
{
    velocity = GridV3(dimensions);
    density = GridF(dimensions);

    velocity.raiseErrorOnBadCoord = false;
    velocity.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    density.raiseErrorOnBadCoord = false;
    density.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
}

void StableFluids::StableFluidsSimulation::step() {
    diffuse();
    project();
    advect();
    project();
    handleCollisions();
}

//void StableFluids::StableFluidsSimulation::addDensity(int x, int y, int z, float amount) {
//    density.at(Vector3(x, y, z)) += amount;
//}

void StableFluids::StableFluidsSimulation::addVelocity(int x, int y, int z, const Vector3& amount) {
    velocity.at(Vector3(x, y, z)) += amount;
}

void StableFluids::StableFluidsSimulation::diffuse() {
    int sizeX = dimensions.x, sizeY = dimensions.y, sizeZ = dimensions.z;
    float a = dt * viscosity * dimensions.x * dimensions.y * dimensions.z * diffusion;
    for (int k = 0; k < solverIterations; ++k) { // 20 iterations is a common choice
        auto oldVelocities = velocity;
        oldVelocities.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
#pragma omp parallel for collapse(3)
        for (int x = 0; x < sizeX; ++x) {
            for (int y = 0; y < sizeY; ++y) {
                for (int z = 0; z < sizeZ; ++z) {
                    Vector3 pos(x, y, z);
                    Vector3 sum = oldVelocities.at(Vector3(x-1, y, z)) + oldVelocities.at(Vector3(x+1, y, z))
                            + oldVelocities.at(Vector3(x, y-1, z)) + oldVelocities.at(Vector3(x, y+1, z))
                            + oldVelocities.at(Vector3(x, y, z-1)) + oldVelocities.at(Vector3(x, y, z+1));
                    velocity.at(pos) = (velocity.at(pos) + a * sum) / (1 + 6 * a);
                }
            }
        }
    }
}

void StableFluids::StableFluidsSimulation::project() {
    GridF divergence(dimensions);
    GridF pressure(dimensions);
    int sizeX = dimensions.x, sizeY = dimensions.y, sizeZ = dimensions.z;

    // Compute divergence of velocity field
#pragma omp parallel for collapse(3)
    for (int x = 1; x < sizeX - 1; ++x) {
        for (int y = 1; y < sizeY - 1; ++y) {
            for (int z = 1; z < sizeZ - 1; ++z) {
                Vector3 pos(x, y, z);
                divergence.at(pos) = -0.5f * (
                            velocity.at(Vector3(x+1, y, z)).x - velocity.at(Vector3(x-1, y, z)).x +
                            velocity.at(Vector3(x, y+1, z)).y - velocity.at(Vector3(x, y-1, z)).y +
                            velocity.at(Vector3(x, y, z+1)).z - velocity.at(Vector3(x, y, z-1)).z
                            );
                pressure.at(pos) = 0;
            }
        }
    }

    // Solve pressure equation
    for (int k = 0; k < solverIterations; ++k) {
        auto oldPressure = pressure;
#pragma omp parallel for collapse(3)
        for (int x = 1; x < sizeX - 1; ++x) {
            for (int y = 1; y < sizeY - 1; ++y) {
                for (int z = 1; z < sizeZ - 1; ++z) {
                    Vector3 pos(x, y, z);
                    pressure.at(pos) = (divergence.at(pos) +
                                        oldPressure.at(Vector3(x+1, y, z)) + oldPressure.at(Vector3(x-1, y, z)) +
                                        oldPressure.at(Vector3(x, y+1, z)) + oldPressure.at(Vector3(x, y-1, z)) +
                                        oldPressure.at(Vector3(x, y, z+1)) + oldPressure.at(Vector3(x, y, z-1))
                                        ) / 6.0f;
                }
            }
        }
    }

    // Subtract pressure gradient from velocity field
#pragma omp parallel for collapse(3)
    for (int x = 1; x < sizeX - 1; ++x) {
        for (int y = 1; y < sizeY - 1; ++y) {
            for (int z = 1; z < sizeZ - 1; ++z) {
                Vector3 pos(x, y, z);
                velocity.at(pos).x -= 0.5f * (pressure.at(Vector3(x+1, y, z)) - pressure.at(Vector3(x-1, y, z)));
                velocity.at(pos).y -= 0.5f * (pressure.at(Vector3(x, y+1, z)) - pressure.at(Vector3(x, y-1, z)));
                velocity.at(pos).z -= 0.5f * (pressure.at(Vector3(x, y, z+1)) - pressure.at(Vector3(x, y, z-1)));
            }
        }
    }
}

void StableFluids::StableFluidsSimulation::advect() {
    GridV3 newVelocity(dimensions);
    int sizeX = dimensions.x, sizeY = dimensions.y, sizeZ = dimensions.z;

#pragma omp parallel for collapse(3)
    for (int x = 0; x < sizeX; ++x) {
        for (int y = 0; y < sizeY; ++y) {
            for (int z = 0; z < sizeZ; ++z) {
                Vector3 pos(x, y, z);

                // Trace back the particle's position in the previous timestep
                Vector3 lastPos = pos - dt * velocity.at(pos);

                // Clamp the position to be within the simulation box
                lastPos.x = std::max(0.0f, std::min((float)sizeX, lastPos.x));
                lastPos.y = std::max(0.0f, std::min((float)sizeY, lastPos.y));
                lastPos.z = std::max(0.0f, std::min((float)sizeZ, lastPos.z));

                // Interpolate the velocity at the new position
                newVelocity.at(pos) = velocity.interpolate(lastPos);
            }
        }
    }

    // Update the velocity field
    velocity = newVelocity;
}

void StableFluids::StableFluidsSimulation::handleCollisions() {
    int sizeX = dimensions.x, sizeY = dimensions.y, sizeZ = dimensions.z;
//#pragma omp parallel for collapse(3)
    for (int x = 0; x < sizeX; ++x) {
        for (int y = 0; y < sizeY; ++y) {
            for (int z = 0; z < sizeZ; ++z) {
                Vector3 pos(x, y, z);

                if (obstacleGrid.at(pos)) {
                    velocity.at(pos) = obstacleGradient.at(pos) * velocity.at(pos).norm();
                } else if (obstacleTrianglesOctree) {
                    // Query Octree for nearby triangles
                    Vector3 endPos = pos + velocity.at(pos).normalized();
                    std::vector<OctreeNodeData> nearbyTriangles = obstacleTrianglesOctree->queryRange(pos, endPos);

                    // Check for intersections with nearby triangles
                    for (auto& triangleData : nearbyTriangles) {
                        auto& triangle = this->triangles[triangleData.index];
                        Vector3 collisionPoint = Collision::segmentToTriangleCollision(pos, endPos, triangle[0], triangle[1], triangle[2]);
                        if (collisionPoint.isValid()) {
                            velocity.at(pos) = (collisionPoint - pos).normalized() * velocity.at(pos).norm();
                            break;
                        }
                    }
                }
            }
        }
    }
}

GridV3 StableFluids::StableFluidsSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ) {
    return (this->velocity * (1.f - this->obstacleGrid)).resize(newSizeX, newSizeY, newSizeZ);
}

void StableFluidsSimulation::setObstacles(GridF new_obstacles)
{
    FluidSimulation::setObstacles(new_obstacles);
    for(size_t i = 0; i < new_obstacles.data.size(); i++) {
        if (new_obstacles[i] > 0.001) {
//            this->velocity[i] = Vector3();
//            this->velocity_old[i] = Vector3();
        }
    }
}
*/

/*
StableFluidsSimulation::StableFluidsSimulation()
{

}

StableFluidsSimulation::StableFluidsSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations)
    : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), dt(dt), diffusionAmount(diffusionAmount), viscosity(viscosity), iterations(iterations)
{
    obstacles = GridF(dimensions);
    density_old = GridF(dimensions);
    density = GridF(dimensions);

    velocity = GridV3(dimensions);
    velocity_old = GridV3(dimensions);

    divergence = GridF(dimensions);
    pressure = GridF(dimensions);

    maxSpeedSquared = maxSpeed * maxSpeed;
}

void StableFluidsSimulation::setObstacles(GridF new_obstacles)
{
    this->obstacles = new_obstacles.resize(dimensions).binarize(0.5);
    for(size_t i = 0; i < this->obstacles.data.size(); i++) {
        if (this->obstacles[i] > 0.001) {
            this->velocity[i] = Vector3();
            this->velocity_old[i] = Vector3();
        }
    }
}

void StableFluidsSimulation::addDensity(int x, int y, int z, float amount)
{
    this->density(x, y, z) += amount;
}

void StableFluidsSimulation::addVelocity(int x, int y, int z, const Vector3& amount)
{
    this->velocity(x, y, z) += amount;
}

void StableFluidsSimulation::setMaxSpeed(float speed)
{
    this->maxSpeed = speed;
    this->maxSpeedSquared = speed * speed;
}

GridV3 StableFluidsSimulation::getVelocities(int rescaleX, int rescaleY, int rescaleZ)
{

    return this->velocity.resize(rescaleX, rescaleY, rescaleZ);
}

void StableFluidsSimulation::step()
{
    this->currentStep ++;
    this->velocityStep();
//    this->densityStep(); // This is actually useless right now. Maybe we could use it to erod the walls depending on the density, but that's not for now.
    this->set_bounds(this->velocity, true, false);
    this->set_bounds(this->density, false, true);
}

void StableFluidsSimulation::velocityStep()
{
    this->velocity += this->velocity_old * this->dt;
//    swapArrays(this->velocity, this->velocity_old); // Originally was on first line
    this->diffuseVelocity(); // Removed for now (considering viscosity = 0)
//    velocity_old = velocity; // Replaced by a simple affectation

//    this->projectVelocity();

//    swapArrays(this->velocity_old, this->velocity);
//    this->advectVelocity();
//    this->projectVelocity();
//    this->setVelocityBounds();
}

void StableFluidsSimulation::diffuseVelocity()
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

void StableFluidsSimulation::advectVelocity()
{
    Vector3 dt0 = Vector3(1.f, 1.f, 1.f) * dt; // = Vector3(dimensions) * dt;

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

void StableFluidsSimulation::densityStep()
{
    swapArrays(this->density, this->density_old);
    this->diffuse(this->density, this->density_old, this->diffusionAmount);
    swapArrays(this->density_old, this->density);
    this->advect(this->density, this->density_old, this->velocity);
}

void StableFluidsSimulation::projectVelocity()
{
    this->pressure = GridF(this->sizeX, this->sizeY, this->sizeZ);

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
    GridF tmp = pressure;
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

//    GridV3 pressureGradient = this->pressure.gradient() / (-h * 2.f);
    // Following ethanjli code :
    GridV3 pressureGradient = this->pressure.gradient() * -h;
    this->velocity += pressureGradient;
    this->setVelocityBounds();
}

void StableFluidsSimulation::setVelocityBounds()
{
    bool nullifyOnBounds = false;
    bool inverseOnBounds = true;
    GridV3 boundariesGradient = this->obstacles.gradient() * (-1.f);

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
*/



StableFluidsSimulation::StableFluidsSimulation()
{

}

StableFluidsSimulation::StableFluidsSimulation(int sizeX, int sizeY, int sizeZ, float dt, float diffusionAmount, float viscosity, int iterations)
    : FluidSimulation(sizeX, sizeY, sizeZ), sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), dt(dt), diffusionAmount(diffusionAmount), viscosity(viscosity), iterations(iterations)
{
//    obstacles = GridF(dimensions);
    density_old = GridF(dimensions);
    density = GridF(dimensions);

    velocity = GridV3(dimensions);
    velocity_old = GridV3(dimensions);

    divergence = GridF(dimensions);
    pressure = GridF(dimensions);

    maxSpeedSquared = maxSpeed * maxSpeed;
}

void StableFluidsSimulation::setObstacles(const GridF &new_obstacles)
{
    /*this->obstacles = new_obstacles.resizeNearest(dimensions);
    for(size_t i = 0; i < this->obstacles.data.size(); i++) {
        if (this->obstacles[i] > 0.001) {
            this->velocity[i] = Vector3();
            this->velocity_old[i] = Vector3();
        }
    }*/

    this->obstacleGrid = new_obstacles.resizeNearest(dimensions);
    for(size_t i = 0; i < this->obstacleGrid.data.size(); i++) {
        if (this->obstacleGrid[i] > 0.001) {
            this->velocity[i] = Vector3();
            this->velocity_old[i] = Vector3();
        }
    }
}

void StableFluidsSimulation::addDensity(int x, int y, int z, float amount)
{
    this->density(x, y, z) += amount;
}

void StableFluidsSimulation::addVelocity(int x, int y, int z, const Vector3 &amount)
{
    this->velocity(x, y, z) += amount;
}

void StableFluidsSimulation::setMaxSpeed(float speed)
{
    this->maxSpeed = speed;
    this->maxSpeedSquared = speed * speed;
}

GridV3 StableFluidsSimulation::getVelocities(int rescaleX, int rescaleY, int rescaleZ)
{
    if (_cachedStep != currentStep) {
        _cachedStep = currentStep;
        std::string pathToOpenFoamSimulation = "OF_Sim_Marcos";
        if (checkPathExists(pathToOpenFoamSimulation)) {
            this->velocity = OpenFoamParser::parseSimulation(pathToOpenFoamSimulation);
            Vector3 aspectRatio = Vector3(rescaleX, rescaleY, rescaleZ).normalize();
            this->velocity *= aspectRatio;
            std::cout << "OpenFOAM simu size: " << velocity.getDimensions() << std::endl;
//            this->velocity = this->velocity.resize(sizeX, sizeY, sizeZ);
        }
        _cachedVelocity = velocity;
//        std::cout << "Recomputed : " << _cachedVelocity.min() << " " << _cachedVelocity.max() << std::endl;
    } else {
//        std::cout << "Not recomputed -> step " << currentStep << std::endl;
    }
    return _cachedVelocity.resize(rescaleX, rescaleY, rescaleZ);
}

Vector3 StableFluidsSimulation::getVelocity(int x, int y, int z)
{
    return FluidSimulation::getVelocity(x, y, z);
}

void StableFluidsSimulation::step()
{
    this->currentStep ++;
    this->velocityStep();
//    this->densityStep(); // This is actually useless right now. Maybe we could use it to erod the walls depending on the density, but that's not for now.
//    this->set_bounds(this->velocity, true, false);
//    this->set_bounds(this->density, false, true);
}

void StableFluidsSimulation::velocityStep()
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

void StableFluidsSimulation::diffuseVelocity()
{
    this->velocity_old.raiseErrorOnBadCoord = false;

    float a = dt * viscosity; // * sizeX * sizeY * sizeZ; // Removing "N^3" following ethanjli code
    for (int i = 0; i < this->iterations; i++) {
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

void StableFluidsSimulation::advectVelocity()
{
    Vector3 dt0 = Vector3(1.f, 1.f, 1.f) * dt; // = Vector3(dimensions) * dt;

    velocity_old.raiseErrorOnBadCoord = false;
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
//    this->setVelocityBounds();
}

void StableFluidsSimulation::densityStep()
{
    swapArrays(this->density, this->density_old);
    this->diffuse(this->density, this->density_old, this->diffusionAmount);
    swapArrays(this->density_old, this->density);
    this->advect(this->density, this->density_old, this->velocity);
}

void StableFluidsSimulation::projectVelocity()
{
    this->pressure = GridF(this->sizeX, this->sizeY, this->sizeZ);

    // What is "h"??? In Josh Stam's code, it's just 1/N
    // So I think that if "h" is big (.5 < h < 1.), there are alternating attraction/repulsion "poles"
    // These poles are having a distance of N/h between them, so if we set h to 1/N, this effect looks canceleld
//    float h = 1.f / (float)(sizeX * sizeY * sizeZ); // (float)std::sqrt(this->sizeX * this->sizeY * this->sizeZ);
//    this->divergence = velocity.divergence() * (h/2.f);
    // In ethanjli's code, the divergence is *= -1, so we try that here
    float h = 1.f/(float)std::max(this->sizeX, std::max(this->sizeY, this->sizeZ));
    this->divergence = velocity.divergence() * -h;

    this->set_bounds(this->divergence);
//    this->set_bounds(this->pressure);

    this->pressure = divergence;
    this->pressure.raiseErrorOnBadCoord = false;
    GridF tmp = pressure;
    for (int i = 0; i < this->iterations; i++) {
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

//    GridV3 pressureGradient = this->pressure.gradient() / (-h * 2.f);
    // Following ethanjli code :
    GridV3 pressureGradient = this->pressure.gradient() * -h;
    this->velocity += pressureGradient;
//    this->setVelocityBounds();
}

void StableFluidsSimulation::setVelocityBounds()
{
    bool nullifyOnBounds = false;
    bool inverseOnBounds = true;
    GridV3 boundariesGradient = this->obstacleGrid.gradient() * (-1.f);

    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                Vector3 origin(x, y, z);
                if (velocity(x, y, z).norm2() > this->maxSpeedSquared) velocity(x, y, z).setMag(this->maxSpeed);
                bool isGoingThroughObstable = (obstacleGrid.at(origin) || obstacleGrid.at(origin + velocity(x, y, z).normalized()));

//                = (velocity.checkCoord(origin + velocity.at(origin) * 1.5f) ? obstacles.at(origin + velocity.at(origin) * 1.5f) > .01 : false);
                if (isGoingThroughObstable) {
                    if (inverseOnBounds) {
                        velocity.at(x, y, z) = boundariesGradient(x, y, z).normalized() * velocity.at(x, y, z).norm();
//                        std::cout << "Inversing at " << origin << "\n";
                    } else if (nullifyOnBounds) {
                        velocity.at(x, y, z) *= 0.f;
                        velocity_old.at(x, y, z) *= 0.f;
//                        std::cout << "Nullify at " << origin << "\n";
                    } else {
//                        std::cout << "WTF?" << "\n";
                    }
                } else {
//                    std::cout << "No obstacle at " << origin << "\n";
                }
            }
        }
    }
}

