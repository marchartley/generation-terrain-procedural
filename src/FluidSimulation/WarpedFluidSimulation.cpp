#include "WarpedFluidSimulation.h"

WarpedFluidSimulation::WarpedFluidSimulation()
{
}

WarpedFluidSimulation::WarpedFluidSimulation(int x, int y, int z)
    : FluidSimulation(x, y, 1)
{
    this->recomputeVelocities();
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
    // No velocity added for now.
//    velocities.at(x, y, z) += amount;
}


void WarpedFluidSimulation::setObstacles(const GridF &obstacles)
{
    this->obstacleGrid = GridF(this->dimensions);
    auto resizedObstacles = obstacles.resize(this->dimensions.x, this->dimensions.y, obstacles.sizeZ);
    for (int x = 0; x < resizedObstacles.sizeX; x++) {
        for (int y = 0; y < resizedObstacles.sizeY; y++) {
            float height = 0;
            while(resizedObstacles.at(x, y, resizedObstacles.sizeZ - height) <= 0 && height < resizedObstacles.sizeZ)
                height++;
            obstacleGrid.at(x, y) = height / float(obstacles.sizeZ);
        }
    }
//    std::cout << obstacleGrid.displayValues() << std::endl;
//    std::cout << obstacleGrid.displayAsPlot() << std::endl;
    this->recomputeVelocities();
}

void WarpedFluidSimulation::recomputeVelocities()
{
    int width = this->dimensions.x;
    int height = this->dimensions.y;
    velocities = GridV3(width, height);

    Vector3 wind = mainDirection;
    GridF augmentedObstacles = obstacleGrid * 1.f;

    std::vector<int> gaussRadii = {3, 5};
    std::vector<float> warpCoefs = {.8f, 1.2f};
    std::vector<float> deviationCoefs = {300.f, 50.f};
    std::vector<GridV3> gradients(gaussRadii.size());
    for (size_t i = 0; i < gradients.size(); i++) {
        auto& grad = gradients[i];
        float radius = gaussRadii[i];
        auto smoothed = augmentedObstacles.meanSmooth(radius, radius, 1, true);
        grad = -smoothed.gradient();

        for (int x = 0; x < grad.sizeX; x++) {
            for (int y = 0; y < grad.sizeY; y++) {
                if (x <= radius || y <= radius || x >= grad.sizeX - radius || y >= grad.sizeY - radius) {
                    grad(x, y) = grad(x + radius, y + radius);
                }
            }
        }
    }

    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            Vector3 vP = wind * (1.f + 5 * float(std::pow(10.f, -3)) * augmentedObstacles(x, y));
            Vector3 sumWivP;

            for (size_t i = 0; i < gradients.size(); i++) {
                auto& gradient1 = gradients[i](x, y);
                float al1 = gradient1.length();
                Vector3 orthGrad1 = (Vector3(0, 0, 1).cross(gradient1)).normalized();
                if (vP.dot(orthGrad1) <= 0)
                {
                    orthGrad1 = orthGrad1 * -1.f;
                }

                sumWivP += warpCoefs[i] * ((1 - al1) * vP + al1 * deviationCoefs[i] * orthGrad1);
            }
            this->velocities(x, y) = sumWivP;
        }
    }
    this->_cachedStep++;
}

void WarpedFluidSimulation::addObstacles(const GridF &obstacles)
{
    GridF addedObstacleGrid;
    addedObstacleGrid = GridF(obstacles.sizeX, obstacles.sizeY);
    for (int x = 0; x < obstacles.sizeX; x++) {
        for (int y = 0; y < obstacles.sizeY; y++) {
            float height = 0;
            while(obstacles.at(x, y, obstacles.sizeZ - height) <= 0 && height < obstacles.sizeZ)
                height++;
            addedObstacleGrid.at(x, y) = height / float(obstacles.sizeZ);
        }
    }
    return this->setObstacles(this->obstacleGrid + addedObstacleGrid);
}
