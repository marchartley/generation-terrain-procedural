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
    std::cout << obstacleGrid.displayValues() << std::endl;
    std::cout << obstacleGrid.displayAsPlot() << std::endl;
    this->recomputeVelocities();
}

void WarpedFluidSimulation::recomputeVelocities()
{
    int width = this->dimensions.x;
    int height = this->dimensions.y;
    velocities = GridV3(width, height);

    Vector3 wind = mainDirection;

//    int rad1 = 3;
//    int rad2 = 5;
    GridF augmentedObstacles = obstacleGrid * 10.f;
//    GridF obstacle1 = augmentedObstacles.meanSmooth(rad1, rad1, 1, false);
//    GridF obstacle2 = augmentedObstacles.meanSmooth(rad2, rad2, 1, false);
//    GridV3 grad1 = -obstacle1.gradient();
//    GridV3 grad2 = -obstacle2.gradient();

//    for (auto& v : grad1)
//        v.normalize();
//    for (auto& v : grad2)
//        v.normalize();

    std::vector<int> gaussRadii = {3, 5};
    std::vector<GridV3> gradients(gaussRadii.size());
    for (size_t i = 0; i < gradients.size(); i++) {
        auto& grad = gradients[i];
        grad = -(augmentedObstacles.meanSmooth(gaussRadii[i], gaussRadii[i], 1, false).gradient());
        for (auto& v : grad)
            v.normalize();
    }

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            Vector3 vP = wind * (1.f + 5 * float(std::pow(10.f, -3)) * augmentedObstacles(i, j));
            Vector3 sumWivP;

            for (size_t iGauss = 0; iGauss < gradients.size(); iGauss++) {
                auto& gradient1 = gradients[iGauss](i, j);
                float al1 = gradient1.length();
                float c1 = 0.8f;
//                float c2 = 0.2f;
                float kt1 = 30.f;
//                float kt2 = 5.f;

                Vector3 orthGrad1;

                float v11 = 1.f;
                if (gradient1 == Vector3(0.f))
                    orthGrad1 = Vector3(0.f);
                else
                {
                    if(gradient1.y == 0.f)
                        orthGrad1 = Vector3(v11, 0.f).normalize() * al1;
                    else
                    {
                        float v21 = (-gradient1.x / gradient1.y);
                        orthGrad1 = Vector3(v11, v21).normalize() * al1;
                    }
                }

                if (vP.dot(orthGrad1) <= 0)
                {
                    orthGrad1 = orthGrad1 * -1.f;
                }

                sumWivP += c1 * ((1 - al1) * vP + al1 * kt1 * orthGrad1);
            }

//            auto& gradient1 = grad1(i, j);
//            auto& gradient2 = grad2(i, j);
//            float al1 = gradient1.length();
//            float al2 = gradient2.length();

//            float c1 = 0.8f;
//            float c2 = 0.2f;
//            float kt1 = 30.f;
//            float kt2 = 5.f;

//            Vector3 orthGrad1, orthGrad2;

//            float v11 = 1.f;
//            if (gradient1 == Vector3(0.f))
//                orthGrad1 = Vector3(0.f);
//            else
//            {
//                if(gradient1.y == 0.f)
//                    orthGrad1 = Vector3(v11, 0.f).normalize() * al1;
//                else
//                {
//                    float v21 = (-gradient1.x / gradient1.y);
//                    orthGrad1 = Vector3(v11, v21).normalize() * al1;
//                }
//            }

//            if (vP.dot(orthGrad1) <= 0)
//            {
//                orthGrad1 = orthGrad1 * -1.f;
//            }

//            float v12 = 1.f;
//            if (gradient2 == Vector3(0.f))
//                orthGrad2 = Vector3(0.f);
//            else
//            {
//                if (gradient2.y == 0.f)
//                    orthGrad2 = Vector3(v12, 0.f).normalize() * al2;
//                else
//                {
//                    float v22 = (-gradient2.x / gradient2.y);
//                    orthGrad2 = Vector3(v12, v22).normalize() * al2;
//                }
//            }

//            if (vP.dot(orthGrad2) <= 0)
//            {
//                orthGrad2 = orthGrad2 * -1.f;
//            }

//            Vector3 sumWivP = c1 * ((1 - al1) * vP + al1 * kt1 * orthGrad1) + c2 * ((1 - al2) * vP + al2 * kt2 * orthGrad2);
            this->velocities(i, j) = sumWivP;
        }
    }
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
