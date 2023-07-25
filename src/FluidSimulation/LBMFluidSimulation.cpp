#include "LBMFluidSimulation.h"

LBMFluidSimulation::LBMFluidSimulation(bool uses3D)
    : FluidSimulation((uses3D ? Vector3(30, 30, 50) : Vector3(100, 100, 1))), uses3D(uses3D)
{
    c = (uses3D ? c3D : c2D);
//    Matrix3<float> heights(dimensions);
    f = Matrix3<Matrix3<float>>(dimensions, Matrix3<float>(c.size(), 1));
    f_next = Matrix3<Matrix3<float>>(dimensions, Matrix3<float>(c.size(), 1));

    f.raiseErrorOnBadCoord = false;
    f.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
    f_next.raiseErrorOnBadCoord = false;
    f_next.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
//    for (int y = 0; y < f.sizeY; y++) {
//        for (int z = 0; z < f.sizeZ; z++) {
//            this->setVelocity(0, y, z, Vector3(.1f, 0, -.1f));
//        }
//    }
}


//LBMFluidSimulation(const Matrix3<float>& heights) :
//    f(heights.sizeX, heights.sizeY, Matrix3<Vector3>(9)),
//    f_next(heights.sizeX, heights.sizeY, Matrix3<Vector3>(9)) {}

Vector3 LBMFluidSimulation::computeMacroscopicVelocity(int x, int y, int z) {
    Vector3 u(0, 0, 0);
    float rho = 0;

    for (size_t i = 0; i < c.size(); i++) {
        rho += f.at(x, y, z)[i];  // Sum up the distribution functions
        u = u + c[i] * f.at(x, y, z)[i];  // Sum up the product of each distribution function and its velocity direction
    }

    if (rho != 0)
        u = u / rho;  // Divide by the sum of the distribution functions

    return u;
}

float LBMFluidSimulation::computeEquilibrium(const Vector3& c_i, const Vector3& u, float rho, float w_i) {
    float u_sq = u.norm2();
    float u_dot_c_i = u.dot(c_i);

    float f_eq_i = w_i * rho * (1 + 3*u_dot_c_i + 4.5*u_dot_c_i*u_dot_c_i - 1.5*u_sq);

    return f_eq_i;
}

void LBMFluidSimulation::handleCollisions() {
    if (uses3D) {
        for (int x = 0; x < f.sizeX; x++) {
            for (int y = 0; y < f.sizeY; y++) {
                for (int z = 0; z < f.sizeZ; z++) {
                    Vector3 segmentStart(x, y, z);
                    if (obstacleGrid.at(segmentStart)) {
                        for (int i = 0; i < c.size(); i++)
                            f.at(segmentStart)[i] = 0;
                    } else if (obstacleTrianglesOctree) {
                        for (int i = 0; i < c.size(); i++) {
                            Vector3 segmentEnd = segmentStart + c[i];
                            Vector3 intersection(false);
                            Vector3 normal(false);

                            Vector3 offset(0, 0, 0);
                            std::vector<OctreeNodeData> nearbyTriangles = obstacleTrianglesOctree->queryRange(segmentStart + offset, segmentEnd + offset);
                            // Check for intersections with nearby triangles
                            for (auto& triangleData : nearbyTriangles) {
                                auto& triangle = this->triangles[triangleData.index];
                                intersection = Collision::segmentToTriangleCollision(segmentStart, segmentEnd, triangle[0], triangle[1], triangle[2]);
                                if (intersection.isValid()) {
                                    normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalize();
        //                            std::cout << segmentStart << " - " << segmentEnd << " intersects at " << intersection << std::endl;
                                    break;
                                }
                            }

                            if (intersection.isValid()) {
        //                        f.at(x, y, z)[i] = 0.f;

                                // Reflect the velocity direction
                                Vector3 reflected = c[i] - 2 * (c[i].dot(normal)) * normal;

                                // Find the closest lattice direction to the reflected direction
                                int closest = 0;
                                float minAngle = std::numeric_limits<float>::max();
                                for (int j = 0; j < c.size(); j++) {
                                    float angle = std::acos(c[j].normalize().dot(reflected.normalize()));
                                    if (angle < minAngle) {
                                        minAngle = angle;
                                        closest = j;
                                    }
                                }

                                // Swap the distribution functions
        //                        std::swap(f.at(x, y, z)[i], f.at(x, y, z)[closest]);
                                f.at(x, y, z)[closest] += f.at(x, y, z)[i];
                                f.at(x, y, z)[i] = 0;
                            }
                        }
                    }
                    for (int i = 0; i < c.size(); i++) {
                        Vector3 segmentEnd = segmentStart + c[i];
                        if (!Vector3::isInBox(segmentEnd, Vector3(), this->dimensions))
                            f.at(segmentStart)[i] = 0;
                    }
                }
            }
        }
    } else {
        for (int x = 0; x < f.sizeX; x++) {
            for (int y = 0; y < f.sizeY; y++) {
                for (int i = 0; i < c.size(); i++) {
                    f.at(x, y)[i] = std::min(f.at(x, y)[i], obstacleGrid.at(Vector3(x, y) + c[i]));
                }
            }
        }
    }
}

void LBMFluidSimulation::stream() {
    // Create a temporary copy of the distribution functions
    Matrix3<Matrix3<float>> f_temp = f;
    Matrix3<Matrix3<float>> counts(f.getDimensions(), f[0].getDimensions());
    counts.raiseErrorOnBadCoord = false;

    for (int i = 0; i < f.size(); i++)
        f[i].reset();

    // Propagate the distribution functions to their neighboring cells
    for (int x = 0; x < f.sizeX; x++) {
        for (int y = 0; y < f.sizeY; y++) {
            for (int z = 0; z < f.sizeZ; z++) {
                for (int i = 0; i < c.size(); i++) {
                    int x_next = x + c[i].x;
                    int y_next = y + c[i].y;
                    int z_next = z + c[i].z;
                    if (x_next >= 0 && x_next < f.sizeX && y && y_next >= 0 && y_next < f.sizeY && z_next >= 0 && z_next < f.sizeZ) {
                        f.at(x_next, y_next, z_next)[i] += f_temp.at(x, y, z)[i];
                        counts.at(x_next, y_next, z_next)[i]++;
                    }
                }
            }
        }
    }
    for (int i = 0; i < f.size(); i++)
        for (int j = 0; j < f[i].size(); j++)
            if (counts[i][j] > 0)
                f[i][j] /= counts[i][j];
}
void LBMFluidSimulation::step() {
    currentStep++;
    float tau = 1000.f;
    // Compute the macroscopic velocity and the equilibrium distribution function at each cell
#pragma omp parallel for collapse(3)
    for (int x = 0; x < f.sizeX; x++) {
        for (int y = 0; y < f.sizeY; y++) {
            for (int z = 0; z < f.sizeZ; z++) {
                Vector3 u = computeMacroscopicVelocity(x, y, z);
                float rho = 0;
                for (int i = 0; i < c.size(); i++) {
                    rho += f.at(x, y, z)[i];
                }
                for (int i = 0; i < c.size(); i++) {
                    float w_i = getWi(i);
                    float f_eq_i = computeEquilibrium(c[i], u, rho, w_i);
                    f_next.at(x, y, z)[i] = f.at(x, y, z)[i] - (1.0 / tau) * (f.at(x, y, z)[i] - f_eq_i);
                }
            }
        }
    }
    f = f_next;

    // Handle interactions with obstacles
    handleCollisions();

    // Perform the streaming step
    stream();

    float minVal = std::numeric_limits<float>::max(), maxVal = std::numeric_limits<float>::min();
    for (int i = 0; i < f.size(); i++) {
        for (int j = 0; j < f[i].size(); j++) {
            minVal = std::min(minVal, std::abs(f[i][j]));
            maxVal = std::max(maxVal, std::abs(f[i][j]));
            if (f[i][j] != f[i][j])
                std::cout << "NaN found!" << std::endl;
        }
    }
    std::cout << minVal << " " << maxVal << std::endl;
}

Matrix3<Vector3> LBMFluidSimulation::getVelocities(int newSizeX, int newSizeY, int newSizeZ)
{
    if (_cachedStep != currentStep) {
        _cachedStep = currentStep;
        Matrix3<Vector3> velocities(this->f.getDimensions());
        for (int x = 0; x < dimensions.x; x++) {
            for (int y = 0; y < dimensions.y; y++) {
                for (int z = 0; z < dimensions.z; z++) {
                    /*Vector3 u(0, 0, 0);
                    float rho = 0;

                    for (int i = 0; i < c.size(); i++) {
                        rho += f.at(x, y)[i];  // Sum up the distribution functions
                        u = u + c[i] * f.at(x, y)[i];  // Sum up the product of each distribution function and its velocity direction
                    }

                    u = u / rho;  // Divide by the sum of the distribution functions
                    velocities.at(x, y) = u;*/
                    velocities.at(x, y, z) = this->getVelocity(x, y, z);
                }
            }
        }
        _cachedVelocity = velocities;
    }
    return _cachedVelocity.resize(newSizeX, newSizeY, newSizeZ);
}

Vector3 LBMFluidSimulation::getVelocity(int x, int y, int z)
{
    return this->computeMacroscopicVelocity(x, y, z);
}

void LBMFluidSimulation::addVelocity(int x, int y, int z, const Vector3& amount)
{
    float rho = 1.0;  // Assume a constant density
    for (int i = 0; i < c.size(); i++) {
        float w_i = getWi(i);
        f.at(x, y, z)[i] += computeEquilibrium(c[i], amount, rho, w_i);
    }
}

void LBMFluidSimulation::setVelocity(int x, int y, int z, const Vector3 &amount)
{
    float rho = 1.0;  // Assume a constant density
    for (int i = 0; i < c.size(); i++) {
        float w_i = getWi(i);
        f.at(x, y, z)[i] = computeEquilibrium(c[i], amount, rho, w_i);
    }
}

void LBMFluidSimulation::setObstacles(const Matrix3<float> &obstacles)
{
    if (uses3D) {
        return FluidSimulation::setObstacles(obstacles);
    } else {
        this->obstacleGrid = Matrix3<float>(obstacles.sizeX, obstacles.sizeY);
        for (int x = 0; x < obstacles.sizeX; x++) {
            for (int y = 0; y < obstacles.sizeY; y++) {
                int height = 0;
                while(obstacles.at(x, y, obstacles.sizeZ - height) <= 0 && height < obstacles.sizeZ)
                    height++;
                obstacleGrid.at(x, y) = height / this->dimensions.z;
            }
        }
    }
}

void LBMFluidSimulation::addObstacles(const Matrix3<float> &obstacles)
{
    Matrix3<float> addedObstacleGrid;
    if (uses3D) {
        addedObstacleGrid = obstacles;
    } else {
        addedObstacleGrid = Matrix3<float>(obstacles.sizeX, obstacles.sizeY);
        for (int x = 0; x < obstacles.sizeX; x++) {
            for (int y = 0; y < obstacles.sizeY; y++) {
                int height = 0;
                while(obstacles.at(x, y, obstacles.sizeZ - height) <= 0 && height < obstacles.sizeZ)
                    height++;
                addedObstacleGrid.at(x, y) = height;
            }
        }
    }
    return this->setObstacles(this->obstacleGrid + addedObstacleGrid);
}

float LBMFluidSimulation::getWi(int i)
{
    if (uses3D)
        return (i == 0) ? 1.f/3.f : ((i < 7) ? 1.f/18.f : 1.f/36.f);
    return (i == 0) ? 4.f/9.f : ((i < 5) ? 1.f/9.f : 1.f/36.f);
}
