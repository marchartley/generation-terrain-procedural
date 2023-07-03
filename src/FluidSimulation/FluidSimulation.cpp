#include "FluidSimulation.h"


void FluidSimulation::setObstacles(const std::vector<std::vector<Vector3> > &triangles) {
    this->triangles = triangles;
    // Create Octree
    Vector3 halfDimension(sizeX / 2.f, sizeY / 2.f, sizeZ / 2.f);
    Vector3 origin = halfDimension; //(0, 0, 0);
    obstacleTrianglesOctree = new Octree(origin, halfDimension);

    // Insert triangles into Octree
    for (size_t iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        const auto& triangle = triangles[iTriangle];
        for (const auto& vertex : triangle) {
            obstacleTrianglesOctree->insert(vertex, iTriangle);
        }
    }
}

void FluidSimulation::setObstacles(Matrix3<float> obstacle) {
    this->obstacleGrid = obstacle.resize(sizeX, sizeY, sizeZ).binarize(0.5);;
    this->obstacleGradient = obstacleGrid.gradient().normalized();
}

void FluidSimulation::addObstacles(const std::vector<std::vector<Vector3> > &triangles)
{
    this->setObstacles(vectorUnion(triangles, this->triangles));
}

void FluidSimulation::addObstacles(Matrix3<float> obstacle)
{
    this->setObstacles(this->obstacleGrid + obstacle);
}
