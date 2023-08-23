#include "FluidSimulation.h"


FluidSimulation::FluidSimulation()
    : FluidSimulation(0, 0, 0)
{

}

FluidSimulation::FluidSimulation(int sizeX, int sizeY, int sizeZ)
    : dimensions(sizeX, sizeY, sizeZ), obstacleGrid(GridI(sizeX, sizeY, sizeZ)), obstacleGradient(GridV3(sizeX, sizeY, sizeZ))
{
//    this->obstacleTrianglesOctree = nullptr;
}

GridV3 FluidSimulation::getVelocities(const Vector3& dimensions)
{
    return this->getVelocities(dimensions.x, dimensions.y, dimensions.z);
}

Vector3 FluidSimulation::getVelocity(const Vector3 &pos)
{
    return this->getVelocity(pos.x, pos.y, pos.z);
}
Vector3 FluidSimulation::getVelocity(int x, int y, int z)
{
    return this->getVelocities(this->dimensions).at(x, y, z);
}

void FluidSimulation::addVelocity(const Vector3 &pos, const Vector3 &amount)
{
    return this->addVelocity(pos.x, pos.y, pos.z, amount);
}

void FluidSimulation::setVelocity(const Vector3 &pos, const Vector3 &amount)
{
    return this->setVelocity(pos.x, pos.y, pos.z, amount);
}

void FluidSimulation::setVelocity(int x, int y, int z, const Vector3 &amount)
{
    this->addVelocity(x, y, z, (amount - this->getVelocity(x, y, z)));
}

void FluidSimulation::setObstacles(const std::vector<std::vector<Vector3> > &triangles) {
    this->triangles = triangles;
//    obstacleTriangleTree = BVHTree(); //(triangles);
    obstacleTriangleTree.build(Triangle::vectorsToTriangles(triangles));
    /*
    // Create Octree
    Vector3 halfDimension(dimensions * .5f);
    Vector3 origin = halfDimension; //(0, 0, 0);
    obstacleTrianglesOctree = new Octree(origin, halfDimension);

    // Insert triangles into Octree
    for (size_t iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        const auto& triangle = triangles[iTriangle];
        obstacleTrianglesOctree->insert(triangle[0], triangle[1], triangle[2], iTriangle);
//        for (const auto& vertex : triangle) {
//            obstacleTrianglesOctree->insert(vertex, iTriangle);
//        }
    }
    */
}

void FluidSimulation::setObstacles(const GridF &obstacle) {
    this->obstacleGrid = obstacle.resize(dimensions).binarize();
    this->obstacleGradient = obstacleGrid.gradient().normalized();
}

void FluidSimulation::addObstacles(const std::vector<std::vector<Vector3> > &triangles)
{
    this->setObstacles(vectorUnion(triangles, this->triangles));
}

void FluidSimulation::addObstacles(const GridF& obstacle)
{
    this->setObstacles(this->obstacleGrid + obstacle);
}
