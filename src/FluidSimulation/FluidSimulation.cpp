#include "FluidSimulation.h"


FluidSimulation::FluidSimulation()
    : FluidSimulation(0, 0, 0)
{

}

FluidSimulation::FluidSimulation(int sizeX, int sizeY, int sizeZ)
    : dimensions(sizeX, sizeY, sizeZ), obstacleGrid(Matrix3<int>(sizeX, sizeY, sizeZ)), obstacleGradient(Matrix3<Vector3>(sizeX, sizeY, sizeZ))
{
    this->obstacleTrianglesOctree = nullptr;
}

Matrix3<Vector3> FluidSimulation::getVelocities(const Vector3& dimensions)
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
}

void FluidSimulation::setObstacles(const Matrix3<float> &obstacle) {
    this->obstacleGrid = obstacle.resize(dimensions).binarize(0.5);;
    this->obstacleGradient = obstacleGrid.gradient().normalized();
}

void FluidSimulation::addObstacles(const std::vector<std::vector<Vector3> > &triangles)
{
    this->setObstacles(vectorUnion(triangles, this->triangles));
}

void FluidSimulation::addObstacles(const Matrix3<float>& obstacle)
{
    this->setObstacles(this->obstacleGrid + obstacle);
}




KDNode::KDNode(size_t pIndex, int a) : particleIndex(pIndex), left(NULL), right(NULL), axis(a) {}


KDTree::KDTree()
{

}

KDTree::~KDTree()
{
    if (this->root)
        delete this->root;
}

KDTree::KDTree(std::vector<Particle> &particles) {
    root = build(particles, 0);
}

KDNode* KDTree::build(std::vector<Particle> particles, int depth) {
    if (particles.empty()) {
        return NULL;
    }

    int axis = depth % 3;
    std::sort(particles.begin(), particles.end(), [axis](Particle a, Particle b) {
        return a.position[axis] < b.position[axis];
    });

    int median = particles.size() / 2;
    KDNode* node = new KDNode(particles[median].index, axis);

    std::vector<Particle> left(particles.begin(), particles.begin() + median);
    std::vector<Particle> right(particles.begin() + median + 1, particles.end());

    node->left = build(left, depth + 1);
    node->right = build(right, depth + 1);

    return node;
}

void KDTree::findNeighbors(std::vector<Particle> &particles, KDNode *node, Vector3 &position, float maxDistance, std::vector<size_t> &neighbors) {
    if (node == NULL) {
        return;
    }

    float d = position[node->axis] - particles[node->particleIndex].position[node->axis];
    KDNode* nearest = d < 0 ? node->left : node->right;
    KDNode* farthest = d < 0 ? node->right : node->left;

    findNeighbors(particles, nearest, position, maxDistance, neighbors);

    if (d * d < maxDistance * maxDistance) {
        if ((particles[node->particleIndex].position - position).magnitude() < maxDistance) {
            neighbors.push_back(node->particleIndex);
        }

        findNeighbors(particles, farthest, position, maxDistance, neighbors);
    }
}
