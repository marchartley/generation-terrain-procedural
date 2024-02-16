#include "KDTree.h"


KDNode::KDNode(size_t pIndex, int a) : particleIndex(pIndex), left(nullptr), right(nullptr), axis(a) {}

KDNode::~KDNode()
{
    if (this->left != nullptr)
        delete this->left;
    if (this->right != nullptr)
        delete this->right;
}


KDTree::KDTree()
{

}

KDTree::~KDTree()
{
    if (this->root != nullptr)
        delete this->root;
}

KDTree::KDTree(std::vector<Particle> &particles) {
    root = build(particles, 0);
}

KDNode* KDTree::build(std::vector<Particle> particles, int depth) {
    if (particles.empty()) {
        return nullptr;
    }
    if (depth == 0) {
        delete this->root;
        for (int i = 0; i < particles.size(); i++)
            particles[i].index = i;
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

    if (depth == 0)
        this->root = node;

    return node;
}

std::vector<size_t> KDTree::findNeighbors(std::vector<Particle> &particles, const Vector3 &position, float maxDistance) const
{
    std::vector<size_t> indices;
    this->findNeighbors(particles, root, position, maxDistance, indices);
    return indices;
}

void KDTree::findNeighbors(std::vector<Particle> &particles, KDNode *node, const Vector3& position, float maxDistance, std::vector<size_t> &neighbors) const
{
    if (node == nullptr) {
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
