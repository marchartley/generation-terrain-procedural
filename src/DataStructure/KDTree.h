#ifndef KDTREE_H
#define KDTREE_H

#include "DataStructure/Particle.h"


class KDNode {
public:
    size_t particleIndex;
    KDNode* left = nullptr;
    KDNode* right = nullptr;
    int axis;

    KDNode(size_t pIndex, int a);
    ~KDNode();
};

class KDTree {
public:
    KDTree();
    ~KDTree();
    KDNode* root = nullptr;

    KDTree(std::vector<Particle>& particles);

    KDNode* build(std::vector<Particle> particles, int depth = 0);

    std::vector<size_t> findNeighbors(std::vector<Particle> &particles, const Vector3 &position, float maxDistance) const;
    void findNeighbors(std::vector<Particle>& particles, KDNode* node, const Vector3& position, float maxDistance, std::vector<size_t>& neighbors) const;
};

#endif // KDTREE_H
