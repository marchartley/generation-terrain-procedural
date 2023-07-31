#ifndef KDTREE_H
#define KDTREE_H

#include "DataStructure/Particle.h"


class KDNode {
public:
    size_t particleIndex;
    KDNode* left;
    KDNode* right;
    int axis;

    KDNode(size_t pIndex, int a);
};

class KDTree {
public:
    KDTree();
    ~KDTree();
    KDNode* root;

    KDTree(std::vector<Particle>& particles);

    KDNode* build(std::vector<Particle> particles, int depth);

    void findNeighbors(std::vector<Particle>& particles, KDNode* node, Vector3& position, float maxDistance, std::vector<size_t>& neighbors);
};

#endif // KDTREE_H
