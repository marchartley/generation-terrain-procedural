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
    KDNode(const KDNode &) = delete;
    KDNode &operator=(const KDNode &) = delete;
};

class KDTree {
public:
    KDTree();
    ~KDTree();
    KDTree(const KDTree &other);
    KDTree &operator=(const KDTree &other);
    KDTree(KDTree &&other) noexcept;
    KDTree &operator=(KDTree &&other) noexcept;
    KDNode* root = nullptr;

    KDTree(std::vector<Particle>& particles);

    // Queries use immutable build-time coordinates. Rebuild after particle motion.
    std::vector<size_t> findNeighbors(const Vector3 &position, float maxDistance) const;
    KDNode* build(std::vector<Particle> particles, int depth = 0);

    std::vector<size_t> findNeighbors(std::vector<Particle> &particles, const Vector3 &position, float maxDistance) const;
    void findNeighbors(std::vector<Particle>& particles, KDNode* node, const Vector3& position, float maxDistance, std::vector<size_t>& neighbors) const;

  private:
    std::vector<Vector3> positions;
    void query(KDNode *node, const Vector3 &p, double radiusSquared, std::vector<size_t> &result) const;
};

#endif // KDTREE_H
