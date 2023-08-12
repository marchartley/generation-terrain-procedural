#ifndef BVH_H
#define BVH_H

//#include <vector>
#include <set>
//#include "DataStructure/Vector3.h"
#include "DataStructure/SpacePartitioning.h"

struct BVHNode {
//    std::vector<std::vector<Vector3>> triangles;
//    std::vector<Triangle> triangles;
    std::vector<size_t> trianglesIndices;
    AABBox box;
    BVHNode* left;
    BVHNode* right;

    BVHNode() : left(nullptr), right(nullptr) {}
    ~BVHNode() {
        if (left != nullptr)
            delete left;
        if (right != nullptr)
            delete right;
    }
};

class BVHTree : public SpacePartitioning {
public:
    BVHTree();
    ~BVHTree();

    virtual SpacePartitioning& build(const std::vector<Triangle>& triangles);
    BVHNode* root = nullptr;

    virtual std::set<size_t> getAllStoredTrianglesIndices() const;

//protected:
    BVHNode* buildBVH(int start, int end);
    virtual std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;
    std::pair<Vector3, size_t> _getIntersectionAndTriangleIndex(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;
    std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, size_t ignoredTriangle) const;
    std::pair<Vector3, size_t> _getIntersectionAndTriangleIndex(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd, size_t ignoredTriangle) const;

    virtual std::vector<std::pair<Vector3, size_t>> getAllIntersectionsAndTrianglesIndices(const Vector3& rayStart, const Vector3& rayEnd) const;
    std::vector<std::pair<Vector3, size_t>> _getAllIntersectionsAndTrianglesIndices(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;

    int findBestSplitSAH(int start, int end, int axis);
};

#endif // BVH_H
