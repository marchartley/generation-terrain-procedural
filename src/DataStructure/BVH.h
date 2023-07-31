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
//    BVHTree(std::vector<std::vector<Vector3>> triangles);
    ~BVHTree();

//    BVHTree& build(const std::vector<std::vector<Vector3>>& triangles);
    virtual SpacePartitioning& build(const std::vector<Triangle>& triangles);
//    std::vector<Vector3> query(const Vector3& rayStart, const Vector3& rayEnd) const;
//    Vector3 getIntersection(const Vector3& rayStart, const Vector3& rayEnd) const;
//    std::pair<Vector3, Vector3> getIntersectionAndNormal(const Vector3& rayStart, const Vector3& rayEnd) const;
    BVHNode* root = nullptr;

    virtual std::set<size_t> getAllStoredTrianglesIndices() const;

protected:
//    BVHNode* buildBVH(std::vector<std::vector<Vector3>>& triangles, int start, int end);
    BVHNode* buildBVH(int start, int end);
//    std::vector<Vector3> BVHQuery(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd);

//    using SpacePartitioning::getIntersection;
//    using SpacePartitioning::getIntersectionAndNormal;
//    Vector3 _getIntersection(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;
//    std::pair<Vector3, Vector3> _getIntersectionAndNormal(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;
    virtual std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd) const;
    std::pair<Vector3, size_t> _getIntersectionAndTriangleIndex(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;

    virtual std::vector<std::pair<Vector3, size_t>> getAllIntersectionsAndTrianglesIndices(const Vector3& rayStart, const Vector3& rayEnd) const;
    std::vector<std::pair<Vector3, size_t>> _getAllIntersectionsAndTrianglesIndices(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;

};

#endif // BVH_H
