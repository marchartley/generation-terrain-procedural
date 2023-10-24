#ifndef BVH_H
#define BVH_H

//#include <vector>
#include <set>
//#include "DataStructure/Vector3.h"
#include "DataStructure/SpacePartitioning.h"
#include "DataStructure/MemoryPool.h"

struct BVHNode {
//    std::vector<std::vector<Vector3>> triangles;
//    std::vector<Triangle> triangles;
    std::vector<size_t> trianglesIndices;
    AABBox box;
    BVHNode* left;
    BVHNode* right;

    BVHNode();
    ~BVHNode();
};

using BVHMemoryPool = MemoryPool<BVHNode>;


class BVHTree : public SpacePartitioning {
public:
    BVHTree();
    ~BVHTree();

    BVHTree(const BVHTree& other);
    BVHTree& operator=(const BVHTree& other);


    virtual SpacePartitioning& build(const std::vector<Triangle>& triangles);
    BVHNode* root = nullptr;

    virtual std::set<size_t> getAllStoredTrianglesIndices() const;

    bool checkIntersection(const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;
    bool _checkIntersection(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;

//protected:
    BVHNode* buildBVH(int start, int end);
    virtual std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;
    std::pair<Vector3, size_t> _getIntersectionAndTriangleIndex(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const;
    std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, size_t ignoredTriangle) const;
    std::pair<Vector3, size_t> _getIntersectionAndTriangleIndex(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd, size_t ignoredTriangle) const;

    virtual std::vector<std::pair<Vector3, size_t>> getAllIntersectionsAndTrianglesIndices(const Vector3& rayStart, const Vector3& rayEnd) const;
    std::vector<std::pair<Vector3, size_t>> _getAllIntersectionsAndTrianglesIndices(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) const;

    // SAH optimisation
    int findBestSplitSAH(int start, int end);

    // QuickSelect optimisation
    int partition(int start, int end, int pivotIdx, int axis);
    int quickSelect(int start, int end, int axis);

    int maxTrianglesPerLeaves = 2;
    // That is often useful
    bool useParallel = false;

    // This has worst building and evaluation times than the other methods, so just ignore it!
    bool useSAH = false;
    // This is very efficient for reducing construction time, but not for evaluation time.
    // Set it to true if the building time is more important than evaluation (< 500.000 eval/build)
    bool useQuickSelect = false;

    BVHMemoryPool* memoryPool;
    BVHNode* allocateNode();
};

#endif // BVH_H
