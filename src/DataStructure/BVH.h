#ifndef BVH_H
#define BVH_H

//#include <vector>
#include <set>
//#include "DataStructure/Vector3.h"
#include "DataStructure/SpacePartitioning.h"
#include "DataStructure/MemoryPool.h"

struct BVHNode {
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

//    static BVHTree aggregate(const BVHTree& A, const BVHTree& B); // Todo one day


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



    void traverseBVH(BVHNode* node, const Vector3& pos, float& minDistance, size_t& closestTriangleIndex, const std::vector<Triangle>& triangles);
    size_t getClosestTriangle(const Vector3& pos/*, const std::vector<Triangle>& triangles*/);
};

#endif // BVH_H
