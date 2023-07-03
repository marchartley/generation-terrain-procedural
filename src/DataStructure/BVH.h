#ifndef BVH_H
#define BVH_H

#include <vector>
#include "DataStructure/Vector3.h"

struct BVHNode {
    std::vector<std::vector<Vector3>> triangles;
    AABBox box;
    BVHNode* left;
    BVHNode* right;

    BVHNode() : left(nullptr), right(nullptr) {}
};

class BVHTree {
public:
    BVHTree();
    BVHTree(std::vector<std::vector<Vector3>> triangles);

    BVHTree& build(std::vector<std::vector<Vector3>>& triangles);
    std::vector<Vector3> query(const Vector3& rayStart, const Vector3& rayEnd);

protected:
    BVHNode* buildBVH(std::vector<std::vector<Vector3>>& triangles, int start, int end);
    std::vector<Vector3> BVHQuery(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd);

    BVHNode* root;
};

#endif // BVH_H
