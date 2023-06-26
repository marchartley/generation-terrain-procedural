#ifndef OCTREE_H
#define OCTREE_H

#include <vector>
#include "DataStructure/Vector3.h"

struct OctreeNodeData {
    Vector3 pos;
    int index;
};

class OctreeNode {
public:
    Vector3 origin; // origin of the node
    Vector3 halfDimension; // half the width/height/depth of the node
    std::vector<OctreeNodeData> data; // triangles contained in this node
    OctreeNode* children[8]; // eight children of the node

    OctreeNode(const Vector3& origin, const Vector3& halfDimension);

    ~OctreeNode();

    bool intersects(const Vector3& start, const Vector3& end) const;
};

class Octree {
public:
    OctreeNode* root;

    Octree();
    Octree(const Vector3& origin, const Vector3& halfDimension);

    ~Octree();

    void insert(const Vector3& point, const int& pointIndex);

    std::vector<OctreeNodeData> queryRange(const Vector3& start, const Vector3& end);

private:
    void insert(OctreeNode* node, const Vector3& point, const int &pointIndex);

    void queryRange(OctreeNode* node, const Vector3& start, const Vector3& end, std::vector<OctreeNodeData> &result);
};

#endif // OCTREE_H
