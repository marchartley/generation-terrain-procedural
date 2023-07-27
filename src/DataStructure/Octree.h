#ifndef OCTREE_H
#define OCTREE_H

#include <vector>
#include "DataStructure/Vector3.h"
#include "DataStructure/Triangle.h"

struct OctreeNodeData {
    Vector3 vertex1;
    Vector3 vertex2;
    Vector3 vertex3;
    int index;

    OctreeNodeData(const Vector3& v1, const Vector3& v2, const Vector3& v3, int idx)
    : vertex1(v1), vertex2(v2), vertex3(v3), index(idx)
    {}
    OctreeNodeData(const Triangle& tri, int idx)
    : OctreeNodeData(tri[0], tri[1], tri[2], idx)
    {}
};

class OctreeNode {
public:
    Vector3 origin; // origin of the node
    Vector3 halfDimension; // half the width/height/depth of the node
    float sphereSqrRadius;
    std::vector<OctreeNodeData> data; // triangles contained in this node
    OctreeNode* children[8]; // eight children of the node
    int maxDataCapacity = 20;

    OctreeNode(const Vector3& origin, const Vector3& halfDimension);

    ~OctreeNode();

    bool intersects(const Vector3& start, const Vector3& end) const;
    bool intersects(const Vector3& v1, const Vector3& v2, const Vector3& v3) const;
    bool intersects(const Triangle& tri) const;

    bool sphereIntersectsTriangle(const Triangle& triangle) const;
};

class Octree {
public:
    OctreeNode* root = nullptr;

    Octree();
    Octree(const Vector3& origin, const Vector3& halfDimension);
    Octree(const std::vector<Triangle>& triangles);
    Octree(const std::vector<std::vector<Vector3>>& triangles);


    ~Octree();

    bool insert(const Vector3& p1, const Vector3& p2, const Vector3& p3, const int& pointIndex);
    bool insert(std::vector<Triangle> triangles);
    bool insert(std::vector<std::vector<Vector3>> triangles);

    std::vector<OctreeNodeData> queryRange(const Vector3& start, const Vector3& end) const;

    std::pair<Vector3, int> _intersectingTriangleIndex(const Vector3& start, const Vector3& end, const std::vector<Triangle > &triangles) const;
    Vector3 getIntersection(const Vector3& start, const Vector3& end, const std::vector<Triangle > &triangles) const;
    std::pair<Vector3, Vector3> getIntersectionAndNormal(const Vector3& start, const Vector3& end, const std::vector<Triangle > &triangles) const;

    Vector3 getIntersection(const Vector3& start, const Vector3& end, const std::vector<std::vector<Vector3>> &triangles) const;
    std::pair<Vector3, Vector3> getIntersectionAndNormal(const Vector3& start, const Vector3& end, const std::vector<std::vector<Vector3>> &triangles) const;

private:
    bool insert(OctreeNode* node, const Vector3& p1, const Vector3& p2, const Vector3& p3, const int &pointIndex);

    void queryRange(OctreeNode* node, const Vector3& start, const Vector3& end, std::vector<OctreeNodeData> &result, int level = 0) const;
};

#endif // OCTREE_H
