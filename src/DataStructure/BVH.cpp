#include "BVH.h"
#include "Utils/Collisions.h"


BVHTree::BVHTree() : root(nullptr)
{

}

BVHTree::BVHTree(std::vector<std::vector<Vector3> > triangles)
{
    this->build(triangles);
}

BVHTree &BVHTree::build(std::vector<std::vector<Vector3> > &triangles)
{
    this->root = this->buildBVH(triangles, 0, triangles.size());
    return *this;
}

std::vector<Vector3> BVHTree::query(const Vector3 &rayStart, const Vector3 &rayEnd)
{
    return this->BVHQuery(this->root, rayStart, rayEnd);
}

BVHNode* BVHTree::buildBVH(std::vector<std::vector<Vector3>>& triangles, int start, int end) {
    BVHNode* node = new BVHNode();

    // Compute bounding box of all triangles from start to end
    Vector3 minPoint = Vector3::max();
    Vector3 maxPoint = Vector3::min();
    for (int i = start; i < end; ++i) {
        for (const auto& vertex : triangles[i]) {
            minPoint = Vector3::min(minPoint, vertex);
            maxPoint = Vector3::max(maxPoint, vertex);
        }
    }
    node->box = AABBox(minPoint, maxPoint);

    // Base case: if 2 or fewer triangles, make a leaf node
    if (end - start <= 2) {
        node->triangles = std::vector<std::vector<Vector3>>(triangles.begin() + start, triangles.begin() + end);
    } else { // Recursive case: partition the triangles and build child nodes
        // Choose an axis and midpoint along that axis to partition the triangles
        Vector3 boxSize = maxPoint - minPoint;
        int axis = boxSize.x > boxSize.y ? (boxSize.x > boxSize.z ? 0 : 2) : (boxSize.y > boxSize.z ? 1 : 2);

        // Sort the triangles based on their midpoint along the chosen axis
        std::sort(triangles.begin() + start, triangles.begin() + end,
                  [axis](const std::vector<Vector3>& triangle1, const std::vector<Vector3>& triangle2) {
                      Vector3 midPoint1 = (triangle1[0] + triangle1[1] + triangle1[2]) / 3;
                      Vector3 midPoint2 = (triangle2[0] + triangle2[1] + triangle2[2]) / 3;
                      return midPoint1[axis] < midPoint2[axis];
                  });

        int mid = start + (end - start) / 2; // Use median splitting
        node->left = buildBVH(triangles, start, mid);
        node->right = buildBVH(triangles, mid, end);
    }

    return node;
}

std::vector<Vector3> BVHTree::BVHQuery(BVHNode* node, const Vector3& rayStart, const Vector3& rayEnd) {
    std::vector<Vector3> intersections;

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return intersections; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->triangles.empty()) { // leaf node
        for (const auto& triangle : node->triangles) {
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid()) {
                intersections.push_back(intersectionPoint);
            }
        }
    } else { // non-leaf node
        std::vector<Vector3> leftIntersections = BVHQuery(node->left, rayStart, rayEnd);
        std::vector<Vector3> rightIntersections = BVHQuery(node->right, rayStart, rayEnd);
        intersections.insert(intersections.end(), leftIntersections.begin(), leftIntersections.end());
        intersections.insert(intersections.end(), rightIntersections.begin(), rightIntersections.end());
    }

    return intersections;
}
