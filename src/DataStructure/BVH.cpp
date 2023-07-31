#include "BVH.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"


BVHTree::BVHTree() : SpacePartitioning(), root(nullptr)
{

}

//BVHTree::BVHTree(std::vector<std::vector<Vector3> > triangles)
//{
//    this->build(triangles);
//}

BVHTree::~BVHTree()
{
    if (root != nullptr)
        delete root;
}

//BVHTree &BVHTree::build(const std::vector<std::vector<Vector3> > &triangles)
//{
//    this->root = this->buildBVH(triangles, 0, triangles.size());
//    return *this;
//}

SpacePartitioning& BVHTree::build(const std::vector<Triangle> &triangles)
{
    this->triangles = triangles;
    this->root = this->buildBVH(0, triangles.size());
    return *this;
}

std::set<size_t> BVHTree::getAllStoredTrianglesIndices() const
{
    std::vector<BVHNode*> queue = {root};
    std::vector<size_t> indices;

    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        indices.insert(indices.end(), current->trianglesIndices.begin(), current->trianglesIndices.end());
        if (current->left)
            queue.push_back(current->left);
        if (current->right)
            queue.push_back(current->right);
    }
    return convertVectorToSet(indices);
}

//std::vector<Vector3> BVHTree::query(const Vector3 &rayStart, const Vector3 &rayEnd) const
//{
//    return this->BVHQuery(this->root, rayStart, rayEnd);
//}
/*
Vector3 BVHTree::getIntersection(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->_getIntersection(root, rayStart, rayEnd);
}

std::pair<Vector3, Vector3> BVHTree::getIntersectionAndNormal(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->_getIntersectionAndNormal(root, rayStart, rayEnd);
}
*/
/*BVHNode* BVHTree::buildBVH(std::vector<std::vector<Vector3>>& triangles, int start, int end) {
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
}*/

BVHNode *BVHTree::buildBVH(/*std::vector<Triangle> &triangles, */int start, int end)
{
    BVHNode* node = new BVHNode();

    // Compute bounding box of all triangles from start to end
    Vector3 minPoint = Vector3::max();
    Vector3 maxPoint = Vector3::min();
    for (int i = start; i < end; ++i) {
        for (int ii = 0; ii < 3; ii++) {
            minPoint = Vector3::min(minPoint, triangles[i][ii]);
            maxPoint = Vector3::max(maxPoint, triangles[i][ii]);
        }
    }
    node->box = AABBox(minPoint, maxPoint);

    // Base case: if 2 or fewer triangles, make a leaf node
    if (end - start <= 2) {
        node->trianglesIndices = std::vector<size_t>(end - start);
        for (int i = 0; i < (end - start); i++)
            node->trianglesIndices[i] = start + i;
    } else { // Recursive case: partition the triangles and build child nodes
        // Choose an axis and midpoint along that axis to partition the triangles
        Vector3 boxSize = maxPoint - minPoint;
        int axis = boxSize.x > boxSize.y ? (boxSize.x > boxSize.z ? 0 : 2) : (boxSize.y > boxSize.z ? 1 : 2);

        // Sort the triangles based on their midpoint along the chosen axis
        std::sort(triangles.begin() + start, triangles.begin() + end,
                  [axis](const Triangle& triangle1, const Triangle& triangle2) {
                      Vector3 midPoint1 = (triangle1[0] + triangle1[1] + triangle1[2]) / 3;
                      Vector3 midPoint2 = (triangle2[0] + triangle2[1] + triangle2[2]) / 3;
                      return midPoint1[axis] < midPoint2[axis];
                  });

        int mid = start + (end - start) / 2; // Use median splitting
        node->left = buildBVH(/*triangles, */start, mid);
        node->right = buildBVH(/*triangles, */mid, end);
    }

    return node;
}
/*
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
}*/

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->_getIntersectionAndTriangleIndex(root, rayStart, rayEnd);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd) const
{
//    Vector3 closestIntersection(false);
//    size_t closestIndex = -1;
    std::pair<Vector3, size_t> result = {Vector3(false), -1};
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                result = {intersectionPoint, triangleIndex};
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getIntersectionAndTriangleIndex(node->left, rayStart, rayEnd);
        auto rightIntersection = _getIntersectionAndTriangleIndex(node->right, rayStart, rayEnd);
        if (!leftIntersection.first.isValid()) {
            result = rightIntersection;
        } else if (!rightIntersection.first.isValid()) {
            result = leftIntersection;
        } else if ((leftIntersection.first - rayStart) < (rightIntersection.first - rayStart)) {
            result = leftIntersection;
        } else {
            result = rightIntersection;
        }
    }

    return result;
}

std::vector<std::pair<Vector3, size_t> > BVHTree::getAllIntersectionsAndTrianglesIndices(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->_getAllIntersectionsAndTrianglesIndices(root, rayStart, rayEnd);
}

std::vector<std::pair<Vector3, size_t> > BVHTree::_getAllIntersectionsAndTrianglesIndices(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    std::vector<std::pair<Vector3, size_t>> result;
//    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid()) {
                result.push_back({intersectionPoint, triangleIndex});
//                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getAllIntersectionsAndTrianglesIndices(node->left, rayStart, rayEnd);
        auto rightIntersection = _getAllIntersectionsAndTrianglesIndices(node->right, rayStart, rayEnd);

        result.insert(result.end(), leftIntersection.begin(), leftIntersection.end());
        result.insert(result.end(), rightIntersection.begin(), rightIntersection.end());
    }

    return result;
}
/*
Vector3 BVHTree::_getIntersection(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    Vector3 closestIntersection(false);
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return closestIntersection; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->triangles.empty()) { // leaf node
        for (const auto& triangle : node->triangles) {
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                closestIntersection = intersectionPoint;
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        Vector3 leftIntersection = _getIntersection(node->left, rayStart, rayEnd);
        Vector3 rightIntersection = _getIntersection(node->right, rayStart, rayEnd);
        if (!leftIntersection.isValid()) {
            closestIntersection = rightIntersection;
        } else if (!rightIntersection.isValid()) {
            closestIntersection = leftIntersection;
        } else if ((leftIntersection - rayStart) < (rightIntersection - rayStart)) {
            closestIntersection = leftIntersection;
        } else {
            closestIntersection = rightIntersection;
        }
    }

    return closestIntersection;
}

std::pair<Vector3, Vector3> BVHTree::_getIntersectionAndNormal(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    Vector3 closestIntersection(false);
    Vector3 normal(false);
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return {closestIntersection, normal}; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->triangles.empty()) { // leaf node
        for (const auto& triangle : node->triangles) {
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                closestIntersection = intersectionPoint;
                normal = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]).normalize();
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto [leftIntersection, leftNormal] = getIntersectionAndNormal(node->left, rayStart, rayEnd);
        auto [rightIntersection, rightNormal] = getIntersectionAndNormal(node->right, rayStart, rayEnd);
        if (!leftIntersection.isValid()) {
            closestIntersection = rightIntersection;
            normal = rightNormal;
        } else if (!rightIntersection.isValid()) {
            closestIntersection = leftIntersection;
            normal = leftNormal;
        } else if ((leftIntersection - rayStart) < (rightIntersection - rayStart)) {
            closestIntersection = leftIntersection;
            normal = leftNormal;
        } else {
            closestIntersection = rightIntersection;
            normal = rightNormal;
        }
    }

    return {closestIntersection, normal};
}
*/
