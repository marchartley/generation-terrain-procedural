#include "BVH.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"


BVHTree::BVHTree() : SpacePartitioning(), root(nullptr)
{

}

BVHTree::~BVHTree()
{
    if (root != nullptr)
        delete root;
}

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

float computeSurfaceArea(const AABBox& box) {
    Vector3 dims = box.dimensions();
    return 2.0f * (dims.x * dims.y + dims.x * dims.z + dims.y * dims.z);
}

int BVHTree::findBestSplitSAH(int start, int end, int axis) {
    // Initially, set split at the midpoint
    int mid = start + (end - start) / 2;
    float bestCost = std::numeric_limits<float>::max();

    for (int i = start; i < end; ++i) {
        AABBox leftBox, rightBox;
        // Compute bounding boxes for triangles on each side of the potential split
        for (int j = start; j <= i; ++j) {
            leftBox.expand(triangles[j].vertices);
        }
        for (int j = i + 1; j < end; ++j) {
            rightBox.expand(triangles[j].vertices);
        }

        float leftSA = computeSurfaceArea(leftBox);
        float rightSA = computeSurfaceArea(rightBox);
        float cost = leftSA * (i - start + 1) + rightSA * (end - i - 1);

        if (cost < bestCost) {
            bestCost = cost;
            mid = i;
        }
    }

    return mid;
}


BVHNode *BVHTree::buildBVH(int start, int end)
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
    node->box = AABBox(minPoint - Vector3(.5f, .5f, .5f), maxPoint + Vector3(.5f, .5f, .5f));

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

//        int mid = start + (end - start) / 2; // Use median splitting
        // Use SAH to determine the best split:
        int mid = findBestSplitSAH(start, end, axis);
        node->left = buildBVH(start, mid);
        node->right = buildBVH(mid, end);
    }

    return node;
}

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    return this->_getIntersectionAndTriangleIndex(root, rayStart, rayEnd, ignoredTriangles);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    std::pair<Vector3, size_t> result = {Vector3(false), -1};
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            if (isIn(triangleIndex, ignoredTriangles)) continue;
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                result = {intersectionPoint, triangleIndex};
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getIntersectionAndTriangleIndex(node->left, rayStart, rayEnd, ignoredTriangles);
        auto rightIntersection = _getIntersectionAndTriangleIndex(node->right, rayStart, rayEnd, ignoredTriangles);
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

std::pair<Vector3, size_t> BVHTree::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, size_t ignoredTriangle) const
{
    return this->_getIntersectionAndTriangleIndex(this->root, rayStart, rayEnd, ignoredTriangle);
}

std::pair<Vector3, size_t> BVHTree::_getIntersectionAndTriangleIndex(BVHNode *node, const Vector3 &rayStart, const Vector3 &rayEnd, size_t ignoredTriangle) const
{
    std::pair<Vector3, size_t> result = {Vector3(false), -1};
    float closestDistance = std::numeric_limits<float>::max();

    if (!node->box.intersects(rayStart, rayEnd).isValid()) {
        return result; // return empty list if ray doesn't intersect bounding box
    }

    if (!node->trianglesIndices.empty()) { // leaf node
        for (const auto& triangleIndex : node->trianglesIndices) {
            if (triangleIndex == ignoredTriangle) continue;
            auto& triangle = this->triangles[triangleIndex];
            Vector3 intersectionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangle[0], triangle[1], triangle[2], true);
            if (intersectionPoint.isValid() && (intersectionPoint - rayStart).norm2() < closestDistance) {
                result = {intersectionPoint, triangleIndex};
                closestDistance = (intersectionPoint - rayStart).norm2();
            }
        }
    } else { // non-leaf node
        auto leftIntersection = _getIntersectionAndTriangleIndex(node->left, rayStart, rayEnd, ignoredTriangle);
        auto rightIntersection = _getIntersectionAndTriangleIndex(node->right, rayStart, rayEnd, ignoredTriangle);
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
