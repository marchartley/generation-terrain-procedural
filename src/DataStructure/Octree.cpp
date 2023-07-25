#include "Octree.h"
#include "Utils/Collisions.h"
#include <set>

OctreeNode::OctreeNode(const Vector3 &origin, const Vector3 &halfDimension)
    : origin(origin), halfDimension(halfDimension) {
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

OctreeNode::~OctreeNode() {
    for (int i = 0; i < 8; ++i)
        if (children[i])
            delete children[i];
}

bool OctreeNode::intersects(const Vector3 &start, const Vector3 &end) const {
    // Check if the two AABBs are disjoint along any axis
    if (start.x > origin.x + halfDimension.x || end.x < origin.x - halfDimension.x) return false;
    if (start.y > origin.y + halfDimension.y || end.y < origin.y - halfDimension.y) return false;
    if (start.z > origin.z + halfDimension.z || end.z < origin.z - halfDimension.z) return false;

    // If the AABBs are not disjoint along any axis, they must intersect
    return true;
}

Octree::Octree() : root(nullptr)
{
}

Octree::Octree(const Vector3 &origin, const Vector3 &halfDimension)
    : root(new OctreeNode(origin, halfDimension))
{
}

Octree::Octree(const std::vector<std::vector<Vector3> > &triangles)
{/*
    Vector3 minPoint = Vector3::max();
    Vector3 maxPoint = Vector3::min();
    for (size_t i = 0; i < triangles.size(); ++i) {
        for (const auto& vertex : triangles[i]) {
            minPoint = Vector3::min(minPoint, vertex);
            maxPoint = Vector3::max(maxPoint, vertex);
        }
    }
    AABBox box(minPoint, maxPoint);
    root = new OctreeNode(box.center(), box.dimensions() * .5f);*/
    this->insert(triangles);
}

Octree::~Octree() {
    if (this->root)
        delete root;
}

bool Octree::insert(const Vector3 &point, const int& pointIndex) {
    if (!Vector3::isInBox(point, this->root->origin - this->root->halfDimension, this->root->origin + this->root->halfDimension))
        return false; // Don't add this point if it's not inside of the octree space
    return insert(root, point, pointIndex);
}

bool Octree::insert(std::vector<std::vector<Vector3> > triangles)
{
    if (!this->root) {
        Vector3 mini = Vector3::max(), maxi = Vector3::min();
        for (const auto& t : triangles) {
            mini = Vector3::min({mini, t[0], t[1], t[2]});
            maxi = Vector3::max({maxi, t[0], t[1], t[2]});
        }
        this->root = new OctreeNode((mini + maxi) * .5f, (maxi - mini) * .5f);
    }
    bool atLeastOneGood = false;
    for (size_t i = 0; i < triangles.size(); i++) {
        atLeastOneGood |= (insert(triangles[i][0], i) || insert(triangles[i][1], i) || insert(triangles[i][2], i));
    }
    return atLeastOneGood;
}

std::vector<OctreeNodeData> Octree::queryRange(const Vector3 &start, const Vector3 &end) const {
    Vector3 _start = Vector3::min(start, end);
    Vector3 _end = Vector3::max(start, end);
    std::vector<OctreeNodeData> result;
    if (this->root)
        queryRange(root, _start, _end, result);
    return result;
}

Vector3 Octree::getIntersection(const Vector3 &start, const Vector3 &end, const std::vector<std::vector<Vector3>>& triangles) const
{
    std::vector<OctreeNodeData> allData = this->queryRange(start, end);
    std::set<int> possibleTriangles;
    for (auto& data : allData)
        possibleTriangles.insert(data.index);

    Vector3 intersectionPoint(false);
    for (int triIndex : possibleTriangles) {
        auto& triangle = triangles[triIndex];
        Vector3 collide = Collision::segmentToTriangleCollision(start, end, triangle[0], triangle[1], triangle[2]);
        if (!intersectionPoint.isValid() || (collide.isValid() && (collide - start) < (intersectionPoint - start))) {
            intersectionPoint = collide;
        }
    }
    return intersectionPoint;
}

std::pair<Vector3, Vector3> Octree::getIntersectionAndNormal(const Vector3 &start, const Vector3 &end, const std::vector<std::vector<Vector3> > &triangles) const
{
    std::vector<OctreeNodeData> allData = this->queryRange(start, end);
    std::set<int> possibleTriangles;
    for (auto& data : allData)
        possibleTriangles.insert(data.index);

    Vector3 intersectionPoint(false);
    Vector3 normal;
    for (int triIndex : possibleTriangles) {
        auto& triangle = triangles[triIndex];
        Vector3 collide = Collision::segmentToTriangleCollision(start, end, triangle[0], triangle[1], triangle[2]);
        if (!intersectionPoint.isValid() || (collide.isValid() && (collide - start) < (intersectionPoint - start))) {
            intersectionPoint = collide;
            normal = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]);
        }
    }
    return { intersectionPoint, normal.normalize() };
}

bool Octree::insert(OctreeNode *node, const Vector3 &point, const int& pointIndex) {
    bool insertValidated = false;
    if (!Vector3::isInBox(point - node->origin, -node->halfDimension, node->halfDimension))
        return false;
    if (node->data.size() < 20 && node->children[0] == nullptr) {
        // If the node has no children and is not full, add the point here
        node->data.push_back(OctreeNodeData({point, pointIndex}));
        insertValidated = true;
    } else {
        // Otherwise, split the node and add the point to the appropriate child
        if (node->children[0] == nullptr) {
            for (int i = 0; i < 8; ++i) {
                Vector3 newOrigin = node->origin;
                newOrigin.x += node->halfDimension.x * ((i & 1) ? 0.5f : -0.5f);
                newOrigin.y += node->halfDimension.y * ((i & 2) ? 0.5f : -0.5f);
                newOrigin.z += node->halfDimension.z * ((i & 4) ? 0.5f : -0.5f);
                node->children[i] = new OctreeNode(newOrigin, node->halfDimension * 0.5f);
            }
        }

        auto dataCopy = node->data;
        node->data.clear();
        for (auto& data : dataCopy) {
            for (size_t child = 0; child < 8; child++)
                if (insert(node->children[child], data.pos, data.index))
                    break;
        }

        // Add the point to the appropriate child
        int octant = (point.x >= node->origin.x) + ((point.y >= node->origin.y) << 1) + ((point.z >= node->origin.z) << 2);
        insertValidated = insert(node->children[octant], point, pointIndex);
    }
    return insertValidated;
}

void Octree::queryRange(OctreeNode *node, const Vector3 &start, const Vector3 &end, std::vector<OctreeNodeData> &result, int level) const {
    // If the node does not intersect with the query range, return
    if (!node->intersects(start, end)) {
        return;
    }

    // If the node is a leaf node, check all points in the node
    if (node->children[0] == nullptr) {
        result.insert(result.end(), node->data.begin(), node->data.end());
//        for (const auto& point : node->data) {
//            if (Vector3::isInBox(point.pos, start, end)) {
//                result.push_back(point);
//            }
//        }
    } else {
        // If the node is not a leaf node, check all children
        for (int i = 0; i < 8; ++i) {
            queryRange(node->children[i], start, end, result, level + 1);
        }
    }
}
