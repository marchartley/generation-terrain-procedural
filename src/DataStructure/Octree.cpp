#include "Octree.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"
#include <set>

OctreeNode::OctreeNode(const Vector3 &origin, const Vector3 &halfDimension)
    : origin(origin), halfDimension(halfDimension), sphereSqrRadius(halfDimension.norm2()) {
    for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

OctreeNode::~OctreeNode() {
    for (int i = 0; i < 8; ++i) delete children[i];
}

bool OctreeNode::intersects(const Vector3 &start, const Vector3 &end) const {
    // Check if the two AABBs are disjoint along any axis
    if (start.x > origin.x + halfDimension.x || end.x < origin.x - halfDimension.x) return false;
    if (start.y > origin.y + halfDimension.y || end.y < origin.y - halfDimension.y) return false;
    if (start.z > origin.z + halfDimension.z || end.z < origin.z - halfDimension.z) return false;

    // If the AABBs are not disjoint along any axis, they must intersect
    return true;
}

bool OctreeNode::intersects(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) const
{
    bool sphereIntersection = this->sphereIntersectsTriangle({v1, v2, v3});
    bool triangleIntersection = Collision::intersectionTriangleAABBox(v1, v2, v3, origin - halfDimension, origin + halfDimension);
    return sphereIntersection && triangleIntersection;
}

bool OctreeNode::intersects(const Triangle &tri) const
{
    return this->intersects(tri[0], tri[1], tri[2]);
}

bool OctreeNode::sphereIntersectsTriangle(const Triangle &triangle) const
{
    Vector3 v1 = triangle[0] - origin;
    Vector3 v2 = triangle[1] - origin;
    Vector3 v3 = triangle[2] - origin;

    // Compute the squared distance from the sphere center to the triangle
    float sqDist = Collision::pointToTriangleDistanceSquared(v1, v2, v3, Vector3());

    // Check if the squared distance is less than the squared sphere radius
    return sqDist <= sphereSqrRadius;
}

Octree::Octree() : root(nullptr)
{
}

Octree::Octree(const Vector3 &origin, const Vector3 &halfDimension)
    : root(new OctreeNode(origin, halfDimension))
{
}

Octree::Octree(const std::vector<Triangle> &triangles)
{
    this->insert(triangles);
}

Octree::Octree(const std::vector<std::vector<Vector3> > &triangles)
{
    this->insert(triangles);

    int maxDepth = 0;
    int nbNodes = 0;
    std::map<int, int> depths = {{0, 0}};
    std::map<int, int> dataStoredPerLevel = {{0, 0}};
    std::vector<int> allElementsWithDuplicates;
    std::vector<std::pair<OctreeNode*, int>> queue = {{root, 1}};
    while (!queue.empty()) {
        nbNodes ++;
        auto [current, depth] = queue.back();
        queue.pop_back();
        maxDepth = std::max(maxDepth, depth);
        depths[depth] = depths[depth] + 1;
        dataStoredPerLevel[depth] = dataStoredPerLevel[depth] + current->data.size();
        for (const auto& data : current->data) {
            allElementsWithDuplicates.push_back(data.index);
        }
        for (auto child : current->children) {
            if (child != nullptr)
                queue.push_back({child, depth + 1});
        }
    }
    std::cout << "Octree : Max depth = " << maxDepth << ", nb nodes = " << nbNodes << ", nb items = " << allElementsWithDuplicates.size() << " without duplicates : " << convertVectorToSet(allElementsWithDuplicates).size() << " for initially " << triangles.size() << " triangles." << std::endl;
    std::cout << "Distribution: \n";
    for (int d = 0; d < depths.size(); d++){
        std::cout << "- level " << d << " : " << depths[d] << " nodes\n";
    }
}

Octree::~Octree() {
    delete root;
}

bool Octree::insert(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const int& pointIndex) {
    if (!this->root->intersects(p1, p2, p3))
        return false; // Don't add this point if it's not inside of the octree space
    return insert(root, p1, p2, p3, pointIndex);
}

bool Octree::insert(std::vector<Triangle> triangles)
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
        auto& triangle = triangles[i];
        bool inserted = insert(triangle[0], triangle[1], triangle[2], i); //(insert(triangles[i][0], i) || insert(triangles[i][1], i) || insert(triangles[i][2], i));
        atLeastOneGood |= inserted;
        if (!inserted) {
            std::cout << "Triangle " << i << "(" << triangle[0] << " " << triangle[1] << " " << triangle[2] << ") ignored" << std::endl;
        }
    }
    return atLeastOneGood;
}

bool Octree::insert(std::vector<std::vector<Vector3> > triangles)
{
    std::vector<Triangle> tris(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++)
        tris[i] = Triangle(triangles[i]);
    return insert(tris);
}

std::vector<OctreeNodeData> Octree::queryRange(const Vector3 &start, const Vector3 &end) const {
    Vector3 _start = Vector3::min(start, end);
    Vector3 _end = Vector3::max(start, end);
    std::vector<OctreeNodeData> result;
    if (this->root)
        queryRange(root, _start, _end, result);
    return result;
}

std::pair<Vector3, int> Octree::_intersectingTriangleIndex(const Vector3 &start, const Vector3 &end, const std::vector<Triangle> &triangles) const
{
    std::vector<OctreeNodeData> allData = this->queryRange(start, end);
    std::set<int> possibleTriangles;
    for (auto& data : allData)
        possibleTriangles.insert(data.index);

    Vector3 intersectionPoint(false);
    int closestIntersectionIndex = -1;
    for (int triIndex : possibleTriangles) {
        auto& triangle = triangles[triIndex];
        Vector3 collide = Collision::segmentToTriangleCollision(start, end, triangle[0], triangle[1], triangle[2]);
        if (!intersectionPoint.isValid() || (collide.isValid() && (collide - start) < (intersectionPoint - start))) {
            intersectionPoint = collide;
            closestIntersectionIndex = triIndex;
        }
    }
    return {intersectionPoint, closestIntersectionIndex};
}

Vector3 Octree::getIntersection(const Vector3 &start, const Vector3 &end, const std::vector<Triangle>& triangles) const
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

std::pair<Vector3, Vector3> Octree::getIntersectionAndNormal(const Vector3 &start, const Vector3 &end, const std::vector<Triangle > &triangles) const
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

Vector3 Octree::getIntersection(const Vector3 &start, const Vector3 &end, const std::vector<std::vector<Vector3> > &triangles) const
{
    std::vector<Triangle> tris(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++)
        tris[i] = Triangle(triangles[i]);
    return this->getIntersection(start, end, tris);
}

std::pair<Vector3, Vector3> Octree::getIntersectionAndNormal(const Vector3 &start, const Vector3 &end, const std::vector<std::vector<Vector3> > &triangles) const
{
    std::vector<Triangle> tris(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++)
        tris[i] = Triangle(triangles[i]);
    return this->getIntersectionAndNormal(start, end, tris);
}

bool Octree::insert(OctreeNode *node, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const int& pointIndex) {
    bool insertValidated = false;
    if (!node->intersects(p1, p2, p3)) {
//        if (node == this->root)
//            std::cout << "Triangle " << pointIndex << "(" << p1 << " " << p2 << " " << p3 << ") rejected by root." << std::endl;
        return false;
    }
    if (node->data.size() < node->maxDataCapacity && node->children[0] == nullptr) {
        // If the node has no children and is not full, add the point here
        node->data.push_back(OctreeNodeData(p1, p2, p3, pointIndex));
        insertValidated = true;
    } else {
        // Otherwise, split the node and add the point to the appropriate child
        if (node->children[0] == nullptr) {
            #pragma omp parallel for
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
            #pragma omp parallel for
            for (size_t child = 0; child < 8; child++)
                insert(node->children[child], data.vertex1, data.vertex2, data.vertex3, data.index);
//                if (insert(node->children[child], data.vertex1, data.vertex2, data.vertex3, data.index))
//                    break;
        }

//        #pragma omp parallel for
        for (size_t child = 0; child < 8; child++)
            insertValidated |= insert(node->children[child], p1, p2, p3, pointIndex);
        // Add the point to the appropriate child
//        int octant = (p1.x >= node->origin.x) + ((p1.y >= node->origin.y) << 1) + ((p1.z >= node->origin.z) << 2);
//        insertValidated = insert(node->children[octant], p1, p2, p3, pointIndex);
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
