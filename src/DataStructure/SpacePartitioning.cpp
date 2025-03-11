#include "SpacePartitioning.h"
#include "Utils/Utils.h"
#include "Utils/Collisions.h"

SpacePartitioning::SpacePartitioning()
{

}

Vector3 SpacePartitioning::getIntersection(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->getIntersectionAndTriangleIndex(rayStart, rayEnd).first;
}

std::pair<Vector3, Vector3> SpacePartitioning::getIntersectionAndNormal(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    if (this->triangles.empty()) return {Vector3(false), Vector3(false)};
    auto [position, index] = this->getIntersectionAndTriangleIndex(rayStart, rayEnd);
    if (position.isValid())
        return {position, triangles[index].normal};
    return {Vector3(false), Vector3(false)};
}

std::vector<Vector3> SpacePartitioning::getAllIntersections(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<Vector3> intersections;
    for (auto& [pos, index] : posAndIndices) {
        if (pos.isValid())
            intersections.push_back(pos);
    }
    return intersections;
}

std::vector<std::pair<Vector3, Vector3> > SpacePartitioning::getAllIntersectionsAndNormals(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<std::pair<Vector3, Vector3>> intersections;
    for (auto& [pos, index] : posAndIndices) {
        if (pos.isValid())
            intersections.push_back({pos, triangles[index].normal});
    }
    return intersections;
}

std::vector<Triangle> SpacePartitioning::getAllStoredTriangles() const
{
    std::vector<size_t> indices = convertSetToVector(getAllStoredTrianglesIndices());
    std::sort(indices.begin(), indices.end());
    std::vector<Triangle> result(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        result[i] = this->triangles[indices[i]];
    }
    return result;
}


NoPartitioning::NoPartitioning() : SpacePartitioning()
{

}

NoPartitioning::NoPartitioning(std::vector<Triangle> triangles)
    : NoPartitioning()
{
    this->triangles = triangles;
}

SpacePartitioning &NoPartitioning::build(const std::vector<Triangle> &triangles)
{
    this->triangles = triangles;
}

std::set<size_t> NoPartitioning::getAllStoredTrianglesIndices() const
{
    std::set<size_t> indices;
    for (size_t i = 0; i < triangles.size(); i++)
        indices.insert(i);
    return indices;
}

std::pair<Vector3, size_t> NoPartitioning::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, std::set<size_t> ignoredTriangles) const
{
    std::map<size_t, Vector3> collisions;
#pragma omp parallel for
    for (size_t i = 0; i < triangles.size(); i++) {
        Vector3 collisionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangles[i][0], triangles[i][1], triangles[i][2]);
        if (collisionPoint.isValid()) {
            collisions[i] = collisionPoint;
        }
    }
    size_t closestIdx;
    Vector3 closestIntersection(rayEnd);
    closestIntersection.setValid(false);
    for (auto& [idx, point] : collisions) {
        if (ignoredTriangles.find(idx) == ignoredTriangles.end() && (rayStart - point).norm2() < (rayStart - closestIntersection).norm2()) {
            closestIdx = idx;
            closestIntersection = point;
        }
    }
    return {closestIntersection, closestIdx};
}

std::pair<Vector3, size_t> NoPartitioning::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd, size_t ignoredTriangle) const
{
    return getIntersectionAndTriangleIndex(rayStart, rayEnd, std::set<size_t>({ignoredTriangle}));
}

std::vector<std::pair<Vector3, size_t> > NoPartitioning::getAllIntersectionsAndTrianglesIndices(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    std::map<size_t, Vector3> collisions;
#pragma omp parallel for
    for (size_t i = 0; i < triangles.size(); i++) {
        Vector3 collisionPoint = Collision::segmentToTriangleCollision(rayStart, rayEnd, triangles[i][0], triangles[i][1], triangles[i][2]);
        if (collisionPoint.isValid()) {
            collisions[i] = collisionPoint;
        }
    }

    std::vector<std::pair<Vector3, size_t> > result(collisions.size());
    int i = 0;
    for (auto& [idx, point] : collisions)
        result[i++] = {point, idx};
    return result;
}
