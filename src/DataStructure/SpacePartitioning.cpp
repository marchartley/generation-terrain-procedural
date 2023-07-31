#include "SpacePartitioning.h"
#include "Utils/Utils.h"

SpacePartitioning::SpacePartitioning()
{

}
/*
SpacePartitioning::SpacePartitioning(const std::vector<Triangle> &triangles)
    : triangles(triangles)
{
}

SpacePartitioning::SpacePartitioning(const std::vector<std::vector<Vector3> > &triangles)
    : SpacePartitioning(Triangle::vectorsToTriangles(triangles))
{
}
*/
//SpacePartitioning &SpacePartitioning::build(const std::vector<Triangle> &triangles)
//{
//    this->triangles = triangles;
//}

//SpacePartitioning &SpacePartitioning::build(const std::vector<std::vector<Vector3> > &triangles) const
//{
//    return this->build(Triangle::vectorsToTriangles(triangles));
//}

Vector3 SpacePartitioning::getIntersection(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    return this->getIntersectionAndTriangleIndex(rayStart, rayEnd).first;
}

std::pair<Vector3, Vector3> SpacePartitioning::getIntersectionAndNormal(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    auto [position, index] = this->getIntersectionAndTriangleIndex(rayStart, rayEnd);
    return {position, triangles[index].normal};
}

std::vector<Vector3> SpacePartitioning::getAllIntersections(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<Vector3> intersections;
    for (auto& [pos, index] : posAndIndices) {
        intersections.push_back(pos);
    }
    return intersections;
}

std::vector<std::pair<Vector3, Vector3> > SpacePartitioning::getAllIntersectionsAndNormals(const Vector3 &rayStart, const Vector3 &rayEnd) const
{
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<std::pair<Vector3, Vector3>> intersections;
    for (auto& [pos, index] : posAndIndices) {
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

//std::pair<Vector3, size_t> SpacePartitioning::getIntersectionAndTriangleIndex(const Vector3 &rayStart, const Vector3 &rayEnd) const
//{

//}
