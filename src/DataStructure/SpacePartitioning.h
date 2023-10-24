#ifndef SPACEPARTITIONING_H
#define SPACEPARTITIONING_H

#include <set>
#include "DataStructure/Triangle.h"

class SpacePartitioning
{
public:
    SpacePartitioning();
//    SpacePartitioning(const std::vector<Triangle>& triangles);
//    SpacePartitioning(const std::vector<std::vector<Vector3>>& triangles);

    virtual SpacePartitioning& build(const std::vector<Triangle>& triangles) = 0;
//    virtual SpacePartitioning& build(const std::vector<std::vector<Vector3>>& triangles) const;

    virtual Vector3 getIntersection(const Vector3& rayStart, const Vector3& rayEnd) const;
    virtual std::pair<Vector3, Vector3> getIntersectionAndNormal(const Vector3& rayStart, const Vector3& rayEnd) const;

    virtual std::vector<Vector3> getAllIntersections(const Vector3& rayStart, const Vector3& rayEnd) const;
    virtual std::vector<std::pair<Vector3, Vector3>> getAllIntersectionsAndNormals(const Vector3& rayStart, const Vector3& rayEnd) const;

    virtual std::set<size_t> getAllStoredTrianglesIndices() const = 0;
    virtual std::vector<Triangle> getAllStoredTriangles() const;

//protected:
    virtual std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, std::set<size_t> ignoredTriangles = {}) const = 0;
    virtual std::pair<Vector3, size_t> getIntersectionAndTriangleIndex(const Vector3& rayStart, const Vector3& rayEnd, size_t ignoredTriangle) const = 0;
    virtual std::vector<std::pair<Vector3, size_t>> getAllIntersectionsAndTrianglesIndices(const Vector3& rayStart, const Vector3& rayEnd) const = 0;
    std::vector<Triangle> triangles;
};

#endif // SPACEPARTITIONING_H
