#include "SpacePartitioning.h"
#include "SpatialQueries.h"
#include "Utils/Utils.h"
#include "Utils/Collisions.h"

namespace {
bool overlapsSegmentBounds(const Triangle &triangle, const Vector3 &lo, const Vector3 &hi) {
    if (triangle.vertices.size() != 3)
        throw std::invalid_argument("Invalid triangle vertex count");
    for (int axis = 0; axis < 3; ++axis) {
        const float a = triangle.vertices[0][axis], b = triangle.vertices[1][axis],
                    c = triangle.vertices[2][axis];
        if (std::min({a, b, c}) > hi[axis] || std::max({a, b, c}) < lo[axis])
            return false;
    }
    return true;
}
} // namespace

SpacePartitioning::SpacePartitioning() {}

Vector3 SpacePartitioning::getIntersection(const Vector3 &rayStart, const Vector3 &rayEnd) const {
    return this->getIntersectionAndTriangleIndex(rayStart, rayEnd).first;
}

std::pair<Vector3, Vector3> SpacePartitioning::getIntersectionAndNormal(const Vector3 &rayStart,
                                                                        const Vector3 &rayEnd) const {
    if (this->triangles.empty())
        return {Vector3::invalid, Vector3::invalid};
    auto [position, index] = this->getIntersectionAndTriangleIndex(rayStart, rayEnd);
    if (position.isValid())
        return {position, triangles.at(index).getNormal()};
    return {Vector3::invalid, Vector3::invalid};
}

std::vector<Vector3> SpacePartitioning::getAllIntersections(const Vector3 &rayStart,
                                                            const Vector3 &rayEnd) const {
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<Vector3> intersections;
    for (auto &[pos, index] : posAndIndices) {
        if (pos.isValid())
            intersections.push_back(pos);
    }
    return intersections;
}

std::vector<std::pair<Vector3, Vector3>>
SpacePartitioning::getAllIntersectionsAndNormals(const Vector3 &rayStart, const Vector3 &rayEnd) const {
    auto posAndIndices = this->getAllIntersectionsAndTrianglesIndices(rayStart, rayEnd);
    std::vector<std::pair<Vector3, Vector3>> intersections;
    for (auto &[pos, index] : posAndIndices) {
        if (pos.isValid())
            intersections.push_back({pos, triangles.at(index).getNormal()});
    }
    return intersections;
}

std::vector<Triangle> SpacePartitioning::getAllStoredTriangles() const {
    std::vector<size_t> indices = convertSetToVector(getAllStoredTrianglesIndices());
    std::sort(indices.begin(), indices.end());
    std::vector<Triangle> result(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        result[i] = this->triangles.at(indices[i]);
    }
    return result;
}

NoPartitioning::NoPartitioning() : SpacePartitioning() {}

NoPartitioning::NoPartitioning(std::vector<Triangle> triangles) : NoPartitioning() {
    spatial_detail::validate(triangles);
    this->triangles = triangles;
}

const SpacePartitioning &NoPartitioning::build(const std::vector<Triangle> &triangles) {
    spatial_detail::validate(triangles);
    this->triangles = triangles;
    return *this;
}

std::set<size_t> NoPartitioning::getAllStoredTrianglesIndices() const {
    std::set<size_t> indices;
    for (size_t i = 0; i < triangles.size(); i++)
        indices.insert(i);
    return indices;
}

std::pair<Vector3, size_t>
NoPartitioning::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                const std::set<size_t> &ignored) const {
    spatial_detail::validate(a, b);
    const auto lo = Vector3::min(a, b), hi = Vector3::max(a, b);
    double best = std::numeric_limits<double>::infinity();
    std::pair<Vector3, size_t> result{Vector3::invalid, size_t(-1)};
    for (size_t i = 0; i < triangles.size(); ++i) {
        if (ignored.count(i))
            continue;
        double t;
        if (overlapsSegmentBounds(triangles[i], lo, hi) && spatial_detail::triangle(triangles[i], a, b, t) &&
            t < best) {
            best = t;
            result = {spatial_detail::point(a, b, t), i};
        }
    }
    return result;
}
std::pair<Vector3, size_t> NoPartitioning::getIntersectionAndTriangleIndex(const Vector3 &a, const Vector3 &b,
                                                                           size_t ignored) const {
    return getIntersectionAndTriangleIndex(a, b, std::set<size_t>{ignored});
}
std::vector<std::pair<Vector3, size_t>>
NoPartitioning::getAllIntersectionsAndTrianglesIndices(const Vector3 &a, const Vector3 &b) const {
    spatial_detail::validate(a, b);
    const auto lo = Vector3::min(a, b), hi = Vector3::max(a, b);
    std::vector<std::pair<Vector3, size_t>> result;
    for (size_t i = 0; i < triangles.size(); ++i) {
        double t;
        if (overlapsSegmentBounds(triangles[i], lo, hi) && spatial_detail::triangle(triangles[i], a, b, t))
            result.emplace_back(spatial_detail::point(a, b, t), i);
    }
    return result;
}
