#ifndef SPATIAL_QUERIES_H
#define SPATIAL_QUERIES_H
#include "DataStructure/Triangle.h"
#include <limits>
#include <array>
#include <stdexcept>

// Shared closed-segment/closed-triangle policy for every partitioning strategy.
namespace spatial_detail {
inline void validate(const Vector3 &a, const Vector3 &b) {
    if (!a.isValid() || !b.isValid())
        throw std::invalid_argument("Nonfinite spatial query");
}
inline void validate(const std::vector<Triangle> &triangles) {
    if (triangles.size() > size_t(std::numeric_limits<int>::max()))
        throw std::length_error("Too many triangles");
    for (const auto &t : triangles) {
        if (t.vertices.size() != 3)
            throw std::invalid_argument("Triangle needs three vertices");
        for (const auto &v : t.vertices)
            if (!v.isValid())
                throw std::invalid_argument("Nonfinite triangle");
    }
}
inline bool box(const AABBox &bounds, const Vector3 &a, const Vector3 &b, double limit, double &near) {
    near = 0;
    double far = std::min(1., limit);
    for (int i = 0; i < 3; ++i) {
        const double lo = bounds.mini[i], hi = bounds.maxi[i], origin = a[i], direction = double(b[i]) - a[i];
        if (lo > hi)
            return false;
        if (direction == 0) {
            if (origin < lo || origin > hi)
                return false;
            continue;
        }
        double x = (lo - origin) / direction, y = (hi - origin) / direction;
        if (x > y)
            std::swap(x, y);
        near = std::max(near, x);
        far = std::min(far, y);
        if (near > far)
            return false;
    }
    return near <= limit;
}
inline bool triangle(const Triangle &tri, const Vector3 &a, const Vector3 &b, double &time) {
    const Vec3<double> edge1 = Vec3<double>(tri[1]) - Vec3<double>(tri[0]),
                       edge2 = Vec3<double>(tri[2]) - Vec3<double>(tri[0]);
    const Vec3<double> direction = Vec3<double>(b) - Vec3<double>(a);
    const auto h = direction.cross(edge2);
    double det = edge1.dot(h);
    const auto s = Vec3<double>(a) - Vec3<double>(tri[0]);
    if (det == 0) {
        // Coplanar segments (including a point query) are clipped in barycentric
        // coordinates. Return the earliest overlap; a degenerate triangle is not a surface.
        const auto n = edge1.cross(edge2);
        const double n2 = n.dot(n);
        if (n2 == 0)
            return false;
        const double plane = s.dot(n);
        const double error = 32 * std::numeric_limits<double>::epsilon() *
                             (std::abs(s.x() * n.x()) + std::abs(s.y() * n.y()) + std::abs(s.z() * n.z()));
        if (std::abs(plane) > error)
            return false;
        const double u0 = s.cross(edge2).dot(n) / n2, v0 = edge1.cross(s).dot(n) / n2;
        const double du = direction.cross(edge2).dot(n) / n2, dv = edge1.cross(direction).dot(n) / n2;
        double enter = 0, exit = 1;
        auto clip = [&](double value, double slope) {
            if (slope == 0)
                return value >= 0;
            const double bound = -value / slope;
            if (slope > 0)
                enter = std::max(enter, bound);
            else
                exit = std::min(exit, bound);
            return enter <= exit;
        };
        if (!clip(u0, du) || !clip(v0, dv) || !clip(1 - u0 - v0, -du - dv))
            return false;
        time = enter;
        return true;
    }
    if (!std::isfinite(det))
        return false;
    const double u = s.dot(h) / det;
    const auto q = s.cross(edge1);
    const double v = direction.dot(q) / det;
    constexpr double tolerance = 32 * std::numeric_limits<double>::epsilon();
    if (u < -tolerance || v < -tolerance || u + v > 1 + tolerance)
        return false;
    time = edge2.dot(q) / det;
    if (time < -tolerance || time > 1 + tolerance)
        return false;
    time = std::clamp(time, 0., 1.);
    return true;
}
inline Vector3 point(const Vector3 &a, const Vector3 &b, double t) {
    return Vector3(double(a.x()) + (double(b.x()) - a.x()) * t, double(a.y()) + (double(b.y()) - a.y()) * t,
                   double(a.z()) + (double(b.z()) - a.z()) * t);
}
inline double triangleDistance2(const Triangle &tri, const Vector3 &point) {
    const Vec3<double> a = tri[0], b = tri[1], c = tri[2], p = point;
    const auto e1 = b - a, e2 = c - a, ap = p - a, n = e1.cross(e2);
    const double n2 = n.dot(n);
    if (n2 > 0) {
        const double u = ap.cross(e2).dot(n) / n2, v = e1.cross(ap).dot(n) / n2;
        if (u >= 0 && v >= 0 && u + v <= 1) {
            const double d = ap.dot(n);
            return d * d / n2;
        }
    }
    double best = std::numeric_limits<double>::infinity();
    for (auto edge : std::array<std::pair<Vec3<double>, Vec3<double>>, 3>{{{a, b}, {b, c}, {c, a}}}) {
        auto d = edge.second - edge.first;
        double length = d.dot(d);
        double t = length ? std::clamp((p - edge.first).dot(d) / length, 0., 1.) : 0;
        auto delta = p - (edge.first + d * t);
        best = std::min(best, delta.dot(delta));
    }
    return best;
}
inline double distance2(const AABBox &box, const Vector3 &p) {
    double sum = 0;
    for (int i = 0; i < 3; ++i) {
        double d = std::max({double(box.mini[i]) - p[i], 0., double(p[i]) - box.maxi[i]});
        sum += d * d;
    }
    return sum;
}
} // namespace spatial_detail
#endif
