#include "DataStructure/Vector3.h"
#include "DataStructure/SpatialQueries.h"

AABBox::AABBox() : AABBox(Vector3(0, 0, 0), Vector3(0, 0, 0)) {}

AABBox::AABBox(const Vector3 &mini, const Vector3 &maxi) : mini(mini), maxi(maxi) {}

AABBox::AABBox(const std::pair<Vector3, Vector3> &minMax)
    : AABBox(std::get<0>(minMax), std::get<1>(minMax)) {}

AABBox::AABBox(const std::vector<Vector3> &points) : mini(Vector3::max()), maxi(Vector3::min()) {
    for (const auto &p : points) {
        if (!p.isValid())
            throw std::invalid_argument("Invalid box point");
        expand(p);
    }
}

Vector3 AABBox::normalize(const Vector3 &p) {
    if (!p.isValid() || !mini.isValid() || !maxi.isValid() || mini.x() > maxi.x() || mini.y() > maxi.y() ||
        mini.z() > maxi.z())
        throw std::invalid_argument("Cannot normalize an invalid box/point");
    Vector3 result;
    for (int i = 0; i < 3; ++i) {
        double extent = double(maxi[i]) - mini[i];
        result[i] = extent ? float((double(p[i]) - mini[i]) / extent) : 0.f;
    }
    return result;
}

AABBox &AABBox::expand(const Vector3 &newPoint) {
    this->mini = Vector3::min(this->mini, newPoint);
    this->maxi = Vector3::max(this->maxi, newPoint);
    return *this;
}

AABBox &AABBox::expand(const std::vector<Vector3> &newPoints) {
    for (const auto &p : newPoints)
        this->expand(p);
    return *this;
}

float AABBox::distanceTo(const Vector3 &p) { return float(std::sqrt(spatial_detail::distance2(*this, p))); }

Vector3 AABBox::random(const Vector3 &mini, const Vector3 &maxi) { return Vector3::random(mini, maxi); }

Vector3 AABBox::intersects(const Vector3 &a, const Vector3 &b) {
    spatial_detail::validate(a, b);
    double entry;
    return spatial_detail::box(*this, a, b, 1, entry) ? spatial_detail::point(a, b, entry) : Vector3::invalid;
}

std::ostream &operator<<(std::ostream &io, const AABBox &bbox) {
    io << bbox.toString();
    return io;
}

std::ostream &operator<<(std::ostream &io, std::shared_ptr<AABBox> bbox) {
    io << bbox->toString();
    return io;
}
