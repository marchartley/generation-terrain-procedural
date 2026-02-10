#include "DataStructure/Vector3.h"


AABBox::AABBox()
    : AABBox(Vector3(0, 0, 0), Vector3(0, 0, 0))
{

}

AABBox::AABBox(const Vector3& mini, const Vector3& maxi) : mini(mini), maxi(maxi)
{

}

AABBox::AABBox(const std::tuple<Vector3, Vector3> &minMax)
    : AABBox(std::get<0>(minMax), std::get<1>(minMax))
{

}

AABBox::AABBox(const std::vector<Vector3> &allPointsToContain)
    : AABBox(Vector3::min(allPointsToContain), Vector3::max(allPointsToContain))
{

}

Vector3 AABBox::normalize(const Vector3& p)
{
    auto mini = this->min();
    auto maxi = this->max();
    Vector3 ret = (p - mini) / (maxi - mini);
    return ret;
}

AABBox &AABBox::expand(const Vector3 &newPoint)
{
    this->mini = Vector3::min(this->mini, newPoint);
    this->maxi = Vector3::max(this->maxi, newPoint);
    return *this;
}

AABBox& AABBox::expand(const std::vector<Vector3>& newPoints)
{
    for (const auto& p : newPoints)
        this->expand(p);
    return *this;
}

float AABBox::distanceTo(const Vector3 &p)
{
    return Vector3::distanceToBoundaries(p, this->min(), this->max());
}

Vector3 AABBox::random(const Vector3& mini, const Vector3& maxi)
{
    return Vector3::random(mini, maxi);
}

Vector3 AABBox::intersects(const Vector3& rayStart, const Vector3& rayEnd) {
    if (this->contains(rayStart)) return rayStart;

    Vector3 direction = (rayEnd - rayStart).normalized();
    Vector3 tMin = (min() - rayStart) / direction;
    Vector3 tMax = (max() - rayStart) / direction;

    if (tMin.x() > tMax.x()) std::swap(tMin.x(), tMax.x());
    if (tMin.y() > tMax.y()) std::swap(tMin.y(), tMax.y());
    if (tMin.z() > tMax.z()) std::swap(tMin.z(), tMax.z());

    float tNear = std::max({tMin.x(), tMin.y(), tMin.z()});
    float tFar = std::min({tMax.x(), tMax.y(), tMax.z()});

    if (tNear > tFar || tFar < 0) {
        return Vector3::invalid();  // return invalid Vector3 if no intersection
    }

    Vector3 intersectionPoint = rayStart + direction * tNear; // calculate intersection point
    return intersectionPoint;
}

std::ostream& operator<<(std::ostream& io, const AABBox& bbox) {
    io << bbox.toString();
    return io;
}

std::ostream& operator<<(std::ostream& io, std::shared_ptr<AABBox> bbox) {
    io << bbox->toString();
    return io;
}
