#include "Curve.h"

float Curve::estimateDistanceFrom(const Vector3 &pos) const {
    return std::sqrt(estimateSqrDistanceFrom(pos));
}

float Curve::estimateSignedDistanceFrom(const Vector3 &pos) const {
    // Only available for 2D paths, otherwise there's no sense
    float t = this->estimateClosestTime(pos);
    Vector3 normal = this->getNormal(t);
    Vector3 posOnCurve = this->getPoint(t);
    float dist = (posOnCurve - pos).norm();
    float sign = (normal.dot(posOnCurve - pos) > 0.f ? 1.f : -1.f);
    return dist * sign;
}

size_t Curve::timeToLowerIndex(float t) const {
    return size_t(std::floor(clamp(t, 0.f, 1.f) * (numSegments()))) - (t >= 1.f ? 1 : 0);
}

Vector3 Curve::getDirection(float x) const {
    return this->getDerivative(x, true);
}

Vector3 Curve::getNormal(float x) const {
    const auto derivative = getDerivative(x);
    return derivative.cross(this->getSecondDerivative(x)).cross(derivative) .normalize();
}

Vector3 Curve::getBinormal(float x) const {
    return this->getDerivative(x).cross(this->getSecondDerivative(x)).normalize();
}

float Curve::getCurvature(float x) const {
    return (getDerivative(x).cross(getSecondDerivative(x))).norm() / (std::pow(getDerivative(x).norm(), 3));
}

Vector3 Curve::containingBoxSize() const {
    auto box = AABBox();
    return box.second - box.first;
}

Curve& Curve::scale(float factor) {
    return scale(Vector3(factor, factor, factor));
}

size_t Curve::size() const {
    return numPoints();
}

size_t Curve::numSegments() const
{
    return numPoints() - 1;
}

size_t Curve::pointIndex(int index) const
{
    return (index + numPoints()) % numPoints();
}

bool Curve::empty() const {
    return numPoints() <= 0;
}

std::ostream& operator<<(std::ostream &io, const Curve& curve) {
    io << curve.toString();
    return io;
}

std::ostream& operator<<(std::ostream &io, std::shared_ptr<Curve> curve) {
    io << curve->toString();
    return io;
}

bool Curve::isClosed() const {
    return closed;
}

Curve& Curve::setClosed(bool close) {
    if (close) this->close();
    this->closed = close;
    return *this;
}

Curve& Curve::close()
{
    this->closed = true;
    return *this;
}
