#include "Curve.h"

#include "Utils/Collisions.h"

Vector3 Curve::linear(float x, const Vector3 &a, const Vector3 &b) const
{
    return a * (1 - x) + b * x;
}

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

Vector3 Curve::center() const
{
    const auto points = getPath();
    if (points.empty()) return Vector3();

    Vector3 previous = Vector3::invalid;
    Vector3 center;
    int count = 0;
    for (size_t i = 0; i < this->size() - (this->closed ? 1 : 0); i++) {
        const auto& p = points[i];
        float diff = (previous == p ? 0.f : 1.f);
        center += p * diff;
        count += int(diff);
    }
    // for (const auto& point : copy.points)
        // center += point;
    return center / (float) count;
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



std::vector<std::pair<size_t, size_t> > Curve::checkAutointersections() const
{
    auto self = this->clone();
    float spacingHeuristic = 1.0f;
    float selfSpacing = self->length() / float(self->size() - 1);
    if (selfSpacing < spacingHeuristic)
        self->scale(spacingHeuristic / selfSpacing);
    std::vector<std::pair<size_t, size_t> > results;
    int s = size();
#pragma omp parallel for collapse(2)
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            int i00 = i;
            int i01 = (i00 + 1) % s;
            int i10 = (j + i + 2) % s;
            int i11 = (i10 + 1) % s;

            if (i10 <= i00 || i11 == i00 || i01 == i10 || (!closed && i11 <= i00)) continue;

            Vector3 intersection = Collision::intersectionBetweenTwoSegments(self->get(i00), self->get(i01), self->get(i10), self->get(i11), 1e-6);

            if (intersection.isValid()) {
#pragma omp critical
                {
                    results.push_back({i00, i10});
                }
            }
        }
    }
    return results;
}
