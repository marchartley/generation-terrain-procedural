#include "Polyline.h"

#include "Utils/Collisions.h"

Polyline::Polyline() : Polyline(std::vector<Vector3>{})
{

}

Polyline::Polyline(const Polyline &s)
{
    *this = s;
}

Polyline::Polyline(const std::vector<Vector3>& points)
    : Curve(), points(points)
{

}

std::vector<Vector3> Polyline::getPath(int numberOfPoints) const
{
    if (numberOfPoints < 0 || size_t(numberOfPoints) == this->points.size()) return this->points;

    Polyline resampled(this->points);
    resampled.resamplePoints(numberOfPoints);
    return resampled.getPath();
}

Vector3 Polyline::getPoint(float x) const
{
    x = clamp(x, 0.f, 1.f);
    size_t startId = timeToLowerIndex(x);
    size_t endId = startId + 1;

    float frac = (x * (points.size() - 1)) - startId;

    return getPoint(frac, points[startId], points[endId]); // Linear
}

Vector3 Polyline::getPoint(float x, const Vector3& a, const Vector3& b) const
{
    x = clamp(x, 0.f, 1.f);
    return a * (1.f - x) + b * x;
}

Vector3 Polyline::getDerivative(float x, bool normalize) const
{
    x = std::max(x, 0.f);
    size_t startId = timeToLowerIndex(x);
    Vector3 derivative;
    if (x >= 1.f) {
        derivative = points[numPoints() - 1] - points[numPoints() - 2];
    } else {
        derivative = points[startId + 1] - points[startId];
    }
    return (normalize ? derivative.normalize() : derivative);
}

Vector3 Polyline::getSecondDerivative(float x, bool normalize) const
{
    return Vector3(0.f, 0.f, 0.f); // Consider no second derivative since linear piece-wise.
}

float Polyline::estimateClosestTime(const Vector3& pos) const
{
    size_t closestStartPoint = 0;
    float minDist = std::numeric_limits<float>::max();

    for (size_t i = 0; i < numPoints() - 1; i++) {
        Vector3 proj = Collision::projectPointOnSegment(pos, points[i], points[i+1]);
        float dist = (proj - pos).norm2();
        if (dist < minDist) {
            minDist = dist;
            closestStartPoint = i;
        }
    }

    return (float(closestStartPoint) + (Collision::projectPointOnSegment(pos, points[closestStartPoint], points[closestStartPoint + 1]) - points[closestStartPoint]).norm() / (points[closestStartPoint + 1] - points[closestStartPoint]).norm()) / float(numPoints() - 1);
}

Vector3 Polyline::estimateClosestPos(const Vector3& pos) const
{
    return getPoint(estimateClosestTime(pos));
}

float Polyline::estimateSqrDistanceFrom(const Vector3& pos) const
{
    float minDist = std::numeric_limits<float>::max();

    for (size_t i = 0; i < numPoints() - 1; i++) {
        Vector3 proj = Collision::projectPointOnSegment(pos, points[i], points[i+1]);
        float dist = (proj - pos).norm2();
        if (dist < minDist) {
            minDist = dist;
        }
    }
    return minDist;
}

float Polyline::length() const
{
    float length = 0;
    if (this->points.empty()) return length;
    for (size_t i = 0; i < this->points.size() - 1; i++) {
        length += (this->points[i] - this->points[i + 1]).norm();
    }
    return length;
}

Polyline& Polyline::setPoint(int i, const Vector3& newPos)
{
    this->points[getIndex(i)] = newPos;
    return *this;
}

Polyline& Polyline::resamplePoints(int newNbPoints)
{
    if (newNbPoints < 0) return *this;

    std::vector<Vector3> newPoints(newNbPoints);
    for (int i = 0; i < newNbPoints; i++) {
        float t = float(i) / (float(numPoints() - 1));
        newPoints[i] = getPoint(t);
    }
    this->points = newPoints;
    return *this;
}

Polyline& Polyline::reverseVertices()
{
    std::reverse(points.begin(), points.end());
    return *this;
}

std::pair<Vector3, Vector3> Polyline::AABBox() const
{
    Vector3 mini = Vector3::max();
    Vector3 maxi = Vector3::min();

    for (auto& p : points) {
        mini.x() = std::min(mini.x(), p.x());
        maxi.x() = std::max(maxi.x(), p.x());
        mini.y() = std::min(mini.y(), p.y());
        maxi.y() = std::max(maxi.y(), p.y());
        mini.z() = std::min(mini.z(), p.z());
        maxi.z() = std::max(maxi.z(), p.z());
    }
    return {mini, maxi};
}

Polyline& Polyline::scale(const Vector3& factor)
{
    for (auto& p : points) {
        p *= factor;
    }
    return *this;
}

Polyline& Polyline::translate(const Vector3& translation)
{
    for (auto& p : points) {
        p += translation;
    }
    return *this;
}

Polyline& Polyline::removeDuplicates()
{
    for (int i = int(numPoints()) - 1; i > 0; i--) {
        if (points[i] == points[i - 1]) {
            points.erase(points.begin() + (i - 1));
        }
    }
    return *this;
}

size_t Polyline::numPoints() const
{
    return points.size();
}

std::string Polyline::toString() const
{
    std::ostringstream out;
    out << "Polyline with " << this->points.size() << " points (" << (closed ? "closed" : "not closed") << ") :\n";
    for (auto& p : this->points)
        out << "- " << p << "\n";
    return out.str();
}

Polyline& Polyline::close()
{
    this->closed = true;
}


Vector3& Polyline::operator[](size_t i) {
    return this->points[i];
}
const Vector3& Polyline::operator[](size_t i) const {
    return this->points[i];
}
