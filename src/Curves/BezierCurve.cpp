#include "BezierCurve.h"

#include "Utils/Collisions.h"

BezierCurve::BezierCurve() : BezierCurve(std::vector<Vector3>())
{}

BezierCurve::BezierCurve(const std::vector<Vector3> &points)
    : BezierCurve(points, std::vector<Vector3>())
{
    handles = std::vector<Vector3>(2 * (points.size() - 1));
    for (size_t i = 1; i < points.size() - 1; i++) {
        Vector3 proj = Collision::projectPointOnSegment(points[i], points[i - 1], points[i + 1]);
        handles[2 * i - 1] = points[i] - (proj - points[i - 1]).normalize() * (points[i] - points[i - 1]).norm() / 3.f;
        handles[2 * i] = points[i] + (points[i + 1] - proj).normalize() * (points[i + 1] - points[i]).norm() / 3.f;
    }
    // Strategy 1 : tangents towards next / previous points
    // handles.front() = points[0] + (points[1] - points[0]) / 3.f;
    // handles.back() = points[numPoints() - 1] + (points[numPoints() - 2] - points[numPoints() - 1]) / 3.f;

    // Strategy 2 : tangent computation use next/previous 2 points
    // Start endpoint: use direction from p0 based on p1 -> p2
    Vector3 t0 = points[2] - points[1];
    handles.front() = points[0] - t0.normalize() * (points[1] - points[0]).norm() / 3.f;
    // End endpoint: use direction based on p[n-3] -> p[n-2]
    size_t n = points.size();
    Vector3 tn = points[n - 2] - points[n - 3];
    handles.back() = points[n - 1] + tn.normalize() * (points[n - 1] - points[n - 2]).norm() / 3.f;
}

BezierCurve::BezierCurve(const std::vector<Vector3>& points, const std::vector<Vector3>& handles)
    : Curve(), points(points), handles(handles)
{
}

std::vector<Vector3> BezierCurve::getPath(int numberOfPoints) const
{
    if (numberOfPoints < 0 || numberOfPoints == numPoints()) return points;
    std::vector<Vector3> points(numberOfPoints);
    for (int i = 0; i < numberOfPoints; i++) {
        float t = float(i) / float(numberOfPoints - 1);
        points[i] = getPoint(t);
    }
    return points;
}

inline Vector3 linearBezier(const Vector3& P0, const Vector3& P1, float t) {
    return P0 * (1.f - t) + P1 * t;
}
inline Vector3 quadraticBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, float t) {
    return linearBezier(P0, P1, t) * (1.f - t) + linearBezier(P1, P2, t) * t;
}
inline Vector3 cubicBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    // return linearBezier(quadraticBezier(P0, P1, P2, t), quadraticBezier(P1, P2, P3, t), t);
    return std::pow(1.f-t, 3) * P0 + 3.f * std::pow(1.f - t, 2) * t * P1 + 3.f * (1.f - t) * t*t * P2 + t*t*t * P3;
}
inline Vector3 cubicBezierDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    return 3.f * std::pow(1.f - t, 2) * (P1 - P0) + 6.f * (1.f - t) * t * (P2 - P1) + 3.f * t * t * (P3 - P2);
}
inline Vector3 cubicBezierSecondDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    return 6.f * (1.f - t) * (P0 + P2 - 2.f * P1) + 6.f * t * (P1 + P3 - 2.f * P2);
}
Vector3 BezierCurve::getPoint(float x) const
{
    int pointId = timeToLowerIndex(x);
    int handleId = pointId * 2;
    const auto& P0 = points[pointId];
    const auto& P1 = handles[handleId];
    const auto& P2 = handles[handleId + 1];
    const auto& P3 = points[pointId + 1];

    float t = (x * float(numPoints() - 1)) - pointId;
    return cubicBezier(P0, P1, P2, P3, t);
}

Vector3 BezierCurve::getPoint(float x, const Vector3 &a, const Vector3 &b) const
{
    return linearBezier(a, b, x);
}

Vector3 BezierCurve::getDerivative(float x, bool normalize) const
{
    int pointId = timeToLowerIndex(x);
    int handleId = pointId * 2;
    const auto& P0 = points[pointId];
    const auto& P1 = handles[handleId];
    const auto& P2 = handles[handleId + 1];
    const auto& P3 = points[pointId + 1];

    float t = (x * float(numPoints() - 1)) - pointId;
    return cubicBezierDerivative(P0, P1, P2, P3, t);
}

Vector3 BezierCurve::getSecondDerivative(float x, bool normalize) const
{
    int pointId = timeToLowerIndex(x);
    int handleId = pointId * 2;
    const auto& P0 = points[pointId];
    const auto& P1 = handles[handleId];
    const auto& P2 = handles[handleId + 1];
    const auto& P3 = points[pointId + 1];

    float t = (x * float(numPoints())) - pointId;
    return cubicBezierSecondDerivative(P0, P1, P2, P3, t);
}

float BezierCurve::estimateClosestTime(const Vector3 &pos) const
{
    int segmentCount = numPoints() - 1;
    int bestSegment = 0;
    float bestT = 0.f;
    float bestDistSq = std::numeric_limits<float>::max();

    constexpr int samplesPerSegment = 20;

    // Coarse search
    for (int segment = 0; segment < segmentCount; ++segment)
    {
        for (int i = 0; i <= samplesPerSegment; ++i)
        {
            float t = float(i) / float(samplesPerSegment);
            Vector3 p = cubicBezier(points[segment], handles[segment * 2], handles[segment * 2 + 1], points[segment + 1], t);

            float dSq = (p - pos).norm2();

            if (dSq < bestDistSq)
            {
                bestDistSq = dSq;
                bestSegment = segment;
                bestT = t;
            }
        }
    }

    // Newton refinement on the closest segment
    constexpr int iterations = 8;

    for (int i = 0; i < iterations; ++i)
    {
        Vector3 p = cubicBezier(points[bestSegment], handles[bestSegment * 2], handles[bestSegment * 2 + 1], points[bestSegment + 1], bestT);
        Vector3 d1 = cubicBezierDerivative(points[bestSegment], handles[bestSegment * 2], handles[bestSegment * 2 + 1], points[bestSegment + 1], bestT);
        Vector3 d2 = cubicBezierSecondDerivative(points[bestSegment], handles[bestSegment * 2], handles[bestSegment * 2 + 1], points[bestSegment + 1], bestT);

        Vector3 r = p - pos;

        float numerator = r.dot(d1);
        float denominator = d1.dot(d1) + r.dot(d2);

        if (std::abs(denominator) < 1e-6f)
            break;

        bestT -= numerator / denominator;
        bestT = std::clamp(bestT, 0.f, 1.f);
    }

    return (float(bestSegment) + bestT) / float(segmentCount);
}

Vector3 BezierCurve::estimateClosestPos(const Vector3 &pos) const
{
    return getPoint(estimateClosestTime(pos));
}

float BezierCurve::estimateSqrDistanceFrom(const Vector3 &pos) const
{
    return (estimateClosestPos(pos) - pos).norm2();
}

float BezierCurve::length() const
{
    auto path = getPath(5 * numPoints());
    float dist = 0;
    for (size_t i = 0; i < path.size() - 1; i++) {
        dist += (path[i + 1] - path[i]).norm();
    }
    return dist;
}

BezierCurve& BezierCurve::setPoint(int i, const Vector3 &newPos)
{
    this->points[i] = newPos;
    return *this;
}

BezierCurve& BezierCurve::resamplePoints(int newNbPoints)
{
    if (newNbPoints < 0 || newNbPoints == numPoints()) return *this;
    std::vector<Vector3> newPoints(newNbPoints);
    std::vector<Vector3> newHandles(newNbPoints);

    for (int i = 0; i < newNbPoints; i++) {
        float t = float(i) / (float(newNbPoints) - 1);
        newPoints[i] = getPoint(t);
    }
    for (int i = 0; i < newNbPoints - 1; i++) {
        float t0 = float(i) / (float(newNbPoints) - 1);
        float t1 = float(i+1) / (float(newNbPoints) - 1);
        newHandles[i * 2] = newPoints[i] + getDerivative(t0);
        newHandles[i * 2 + 1] = newPoints[i+1] - getDerivative(t1);
    }
    this->points = newPoints;
    this->handles = newHandles;
    return *this;
}

BezierCurve& BezierCurve::reverseVertices()
{
    std::reverse(points.begin(), points.end());
    std::reverse(handles.begin(), handles.end());
    return *this;
}


std::pair<Vector3, Vector3> BezierCurve::AABBox() const
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
    for (auto& p : handles) {
        mini.x() = std::min(mini.x(), p.x());
        maxi.x() = std::max(maxi.x(), p.x());
        mini.y() = std::min(mini.y(), p.y());
        maxi.y() = std::max(maxi.y(), p.y());
        mini.z() = std::min(mini.z(), p.z());
        maxi.z() = std::max(maxi.z(), p.z());
    }
    return {mini, maxi};
}

BezierCurve& BezierCurve::scale(const Vector3 &factor)
{
    for (auto& p : points)
        p *= factor;
    for (auto& p : handles)
        p *= factor;
    return *this;
}

BezierCurve& BezierCurve::translate(const Vector3 &translation)
{
    for (auto& p : points)
        p += translation;
    for (auto& p : handles)
        p += translation;
    return *this;
}

BezierCurve& BezierCurve::removeDuplicates()
{
    for (int i = int(numPoints() - 1); i > 0; i--) {
        if (points[i] == points[i - 1]) {
            points.erase(points.begin() + (i - 1));
            handles[2 * i - 2] = (handles[2 * i] + handles[2 * i - 2]) / 2.f;
            handles[2 * i - 1] = (handles[2 * i + 1] + handles[2 * i - 1]) / 2.f;
        }
    }
    return *this;
}

size_t BezierCurve::numPoints() const
{
    return this->points.size();
}

std::string BezierCurve::toString() const
{
    std::ostringstream out;
    out << "Bezier with " << this->points.size() << " points (" << (closed ? "closed" : "not closed") << ") :\n";
    for (auto& p : this->points)
        out << "- " << p << "\n";
    return out.str();
}

BezierCurve& BezierCurve::close()
{
    this->closed = true;
    return *this;
}



Vector3& BezierCurve::operator[](size_t i) {
    return this->points[i];
}
const Vector3& BezierCurve::operator[](size_t i) const {
    return this->points[i];
}
