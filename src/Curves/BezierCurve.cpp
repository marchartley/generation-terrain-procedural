#include "BezierCurve.h"

#include "Utils/Collisions.h"

#include "GUIElements/ImageViewer.h"

BezierCurve::BezierCurve() : BezierCurve(std::vector<Vector3>())
{}

BezierCurve::BezierCurve(const std::vector<Vector3> &points)
    : BezierCurve(points, std::vector<Vector3>())
{
    handles = std::vector<Vector3>(2 * (points.size()));
    for (size_t i = 0; i < points.size(); i++) {
        // Vector3 proj = Collision::projectPointOnSegment(points[i], points[i - 1], points[i + 1]);
        // handles[2 * i - 1] = points[i] - (proj - points[i - 1]).normalize() * (points[i] - points[i - 1]).norm() / 3.f;
        // handles[2 * i] = points[i] + (points[i + 1] - proj).normalize() * (points[i + 1] - points[i]).norm() / 3.f;
        this->autosmooth(i);
    }
    // Strategy 1 : tangents towards next / previous points
    // handles.front() = points[0] + (points[1] - points[0]) / 3.f;
    // handles.back() = points[numSegments()] + (points[numPoints() - 2] - points[numSegments()]) / 3.f;

    // Strategy 2 : tangent computation use next/previous 2 points
    // Start endpoint: use direction from p0 based on p1 -> p2
    // Vector3 t0 = points[2] - points[1];
    // handles.front() = points[0] - t0.normalize() * (points[1] - points[0]).norm() / 3.f;
    // End endpoint: use direction based on p[n-3] -> p[n-2]
    // size_t n = points.size();
    // Vector3 tn = points[n - 2] - points[n - 3];
    // handles.back() = points[n - 1] + tn.normalize() * (points[n - 1] - points[n - 2]).norm() / 3.f;
}

BezierCurve::BezierCurve(const std::vector<Vector3>& points, const std::vector<Vector3>& handles)
    : Curve(), points(points), handles(handles)
{
}

// BezierCurve::BezierCurve(const BezierCurve &curve)
//     : BezierCurve(curve.points, curve.handles)
// {}

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

inline Vector3 BezierCurve::linearBezier(const Vector3& P0, const Vector3& P1, float t) {
    return P0 * (1.f - t) + P1 * t;
}
inline Vector3 BezierCurve::quadraticBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, float t) {
    return linearBezier(P0, P1, t) * (1.f - t) + linearBezier(P1, P2, t) * t;
}
inline Vector3 BezierCurve::cubicBezier(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    // return linearBezier(quadraticBezier(P0, P1, P2, t), quadraticBezier(P1, P2, P3, t), t);
    return std::pow(1.f-t, 3) * P0 + 3.f * std::pow(1.f - t, 2) * t * P1 + 3.f * (1.f - t) * t*t * P2 + t*t*t * P3;
}
inline Vector3 BezierCurve::cubicBezierDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    return 3.f * std::pow(1.f - t, 2) * (P1 - P0) + 6.f * (1.f - t) * t * (P2 - P1) + 3.f * t * t * (P3 - P2);
}
inline Vector3 BezierCurve::cubicBezierSecondDerivative(const Vector3& P0, const Vector3& P1, const Vector3& P2, const Vector3& P3, float t) {
    return 6.f * (1.f - t) * (P0 + P2 - P1 * 2.f) + 6.f * t * (P1 + P3 - P2 * 2.f);
}
Vector3 BezierCurve::getPoint(float x) const
{
    int pointId = timeToLowerIndex(x);
    int nextPointId = pointIndex(pointId + 1);
    const auto& P0 = points[pointId];
    const auto& P1 = handlePos(handleOut(pointId));
    const auto& P2 = handlePos(handleIn(nextPointId));
    const auto& P3 = points[nextPointId];

    float t = (x * float(numSegments())) - pointId;
    return cubicBezier(P0, P1, P2, P3, t);
}

Vector3 falseDerivativeOnTime0(const BezierCurve& curve) {
    Vector3 possibleDerivative = -curve.handles[curve.handleIn(0)];
    if (possibleDerivative.norm2() == 0) possibleDerivative = curve.handles[curve.handleOut(0)];
    if (possibleDerivative.norm2() == 0) possibleDerivative = (curve.getPoint(.01f) - curve.get(0)) / .01f;
    return possibleDerivative;
}
Vector3 falseDerivativeOnTime1(const BezierCurve& curve) {
    Vector3 possibleDerivative = -curve.handles[curve.handleOut(-1)];
    if (possibleDerivative.norm2() == 0) possibleDerivative = curve.handles[curve.handleIn(-1)];
    if (possibleDerivative.norm2() == 0) possibleDerivative = (curve.get(-1) - curve.getPoint(.99f)) / .01f;
    return possibleDerivative;
}
Vector3 BezierCurve::getDerivative(float x, bool normalize) const
{
    if (x <= 0.f && !isClosed()) {
        Vector3 res = falseDerivativeOnTime0(*this);
        return (normalize ? res.normalize() : res);
    }
    else if (x >= 1.f && !isClosed()) {
        Vector3 res = falseDerivativeOnTime1(*this);
        return (normalize ? res.normalize() : res);
    }
    int pointId = timeToLowerIndex(x);
    int nextPointId = pointIndex(pointId + 1);
    const auto& P0 = points[pointId];
    const auto& P1 = handlePos(handleOut(pointId));
    const auto& P2 = handlePos(handleIn(nextPointId));
    const auto& P3 = points[nextPointId];

    float t = (x * float(numSegments())) - pointId;
    Vector3 res = cubicBezierDerivative(P0, P1, P2, P3, t);
    return (normalize ? res.normalize() : res);
}

Vector3 BezierCurve::getSecondDerivative(float x, bool normalize) const
{
    int pointId = timeToLowerIndex(x);
    int nextPointId = pointIndex(pointId + 1);
    const auto& P0 = points[pointId];
    const auto& P1 = handlePos(handleOut(pointId));
    const auto& P2 = handlePos(handleIn(nextPointId));
    const auto& P3 = points[nextPointId];

    float t = (x * float(numSegments())) - pointId;
    Vector3 res = cubicBezierSecondDerivative(P0, P1, P2, P3, t);
    return (normalize ? res.normalize() : res);
}

float BezierCurve::estimateClosestTime(const Vector3 &pos) const
{
    int segmentCount = numSegments();
    int bestSegment = 0;
    float bestT = 0.f;
    float bestDistSq = std::numeric_limits<float>::max();

    constexpr int samplesPerSegment = 20;

    // Coarse search
    for (int segment = 0; segment < segmentCount; ++segment)
    {
        const auto point0 = get(segment);
        const auto point1 = handlePos(handleOut(segment));
        const auto point2 = handlePos(handleIn(segment + 1));
        const auto point3 = get(segment + 1);
        for (int i = 0; i <= samplesPerSegment; ++i)
        {
            float t = float(i) / float(samplesPerSegment);
            const Vector3 p = cubicBezier(point0, point1, point2, point3, t);

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

    const auto bestPoint0 = get(bestSegment);
    const auto bestPoint1 = handlePos(handleOut(bestSegment));
    const auto bestPoint2 = handlePos(handleIn(bestSegment + 1));
    const auto bestPoint3 = get(bestSegment + 1);
    for (int i = 0; i < iterations; ++i)
    {
        const Vector3 p = cubicBezier(bestPoint0, bestPoint1, bestPoint2, bestPoint3, bestT);
        const Vector3 d1 = cubicBezierDerivative(bestPoint0, bestPoint1, bestPoint2, bestPoint3, bestT);
        const Vector3 d2 = cubicBezierSecondDerivative(bestPoint0, bestPoint1, bestPoint2, bestPoint3, bestT);

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

Vector3 &BezierCurve::get(int i) {
    return this->points[pointIndex(i)];
}

Vector3 BezierCurve::get(int i) const {
    return this->points[pointIndex(i)];
}

/*size_t BezierCurve::pointIndex(int i) const {
    return (i + numPoints()) % numPoints();
}*/

BezierCurve& BezierCurve::setPoint(int _i, const Vector3 &newPos)
{
    size_t i = pointIndex(_i);
    Vector3 translation = newPos - points[i];
    this->points[i] = newPos;
    // this->handles[handleIn(i)] += translation;
    // this->handles[handleOut(i)] += translation;
    return *this;
}

int computeNewNbMultiplicatorForResampling(int initialNb, int wantedNb) {
    int k = (wantedNb - 1) / (initialNb - 1);
    return k;
}
int computeNewNbPointsForResampling(int initialNb, int wantedNb) {
    return computeNewNbMultiplicatorForResampling(initialNb, wantedNb) * (initialNb - 1) + 1;
}
BezierCurve& BezierCurve::resamplePoints(int newNbPoints)
{
    if (newNbPoints < 0 || newNbPoints == numPoints()) return *this;
    int slicesPerSegment = (newNbPoints - (isClosed() ? 0 : 1)) / numSegments();
    int realNewNbPoints = (slicesPerSegment - 1) * numSegments() + numPoints();
    newNbPoints = realNewNbPoints;
    int newNbSegments = newNbPoints - (isClosed() ? 0 : 1);
    std::vector<Vector3> newPoints(newNbPoints);
    std::vector<Vector3> newHandles(newNbPoints * 2);

    for (int i = 0; i < newNbPoints; i++) {
        float t = float(i) / (float(newNbSegments));
        newPoints[i] = getPoint(t);
    }
    for (int i = 0; i < numPoints(); i++) {
        newHandles[BezierCurve::handleIn(i * slicesPerSegment, newNbPoints, isClosed())] = handlePos(handleIn(i));
        newHandles[BezierCurve::handleOut(i * slicesPerSegment, newNbPoints, isClosed())] = handlePos(handleOut(i));
    }
    for (int i = 1; i < newNbPoints; i++) {
        float t = float(i) / (float(newNbSegments));
        int currentIdx = timeToLowerIndex(t);
        int nextIdx = pointIndex(currentIdx + 1);
        if (i % slicesPerSegment == 0) {
            newHandles[BezierCurve::handleIn(i, newNbPoints, isClosed())] = handlePos(handleIn(i / slicesPerSegment));
            newHandles[BezierCurve::handleOut(i, newNbPoints, isClosed())] = handlePos(handleOut(i / slicesPerSegment));
            continue;
        }
        float u = 1.f / float((slicesPerSegment) - (i % slicesPerSegment) + 1);

        const Vector3 P0 = newPoints[i - 1];
        const Vector3 P1 = newHandles[BezierCurve::handleOut(i - 1, newNbPoints, isClosed())];
        const Vector3 P2 = handlePos(handleIn(nextIdx));
        const Vector3 P3 = points[nextIdx];
        const Vector3 Q0 = P0 + (P1 - P0) * u;
        const Vector3 Q1 = P1 + (P2 - P1) * u;
        const Vector3 Q2 = P2 + (P3 - P2) * u;
        const Vector3 R0 = Q0 + (Q1 - Q0) * u;
        const Vector3 R1 = Q1 + (Q2 - Q1) * u;

        newHandles[BezierCurve::handleOut(i - 1, newNbPoints, isClosed())] = Q0;
        newHandles[BezierCurve::handleIn(i, newNbPoints, isClosed())] = R0;
        newHandles[BezierCurve::handleOut(i, newNbPoints, isClosed())] = R1;
        handles[handleIn(nextIdx)] = Q2 - get(nextIdx);
    }
    if (!isClosed()) {
        newHandles[BezierCurve::handleIn(-1, newNbPoints, isClosed())] = handlePos(handleIn(-1));
        newHandles[BezierCurve::handleOut(-1, newNbPoints, isClosed())] = handlePos(handleOut(-1));
    } else {
        newHandles[0] = handlePos(0);
    }
    for (size_t i = 0; i < newNbPoints; i++) {
        newHandles[BezierCurve::handleIn(i, newNbPoints, isClosed())] -= newPoints[i];
        newHandles[BezierCurve::handleOut(i, newNbPoints, isClosed())] -= newPoints[i];
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


AABBox BezierCurve::getAABBox() const
{
    std::vector<Vector3> handlePositions(handles.size());
    for (size_t i = 0; i < handlePositions.size(); i++) {
        handlePositions[i] = handlePos(i);
    }
    return AABBox(points).expand(handlePositions);
    /*
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
    */
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
    return *this;
}

BezierCurve& BezierCurve::removeDuplicates()
{
    for (int i = int(numSegments()); i > 0; i--) {
        if (points[pointIndex(i)] == points[i - 1]) {
            points.erase(points.begin() + (i - 1));
            handles[2 * i - 2] = (handles[handleIndex(2 * i)] + handles[2 * i - 2]) / 2.f;
            handles[2 * i - 1] = (handles[handleIndex(2 * i + 1)] + handles[2 * i - 1]) / 2.f;
        }
    }
    return *this;
}

size_t BezierCurve::size() const {
    return numPoints();
}

size_t BezierCurve::numPoints() const
{
    return this->points.size();
}

std::vector<Vector3>::const_iterator BezierCurve::begin() const {
    return points.begin();
}

std::vector<Vector3>::const_iterator BezierCurve::end() const {
    return points.end();
}

std::vector<Vector3>::iterator BezierCurve::begin() {
    return points.begin();
}

std::vector<Vector3>::iterator BezierCurve::end() {
    return points.end();
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
    Curve::close();
    // if (this->points.size() > 0 && this->points.front() != this->points.back()) {
        // handles.resize(handles.size() + 2);

        /*
        // Handle connected to first (newly last) point
        auto lastHandle = points.front() - handles.front();
        // scaling
        lastHandle *= (points.back() - points.front()).norm() / (points[1] - points.front()).norm();
        this->handles[handles.size() - 1] = points.front() + lastHandle;
        */

        // Or maybe apply autosmooth to the first point...
        // this->autosmooth(0);



        // recompute handles for previously last point
        // this->autosmooth(-1);
        // this->autosmooth(-2);
        // Add the point
        // this->points.push_back(points.front());

    // }
    return *this;
}

BezierCurve &BezierCurve::reset() {
    points.clear();
    handles.clear();
    return *this;
}

size_t BezierCurve::pointIndex(int index, int nbPoints, bool closed) {
    return (index + nbPoints) % nbPoints;
}
size_t BezierCurve::handleIndex(int index, int nbPoints, bool closed)
{
    return (index + (2 * nbPoints)) % (2 * nbPoints);
}

size_t BezierCurve::handleIn(int pointIdx, int nbPoints, bool closed)
{
    return handleIndex(2 * pointIdx, nbPoints, closed);
}

size_t BezierCurve::handleOut(int pointIdx, int nbPoints, bool closed)
{
    return handleIndex(2 * pointIdx + 1, nbPoints, closed);
}

size_t BezierCurve::handleIndex(int index) const
{
    return BezierCurve::handleIndex(index, numPoints(), closed);
}

size_t BezierCurve::handleIn(int pointIdx) const
{
    return handleIn(pointIdx, numPoints(), closed);
}

size_t BezierCurve::handleOut(int pointIdx) const
{
    return handleOut(pointIdx, numPoints(), closed);
}

size_t BezierCurve::pointIndexFromHandleIndex(int handleIdx) const
{
    return handleIndex(handleIdx) / 2;
}

Vector3 BezierCurve::handlePos(int handleIdx) const
{
    return points[pointIndexFromHandleIndex(handleIdx)] + handles[handleIndex(handleIdx)];
}

std::vector<Vector3> BezierCurve::getPoints() const {
    return points;
}

std::vector<Vector3> BezierCurve::getHandles() const {
    return handles;
}

BezierCurve& BezierCurve::autosmooth(int pointIdx)
{
    const size_t currentIdx = pointIndex(pointIdx);
    const size_t nextIdx = pointIndex(pointIdx + 1);
    const size_t prevIdx = pointIndex(pointIdx - 1);
    const auto& P = points[currentIdx];
    const auto& P_next = points[nextIdx]; //(nextIdx == 0 ? nextIdx + 1 : nextIdx)];
    const auto& P_prev = points[prevIdx];
    // Vector3 proj = Collision::projectPointOnLine(P, P_prev, P_next);

    Vector3 startToPoint = P - P_prev;
    Vector3 segment = P_next - P_prev;
    float t = (segment.norm2() > 0 ? startToPoint.dot(segment) / segment.dot(segment) : 0.f);
    auto proj = P_prev + segment * t;
    if ((currentIdx != 0 && nextIdx != 0) || isClosed())
        handles[handleIn(pointIdx)] = -(proj - P_prev).normalize() * (P - P_prev).norm() / 3.f * (t < 0 ? -1.f : 1.f);
    else
        handles[handleIn(pointIdx)] = Vector3::origin;

    if ((currentIdx != 0 && nextIdx != 0) || isClosed())
        handles[handleOut(pointIdx)] = (P_next - proj).normalize() * (P_next - P).norm() / 3.f * (t > 1 ? -1.f : 1.f);
    else
        handles[handleOut(pointIdx)] = Vector3::origin;

    return *this;
}

BezierCurve& BezierCurve::autosmooth()
{
    if (numPoints() <= 2) return *this;

    for (size_t i = 0; i < numPoints(); i++) {
        autosmooth(i);
    }
    return *this;
}

std::vector<std::shared_ptr<Curve>> BezierCurve::slice(float t) const
{
    int lowerIndex = timeToLowerIndex(t);
    int nextIdx = pointIndex(lowerIndex + 1);
    int firstHalfNbPoints = lowerIndex + 2;
    int secondHalfNbPoints = (numPoints() - (lowerIndex + 1)) + 1;


    std::vector<Vector3> firstHalfPoints(firstHalfNbPoints);
    std::vector<Vector3> firstHalfHandles(firstHalfNbPoints * 2);

    std::vector<Vector3> secondHalfPoints(secondHalfNbPoints);
    std::vector<Vector3> secondHalfHandles(secondHalfNbPoints * 2);

    for (int i = 0; i < lowerIndex + 1; i++) {
        firstHalfPoints[i] = points[pointIndex(i)];
        firstHalfHandles[BezierCurve::handleIn(i, firstHalfNbPoints, false)] = handles[handleIn(i)];
        firstHalfHandles[BezierCurve::handleOut(i, firstHalfNbPoints, false)] = handles[handleOut(i)];
    }
    for (int i = 0; i < secondHalfNbPoints - 1; i++) {
        secondHalfPoints[i + 1] = points[pointIndex(i + lowerIndex + 1)];
        secondHalfHandles[BezierCurve::handleIn(i + 1, secondHalfNbPoints, false)] = handles[handleIn(i + lowerIndex + 1)];
        secondHalfHandles[BezierCurve::handleOut(i + 1, secondHalfNbPoints, false)] = handles[handleOut(i + lowerIndex + 1)];
    }

    float u = ::map(t, float(lowerIndex) / float(numSegments()), float(nextIdx) / float(numSegments()), 0.f, 1.f);

    const Vector3 P0 = points[lowerIndex];
    const Vector3 P1 = handlePos(handleOut(lowerIndex));
    const Vector3 P2 = handlePos(handleIn(nextIdx));
    const Vector3 P3 = points[nextIdx];
    const Vector3 Q0 = P0 + (P1 - P0) * u;
    const Vector3 Q1 = P1 + (P2 - P1) * u;
    const Vector3 Q2 = P2 + (P3 - P2) * u;
    const Vector3 R0 = Q0 + (Q1 - Q0) * u;
    const Vector3 R1 = Q1 + (Q2 - Q1) * u;
    const Vector3 S  = R0 + (R1 - R0) * u;

    firstHalfPoints[BezierCurve::pointIndex(-1, firstHalfNbPoints, false)] = S;
    firstHalfHandles[BezierCurve::handleOut(-2, firstHalfNbPoints, false)] = Q0 - P0;
    firstHalfHandles[BezierCurve::handleIn(-1, firstHalfNbPoints, false)] = R0 - S;
    firstHalfHandles[BezierCurve::handleOut(-1, firstHalfNbPoints, false)] = R1 - S;

    secondHalfPoints[0] = S;
    secondHalfHandles[handleIn(0)] = R0 - S;
    secondHalfHandles[handleOut(0)] = R1 - S;
    if (this->isClosed()) {
        secondHalfPoints.push_back(points[0]);
        secondHalfHandles.push_back(handles[0]);
        secondHalfHandles.push_back(handles[1]);
    }
    if (secondHalfPoints.size() > 1) {
        secondHalfHandles[handleIn(1)] = Q2 - P3;
    }

    // ImageViewer::get().addScatter({P0, P1, P2, P3, Q0, Q1, Q2, R0, R1, S}, "", {}, Vector3::black);

    BezierCurve firstCurve(firstHalfPoints, firstHalfHandles);
    BezierCurve secondCurve(secondHalfPoints, secondHalfHandles);
    return {std::make_shared<BezierCurve>(firstCurve), std::make_shared<BezierCurve>(secondCurve)};
}

std::vector<std::shared_ptr<Curve>> BezierCurve::slice(const std::vector<float>& _ts) const
{
    auto ts = _ts;
    std::sort(ts.begin(), ts.end(), std::greater<float>()); // sort in descending order, just to optimize the poping

    std::vector<std::shared_ptr<Curve>> subCurves;
    // auto original = *this;
    std::shared_ptr<BezierCurve> remaining(this->clone());
    if (remaining->isClosed()) {
        remaining->setClosed(false);
        remaining->points.push_back(points[0]);
        remaining->handles.push_back(handles[0]);
        remaining->handles.push_back(handles[1]);
    }
    float previousT = -1.0f;
    int nbPointsPassed = 0;
    while (!ts.empty()) {
        auto t = ts.back();
        ts.pop_back();

        int prevIndex = this->timeToLowerIndex(t);
        bool previousIsSlice = (this->indexToTime(prevIndex) < previousT);
        float oldTime = (previousIsSlice ? previousT : this->indexToTime(prevIndex));

        float remapedU = (t - oldTime) / (this->indexToTime(prevIndex + 1) - oldTime);
        float localT = ::map(remapedU, 0.f, 1.f, remaining->indexToTime(prevIndex - nbPointsPassed), remaining->indexToTime((prevIndex - nbPointsPassed) + 1));

        nbPointsPassed += remaining->timeToLowerIndex(localT);

        auto subs = remaining->slice(localT);
        subCurves.push_back(subs[0]);
        remaining = std::dynamic_pointer_cast<BezierCurve>(subs[1]);

        previousT = t;
    }
    subCurves.push_back(remaining);
    return subCurves;
}



Vector3& BezierCurve::operator[](size_t i) {
    return this->points[i];
}

void BezierCurve::addPoint(const Vector3 &newPoint) {
    this->points.push_back(newPoint);
    this->handles.resize(handles.size() + 2);
    this->autosmooth(-2);
    this->autosmooth(-1);
}

BezierCurve &BezierCurve::insertPoint(int i, const Vector3 &newPos) {
    int newIdx = pointIndex(i);
    this->points.insert(points.begin() + newIdx, newPos);
    this->handles.insert(handles.begin() + handleIn(newIdx + 1), {Vector3(), Vector3()});
    autosmooth(newIdx);
    return *this;
}

BezierCurve &BezierCurve::removePoint(int i) {
    int newIdx = pointIndex(i);
    int handleIdx = handleIn(newIdx);
    this->handles.erase(handles.begin() + handleIdx, handles.begin() + handleIdx + 1); // ???? OR 2
    this->points.erase(points.begin() + newIdx);
    return *this;
}

const Vector3& BezierCurve::operator[](size_t i) const {
    return this->points[i];
}
