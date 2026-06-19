#include "Curves/CatmullRomSpline.h"

#include "Utils/Utils.h"
#include "Utils/Collisions.h"

#include <exception>
#include <sstream>

CatmullRomSpline::CatmullRomSpline()
{

}
CatmullRomSpline::CatmullRomSpline(const std::vector<Vector3>& points)
    : points(points)
{
}
CatmullRomSpline::CatmullRomSpline(std::vector<CatmullRomSpline> subsplines)
{
    for (CatmullRomSpline& spline : subsplines) {
        bool ignoreFirstPoint = (this->points.empty() ? false : this->points.back() == spline.front());
        this->points.insert(this->points.end(), spline.points.begin() + (ignoreFirstPoint ? 1 : 0), spline.points.end());
    }
}
CatmullRomSpline::CatmullRomSpline(const Curve& curve)
    : CatmullRomSpline(curve.getPath())
{}


CatmullRomSpline& CatmullRomSpline::reverseVertices()
{
    std::reverse(this->points.begin(), this->points.end());
    return *this;
}

CatmullRomSpline::operator bool() const {
    return (this->points.size() > 0);
}

std::vector<Vector3> CatmullRomSpline::getPath(int numberOfPoints) const
{
    if (numberOfPoints < 0 || numberOfPoints == numPoints()) return this->points;
    /// I'm really not sure this is the best solution, but an easy fix :
    /// forcing user to have at least 2 points
    numberOfPoints = std::max(numberOfPoints, 2);

    std::vector<Vector3> path;
    for (int i = 0; i < numberOfPoints; i ++) {
        float t = float(i) / (float)(numberOfPoints - 1);
        path.push_back(this->getPoint(t));
    }
    return path;
}

Vector3 CatmullRomSpline::getPoint(float x) const
{
    if (this->closed) {
        x = x - std::floor(x); // Warp around if x < 0 or x > 1
    } else {
        x = std::clamp(x, 0.f, 1.f);
    }
    if (points.size() > 2)
        return this->getCatmullPoint(x);

    if(points.size() == 0)
        return Vector3();

    std::vector<Vector3> controls = points;
    while (controls.size() > 1)
    {
        std::vector<Vector3> newCtrls;
        for (size_t i = 0; i < controls.size() - 1; i++) {
            newCtrls.push_back(Curve::linear(x, controls[i], controls[i+1]));
        }
        controls = newCtrls;
    }
    return controls[0];
}

Vector3 CatmullRomSpline::getDerivative(float x, bool normalize) const
{
    /*float previousTime = std::clamp(x - 0.01f, 0.f, 1.f);
    float nextTime = std::clamp(x + 0.01f, 0.f, 1.f);
    float e = nextTime - previousTime; // Case for start/end of curve

    Vector3 v = (getPoint(nextTime) - getPoint(previousTime));
    return (normalize ? v.normalized() : v / e);
    */
    Vector3 v = std::get<1>(this->pointAndDerivativeAndSecondDerivative(x));
    return (normalize ? v.normalized() : v);
}

Vector3 CatmullRomSpline::getSecondDerivative(float x, bool normalize) const
{
    /*float previousTime = std::clamp(x - 0.01f, 0.f, 1.f);
    float nextTime = std::clamp(x + 0.01f, 0.f, 1.f);
    float e = nextTime - previousTime; // Case for start/end of curve

    Vector3 v = (getDerivative(nextTime) - getDerivative(previousTime));
    return (normalize ? v.normalized() : v / e);
    */
    Vector3 v = std::get<2>(this->pointAndDerivativeAndSecondDerivative(x));
    return (normalize ? v.normalized() : v);
}

/*
 * Everything here is useless...
float getSegmentClosestTime(const BSpline& curve, float t1, float t2, const Vector3& pos) {
    // Returns the parameter t of the closest point on the curve segment between t1 and t2 to the point pos.
    Vector3 pointAtT1 = curve.getPoint(t1);
    Vector3 pointAtT2 = curve.getPoint(t2);
    Vector3 line = pointAtT2 - pointAtT1;
    float tClosest = ((pos - pointAtT1).dot(line)) / line.norm2();
    return std::clamp(tClosest, 0.0f, 1.0f);  // Ensure tClosest is within [0, 1] relative to t1 and t2
}

float closestTimeSubdivision(const BSpline& curve, const Vector3& pos, float tStart = 0.0f, float tEnd = 1.0f, int depth = 0, const int maxDepth = 10, float epsilon = 0.001f) {
    // Calculate the midpoint and its position
    float tMid = (tStart + tEnd) / 2.0f;
    Vector3 midPoint = curve.getPoint(tMid);

    // Check if the distance at tMid is within the acceptable epsilon range
    float midDistance = (midPoint - pos).norm2();
    if (midDistance < epsilon || depth >= maxDepth) {
        // If the closest distance is less than epsilon or the maximum depth is reached, return tMid
        return tMid;
    }

    // Subdivide the curve and check the distance at the midpoints
    float tLeftMid = (tStart + tMid) / 2.0f;
    Vector3 leftMidPoint = curve.getPoint(tLeftMid);

    float tRightMid = (tMid + tEnd) / 2.0f;
    Vector3 rightMidPoint = curve.getPoint(tRightMid);

    // Compute distances to pos
    float startDistance = (curve.getPoint(tStart) - pos).norm2();
    float leftMidDistance = (leftMidPoint - pos).norm2();
    float rightMidDistance = (rightMidPoint - pos).norm2();
    float endDistance = (curve.getPoint(tEnd) - pos).norm2();

    // Use a lambda function to find the minimum distance and its corresponding parameter
    auto [minDistance, minT] = std::min({std::make_pair(startDistance, tStart),
                                         std::make_pair(leftMidDistance, tLeftMid),
                                         std::make_pair(midDistance, tMid),
                                         std::make_pair(rightMidDistance, tRightMid),
                                         std::make_pair(endDistance, tEnd)},
                                        [](const auto& a, const auto& b) { return a.first < b.first; });

    // If the closest distance is within the epsilon threshold, return the corresponding t value
    if (minDistance < epsilon) {
        return minT;
    }

    // Recur on the appropriate half of the curve
    if (minT <= tLeftMid) {
        return closestTimeSubdivision(curve, pos, tStart, tMid, depth + 1, maxDepth, epsilon);
    } else if (minT <= tRightMid) {
        return closestTimeSubdivision(curve, pos, tLeftMid, tRightMid, depth + 1, maxDepth, epsilon);
    } else {
        return closestTimeSubdivision(curve, pos, tMid, tEnd, depth + 1, maxDepth, epsilon);
    }
}


float closestTimeIntervalAnalysis(const BSpline& curve, const Vector3& pos, float tolerance = 0.001f) {
    float lowerBound = 0.0f;
    float upperBound = 1.0f;
    float closestTime = 0.0f;
    float closestDistance = std::numeric_limits<float>::max();

    while ((upperBound - lowerBound) > tolerance) {
        float mid = (lowerBound + upperBound) / 2.0f;
        float interval = (upperBound - lowerBound) / 4.0f;

        // Sample the curve at the mid-point and quarter points within the current interval
        std::vector<float> sampleTimes = {lowerBound, lowerBound + interval, mid, upperBound - interval, upperBound};
        for (float t : sampleTimes) {
            Vector3 samplePoint = curve.getPoint(t);
            float distance = (samplePoint - pos).norm2(); // Using squared distance to avoid sqrt computation

            if (distance < closestDistance) {
                closestDistance = distance;
                closestTime = t;
            }
        }

        // Narrow down the interval to the one containing the closestTime
        if (closestTime <= mid) {
            upperBound = mid;
        } else {
            lowerBound = mid;
        }
    }

    return closestTime;
}


struct KDNode {
    Vector3 point;
    float t; // parameter t corresponding to the point on the spline
    std::shared_ptr<KDNode> left;
    std::shared_ptr<KDNode> right;
    int axis;

    KDNode(const Vector3& pt, float param, int ax)
        : point(pt), t(param), axis(ax), left(nullptr), right(nullptr) {}
};

class KDTree {
private:
    std::shared_ptr<KDNode> root;

    std::shared_ptr<KDNode> insert(std::shared_ptr<KDNode> node, const Vector3& point, float t, int depth) {
        if (!node) {
            return std::make_shared<KDNode>(point, t, depth % 3);
        }

        int axis = depth % 3;
        if (point[axis] < node->point[axis]) {
            node->left = insert(std::move(node->left), point, t, depth + 1);
        } else {
            node->right = insert(std::move(node->right), point, t, depth + 1);
        }

        return node;
    }

    void nearest(const KDNode* node, const Vector3& point, float& bestDist, float& bestT, int depth) const {
        if (!node) return;

        float d = (node->point - point).norm2();
        if (d < bestDist) {
            bestDist = d;
            bestT = node->t;
        }

        int axis = depth % 3;
        float delta = point[axis] - node->point[axis];
        float delta2 = delta * delta;

        const KDNode* near = delta < 0 ? node->left.get() : node->right.get();
        const KDNode* far = delta < 0 ? node->right.get() : node->left.get();

        nearest(near, point, bestDist, bestT, depth + 1);

        // if there might be a closer point on the other side of the splitting plane
        if (delta2 < bestDist) {
            nearest(far, point, bestDist, bestT, depth + 1);
        }
    }

public:
    void insert(const Vector3& point, float t) {
        root = insert(std::move(root), point, t, 0);
    }

    float nearest(const Vector3& point) const {
        float bestDist = std::numeric_limits<float>::max();
        float bestT = -1;
        nearest(root.get(), point, bestDist, bestT, 0);
        return bestT;
    }
};

float closestTimeSpatialIndex(const BSpline& curve, const Vector3& pos, int initialSamples = 100, float refinementThreshold = 0.001f) {
    // Create the spatial index
    KDTree index;

    // Sample the curve at regular intervals and add the points to the index
    for (int i = 0; i <= initialSamples; ++i) {
        float t = i / static_cast<float>(initialSamples);
        Vector3 point = curve.getPoint(t);
        index.insert(point, t);
    }

    return index.nearest(pos);
}*/

float CatmullRomSpline::estimateClosestTime(const Vector3& pos) const
{
    if (this->points.size() == 0) {
        return 0;
    } else if (this->points.size() == 1) {
        return 0;
    } else if (this->points.size() == 2) {
        Vector3 line = (this->points[1] - this->points[0]);
        float time = clamp((pos - this->points[0]).dot(line) / line.dot(line), 0.f, 1.f);
        return time;
    }

    const float initialEpsilon = 1e-4;
    const float nbChecksFactor = 2.f;

    int numberOfChecks = std::max(8, int(points.size() * nbChecksFactor));
    float step0 = 2.0f / numberOfChecks;

    float closestTime = 0.0f;

    float step = step0;
    float center = 0.5f;

    while (step > initialEpsilon) {
        float minD2 = std::numeric_limits<float>::max();
        float halfSpan = 0.5f * numberOfChecks * step;
        float tMin = std::max(0.0f, center - halfSpan);
        float tMax = std::min(1.0f, center + halfSpan);

        Vector3 previousPoint = getPoint(tMin);
        // Vector3 nextPoint;

        // sample across [tMin, tMax]
        for (int i = 0; i <= numberOfChecks; ++i) {
            float a = float(i) / float(numberOfChecks);
            float t = tMin + a * (tMax - tMin);
            // float a1 = float(i-1) / float(numberOfChecks);
            // float t1 = tMin + a1 * (tMax - tMin);
            float a2 = float(i+1) / float(numberOfChecks);
            float t2 = tMin + a2 * (tMax - tMin);

            const Vector3 nextPoint = getPoint(t2);
            auto projection = Collision::projectPointOnSegment(pos, previousPoint, nextPoint);
            // auto projection = getPoint(t);
            float d2 = (projection - pos).norm2();
            if (d2 < minD2) {
                minD2 = d2;
                closestTime = t;
            }

            previousPoint = nextPoint;
        }

        center = closestTime;
        step *= step0; //step0; // or step *= 0.5f; (often more stable)
    }

    // Second method: get distance to bspline points directly, return closest
    // This works only because the BSpline go through all the points.
    /*for (int i = 0; i < this->size(); i++) {
        float dist = (this->points[i] - pos).norm2();
        if (dist < minD2) {
            minD2 = dist;
            closestTime = float(i) / float(this->size() - 1);
        }
    }*/
    return closestTime;
}

Vector3 CatmullRomSpline::estimateClosestPos(const Vector3& pos) const
{
    return this->getPoint(this->estimateClosestTime(pos));
}

float CatmullRomSpline::estimateSqrDistanceFrom(const Vector3& pos) const
{
    return (this->estimateClosestPos(pos) - pos).norm2();
}

float CatmullRomSpline::length() const
{
    // Should be using approximation of arclength instead
    float length = 0;
    if (this->points.empty()) return length;
    for (size_t i = 0; i < this->points.size() - 1; i++) {
        length += (this->points[i] - this->points[i + 1]).norm();
    }
    return length;
}

std::vector<Vector3> CatmullRomSpline::getPoints() const {
    return this->points;
}

Vector3 &CatmullRomSpline::get(int i) {
    return this->points[pointIndex(i)];
}

Vector3 CatmullRomSpline::get(int i) const {
    return this->points[pointIndex(i)];
}

CatmullRomSpline CatmullRomSpline::smooth(float factor) const
{
    CatmullRomSpline newCurve = *this;
    for (size_t i = 0; i < this->size(); i++) {
        if (i == 0 || i == this->size() - 1) continue;

        newCurve[i] = (*this)[i] - factor * ((*this)[i] - ((*this)[i+1] + (*this)[i-1]) * .5f);
    }
    return newCurve;
}

CatmullRomSpline CatmullRomSpline::taubinSmooth(float factor) const
{
    auto initCurve = *this;
    auto newCurve = *this;
    for (size_t i = 1; i < initCurve.size() - 1; i++) {
        newCurve[i] = initCurve[i] - factor * (initCurve[i] - (initCurve[i+1] + initCurve[i-1]) * .5f).normalized();
    }
    initCurve = newCurve;
    for (size_t i = 1; i < newCurve.size() - 1; i++) {
        initCurve[i] = newCurve[i] + factor * (newCurve[i] - (newCurve[i+1] + newCurve[i-1]) * .5f).normalized();
    }
    return initCurve;
}

CatmullRomSpline& CatmullRomSpline::setPoint(int i, const Vector3 &newPos)
{
    this->points[pointIndex(i)] = newPos;
    return *this;
}

CatmullRomSpline& CatmullRomSpline::resamplePoints(int newNbPoints)
{
    this->cleanPoints();
    if (this->size() == 0) return *this;
    if (newNbPoints < 0)
        newNbPoints = this->points.size();

    std::vector<Vector3> newPoints; //(newNbPoints);

    float totalLength = this->length();

    if (totalLength == 0) {
        this->points = std::vector<Vector3>(newNbPoints, this->points[0]);
        return *this;
    }
    float currentDistance = 0.f;
    int currentPointIndex = 0;
    Vector3 currentPos = points[currentPointIndex];
    float nextObjectiveDistance = totalLength / float(newNbPoints - 1);

    newPoints.push_back(this->points.front());

    while (int(newPoints.size()) < newNbPoints - 1) {
        auto& p1 = points[currentPointIndex + 1];

        float edgeDist = (p1 - currentPos).norm();
        if (currentDistance + edgeDist > nextObjectiveDistance) {
            float t = (nextObjectiveDistance - currentDistance) / edgeDist;
            newPoints.push_back(lerp(t, currentPos, p1));
//            currentDistance += t * edgeDist;
            currentDistance = 0;
            currentPos = newPoints.back();
        } else {
            currentDistance += edgeDist;
            currentPointIndex ++;
            currentPos = p1;
        }
    }
    newPoints.push_back(this->points.back());
    this->points = newPoints;
    return *this;
}

std::tuple<Vector3, Vector3, Vector3> CatmullRomSpline::getFrenetFrame(float x) const
{
    return {getFrenetDirection(x), getFrenetNormal(x), getFrenetBinormal(x)};
}

Vector3 CatmullRomSpline::getFrenetDirection(float x) const
{
    return getDirection(x);
}

Vector3 CatmullRomSpline::getFrenetNormal(float x) const
{
    Vector3 new_dir = this->getFrenetDirection(x);
    Vector3 forward(0, 1, 0);
    Vector3 up(0, 0, 1);
    Vector3 right(1, 0, 0);
    if (!new_dir.isAlmostVertical())
        right = Vector3(0, 0, 1).cross(new_dir);
    else
        right = Vector3(0, 0.01, 1).cross(new_dir);
    return right.normalize();
}

Vector3 CatmullRomSpline::getFrenetBinormal(float x) const
{
    return this->getFrenetDirection(x).cross(this->getFrenetNormal(x)).normalize();
}

Vector3 CatmullRomSpline::getCenterCircle(float x) const
{
    return this->getPoint(x) + this->getNormal(x) * (1 / this->getCurvature(x));
}

CatmullRomSpline& CatmullRomSpline::close()
{
    Curve::close();
    // if (this->points.size() > 0 && this->points.front() != this->points.back())
        // this->points.push_back(points.front());
    return *this;
}

CatmullRomSpline& CatmullRomSpline::cleanPoints()
{
    for (int i = this->size() - 1; i >= 0; i--) {
        if (!this->points[i].isValid()) this->points.erase(this->points.begin() + i);
    }
    return *this;
}

float CatmullRomSpline::CatmullNextT(const Vector3& P0, const Vector3& P1, float t_prev, float alpha)
{
    float norm = std::max(1e-5f, (P0 - P1).norm2());
    return std::pow(norm, alpha * 0.5f) + t_prev;
}

Vector3 CatmullRomSpline::getCatmullPoint(float x) const
{
    Vector3 v = std::get<0>(this->pointAndDerivativeAndSecondDerivative(x));
    return v;
    /*
    alpha /= 2.f;

    std::vector<Vector3> displayedPoints = this->points;
    if (this->closed)
        displayedPoints.push_back(displayedPoints.front());

    size_t lastPointIndex = displayedPoints.size() - 1;
    size_t nbPoints = displayedPoints.size(); // + (this->closed ? 1 : 0);

    if (x == 0.f) return displayedPoints[0];
    if (x == 1.f) return displayedPoints[lastPointIndex];


    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    Vector3 P0 = displayedPoints[(iFloor == 0 ? (this->closed ? int(nbPoints-2) : 1) : iFloor - 1)];
    Vector3 P1 = displayedPoints[iFloor - 0];
    Vector3 P2 = displayedPoints[iCeil + 0];
    Vector3 P3 = displayedPoints[(iCeil >= int(nbPoints-1) ? (this->closed ? 1 : displayedPoints.size()-2) : iCeil + 1)];

    float t0 = 0;
    float t1 = CatmullNextT(P0, P1, t0, alpha);
    float t2 = CatmullNextT(P1, P2, t1, alpha);
    float t3 = CatmullNextT(P2, P3, t2, alpha);

    float t = map(x_prime, 0.f, 1.f, t1, t2);

    Vector3 A1 = P0 * (t1 - t) / (t1 - t0) + P1 * (t - t0) / (t1 - t0);
    Vector3 A2 = P1 * (t2 - t) / (t2 - t1) + P2 * (t - t1) / (t2 - t1);
    Vector3 A3 = P2 * (t3 - t) / (t3 - t2) + P3 * (t - t2) / (t3 - t2);

    Vector3 B1 = A1 * (t2 - t) / (t2 - t0) + A2 * (t - t0) / (t2 - t0);
    Vector3 B2 = A2 * (t3 - t) / (t3 - t1) + A3 * (t - t1) / (t3 - t1);

    Vector3 C  = B1 * (t2 - t) / (t2 - t1) + B2 * (t - t1) / (t2 - t1);
    return C;
    */
}

CatmullRomSpline CatmullRomSpline::simplifyByRamerDouglasPeucker(float epsilon, CatmullRomSpline subspline)
{
    if (subspline.points.empty()) {
        if (this->points.empty()) return *this; // We are just trying to do a simplification from an empty curve, that's pointless
        else subspline = *this;
    }
    if (subspline.points.size() == 1) return subspline;

    // Find farest point from the line going from start to end of this (sub)spline
    Vector3 vecAB = (subspline.front() - subspline.back()).normalized();
    float maxDist = 0;
    int index = -1;
    for (size_t i = 1; i < subspline.points.size() - 1; i++) {
        float d = vecAB.cross((subspline.points[i] - subspline.front())).norm2();
        if (d > maxDist) {
            maxDist = d;
            index = i;
        }
    }
    // Now we split the spline in two subsplines, and apply recursively the algorithm until all points are "close" enough ( dist < epsilon)
    CatmullRomSpline returningSpline;
    if (maxDist > epsilon * epsilon) {
        CatmullRomSpline sub1 = this->simplifyByRamerDouglasPeucker(epsilon, std::vector<Vector3>(subspline.points.begin(), subspline.points.begin() + index));
        CatmullRomSpline sub2 = this->simplifyByRamerDouglasPeucker(epsilon, std::vector<Vector3>(subspline.points.begin() + index, subspline.points.end()));
        returningSpline = CatmullRomSpline({sub1, sub2});
    } else {
        returningSpline.points = {subspline.front(), subspline.back()};
    }
    return returningSpline;
}

AABBox CatmullRomSpline::getAABBox() const
{
    std::vector<Vector3> samples(numSegments());
    for (size_t i = 0; i < samples.size(); i++) {
        float t = float(i + 0.5f) / float(samples.size());
        samples[i] = getPoint(t);
    }
    return AABBox(points).expand(samples);  // I guess good enough approximation, using control points and middle of each subcurve.
    /*
    if (this->points.empty()) return {Vector3::invalid, Vector3::invalid};
    if (this->points.size() == 1) return {points[0], points[0]};

    float minDim = std::numeric_limits<float>::lowest();
    float maxDim = std::numeric_limits<float>::max();
    Vector3 minVec = Vector3(maxDim, maxDim, maxDim),
            maxVec = Vector3(minDim, minDim, minDim);
    for (const auto& point : points) {
        minVec.x() = std::min(point.x(), minVec.x());
        minVec.y() = std::min(point.y(), minVec.y());
        minVec.z() = std::min(point.z(), minVec.z());
        maxVec.x() = std::max(point.x(), maxVec.x());
        maxVec.y() = std::max(point.y(), maxVec.y());
        maxVec.z() = std::max(point.z(), maxVec.z());
    }
    return {minVec, maxVec};
    */
}

CatmullRomSpline& CatmullRomSpline::scale(const Vector3 &factor)
{
    for (auto& vert : this->points)
        vert *= factor;
    return *this;
}
CatmullRomSpline CatmullRomSpline::scaled(float factor)
{
    return this->scaled(Vector3(factor, factor, factor));
}

CatmullRomSpline CatmullRomSpline::scaled(const Vector3& factor)
{
    CatmullRomSpline copy = *this;
    return copy.scale(factor);
}


CatmullRomSpline &CatmullRomSpline::translate(const Vector3& translation)
{
    for (auto& point : points)
        point += translation;
    return *this;
}

CatmullRomSpline &CatmullRomSpline::displacePointsRandomly(float maxDistance)
{
    /*std::vector<Vector3> newPoints(this->size());
    for (int i = 0; i < this->size(); i++) {
        newPoints[i] = points[i] + this->getNormal(float(i) / float(this->size())) * random_gen::generate(-maxDistance, maxDistance);
    }
    this->points = newPoints;
    return *this;*/
    return this->displacePointsRandomly(Vector3(maxDistance, maxDistance, maxDistance));
}

CatmullRomSpline &CatmullRomSpline::displacePointsRandomly(const Vector3 &maxDistance)
{
    std::vector<Vector3> newPoints(this->size());
    for (size_t i = 0; i < this->size(); i++) {
        float x = float(i) / float(this->size());
        auto [pos, dir, normal] = this->pointAndDerivativeAndSecondDerivative(x);
        Vector3 binormal = this->getBinormal(x).normalized();
        dir.normalize();
        normal.normalize();
        newPoints[i] = points[i] + (dir * random_gen::generate(-1, 1) * maxDistance.x()) + (normal * random_gen::generate(-1, 1) * maxDistance.y()) + (binormal * random_gen::generate(-1, 1) * maxDistance.z());
    }
    this->points = newPoints;
    return *this;
}

CatmullRomSpline &CatmullRomSpline::displacePointsRandomlyPerlin(float maxDistance, float scale, bool loop)
{
    return this->displacePointsRandomlyPerlin(Vector3(maxDistance, maxDistance, maxDistance), scale, loop);
}

CatmullRomSpline &CatmullRomSpline::displacePointsRandomlyPerlin(const Vector3& maxDistance, float scale, bool loop)
{
    std::vector<Vector3> newPoints(this->size());
    for (size_t i = 0; i < this->size(); i++) {
        float x = float(i) / float(this->size());
        auto [pos, dir, normal] = this->pointAndDerivativeAndSecondDerivative(x);
        Vector3 binormal = this->getBinormal(x).normalized();
        dir.normalize();
        normal.normalize();
        newPoints[i] = points[i] + (dir * random_gen::generate_perlin(points[i].x() * scale, points[i].y() * scale, points[i].z() * scale + (!loop ? 10 * i : 0)) * maxDistance.x()) + (normal * random_gen::generate_perlin(100 + points[i].x() * scale, points[i].y() * scale, points[i].z() * scale + (!loop ? 10 * i : 0)) * maxDistance.y()) + (binormal * random_gen::generate_perlin(points[i].x() * scale, 100 + points[i].y() * scale, points[i].z() * scale + (!loop ? 10 * i : 0)) * maxDistance.z());
    }
    this->points = newPoints;
    return *this;
}

CatmullRomSpline& CatmullRomSpline::removeDuplicates()
{
    std::vector<Vector3> newPoints;
    for (const auto& point : this->points) {
        if (newPoints.empty() || (point - newPoints.back()).norm() > 0.001f)
            newPoints.push_back(point);
    }
    this->points = newPoints;
    return *this;
}

std::string CatmullRomSpline::toString() const
{
    std::ostringstream out;
    out << "BSpline with " << this->points.size() << " points (" << (closed ? "closed" : "not closed") << ") :\n";
    for (auto& p : this->points)
        out << "- " << p << "\n";
    return out.str();
}

std::vector<Vector3>::const_iterator CatmullRomSpline::begin() const {
    return points.begin();
}

std::vector<Vector3>::const_iterator CatmullRomSpline::end() const {
    return points.end();
}

std::vector<Vector3>::iterator CatmullRomSpline::begin() {
    return points.begin();
}

std::vector<Vector3>::iterator CatmullRomSpline::end() {
    return points.end();
}

std::size_t CatmullRomSpline::size() const {
    return end() - begin();
}

std::size_t CatmullRomSpline::numPoints() const {
    return size();
}

std::size_t CatmullRomSpline::numVertices() const {
    return size();
}

bool CatmullRomSpline::empty() const {
    return begin() == end();
}

Vector3 &CatmullRomSpline::operator[](size_t i)
{
    return this->points[(i + size()) % size()];
}

std::string CatmullRomSpline::display1DPlot(int sizeX, int sizeY) const
{
    std::ostringstream oss;
    CatmullRomSpline translated = *this;
    std::vector<Vector3> roundedPositions(translated.size());
    auto [mini, maxi] = getAABBox().minMax();
    for (auto& p : translated) {
        p = (p - mini) / (maxi - mini); // between (0, 0) and (1, 1)
        p.y() = 1.f - p.y();
        p *= Vector3(sizeX, sizeY);
        roundedPositions.push_back(p.roundedDown());
    }

    for (int y = 0; y < sizeY; y++) {
        if (y == 0) oss << std::fixed << std::setprecision(2) << maxi.y() << "|";
        else if (y == sizeY - 1) oss << std::fixed << std::setprecision(2) << mini.y() << "|";
        else oss << "    |";
        for (int x = 0; x < sizeX; x++) {
            Vector3 pos(x, y);
            if (isIn(pos, roundedPositions)) {
                oss << "X";
            } else {
                float dist = translated.estimateSqrDistanceFrom(pos);
                oss << (dist < 2 * 2 ? "#" : (dist < 4 * 4 ? "+" : "-"));
            }
        }
        oss << std::endl;
    }
    oss << std::fixed << std::setprecision(2) << "   " << mini.x();
    for (int x = 0; x < sizeX - 4; x++) oss << " ";
    oss << std::fixed << std::setprecision(2) << maxi.x();

    return oss.str();
}

Vector3 CatmullRomSpline::computeDerivative(float x) const
{
    std::vector<Vector3> displayedPoints = this->points;
    if (this->closed)
        displayedPoints.push_back(displayedPoints.front());

    size_t lastPointIndex = displayedPoints.size() - 1;
    size_t nbPoints = displayedPoints.size(); // + (this->closed ? 1 : 0);

    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    Vector3 P0 = displayedPoints[(iFloor == 0 ? (this->closed ? int(nbPoints-2) : 1) : iFloor - 1)];
    Vector3 P1 = displayedPoints[iFloor - 0];
    Vector3 P2 = displayedPoints[iCeil + 0];
    Vector3 P3 = displayedPoints[(iCeil >= int(nbPoints-1) ? (this->closed ? 1 : displayedPoints.size()-2) : iCeil + 1)];

    float t0 = 0;
    float t1 = CatmullNextT(P0, P1, t0, alpha);
    float t2 = CatmullNextT(P1, P2, t1, alpha);
    float t3 = CatmullNextT(P2, P3, t2, alpha);

    float t = map(x_prime, 0.f, 1.f, t1, t2);


    Vector3 A1 = P0 * (t1 - t) / (t1 - t0) + P1 * (t - t0) / (t1 - t0);
    Vector3 A2 = P1 * (t2 - t) / (t2 - t1) + P2 * (t - t1) / (t2 - t1);
    Vector3 A3 = P2 * (t3 - t) / (t3 - t2) + P3 * (t - t2) / (t3 - t2);

    Vector3 B1 = A1 * (t2 - t) / (t2 - t0) + A2 * (t - t0) / (t2 - t0);
    Vector3 B2 = A2 * (t3 - t) / (t3 - t1) + A3 * (t - t1) / (t3 - t1);

    Vector3 A1_prim = (P1 - P0) / (t1 - t0);
    Vector3 A2_prim = (P2 - P1) / (t2 - t1);
    Vector3 A3_prim = (P3 - P2) / (t3 - t2);

    Vector3 B1_prim = (A2 - A1) / (t2 - t0) + ((t2 - t) / (t2 - t0)) * A1_prim + ((t - t0) / (t2 - t0)) * A2_prim;
    Vector3 B2_prim = (A3 - A2) / (t3 - t1) + ((t3 - t) / (t3 - t1)) * A2_prim + ((t - t1) / (t3 - t1)) * A3_prim;

    Vector3 C_prim  = (B2 - B1) / (t2 - t1) + ((t2 - t) / (t2 - t1)) * B1_prim + ((t - t1) / (t2 - t1)) * B2_prim;
    return C_prim;
}

std::pair<Vector3, Vector3> CatmullRomSpline::pointAndDerivative(float x) const
{
    std::vector<Vector3> displayedPoints = this->points;
    if (this->closed)
        displayedPoints.push_back(displayedPoints.front());

    size_t lastPointIndex = displayedPoints.size() - 1;
    size_t nbPoints = displayedPoints.size(); // + (this->closed ? 1 : 0);

    if (nbPoints == 0) return {Vector3::invalid, Vector3::invalid};
    else if (nbPoints == 1) return {displayedPoints[0], Vector3::invalid};
    else if (nbPoints == 2) return {displayedPoints[0] * (1.f - x) + displayedPoints[1] * x, displayedPoints[1] - displayedPoints[0]};

    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    Vector3 P0 = displayedPoints[(iFloor == 0 ? (this->closed ? int(nbPoints-2) : 1) : iFloor - 1)];
    Vector3 P1 = displayedPoints[iFloor - 0];
    Vector3 P2 = displayedPoints[iCeil + 0];
    Vector3 P3 = displayedPoints[(iCeil >= int(nbPoints-1) ? (this->closed ? 1 : displayedPoints.size()-2) : iCeil + 1)];

    float t0 = 0;
    float t1 = CatmullNextT(P0, P1, t0, alpha);
    float t2 = CatmullNextT(P1, P2, t1, alpha);
    float t3 = CatmullNextT(P2, P3, t2, alpha);

    float t = map(x_prime, 0.f, 1.f, t1, t2);


    Vector3 A1 = P0 * (t1 - t) / (t1 - t0) + P1 * (t - t0) / (t1 - t0);
    Vector3 A2 = P1 * (t2 - t) / (t2 - t1) + P2 * (t - t1) / (t2 - t1);
    Vector3 A3 = P2 * (t3 - t) / (t3 - t2) + P3 * (t - t2) / (t3 - t2);

    Vector3 B1 = A1 * (t2 - t) / (t2 - t0) + A2 * (t - t0) / (t2 - t0);
    Vector3 B2 = A2 * (t3 - t) / (t3 - t1) + A3 * (t - t1) / (t3 - t1);

    Vector3 C  = B1 * (t2 - t) / (t2 - t1) + B2 * (t - t1) / (t2 - t1);

    Vector3 A1_prim = (P1 - P0) / (t1 - t0);
    Vector3 A2_prim = (P2 - P1) / (t2 - t1);
    Vector3 A3_prim = (P3 - P2) / (t3 - t2);

    Vector3 B1_prim = (A2 - A1) / (t2 - t0) + ((t2 - t) / (t2 - t0)) * A1_prim + ((t - t0) / (t2 - t0)) * A2_prim;
    Vector3 B2_prim = (A3 - A2) / (t3 - t1) + ((t3 - t) / (t3 - t1)) * A2_prim + ((t - t1) / (t3 - t1)) * A3_prim;

    Vector3 C_prim  = (B2 - B1) / (t2 - t1) + ((t2 - t) / (t2 - t1)) * B1_prim + ((t - t1) / (t2 - t1)) * B2_prim;
    return {C, C_prim};
}

std::tuple<Vector3, Vector3, Vector3> CatmullRomSpline::pointAndDerivativeAndSecondDerivative(float _x) const
{
    float x = std::clamp(_x, 0.f, 1.f);

    std::vector<Vector3> displayedPoints = this->points;
    if (closed && points.size() > 0) displayedPoints.push_back(points[0]);
    size_t nbPoints = displayedPoints.size();

    if (nbPoints == 0) return {Vector3::invalid, Vector3::invalid, Vector3::invalid};
    else if (nbPoints == 1) return {displayedPoints[0], Vector3::invalid, Vector3::invalid};
    else if (nbPoints == 2) return {displayedPoints[0] * (1.f - x) + displayedPoints[1] * x, displayedPoints[1] - displayedPoints[0], Vector3(0, 0, 0)};

    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);
    if (_x >= 1.f) {
        iFloor--;
        iCeil--;
        x_prime = 1.f;
    }

    Vector3 P0;
    Vector3 P1;
    Vector3 P2;
    Vector3 P3;

    // No problem for these two points
    P1 = displayedPoints[iFloor - 0];
    P2 = displayedPoints[iCeil + 0];

    const float artificialCurvature = 0.1f;
    if (!isClosed() && iFloor == 0) {
        P3 = displayedPoints[iCeil + 1];
        P0 = (2.f * P1 - P2) - (Collision::projectPointOnSegment(P2, P1, P3) - P2) * artificialCurvature; // Introduce a bit of curvature
    } else if (!isClosed() && iCeil == nbPoints - 1) {
        P0 = displayedPoints[iFloor - 1];
        P3 = (2.f * P2 - P1) + (Collision::projectPointOnSegment(P1, P0, P2) - P1) * artificialCurvature; // Introduce a bit of curvature
    } else {
        P0 = get(iFloor - 1);
        P3 = get(iCeil + 1);
    }

    float t0 = 0;
    float t1 = CatmullNextT(P0, P1, t0, alpha);
    float t2 = CatmullNextT(P1, P2, t1, alpha);
    float t3 = CatmullNextT(P2, P3, t2, alpha);

    float t = map(x_prime, 0.f, 1.f, t1, t2);

    float dtdx = (t2 - t1) * float(nbPoints - 1);


    const Vector3 A1 = P0 * (t1 - t) / (t1 - t0) + P1 * (t - t0) / (t1 - t0);
    const Vector3 A2 = P1 * (t2 - t) / (t2 - t1) + P2 * (t - t1) / (t2 - t1);
    const Vector3 A3 = P2 * (t3 - t) / (t3 - t2) + P3 * (t - t2) / (t3 - t2);

    const Vector3 B1 = A1 * (t2 - t) / (t2 - t0) + A2 * (t - t0) / (t2 - t0);
    const Vector3 B2 = A2 * (t3 - t) / (t3 - t1) + A3 * (t - t1) / (t3 - t1);

    const Vector3 C  = B1 * (t2 - t) / (t2 - t1) + B2 * (t - t1) / (t2 - t1);

    const Vector3 A1p = (P1 - P0) / (t1 - t0);
    const Vector3 A2p = (P2 - P1) / (t2 - t1);
    const Vector3 A3p = (P3 - P2) / (t3 - t2);

    const Vector3 B1p = (A2 - A1) / (t2 - t0) + ((t2 - t) / (t2 - t0)) * A1p + ((t - t0) / (t2 - t0)) * A2p;
    const Vector3 B2p = (A3 - A2) / (t3 - t1) + ((t3 - t) / (t3 - t1)) * A2p + ((t - t1) / (t3 - t1)) * A3p;

    const Vector3 Cp  = (B2 - B1) / (t2 - t1) + ((t2 - t) / (t2 - t1)) * B1p + ((t - t1) / (t2 - t1)) * B2p;

    const Vector3 B1pp = 2.f * (A2p - A1p) / (t2 - t0);
    const Vector3 B2pp = 2.f * (A3p - A2p) / (t3 - t1);

    const Vector3 Cpp = (B1pp * (t2 - t) + B2pp * (t - t1) + 2.f * (B2p - B1p)) / (t2 - t1);

    return {C, Cp * dtdx, Cpp * (dtdx * dtdx)};
}
const Vector3& CatmullRomSpline::operator[](size_t i) const
{
    return this->points[(i + size()) % size()];
}

/*std::ostream& operator<<(std::ostream& io, const BSpline& s) {
    io << s.toString();
    return io;
}

std::ostream& operator<<(std::ostream& io, std::shared_ptr<BSpline> s) {
    io << s->toString();
    return io;
}*/



CatmullRomSpline CatmullRomSpline::random(int numberOfPoints) {
    CatmullRomSpline curve;
    for(int i = 0; i < numberOfPoints; i++) {
        curve.addPoint(Vector3::random());
    }
    return curve;
}

void CatmullRomSpline::setAlpha(float newAlpha) {
    this->alpha = newAlpha;
}

float CatmullRomSpline::getAlpha() const {
    return alpha;
}

void CatmullRomSpline::addPoint(const Vector3 &newPoint) {
    this->points.push_back(newPoint);
}

CatmullRomSpline &CatmullRomSpline::insertPoint(int i, const Vector3 &newPos) {
    this->points.insert(points.begin() + i, newPos); return *this;
}

CatmullRomSpline &CatmullRomSpline::removePoint(int i) {
    this->points.erase(points.begin() + i); return *this;
}

CatmullRomSpline &CatmullRomSpline::reset() {
    this->points.clear(); return *this;
}

std::vector<std::shared_ptr<Curve> > CatmullRomSpline::slice(float t) const
{
    size_t lowerIndex = this->timeToLowerIndex(t);
    std::vector<Vector3> firstPoints(lowerIndex + 2);
    std::vector<Vector3> lastPoints(numPoints() - lowerIndex);

    for (size_t i = 0; i < numPoints(); i++) {
        if (i <= lowerIndex) {
            firstPoints[i] = points[i];
        } else {
            lastPoints[i - lowerIndex] = points[pointIndex(i)];
        }
    }
    const auto p = getPoint(t);
    firstPoints.back() = p;
    lastPoints.front() = p;

    if (closed) lastPoints.push_back(points[0]);

    CatmullRomSpline p1(firstPoints);
    CatmullRomSpline p2(lastPoints);
    return {std::make_shared<CatmullRomSpline>(p1), std::make_shared<CatmullRomSpline>(p2)};
}


