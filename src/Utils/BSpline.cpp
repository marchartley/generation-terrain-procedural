#include "Utils/BSpline.h"

BSpline::BSpline()
{

}
BSpline::BSpline(int numberOfPoints) {
    for(int i = 0; i < numberOfPoints; i++) {
        this->points.push_back(Vector3::random());
    }
}

BSpline::BSpline(std::vector<Vector3> points)
    : points(points)
{
}
BSpline::BSpline(std::vector<BSpline> subsplines)
{
    for (BSpline& spline : subsplines) {
        bool ignoreFirstPoint = (this->points.empty() ? false : this->points.back() == spline.points.front());
        this->points.insert(this->points.end(), spline.points.begin() + (ignoreFirstPoint ? 1 : 0), spline.points.end());
    }
}

std::vector<Vector3> BSpline::getPath(float resolution)
{
    std::vector<Vector3> path;
    for (float x = 0.0; x < 1.0; x += resolution)
        path.push_back(this->getPoint(x));
    path.push_back(this->getPoint(1.0));
    return path;
}

Vector3 BSpline::getPoint(float x)
{
    if (points.size() >= 4)
        return this->getCatmullPoint(x);

    if(points.size() == 0)
        return Vector3();

    std::vector<Vector3> controls = points;
    while (controls.size() > 1)
    {
        std::vector<Vector3> newCtrls;
        for (size_t i = 0; i < controls.size() - 1; i++) {
            newCtrls.push_back(this->getPoint(x, controls[i], controls[i+1]));
        }
        controls = newCtrls;
    }
    return controls[0];
}
Vector3 BSpline::getPoint(float x, Vector3 a, Vector3 b)
{
    return a * (1 - x) + b * x;
}

Vector3 BSpline::getDerivative(float x)
{
    float previousTime = std::max(0.f, x - 0.001f);
    float nextTime = std::min(1.f, x + 0.001f);
    return (getPoint(nextTime) - getPoint(previousTime)).normalized();
}

float BSpline::estimateClosestTime(Vector3 pos)
{
    int checks = 0;
    float distDiff = std::numeric_limits<float>::max();
//    Vector3 bestPosPrediction = pos;
    size_t closestPoint;
    for (size_t i = 0; i < this->points.size(); i++) {
        float dist = (pos - this->points[i]).norm2();
        if (dist < distDiff) {
            distDiff = dist;
            closestPoint = i;
        }
    }
    // We know our point is between t(closest - 1) and t(closest + 1)
    float halfSpace = 1 / (this->points.size() - 1.f);
    float currentTime = closestPoint * halfSpace;
    while (halfSpace > 0.00001) {
        halfSpace *= .5f;
        Vector3 low = getPoint(currentTime - halfSpace), cur = getPoint(currentTime), hig = getPoint(currentTime + halfSpace);
        float lowDist = (low - pos).norm2(), curDist = (cur - pos).norm2(), higDist = (hig - pos).norm2();
        if (lowDist < curDist && lowDist < higDist) { // Lowest is closest
            currentTime -= halfSpace * .5f;
        }
        else if (curDist < lowDist && curDist < higDist) { // Current is closest
            continue;
        }
        else if (higDist < lowDist && higDist < curDist) { // Highest is closest
            currentTime += halfSpace * .5f;
        }
        checks++;
    }
    return std::clamp(currentTime, 0.f, 1.f);
}

Vector3 BSpline::estimateClosestPos(Vector3 pos)
{
    return this->getPoint(this->estimateClosestTime(pos));
}

float BSpline::length()
{
    float length = 0;
    for (size_t i = 0; i < this->points.size() - 1; i++) {
        length += (this->points[i] - this->points[i + 1]).norm();
    }
    return length;
}

BSpline& BSpline::close()
{
    if (this->points.size() > 0 && this->points.front() != this->points.back())
        this->points.push_back(this->points.front());
    return *this;
}

float CatmullNextT(Vector3 P0, Vector3 P1, float t_prev, float alpha)
{
    return std::pow((P0 - P1).norm2(), alpha) + t_prev;
}
template <class T>
T map(T x, T prev_min, T prev_max, T new_min, T new_max)
{
    return ((x - prev_min) / (prev_max - prev_min)) * (new_max - new_min) + new_min;
}

Vector3 BSpline::getCatmullPoint(float x)
{
    float alpha = 0.5f;

    //alpha /= 2.f;

    float res = 1 / (float)(this->points.size() - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);

    if (iFloor <= 0) {
        return this->getPoint(x_prime, this->points[0], this->points[1]);
    } else if (iCeil >= int(this->points.size()) - 1) {
        return this->getPoint(x_prime, this->points[iFloor], this->points[iCeil]);
    }

    Vector3 P0 = points[iFloor - 1];
    Vector3 P1 = points[iFloor - 0];
    Vector3 P2 = points[iCeil + 0];
    Vector3 P3 = points[iCeil + 1];

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
}

BSpline BSpline::simplifyByRamerDouglasPeucker(float epsilon, BSpline subspline)
{
    if (subspline.points.empty()) {
        if (this->points.empty()) return *this; // We are just trying to do a simplification from an empty curve, that's pointless
        else subspline = *this;
    }
    if (subspline.points.size() == 1) return subspline;

    // Find farest point from the line going from start to end of this (sub)spline
    Vector3 vecAB = (subspline.points.front() - subspline.points.back()).normalized();
    float maxDist = 0;
    int index = -1;
    for (size_t i = 1; i < subspline.points.size() - 1; i++) {
        float d = vecAB.cross((subspline.points[i] - subspline.points.front())).norm2();
        if (d > maxDist) {
            maxDist = d;
            index = i;
        }
    }
    // Now we split the spline in two subsplines, and apply recursively the algorithm until all points are "close" enough ( dist < epsilon)
    BSpline returningSpline;
    if (maxDist > epsilon * epsilon) {
        BSpline sub1 = this->simplifyByRamerDouglasPeucker(epsilon, std::vector<Vector3>(subspline.points.begin(), subspline.points.begin() + index));
        BSpline sub2 = this->simplifyByRamerDouglasPeucker(epsilon, std::vector<Vector3>(subspline.points.begin() + index, subspline.points.end()));
        returningSpline = BSpline({sub1, sub2});
    } else {
        returningSpline.points = {subspline.points.front(), subspline.points.back()};
    }
    return returningSpline;
}
