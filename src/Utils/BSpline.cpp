#include "Utils/BSpline.h"

#include "Utils/Utils.h"
#include "Utils/Collisions.h"

#include <exception>

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

std::vector<Vector3> BSpline::getPath(int numberOfPoints)
{
    /// I'm really not sure this is the best solution, but an easy fix :
    /// forcing user to have at least 2 points
    numberOfPoints = std::max(numberOfPoints, 2);
//    if (numberOfPoints < 2)
//        throw std::invalid_argument("At least 2 points needed to get the path");
    std::vector<Vector3> path;
    float resolution = 1.f / (float)(numberOfPoints - 1);
    for (int i = 0; i < numberOfPoints; i ++)
        path.push_back(this->getPoint(i * resolution));
//    path.push_back(this->getPoint(1.0));
    return path;
}

Vector3 BSpline::getPoint(float x)
{
    if (points.size() >= 2)
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

Vector3 BSpline::getDerivative(float x, bool verbose)
{
    float previousTime = std::max(0.f, x - 0.001f);
    float nextTime = std::min(1.f, x + 0.001f);
    if (verbose) {
        std::cout << "Getting the derivative of a path with " << this->points.size() << " points : \n";
        for (const auto& p : this->points)
            std::cout << " - " << p << "\n";
        std::cout << "Position used : " << previousTime << " and " << nextTime << " -> " << getPoint(nextTime) << " - " << getPoint(previousTime) << " = " << getPoint(nextTime) - getPoint(previousTime) << "\n";
    }

    return (getPoint(nextTime) - getPoint(previousTime)).normalized();
}

Vector3 BSpline::getSecondDerivative(float x)
{
    float previousTime = std::max(0.f, x - 0.001f);
    float nextTime = std::min(1.f, x + 0.001f);
    return (getDerivative(nextTime) - getDerivative(previousTime)).normalized();
}

float BSpline::estimateClosestTime(Vector3 pos)
{
    if (this->points.size() == 0) {
        return 0;
    } else if (this->points.size() == 1) {
        return 0;
    } else if (this->points.size() == 2) {
        Vector3 line = (this->points[1] - this->points[0]);
        return clamp((pos - this->points[0]).dot(line) / line.dot(line), 0.f, 1.f);
    }
    float closestTime = 0;
    float minDistance = std::numeric_limits<float>::max();
    int numberOfChecks = int(this->points.size());
    float precisionFactor = 1.f / numberOfChecks;
    float precision = precisionFactor;
    float searchRangeMin = 0.f;

    while (precision > 1e-6) {
        for (int i = 0; i < numberOfChecks + 1; i++) {
            float time = searchRangeMin + (i * precision);
            float distance = (getPoint(time) - pos).norm2();
            if (distance < minDistance) {
                minDistance = distance;
                closestTime = time;
            }
        }
        searchRangeMin = clamp(closestTime - precision/2.f, 0.f, 1.f);
        precision *= precisionFactor;
    }
    return closestTime;
    /*
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
    while (halfSpace > 1e-6) {
        Vector3 low = getPoint(clamp(currentTime - halfSpace, 0.f, 1.f)),
                cur = getPoint(currentTime),
                hig = getPoint(clamp(currentTime + halfSpace, 0.f, 1.f));
        float lowDist = (low - pos).norm2(),
                curDist = (cur - pos).norm2(),
                higDist = (hig - pos).norm2();
        std::cout << "-Checking between " << low << " and " << hig << " (current time = " << currentTime << " at pos " << cur << "\n";

        currentTime = clamp(currentTime, 0.f, 1.f);
        halfSpace *= .5f;

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
    return std::clamp(currentTime, 0.f, 1.f);*/
}

Vector3 BSpline::estimateClosestPos(Vector3 pos)
{
    return this->getPoint(this->estimateClosestTime(pos));
}
float BSpline::estimateDistanceFrom(Vector3 pos)
{
    return (this->estimateClosestPos(pos) - pos).norm();
}

float BSpline::length()
{
    float length = 0;
    for (size_t i = 0; i < this->points.size() - 1; i++) {
        length += (this->points[i] - this->points[i + 1]).norm();
    }
    return length;
}

Vector3 BSpline::getFrenetDirection(float x)
{
    return getDirection(x);
}

Vector3 BSpline::getFrenetNormal(float x)
{
    return this->getFrenetDirection(x).cross(this->getFrenetBinormal(x)).normalize();
    /*
    Vector3 new_dir = this->getFrenetDirection(x);
    Vector3 forward(0, 1, 0);
    Vector3 up(0, 0, 1);
    Vector3 right(1, 0, 0);
    if (new_dir.isAlmostVertical())
        right = Vector3(0, 1, 1).cross(new_dir);
    else
        right = Vector3(0, 0, 1).cross(new_dir);
    return new_dir.cross(right).normalize();*/
}

Vector3 BSpline::getFrenetBinormal(float x)
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

Vector3 BSpline::getCenterCircle(float x)
{
    return this->getPoint(x) + this->getNormal(x) * (1 / this->getCurvature(x));
}

Vector3 BSpline::getDirection(float x)
{
    return this->getDerivative(x).normalize();
}

Vector3 BSpline::getNormal(float x)
{
    return this->getSecondDerivative(x).normalize();
}

Vector3 BSpline::getBinormal(float x)
{
    return this->getDirection(x).cross(this->getNormal(x)).normalize();
}

float BSpline::getCurvature(float x)
{
    return (getDerivative(x).cross(getSecondDerivative(x))).norm() / (std::pow(getDerivative(x).norm(), 3));
}

BSpline& BSpline::close()
{
    if (this->points.size() > 1 && !this->closed) { // && this->points.front() != this->points.back()) {
        this->closed = true;
    }
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
    std::vector<Vector3> displayedPoints = this->points;
    if (this->closed)
        displayedPoints.push_back(displayedPoints.front());

    size_t lastPointIndex = (this->closed ? displayedPoints.size() - 1 : displayedPoints.size() - 1);
    size_t nbPoints = displayedPoints.size(); // + (this->closed ? 1 : 0);

    if (x == 0.f) return displayedPoints[0];
    if (x == 1.f) return displayedPoints[lastPointIndex];
    float alpha = 2.f; // 2 = very round, 1 = quite normal, 0.5 = almost linear

    //alpha /= 2.f;

    float res = 1 / (float)(nbPoints - 1);
    int iFloor = int(x / res);
    int iCeil = int(x / res) + 1;
    float resFloor = iFloor * res;
    float resCeil = iCeil * res;
    float x_prime = map(x, resFloor, resCeil, 0.f, 1.f);
/*
    if (iFloor <= 0) {
        return this->getPoint(x_prime, displayedPoints[0], displayedPoints[1]);
    } else if (iCeil >= int(displayedPoints.size()) - 1) {
        return this->getPoint(x_prime, displayedPoints[iFloor], displayedPoints[iCeil]);
    }*/

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

std::tuple<Vector3, Vector3> BSpline::AABBox()
{
    if (this->points.empty()) return std::make_tuple(Vector3(), Vector3());
    if (this->points.size() == 1) return std::make_tuple(points[0], points[0]);

    float minDim = std::numeric_limits<float>::min();
    float maxDim = std::numeric_limits<float>::max();
    Vector3 minVec = Vector3(maxDim, maxDim, maxDim),
            maxVec = Vector3(minDim, minDim, minDim);
    for (const auto& point : points) {
        minVec.x = std::min(point.x, minVec.x);
        minVec.y = std::min(point.y, minVec.y);
        minVec.z = std::min(point.z, minVec.z);
        maxVec.x = std::max(point.x, maxVec.x);
        maxVec.y = std::max(point.y, maxVec.y);
        maxVec.z = std::max(point.z, maxVec.z);
    }
    return std::make_tuple(minVec, maxVec);
}
