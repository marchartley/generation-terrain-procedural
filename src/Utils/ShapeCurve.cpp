#include "ShapeCurve.h"
#include "Utils/Collisions.h"

ShapeCurve::ShapeCurve()
    : ShapeCurve(std::vector<Vector3>())
{

}

ShapeCurve::ShapeCurve(std::vector<Vector3> points)
    : ShapeCurve(BSpline(points))
{

}

ShapeCurve::ShapeCurve(BSpline path)
    : BSpline(path)
{
    if (!this->closed) {
        if (!this->points.empty() && this->points.front() != this->points.back())
            this->close();
    }
}

bool ShapeCurve::inside(Vector3 pos, bool useNativeShape)
{
    std::vector<Vector3> pointsUsed;
    if (useNativeShape) {
        pointsUsed = this->points;
        if (!this->points.empty())
            pointsUsed.insert(pointsUsed.end(), this->points.front());
    } else {
        pointsUsed = this->getPath(10);
    }
    return Collision::pointInPolygon(pos, pointsUsed);
    /*
    if (this->points.size() < 2) return false;

    size_t firstIndex = 0, secondIndex = 0;

    Vector3 firstVertex = this->points[firstIndex];
    Vector3 secondVertex = this->points[secondIndex];
    Vector3 center;
    for (const auto& p : this->points)
        center += p;
    center /= (float)this->points.size();
    // First, lets find a good plane to define our shape
    while (std::abs((firstVertex - center).normalized().dot((secondVertex - center).normalized())) < (1 - 1e-5)) {
        secondIndex ++;
        if (secondIndex >= this->points.size()) {
            secondIndex = 0;
            firstIndex++;
            if (firstIndex >= this->points.size()) return false;
            firstVertex = this->points[firstIndex];
        }
        secondVertex = this->points[secondIndex];
    }
    // At this point, we can check if the "pos" is in the same plane as the shape

    std::vector<Vector3> path = this->getPath(10); // Should be enough for now
    Vector3 ray = center + (firstVertex - center).normalized() * this->length(); // Same, should be outside the shape

    // Check the intersection of the ray with all the segments of the shape
    int nb_intersections = 0;
    for (size_t i = 0; i < path.size() - 1; i++) {
        if (Collision::intersectionBetweenTwoSegments(pos, ray, path[i], path[i+1]).isValid())
            nb_intersections++;
    }
    // If there's a odd number of intersections, the point is inside
    return (nb_intersections % 2) == 1;*/
}

float ShapeCurve::estimateDistanceFrom(Vector3 pos)
{
    float dist = BSpline::estimateDistanceFrom(pos);
    return dist * (inside(pos) ? -1.f : 1.f); // Negative distance if it's currently inside
}

float ShapeCurve::computeArea()
{
    float area = 0;
    for (size_t i = 1; i < this->points.size() + 1; i++){
        area += points[i % points.size()].x * (points[(i+1) % points.size()].y - points[(i-1) % points.size()].y);
    }
    return std::abs(area) / 2.f;
}

Vector3 ShapeCurve::planeNormal()
{
    Vector3 normal;
    int numberOfSamples = 10;
    for (int i = 0; i < numberOfSamples; i++) {
        float t = i / (float)(numberOfSamples);
        Vector3 binormal = this->getBinormal(t);
        if (binormal.isValid())
            normal += binormal;
    }
    return normal.normalize();
}

std::vector<Vector3> ShapeCurve::randomPointsInside(int numberOfPoints)
{
    if (this->points.empty()) return std::vector<Vector3>();
    int maxFailures = 10000 * numberOfPoints;
    std::vector<Vector3> returnedPoints;
    Vector3 minVec, maxVec;
    std::tie(minVec, maxVec) = this->AABBox();
    minVec.z = -1;
    maxVec.z =  1;
    Vector3 normalRay = this->planeNormal() * (maxVec - minVec).norm();

    for (int i = 0; i < numberOfPoints; i++) {
        // Check the collision from a point below and a point above the plane
        Vector3 randomPoint = Vector3::random(minVec, maxVec) - normalRay;
        Vector3 intersectionPoint = Collision::intersectionRayPlane(randomPoint, normalRay * 2.f, this->points[0], normalRay);
//        if (!intersectionPoint.isValid()) {
//            intersectionPoint = Collision::intersectionRayPlane(randomPoint, normalRay * -1.f, this->points[0], normalRay);
//        }
        if (intersectionPoint.isValid()) {
            if (Collision::pointInPolygon(intersectionPoint, points)) { //getPath(points.size()))) {
                returnedPoints.push_back(intersectionPoint);
                continue;
            }
        }
        i --;
        maxFailures --;
        if (maxFailures < 0)
            break;
    }
    return returnedPoints;
}


ShapeCurve ShapeCurve::grow(float increase)
{
    Vector3 normal = this->planeNormal();
    std::vector<Vector3> newPoints = this->points;
    for (size_t i = 0; i < this->points.size(); i++) {
//        Vector3 point = this->points[i];
        Vector3 next_point = this->points[i + 1];
        Vector3 prev_point = this->points[(this->points.size() + i - 1) % this->points.size()];
        Vector3 dir = (next_point - prev_point).normalize();
        /*
        float time = estimateClosestTime(this->points[i]);
        Vector3 normal = getNormal(estimateClosestTime(this->points[i]));*/
        newPoints[i] = this->points[i] + dir.cross(normal) * increase; //+ getNormal(estimateClosestTime(this->points[i])) * increase;
    }
    this->points = newPoints;
    return *this;
}

ShapeCurve ShapeCurve::shrink(float decrease)
{
    return this->grow(-decrease);
}
