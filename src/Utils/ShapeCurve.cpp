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

bool ShapeCurve::inside(Vector3 pos)
{
    return Collision::pointInPolygon(pos, this->getPath(10));
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
