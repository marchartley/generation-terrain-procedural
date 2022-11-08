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

ShapeCurve &ShapeCurve::translate(Vector3 translation)
{
    for (auto& point : points)
        point += translation;
    return *this;
}

bool ShapeCurve::inside(Vector3 pos, bool useNativeShape)
{
    std::vector<Vector3> pointsUsed;
    if (useNativeShape) {
        pointsUsed = this->points;
//        if (!this->points.empty())
//            pointsUsed.insert(pointsUsed.end(), this->points.front());
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

struct ClipVertex {
    Vector3 coord;
    bool isEntry;
    bool isExit;
    bool isInside;
    ClipVertex* neighbor = nullptr;
    ClipVertex* prev = nullptr;
    ClipVertex* next = nullptr;
    int shapeID = -1;
    int index = -1;

    ClipVertex(Vector3 p) : ClipVertex(p, false, false, false) {}
    ClipVertex(Vector3 p, bool entry, bool exit, bool inside)
        : coord(p), isEntry(entry), isExit(exit), isInside(inside)
    {}
};
int getIndex(int index, size_t size) {
    return (index + size) % size;
}
int markEntriesExits(std::vector<ClipVertex>& poly, bool currentlyInside, int shapeID) {
    int firstIntersectionIndex = -1;
    for (size_t i = 0; i < poly.size(); i++) {
        poly[i].shapeID = shapeID;
        poly[i].index = i;
        poly[i].prev = &poly[getIndex(i - 1, poly.size())];
        poly[i].next = &poly[getIndex(i + 1, poly.size())];

        int id0 = i;
        int id1 = getIndex(i + 1, poly.size());
        auto& p0 = poly[id0];
        auto& p1 = poly[id1];

        if (p0.neighbor && firstIntersectionIndex < 0) {
            firstIntersectionIndex = p0.index;
        }

        if (p0.neighbor && p1.neighbor) {
            // Both are on the intersection, one is entry, the other is exit
            p0.isEntry = !currentlyInside;
            p0.isExit = currentlyInside;
            currentlyInside = !currentlyInside;
        } else if (p0.neighbor && !p1.neighbor) {
            // One intersection, change "isInside"
            p0.isEntry = !currentlyInside;
            p0.isExit = currentlyInside;
            currentlyInside = !currentlyInside;
        } else if (!p0.neighbor && p1.neighbor) {
            // One intersection, change "isInside"
            p0.isInside = currentlyInside;
//            currentlyInside = !currentlyInside;
        } else {
            // Completly inside or outside
            p0.isInside = currentlyInside;
        }
    }
    return firstIntersectionIndex;
}

ShapeCurve ShapeCurve::intersect(ShapeCurve other)
{
    std::vector<std::vector<Vector3>> result;

    ShapeCurve polyShape = *this;
    polyShape = polyShape.removeDuplicates();
    ShapeCurve clipShape = other;
    other = other.removeDuplicates();

    std::vector<ClipVertex> poly, clip;
    for (size_t i = 0; i < polyShape.points.size(); i++) poly.push_back(polyShape.points[i]);
    for (size_t i = 0; i < clipShape.points.size(); i++) clip.push_back(clipShape.points[i]);

    bool foundIntersection = false;

    for (size_t i = 0; i < poly.size(); i++) {
        auto& A = poly[i];
        auto& B = poly[getIndex(i + 1, poly.size())];

        for (size_t j = 0; j < clip.size(); j++) {
            auto& C = clip[j];
            auto& D = clip[getIndex(j + 1, clip.size())];

            Vector3 intersection = Collision::intersectionBetweenTwoSegments(A.coord, B.coord, C.coord, D.coord);
            if (intersection.isValid()) {
                foundIntersection = true;

                poly.insert(poly.begin() + i + 1, intersection);
                clip.insert(clip.begin() + j + 1, intersection);

                poly[i + 1].neighbor = &clip[j + 1];
                clip[j + 1].neighbor = &poly[i + 1];

                j++;
                i++;
            }
        }
    }

    /// TODO : There is a special case where P0 is on an intersection...
    bool currentlyInside = clipShape.inside(poly[0].coord);

    if (!foundIntersection) {
        // Shape is completely inside or outside
        if (currentlyInside) {
            // Poly is inside
            return polyShape;
        } else {
            if (polyShape.inside(clip[0].coord)) {
                // Clipping shape is inside
                return clipShape;
            } else {
                // No intersection
                return ShapeCurve();
            }
        }
    } else {
        int firstIntersectionIndex = markEntriesExits(poly, clipShape.inside(poly[0].coord), 0);
        markEntriesExits(clip, polyShape.inside(clip[0].coord), 1);

        std::vector<Vector3> resultingShape;
        auto& firstVertex = poly[firstIntersectionIndex];
        auto current = firstVertex;
//        resultingShape.push_back(firstVertex.coord);
        while (true) {
            resultingShape.push_back(current.coord);
            if ((current.shapeID == 0 && current.isEntry) || (current.shapeID == 1 && current.isExit)) {
                current = *current.next;
            } else if ((current.shapeID == 0 && current.isExit) || (current.shapeID == 1 && current.isEntry)) {
                current = *current.neighbor->next;
            } else {
                current = *current.next;
            }
            // If we're back to the start, end the loop
            if ((current.shapeID == 0 && current.index == firstIntersectionIndex) ||
                    (current.shapeID == 1 && current.neighbor && current.neighbor->index == firstIntersectionIndex))
                break;
        }

        return ShapeCurve(resultingShape);
    }
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

ShapeCurve &ShapeCurve::scale(float factor)
{
    for (auto& vert : this->points)
        vert *= factor;
    return *this;
}


ShapeCurve ShapeCurve::grow(float increase)
{
    ShapeCurve copy = *this;
    std::vector<Vector3> newPoints = copy.removeDuplicates().points;
    Vector3 normal = copy.planeNormal();
    for (size_t i = 0; i < copy.points.size(); i++) {
//        Vector3 point = this->points[i];
        Vector3 next_point = copy.points[(i + 1) % copy.points.size()];
        Vector3 prev_point = copy.points[(copy.points.size() + i - 1) % copy.points.size()];
        Vector3 dir = (next_point - prev_point).normalize();
        /*
        float time = estimateClosestTime(this->points[i]);
        Vector3 normal = getNormal(estimateClosestTime(this->points[i]));*/
        newPoints[i] = copy.points[i] + dir.cross(normal) * increase; //+ getNormal(estimateClosestTime(this->points[i])) * increase;
    }
    copy.points = newPoints;
    return copy;
}

ShapeCurve ShapeCurve::shrink(float decrease)
{
    return this->grow(-decrease);
}

ShapeCurve &ShapeCurve::removeDuplicates()
{
    BSpline::removeDuplicates();
    std::vector<Vector3> foundPattern;
    for (size_t i = 0; i < points.size(); i++) {
        if (i < 3)
            foundPattern.push_back(points[i]);
        else {
            if (points[i] == foundPattern[0]) {
                size_t ii;
                for (ii = i; ii < points.size(); ii++) {
                    int patternIndex = (ii - i) % foundPattern.size();
                    if (points[ii] != foundPattern[patternIndex]) {
                        break;
                    }
                }
                if (ii - i > 1) {
                    // This will ignore the last checked vertex. If you want to continue exactly, substract 1 to the distance.
                    i = std ::distance(points.begin(), points.erase(points.begin() + i, points.begin() + ii));
                }
            }
        }
    }
    return *this;
}
