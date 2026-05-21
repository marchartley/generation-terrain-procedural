#include "Contour.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"

Contour::Contour()
    : Contour(std::vector<Vector3>())
{

}

Contour::Contour(std::vector<Vector3> points)
    : Contour(Polyline(points))
{

}

Contour::Contour(const Curve& path)
{
    this->curve = std::shared_ptr<Curve>(path.clone());
    curve->close();
}
Contour::Contour(Curve* path)
{
    this->curve = std::shared_ptr<Curve>(path);
    curve->close();
}
Contour::Contour(std::shared_ptr<Curve> path)
{
    this->curve = path;
    curve->close();
}

Contour::~Contour()
{
    // std::cout << "~Contour" << std::endl;
}

Contour& Contour::translate(const Vector3& translation)
{
    curve->translate(translation);
    return *this;
}

bool Contour::contains(const Vector3& pos, bool useNativeShape) const
{
    auto points = this->curve->getPath(200);
    return Collision::pointInPolygon(pos, points);
    /*
    std::vector<Vector3> pointsUsed;
    if (useNativeShape) {
        pointsUsed = points;
//        if (!points.empty())
//            pointsUsed.insert(pointsUsed.end(), points.front());
    } else {
        pointsUsed = this->getPath(200);
    }
    return Collision::pointInPolygon(pos, pointsUsed);
    */
}

bool Contour::containsXY(const Vector3 &pos, bool useNativeShape, int increaseAccuracy) const
{
    return this->contains(pos);
    /*
    if (this->size() == 0) return false;
    std::vector<Vector3> pointsUsed;
    if (useNativeShape) {
        pointsUsed = points;
    } else {
        pointsUsed = this->getPath(20);
    }
    for (auto& p :  pointsUsed)
        p = p.xy();
    if (increaseAccuracy <= 0)
        return Collision::pointInPolygon(pos.xy(), pointsUsed);

    int nbTests = increaseAccuracy * 2 + 1;
    int positiveCounts = 0;
    for (int i = 0; i < nbTests; i++)
        positiveCounts += (Collision::pointInPolygon(pos.xy() + Vector3::random().xy() * .1f, pointsUsed) ? 1 : 0);
    return positiveCounts > increaseAccuracy;
    */
}

float Contour::estimateDistanceFrom(const Vector3& pos) const
{
    float dist = curve->estimateDistanceFrom(pos);
    // float dist = /*BSpline(this->closedPath()).*/CatmullRomSpline::estimateDistanceFrom(pos);
    return dist * (contains(pos, false) ? -1.f : 1.f); // Negative distance if it's currently inside
}

float Contour::estimateSignedDistanceFrom(const Vector3& pos) const
{
    return this->estimateDistanceFrom(pos);
    // return /*BSpline(this->closedPath()).*/CatmullRomSpline::estimateSignedDistanceFrom(pos);
}

float Contour::computeArea() const
{
    return std::abs(this->computeSignedArea());
}

float Contour::computeSignedArea() const
{
    float area = 0;
    auto points = this->curve->getPath(this->curve->numPoints() * 3);
    for (size_t i = 1; i < points.size() + 1; i++){
        area += points[i % points.size()].x() * (points[(i+1) % points.size()].y() - points[(i-1) % points.size()].y());
    }
    return area / 2.f;
}

Vector3 Contour::centroid() const
{
    /*Vector3 centroid;
    float totalArea = 0;
    const auto points = curve->getPath();
    for (int i = 0; i < points.size(); i++) {
        int j = (i + 1) % (points.size());
        int k = (i - 1 + points.size()) % (points.size());

        float triangleArea = 0.5f * (points[i] - points[j]).cross(points[i] - points[k]).norm();
        centroid += triangleArea * points[i];
        totalArea += triangleArea;
    }
    return centroid / totalArea;*/
    Vector3 centroid;
    const auto points = curve->getPath();
    for (size_t i = 0; i < points.size(); i++) {
        const auto& p0 = curve->get(i);
        const auto& p1 = curve->get(i+1);
        centroid += Vector3(p0.x() + p1.x(), p0.y() + p1.y()) * (p0.x() * p1.y() - p1.x() * p0.y());
    }
    return (centroid / (6.f * computeArea())) + Vector3(0, 0, points[0].z());
}

int getIndex(int index, size_t size) {
    return (index + size) % size;
}
int markEntriesExits(std::vector<ClipVertex*>& poly, bool currentlyInside, int shapeID) {
    int firstIntersectionIndex = -1;
    for (size_t i = 0; i < poly.size(); i++) {
        poly[i]->shapeID = shapeID;
        poly[i]->index = i;
        poly[i]->prev = poly[getIndex(i - 1, poly.size())];
        poly[i]->next = poly[getIndex(i + 1, poly.size())];

        int id0 = i;
        int id1 = getIndex(i + 1, poly.size());
        auto& p0 = poly[id0];
        auto& p1 = poly[id1];

        if (p0->neighbor && firstIntersectionIndex < 0) {
            firstIntersectionIndex = p0->index;
        }

        if (p0->neighbor && p1->neighbor) {
            // Both are on the intersection, one is entry, the other is exit
            p0->isEntry = !currentlyInside;
            p0->isExit = currentlyInside;
            currentlyInside = !currentlyInside;
        } else if (p0->neighbor && !p1->neighbor) {
            // One intersection, change "isInside"
            p0->isEntry = !currentlyInside;
            p0->isExit = currentlyInside;
            currentlyInside = !currentlyInside;
        } else if (!p0->neighbor && p1->neighbor) {
            // One intersection, change "isInside"
            p0->isInside = currentlyInside;
//            currentlyInside = !currentlyInside;
        } else {
            // Completly inside or outside
            p0->isInside = currentlyInside;
        }
    }
    return firstIntersectionIndex;
}
/*
Contour Contour::intersect(Contour other)
{
    std::vector<std::vector<Vector3>> result;

    Contour polyShape = *this;
    polyShape = polyShape.removeDuplicates();
    Contour clipShape = other;
    other = other.removeDuplicates();

    std::vector<ClipVertex*> poly, clip;
    for (size_t i = 0; i < polyShape.points.size(); i++) poly.push_back(new ClipVertex(polyShape.points[i]));
    for (size_t i = 0; i < clipShape.points.size(); i++) clip.push_back(new ClipVertex(clipShape.points[i]));

    bool foundIntersection = false;

    for (size_t i = 0; i < poly.size(); i++) {
        auto& A = poly[i];
        auto& B = poly[getIndex(i + 1, poly.size())];

        for (size_t j = 0; j < clip.size(); j++) {
            auto& C = clip[j];
            auto& D = clip[getIndex(j + 1, clip.size())];

            Vector3 intersection = Collision::intersectionBetweenTwoSegments(A->coord, B->coord, C->coord, D->coord);
            if (intersection.isValid()) {
                foundIntersection = true;

                poly.insert(poly.begin() + i + 1, new ClipVertex(intersection));
                clip.insert(clip.begin() + j + 1, new ClipVertex(intersection));

                poly[i + 1]->neighbor = clip[j + 1];
                clip[j + 1]->neighbor = poly[i + 1];

                j++;
                i++;
            }
        }
    }

    /// TODO : There is a special case where P0 is on an intersection...
    bool currentlyInside = clipShape.contains(poly[0]->coord);

    if (!foundIntersection) {
        // Shape is completely inside or outside
        if (currentlyInside) {
            // Poly is inside
            return polyShape;
        } else {
            if (polyShape.contains(clip[0]->coord)) {
                // Clipping shape is inside
                return clipShape;
            } else {
                // No intersection
                return Contour();
            }
        }
    } else {
        int firstIntersectionIndex = markEntriesExits(poly, clipShape.contains(poly[0]->coord), 0);
        markEntriesExits(clip, polyShape.contains(clip[0]->coord), 1);

        std::vector<Vector3> resultingShape;
        auto& firstVertex = poly[firstIntersectionIndex];
        auto current = firstVertex;
//        resultingShape.push_back(firstVertex.coord);
        while (true) {
            resultingShape.push_back(current->coord);
            if ((current->shapeID == 0 && current->isEntry) || (current->shapeID == 1 && current->isExit)) {
                current = current->next;
            } else if ((current->shapeID == 0 && current->isExit) || (current->shapeID == 1 && current->isEntry)) {
                current = current->neighbor->next;
            } else {
                current = current->next;
            }
            // If we're back to the start, end the loop
            if ((current->shapeID == 0 && current->index == firstIntersectionIndex) ||
                    (current->shapeID == 1 && current->neighbor && current->neighbor->index == firstIntersectionIndex))
                break;
        }

        return Contour(resultingShape);
    }
}
*/
Vector3 Contour::planeNormal() const
{
    Vector3 normal;
    /*int numberOfSamples = 10;
    for (int i = 0; i < numberOfSamples; i++) {
        float t = i / (float)(numberOfSamples);
        Vector3 binormal = this->curve->getBinormal(t);
        if (binormal.isValid())
            normal += binormal;
    }*/
    const auto& points = curve->getPath();
    for (int i = 0; i < points.size() - 2; i++) {
        Vector3 n = (points[i + 1] - points[i]).cross(points[i + 2] - points[i]);
        normal += n * sign(n.dot(normal));
    }
    return normal.normalize();
}

std::vector<Vector3> Contour::randomPointsInside(int numberOfPoints)
{
    std::vector<Vector3> returnedPoints;
    const auto points = curve->getPath();
    if (points.size() < 3) return std::vector<Vector3>();

    if (points.size() == 3) {
        for (int i = 0; i < numberOfPoints; i++) {
            returnedPoints.push_back(randomPointInTriangle(points[0], points[1], points[2]));
        }
        return returnedPoints;
    }

    int maxFailures = 10000 * numberOfPoints;
    Vector3 minVec, maxVec;
    std::tie(minVec, maxVec) = this->curve->AABBox();
    minVec.z() = -1;
    maxVec.z() =  1;
    Vector3 normalRay = this->planeNormal() * (maxVec - minVec).norm();

    for (int i = 0; i < numberOfPoints; i++) {
        // Check the collision from a point below and a point above the plane
        Vector3 randomPoint = Vector3::random(minVec, maxVec) - normalRay;
        Vector3 intersectionPoint = Collision::intersectionRayPlane(randomPoint, normalRay * 2.f, points[0], normalRay);
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

Contour __sub_merge(Contour self, Contour other)
{
    std::vector<Vector3> resultingShape;

    Contour polyShape = self;
    polyShape = polyShape.removeDuplicates();
    Contour clipShape = other;
    clipShape = clipShape.removeDuplicates();

    std::vector<ClipVertex*> poly, clip;

    const auto polyPoints = polyShape.curve->getPath();
    const auto clipPoints = clipShape.curve->getPath();

    poly.reserve((polyPoints.size() + clipPoints.size()) * 4);
    clip.reserve((polyPoints.size() + clipPoints.size()) * 4);
    for (size_t i = 0; i < polyPoints.size(); i++) poly.push_back(new ClipVertex(polyPoints[i]));
    for (size_t i = 0; i < clipPoints.size(); i++) clip.push_back(new ClipVertex(clipPoints[i]));

    bool foundIntersection = false;

    for (size_t i = 0; i < poly.size(); i++) {
        auto& A = poly[i];
        auto& B = poly[getIndex(i + 1, poly.size())];

        for (size_t j = 0; j < clip.size(); j++) {
            auto& C = clip[j];
            auto& D = clip[getIndex(j + 1, clip.size())];
            if ((D->coord - C->coord).norm2() < 0.00001)
                continue;

            if (A->coord == C->coord) {
                poly[getIndex(i, poly.size())]->neighbor = clip[getIndex(j, clip.size())];
                clip[getIndex(j, clip.size())]->neighbor = poly[getIndex(i, poly.size())];
                foundIntersection = true;
            } else if (A->coord == D->coord) {
                poly[getIndex(i, poly.size())]->neighbor = clip[getIndex(j+1, clip.size())];
                clip[getIndex(j+1, clip.size())]->neighbor = poly[getIndex(i, poly.size())];
                foundIntersection = true;
            } else if (B->coord == C->coord) {
                poly[getIndex(i+1, poly.size())]->neighbor = clip[getIndex(j, clip.size())];
                clip[getIndex(j, clip.size())]->neighbor = poly[getIndex(i+1, poly.size())];
                foundIntersection = true;
            } else if (B->coord == D->coord) {
                poly[getIndex(i+1, poly.size())]->neighbor = clip[getIndex(j+1, clip.size())];
                clip[getIndex(j+1, clip.size())]->neighbor = poly[getIndex(i+1, poly.size())];
                foundIntersection = true;
            } else {

                Vector3 intersection = Collision::intersectionBetweenTwoSegments(A->coord, B->coord, C->coord, D->coord);
                if (intersection.isValid()) {
                    foundIntersection = true;

                    int ii = i;
                    int jj = j;
                    if (intersection == A->coord) {}
                    else if (intersection == B->coord) { ii++; }
                    else {
                        poly.insert(poly.begin() + i + 1, new ClipVertex(intersection));
                        ii ++;
                    }
                    if (intersection == C->coord) {}
                    else if (intersection == D->coord) { jj++; }
                    else {
                        clip.insert(clip.begin() + j + 1, new ClipVertex(intersection));
                        jj ++;
                    }

                    poly[ii]->neighbor = clip[jj];
                    clip[jj]->neighbor = poly[ii];

                    j++;
                    i++;
                }
            }
        }
    }

    /// TODO : There is a special case where P0 is on an intersection...
    bool currentlyInside = clipShape.contains(poly[0]->coord);

    if (!foundIntersection) {
        // Shape is completely inside or outside
        if (currentlyInside) {
            // Poly is inside
            return polyShape;
        } else {
            if (polyShape.contains(clip[0]->coord)) {
                // Clipping shape is inside
                return clipShape;
            } else {
                // No intersection
                return Contour();
            }
        }
    } else {
        int firstIntersectionIndex = getIndex(markEntriesExits(poly, clipShape.contains(poly[0]->coord), 0) - 1, poly.size());
        markEntriesExits(clip, polyShape.contains(clip[0]->coord), 1);

        auto& firstVertex = poly[firstIntersectionIndex];
        auto current = firstVertex;
//        resultingShape.push_back(firstVertex.coord);
        while (true) {
            resultingShape.push_back(current->coord);
            if (!current->neighbor) {
                current = current->next;
            } else {
                current = current->neighbor->next;
            }

            // If we're back to the start, end the loop
            if ((current->shapeID == 0 && current->index == firstIntersectionIndex) ||
                    (current->shapeID == 1 && current->neighbor && current->neighbor->index == firstIntersectionIndex))
                break;
        }

        return Contour(resultingShape).removeDuplicates();
    }
}

Contour Contour::merge(Contour other) {
    if (other.curve->empty()) return *this;
    else if (curve->empty()) return other;
    auto res1 = __sub_merge(*this, other);
    auto res2 = __sub_merge(*this, other.curve->reverseVertices());
    return (res1.computeArea() > res2.computeArea() ? res1 : res2);
}

Contour& Contour::resamplePoints(int newNbPoints)
{
    /*CatmullRomSpline::resamplePoints(newNbPoints);
    if (points.size() > 0) {
        points.pop_back(); // Remove last point (as it should be also the first one)
    }*/
    this->curve->resamplePoints();
    return *this;
}

Contour& Contour::setPoint(int i, const Vector3 &newPos)
{
    // CatmullRomSpline::setPoint(i, newPos);
    curve->setPoint(i, newPos);
    return *this;
}

Polyline Contour::computeConvexHull() const
{
    auto points = curve->getPath();
    if (points.empty()) return Polyline();
    // Graham scan's algorithm
    std::vector<Vector3> stack;
    Vector3 start = points[0];
    // Get point with lowest Y (and lowest X in case of tie)
    for (size_t i = 0; i < points.size(); i++) {
        Vector3 p = points[i];
        if (p.y() < start.y() || (p.y() == start.y() && p.x() < start.x())) {
            start = p;
        }
    }
    // Sort points by the minimum angle from the "starting point"
    std::map<float, Vector3> points_angle;
    for (size_t i = 0; i < points.size(); i++) {
        Vector3 dir = (points[i] - start).normalize();
        if (dir == Vector3()) continue; // Ignore the starting point
        float angle = -dir.x(); // Sort from "most right" to "more left"
        if (points_angle.count(angle) == 0 || ((points_angle[angle] - start).norm2() < (points[i] - start).norm2())) {
            points_angle[angle] = points[i];
        }
    }
    // Add the starting point on the stack
    stack.push_back(start);
    // Iterate over all the points:
    while (points_angle.begin() != points_angle.end()) {
        // Remove the points from the stack if they create a "left turn"
        // This can be checked if the Z component of (P1-P0).cross(P2-P0) <= 0
        // With P0 the current point, P1 the top of stack and P2 the second top
        while(stack.size() > 1 && ((stack[stack.size() - 1] - stack[stack.size() - 2]).cross((points_angle.begin()->second - stack[stack.size() - 2])).z() <= 0)) {
            stack.pop_back();
        }
        // Add the point at the end of the stack
        stack.push_back(points_angle.begin()->second);
        points_angle.erase(points_angle.begin());
    }
    return Polyline(stack);
}

std::vector<Vector3> Contour::circle(float radius, const Vector3 &center, int nbPoints)
{
    std::vector<Vector3> points;
    for (int i = 0; i < nbPoints; i++) {
        float angle = (float(i) * 2.f * M_PI) / float(nbPoints);
        points.push_back(Vector3(std::cos(angle) * radius, std::sin(angle) * radius, 0) + center);
    }
    return points;
}


Contour Contour::grow(float increase)
{
    Contour copy(curve->clone());
    std::vector<Vector3> newPoints = copy.removeDuplicates().curve->getPath();
    float previousArea = copy.computeArea();
    const auto points = copy.curve->getPath();
    Vector3 normal = copy.planeNormal();
    std::vector<Vector3> displacement(points.size());

    for (size_t i = 0; i < points.size(); i++) {
        float t = float(i) / float(curve->numSegments());
        displacement[i] = curve->getNormal(t, normal) * increase;
    }
    for (size_t i = 0; i < points.size(); i++) {
        std::cout << i << " " << displacement[i] << std::endl;
        copy.curve->setPoint(i, copy.curve->get(i) + displacement[i]);
    }
    if (sign(copy.computeArea() - previousArea) != sign(increase)) {
        // Made an error in the direction, need to rectify
        for (size_t i = 0; i < points.size(); i++) {
            copy.curve->setPoint(i, copy.curve->get(i) - 2.f * displacement[i]);
        }
    }
    // Vector3 mid = (copy.curve->get(0) + copy.curve->get(-1)) * .5f;
    // copy.curve->setPoint(0, mid);
    // copy.curve->setPoint(-1, mid);
    return copy;
}

Contour Contour::shrink(float decrease)
{
    return this->grow(-decrease);
}

Contour& Contour::removeDuplicates()
{
    curve->removeDuplicates();
    return *this;
    /*
    CatmullRomSpline::removeDuplicates();
    std::vector<Vector3> foundPattern;
    for (size_t i = 0; i < points.size(); i++) {
        for (size_t j = i + 1; j < points.size(); j++) {
            if (points[i] == points[j]) {
                points.erase(points.begin() + j, points.end());
                break;
            }
        }
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
    if (points.size() > 1 && (points[0] - points.back()).norm2() < 0.01) {
        points.pop_back();
        this->closed = true;
    }
    return *this;
    */
}

std::vector<Vector3> Contour::closedPath() const
{
    /*std::vector<Vector3> res = points;
    if (!points.empty())
        res.push_back(points[0]);
    return res;
    */
    return curve->getPath();
}




std::vector<float> computeGreenCoordinates(const Vector3& p, const Contour& polygon) {
    auto vertices = polygon.curve->getPath();
    size_t n = vertices.size();
    std::vector<float> weights(n, 0.0);
    std::vector<float> tans(n, 0.0);

    for (size_t i = 0; i < n; i++) {
        Vector3 v1 = vertices[i] - p;
        Vector3 v2 = vertices[(i + 1) % n] - p;

        float dist1 = v1.magnitude();
        float dist2 = v2.magnitude();

        // Avoid division by zero by ensuring no zero-length vectors
        if (dist1 == 0 || dist2 == 0) return {};

        float cosine = v1.dot(v2) / (dist1 * dist2);
        float sine = v1.cross(v2).norm() / (dist1 * dist2);

        // Calculate tan of half the angle between v1 and v2
        float angle = atan2(sine, cosine);
        tans[i] = tan(angle / 2);
    }

    float totalWeight = 0.0;
    for (size_t i = 0; i < n; i++) {
        size_t prev = (i + n - 1) % n;
        float weight = (tans[prev] + tans[i]) / (vertices[i] - p).magnitude();
        weights[i] = weight;
        totalWeight += weight;
    }

    // Normalize weights
    for (float& weight : weights) {
        weight /= totalWeight;
    }

    return weights;

    /*std::vector<float> greenCoords;

    for (size_t i = 0; i < polygon.size(); i ++) {
        const Vector3& a = polygon[i];
        const Vector3& b = polygon[i + 1];
        const Vector3& c = polygon[i + 2];

        // Check if the point is inside the triangle formed by a, b, and c
        if (Collision::pointInPolygon(p, {a, b, c})) {
            // Compute the barycentric coordinates of the point with respect to triangle abc
            Vector3 v0 = c - a;
            Vector3 v1 = b - a;
            Vector3 v2 = p - a;

            float dot00 = v0.dot(v0);
            float dot01 = v0.dot(v1);
            float dot02 = v0.dot(v2);
            float dot11 = v1.dot(v1);
            float dot12 = v1.dot(v2);

            float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
            float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
            float w = 1.0f - u - v;

            // Store the barycentric coordinates
            greenCoords.push_back(w);
            greenCoords.push_back(v);
            greenCoords.push_back(u);

            // Assuming the polygon is simple, return the green coordinates once found
            // return greenCoords;
        } else {
            greenCoords.push_back(0);
            greenCoords.push_back(0);
            greenCoords.push_back(0);
        }
    }

    // Point is outside the polygon
    return greenCoords;*/
}

Vector3 computePointFromGreenCoordinates(const std::vector<float>& greenCoords, const Contour& polygon) {
    const auto& vertices = polygon.curve->getPath();
    const auto& weights = greenCoords;
    Vector3 p;
    for (size_t i = 0; i < vertices.size(); i++) {
        p.x() += vertices[i].x() * weights[i];
        p.y() += vertices[i].y() * weights[i];
    }
    return p;

    /*
    Vector3 p(0.0, 0.0, 0.0); // Initialize point P

    // Interpolate the position of the point based on the barycentric coordinates
    for (size_t i = 0; i < polygon.size(); i ++) {
        const Vector3& a = polygon[i];
        const Vector3& b = polygon[i + 1];
        const Vector3& c = polygon[i + 2];

        // Extract barycentric coordinates for triangle abc
        float u = greenCoords[i];
        float v = greenCoords[i + 1];
        float w = greenCoords[i + 2];

        // Compute the interpolated position of the point P using barycentric coordinates
        p += (u * a + v * b + w * c);
    }

    return p;*/
}

Vector3 randomPointInTriangle(const Vector3 &A, const Vector3 &B, const Vector3 &C)
{
    float s = random_gen::generate();
    float t = random_gen::generate();

    return A + (s + t <= 1 ? s * (B - A) + t * (C - A) : (1 - s) * (B - A) + (1 - t) * (C - A));
}
