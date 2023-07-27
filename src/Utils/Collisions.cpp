#include "Collisions.h"


// Source : http://paulbourke.net/geometry/pointlineplane/
Vector3 Collision::intersectionBetweenTwoSegments(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4)
{
    Vector3 l21 = (p1 - p2); // .normalized();
    Vector3 l13 = (p3 - p1); // .normalized();
    Vector3 l43 = (p3 - p4); // .normalized();

    float d1321 = l13.dot(l21);
    float d1343 = l13.dot(l43);
    float d4321 = l43.dot(l21);
    float d4343 = l43.dot(l43);
    float d2121 = l21.dot(l21);
    float d2143 = l21.normalized().dot(l43.normalized());

    if (std::abs(std::abs(d2143) - 1) < 0.001 || std::abs((d2121*d4343 - d4321*d4321)) < 0.001) return Vector3(false); // Parallel lines?
    float mu_a = (d1343*d4321 - d1321*d4343) / (d2121*d4343 - d4321*d4321);
    float mu_b = (d1343 + mu_a*d4321) / d4343;

    if (mu_a < 0.0 || 1.0 < mu_a || mu_b < 0.0 || 1.0 < mu_b)
        return Vector3(false);

    return p1 + (p2 - p1) * mu_a;
}

float Collision::shortestDistanceBetweenSegments(const Vector3& p11, const Vector3& p12, const Vector3& p21, const Vector3& p22)
{
    Vector3 d1 = p12 - p11;
    Vector3 d2 = p22 - p21;

    Vector3 c1, c2;
    // Extreme case, same first point
    if (p11 == p21)
        return 0;
    // Case of parallel lines
    else if (std::abs(d1.dot(d2)) == (d1.norm() * d2.norm())) {
        Vector3 v = d1.normalized();
        c1 = p11 - v * (p11.dot(v));
        c2 = p21 - v * (p21.dot(v));
        //return (p21 - p11).cross(d1).norm();
    } else {
        Vector3 n = d1.cross(d2);
        Vector3 n1 = d1.cross(n);
        Vector3 n2 = d2.cross(n);

        // Closest point on line p11-p12
        float t1 = (p21 - p11).dot(n2)/d1.dot(n2);
        t1 = std::max(std::min(t1, 1.f), 0.f);
        c1 = p11 + d1 * t1;
        // Closest point on line p21-p22
        float t2 = (p11 - p21).dot(n1)/d2.dot(n1);
        t2 = std::max(std::min(t2, 1.f), 0.f);
        c2 = p21 + d2 * t2;

    }
    return (c2 - c1).norm();
}
float Collision::tetrahedronSignedVolume(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    return (1.f/6.f) * (c-a).cross((b-a)).dot((d-a));
}
int sign(float a){return (a < 0 ? -1 : (a > 0 ? 1 : 0)); }
Vector3 Collision::segmentToTriangleCollision(const Vector3& s1, const Vector3& s2, const Vector3& t1, const Vector3& t2, const Vector3& t3, bool strict)
{
    // Möller–Trumbore intersection algorithm
//    if (int(s1.x) == 10 && int(s1.z) == 10) {
//        int a = 0;
//    }
    Vector3 rayOrigin = s1;
    Vector3 rayDir = (s2 - s1);

    Vector3 triEdge1 = (t2 - t1);
    Vector3 triEdge2 = (t3 - t1);

    Vector3 h = rayDir.cross(triEdge2);
    float dot = triEdge1.dot(h);
    if (std::abs(dot) <  1.e-8)
        return Vector3(false); // Ray parallel to the triangle
    float f = 1.f/dot;
    Vector3 s = (rayOrigin - t1);
    float u = f * s.dot(h);
    if (u < 0.f || 1.f < u)
        return Vector3(false); // Ray did not reach the triangle
    Vector3 q = s.cross(triEdge1);
    float v = f * rayDir.dot(q);
    if (v < 0.f || 1.f < (u + v))
        return Vector3(false);

    float t = f * triEdge2.dot(q);

    if (t < 0.f || 1.f < t)
        return Vector3(false); // Intersection before or after the ray

    Vector3 intersectionPoint = rayOrigin + rayDir * t;
    if (!strict && (intersectionPoint == t1 || intersectionPoint == t2 || intersectionPoint == t3)) // Edge case of collision exactly on the corner
        return Vector3(false);

    return intersectionPoint;
}

Vector3 Collision::intersectionRayPlane(const Vector3& rayOrigin, const Vector3& rayDir, const Vector3& planeCenter, const Vector3& planeNormal, bool limitRayLength)
{
    // assuming vectors are all normalized
    float denom = planeNormal.dot(rayDir);
    if (std::abs(denom) > 1e-6) {
        Vector3 planeToRay = rayOrigin - planeCenter;
        if (planeToRay.dot(rayDir) > 0) {
            return Vector3(false); // Direction pointing the wrong side = no intersection
        }
        float t = std::abs(planeToRay.dot(planeNormal) / denom);
        if (limitRayLength && (t < 0 || 1 < t))
            return Vector3(false);
        return rayOrigin + rayDir * t;
    }

    return Vector3(false);
}

Vector3 Collision::intersectionRaySphere(const Vector3& rayOrigin, const Vector3& _rayDir, const Vector3& sphereCenter, float sphereRadius, bool returnClosestPoint)
{
    // Using Joachimsthal's Equation
    // Found from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // ad^2 + bd + c = 0 with a = ||u||^2, b = 2 * (u dot rayToCenter) and c = ||rayToCenter||^2 - r^2
    // Here ||u|| = 1
    Vector3 rayDir = _rayDir.normalized();
    Vector3 centerToRay = rayOrigin - sphereCenter;
    //float a = 1.f;
    //float b = 2 * rayDir.dot(centerToRay);
    float c = centerToRay.norm2() - sphereRadius * sphereRadius;

    // Solutions are d = -u dot centerToRay +/- sqrt(root)
    // with root = (u dot centerToRay)^2 - (||centerToRay||^2 - radius^2)
    float root = rayDir.dot(centerToRay) * rayDir.dot(centerToRay) - c;

    float d1 = -(rayDir.dot(centerToRay));
    float d2 = -(rayDir.dot(centerToRay));
    if (root < 0) {
        // No solution
        return Vector3(false);
    } else if (root == 0) {
        // One unique solution
        // Use d1
    } else {
        // Two solutions, choose the closest or furthest
        float sqrtRoot = std::sqrt(root);
        d1 -= sqrtRoot;
        d2 += sqrtRoot;
    }
    if (d1 < 0 && d2 >= 0) // d1 is behind
        return rayOrigin + rayDir * d2;
    if (d1 >= 0 && d2 < 0) // d2 is behind
        return rayOrigin + rayDir * d1;
    if (d1 < 0 && d2 < 0) // All behind
        return Vector3(false);

    if (returnClosestPoint) {
        // The closest must be d1
        return rayOrigin + rayDir * d1;
    } else {
        return rayOrigin + rayDir * d2;
    }
}

bool Collision::intersectionTriangleAABBox(const Vector3& t0, const Vector3& t1, const Vector3& t2, const Vector3& minAABBox, const Vector3& maxAABBox)
{
    Vector3 box = maxAABBox - minAABBox;
    bool planeIntersection = Collision::intersectionAABBoxPlane(minAABBox, maxAABBox, {t0, t1, t2});
    bool triangleIntersection = Collision::intersectionTriangleAABBox(t0, t1, t2, (minAABBox + maxAABBox)/2, Vector3(box.x * .5f, 0, 0), Vector3(0, box.y * .5f, 0), Vector3(0, 0, box.z * .5f));
    return planeIntersection && triangleIntersection;
}

std::tuple<float, float> project(std::vector<Vector3> vertices, const Vector3& _axis) {
    float proj_min = std::numeric_limits<float>::max(), proj_max = std::numeric_limits<float>::min();
    float projection;
    Vector3 axis = _axis.normalized();

    for (auto& v : vertices) {
        projection = axis.dot(v);
        proj_min = std::min(proj_min, projection);
        proj_max = std::max(proj_max, projection);
    }
    return std::make_tuple(proj_min, proj_max);
}
bool Collision::intersectionTriangleAABBox(const Vector3& t0, const Vector3& t1, const Vector3& t2, const Vector3& boxCenter, const Vector3& halfSizeX, const Vector3& halfSizeY, const Vector3& halfSizeZ)
{
    // Taken from https://stackoverflow.com/a/17503268/9863298

    std::vector<Vector3> boxNormals = {halfSizeX.normalized(), halfSizeY.normalized(), halfSizeZ.normalized()};
    std::vector<Vector3> triPoints = {t0, t1, t2};
    std::vector<Vector3> triEdges = {t1 - t0, t2 - t1, t0 - t2};
    std::vector<Vector3> cubeVertices = {boxCenter - halfSizeX - halfSizeY - halfSizeZ,
                                         boxCenter - halfSizeX - halfSizeY + halfSizeZ,
                                         boxCenter - halfSizeX + halfSizeY - halfSizeZ,
                                         boxCenter - halfSizeX + halfSizeY + halfSizeZ,
                                         boxCenter + halfSizeX - halfSizeY - halfSizeZ,
                                         boxCenter + halfSizeX - halfSizeY + halfSizeZ,
                                         boxCenter + halfSizeX + halfSizeY - halfSizeZ,
                                         boxCenter + halfSizeX + halfSizeY + halfSizeZ,
                                         };

    // Check that the projected triangle intersect the projected box on all the axis
    float proj_min, proj_max, proj_min2, proj_max2;
    for (int i = 0; i < 3; i++) {
        std::tie(proj_min, proj_max) = project({t0, t1, t2}, boxNormals[i]);
        std::tie(proj_min2, proj_max2) = project(cubeVertices, boxNormals[i]);
        if (proj_min > proj_max2 || proj_max < proj_min2)
            return false;
    }

    Vector3 normal = (t0 - t1).cross(t2 - t1).normalize();
    float triOffset = normal.dot(t0);
    std::tie(proj_min, proj_max) = project(cubeVertices, normal);
    if (proj_min > triOffset || proj_max < triOffset)
        return false;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Vector3 axis = triEdges[i].cross(boxNormals[j]);
            std::tie(proj_min, proj_max) = project(cubeVertices, axis);
            std::tie(proj_min2, proj_max2) = project(triPoints, axis);
            if(proj_max <= proj_min2 || proj_min >= proj_max2)
                return false;
        }
    }
    return true;
}

Vector3 Collision::intersectionRayAABBox(const Vector3& orig, const Vector3& dir, const Vector3& boxMin, const Vector3& boxMax)
{
//    Vector3 box = (boxMax - boxMin);
//    orig -= boxMin;

    float tmin = (boxMin.x - orig.x) / dir.x;
    float tmax = (boxMax.x - orig.x) / dir.x;

    if (tmin > tmax) std::swap(tmin, tmax);

    float tymin = (boxMin.y - orig.y) / dir.y;
    float tymax = (boxMax.y - orig.y) / dir.y;

    if (tymin > tymax) std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
        return Vector3(false);

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (boxMin.z - orig.z) / dir.z;
    float tzmax = (boxMax.z - orig.z) / dir.z;

    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
        return Vector3(false);

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    if (tmin < 0) {
        if (tmax < 0) {
            return Vector3(false);
        } else {
            return orig + dir * tmax;
        }
    } else {
        return orig + dir * tmin;
    }

}

bool Collision::pointInPolygon(const Vector3& point, std::vector<Vector3> _polygon)
{
    std::vector<Vector3> polygon;
    for (auto& p : _polygon)
        if (polygon.empty() || (p - polygon.back()).norm2() > 0.01)
            polygon.push_back(p);

    if (polygon.size() < 2) return false;

    size_t firstIndex = 0, secondIndex = 0;

    Vector3 firstVertex = polygon[firstIndex];
    Vector3 secondVertex = polygon[secondIndex];
    Vector3 center;
    float polygonLength = 0;
    for (size_t i = 0; i < polygon.size(); i++) {
        center += polygon[i];
        if (i > 0)
            polygonLength += (polygon[i] - polygon[i - 1]).norm();
    }
    center /= (float)polygon.size();
    // First, lets find a good plane to define our shape
    while (std::abs((firstVertex - center).normalized().dot((secondVertex - center).normalized())) < (1 - 1e-5)) {
        secondIndex ++;
        if (secondIndex >= polygon.size()) {
            secondIndex = 0;
            firstIndex++;
            if (firstIndex >= polygon.size()) return false;
            firstVertex = polygon[firstIndex];
        }
        secondVertex = polygon[secondIndex];
    }
    // At this point, we can check if the "pos" is in the same plane as the shape

    Vector3 ray = center + (firstVertex - center).normalized() * polygonLength; // Same, should be outside the shape
    polygon.push_back(polygon[0]);
    // Check the intersection of the ray with all the segments of the shape
    int nb_intersections = 0;
    for (size_t i = 0; i < polygon.size() - 1; i++) {
        if (Collision::intersectionBetweenTwoSegments(point, ray, polygon[i], polygon[i+1]).isValid())
            nb_intersections++;
    }
    // If there's a odd number of intersections, the point is inside
    return (nb_intersections % 2) == 1;
}

Vector3 Collision::projectPointOnSegment(const Vector3& point, const Vector3& segmentStart, const Vector3& segmentEnd)
{
    Vector3 startToPoint = point - segmentStart;
    Vector3 segment = segmentEnd - segmentStart;
    float t = startToPoint.dot(segment) / segment.norm2();
    t = std::max(0.f, std::min(1.f, t)); // Stay on segment
    return segmentStart + t * segment;
}

Vector3 Collision::projectPointOnSphere(const Vector3& point, const Vector3& sphereCenter, float sphereRadius)
{
    return sphereCenter + (point - sphereCenter).normalize() * sphereRadius;
}

float Collision::pointToTriangleDistanceSquared(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& point) {
    Vector3 edge0 = p2 - p1;
    Vector3 edge1 = p3 - p1;
    Vector3 v0 = p1 - point;

    float a = edge0.dot(edge0);
    float b = edge0.dot(edge1);
    float c = edge1.dot(edge1);
    float d = edge0.dot(v0);
    float e = edge1.dot(v0);

    float det = a*c - b*b;
    float s = b*e - c*d;
    float t = b*d - a*e;

    if (s + t < det) {
        if (s < 0.0f) {
            if (t < 0.0f) {
                if (d < 0.0f) {
                    s = std::clamp(-d/a, 0.0f, 1.0f);
                    t = 0.0f;
                } else {
                    s = 0.0f;
                    t = std::clamp(-e/c, 0.0f, 1.0f);
                }
            } else {
                s = 0.0f;
                t = std::clamp(-e/c, 0.0f, 1.0f);
            }
        } else if (t < 0.0f) {
            s = std::clamp(-d/a, 0.0f, 1.0f);
            t = 0.0f;
        } else {
            float invDet = 1.0f / det;
            s *= invDet;
            t *= invDet;
        }
    } else {
        if (s < 0.0f) {
            float tmp0 = b + d;
            float tmp1 = c + e;
            if (tmp1 > tmp0) {
                float numer = tmp1 - tmp0;
                float denom = a - 2*b + c;
                s = std::clamp(numer/denom, 0.0f, 1.0f);
                t = 1 - s;
            } else {
                t = std::clamp(-e/c, 0.0f, 1.0f);
                s = 0.0f;
            }
        } else if (t < 0.0f) {
            if (a + d > b + e) {
                float numer = c + e - b - d;
                float denom = a - 2*b + c;
                s = std::clamp(numer/denom, 0.0f, 1.0f);
                t = 1 - s;
            } else {
                s = std::clamp(-e/c, 0.0f, 1.0f);
                t = 0.0f;
            }
        } else {
            float numer = c + e - b - d;
            float denom = a - 2*b + c;
            s = std::clamp(numer/denom, 0.0f, 1.0f);
            t = 1 - s;
        }
    }

    Vector3 closestPointOnTriangle = p1 + s * edge0 + t * edge1;
    Vector3 diff = closestPointOnTriangle - point;

    return diff.dot(diff); // Return the squared distance
}

bool Collision::intersectionAABBoxPlane(const Vector3 &boxMin, const Vector3 &boxMax, const Triangle &triangle)
{
    // Get box center and half-diagonal (vector from center to corner)
    Vector3 boxCenter = (boxMin + boxMax) / 2;
    Vector3 boxHalfDiagonal = boxMax - boxCenter;

    // Compute dot product of box half-diagonal with plane normal
    float dot = boxHalfDiagonal.dot(triangle.normal);

    // Compute signed distance from box center to plane
    float dist = triangle.normal.dot(boxCenter) + triangle.d;

    // Box intersects plane if distance to plane is within box's extent along plane normal
    return std::abs(dist) <= std::abs(dot * 1.5f); // (very big epsilon...)
}
