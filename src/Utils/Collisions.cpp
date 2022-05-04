#include "Collisions.h"


// Source : http://paulbourke.net/geometry/pointlineplane/
Vector3 Collision::intersectionBetweenTwoSegments(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4)
{
    Vector3 l21 = (p1 - p2);
    Vector3 l13 = (p3 - p1);
    Vector3 l43 = (p3 - p4);

    float d1321 = l13.dot(l21);
    float d1343 = l13.dot(l43);
    float d4321 = l43.dot(l21);
    float d4343 = l43.dot(l43);
    float d2121 = l21.dot(l21);

    if (std::abs((d2121*d4343 - d4321*d4321)) < 0.001) return Vector3(false); // Parallel lines?
    float mu_a = (d1343*d4321 - d1321*d4343) / (d2121*d4343 - d4321*d4321);
    float mu_b = (d1343 + mu_a*d4321) / d4343;

    if (mu_a < 0.0 || 1.0 < mu_a || mu_b < 0.0 || 1.0 < mu_b)
        return Vector3(false);

    return p1 + (p2 - p1) * mu_a;
}/*
bool intersectionBetweenTwoSegments(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4)
{
    Vector3 l21 = (p1 - p2);
    Vector3 l13 = (p3 - p1);
    Vector3 l43 = (p3 - p4);

    float d1321 = l13.dot(l21);
    float d1343 = l13.dot(l43);
    float d4321 = l43.dot(l21);
    float d4343 = l43.dot(l43);
    float d2121 = l21.dot(l21);

    if (std::abs((d2121*d4343 - d4321*d4321)) < 0.001) return false; // Parallel lines?
    float mu_a = (d1343*d4321 - d1321*d4343) / (d2121*d4343 - d4321*d4321);
    float mu_b = (d1343 + mu_a*d4321) / d4343;
//    std::cout << "(mu_a = " << mu_a << " and mu_b = " << mu_b << ")";
    return (0 <= mu_a) && (mu_a <= 1.0) && (0.0 <= mu_b) && (mu_b <= 1.0);
}*/

float Collision::shortestDistanceBetweenSegments(Vector3 p11, Vector3 p12, Vector3 p21, Vector3 p22)
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
float Collision::tetrahedronSignedVolume(Vector3 a, Vector3 b, Vector3 c, Vector3 d) {
    return (1.f/6.f) * (c-a).cross((b-a)).dot((d-a));
}
int sign(float a){return (a < 0 ? -1 : (a > 0 ? 1 : 0)); }
Vector3 Collision::segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3)
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
    if (std::abs(dot) <  1.e-6)
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
    return rayOrigin + rayDir * t;
}

Vector3 Collision::intersectionRayPlane(Vector3 rayOrigin, Vector3 rayDir, Vector3 planeCenter, Vector3 planeNormal)
{
    // assuming vectors are all normalized
    float denom = planeNormal.dot(rayDir);
    if (std::abs(denom) > 1e-6) {
        Vector3 planeToRay = rayOrigin - planeCenter;
        if (planeToRay.dot(rayDir) > 0) {
            return Vector3(false); // Direction pointing the wrong side = no intersection
        }
        float t = std::abs(planeToRay.dot(planeNormal) / denom);
        return rayOrigin + rayDir * t;
    }

    return Vector3(false);
}

Vector3 Collision::intersectionRaySphere(Vector3 rayOrigin, Vector3 rayDir, Vector3 sphereCenter, float sphereRadius, bool returnClosestPoint)
{
    // Using Joachimsthal's Equation
    // Found from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // ad^2 + bd + c = 0 with a = ||u||^2, b = 2 * (u dot rayToCenter) and c = ||rayToCenter||^2 - r^2
    // Here ||u|| = 1
    rayDir.normalize();
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
