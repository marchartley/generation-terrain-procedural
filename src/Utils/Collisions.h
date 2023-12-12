#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "DataStructure/Vector3.h"
#include "DataStructure/Triangle.h"

namespace Collision {
Vector3 intersectionBetweenTwoSegments(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);
float shortestDistanceBetweenSegments(const Vector3& p11, const Vector3& p12, const Vector3& p21, const Vector3& p22);
float tetrahedronSignedVolume(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d);
int sign(float a);
Vector3 segmentToTriangleCollision(const Vector3& s1, const Vector3& s2, const Vector3& t1, const Vector3& t2, const Vector3& t3, bool strict = false);
Vector3 intersectionRayPlane(const Vector3& rayOrigin, const Vector3& rayDir, const Vector3& planeCenter, const Vector3& planeNormal, bool limitRayLength = false);
Vector3 intersectionRaySphere(const Vector3& rayOrigin, const Vector3& rayDir, const Vector3& sphereCenter, float sphereRadius, bool returnClosestPoint = true);
bool intersectionTriangleAABBox(const Vector3& t0, const Vector3& t1, const Vector3& t2, const Vector3& minAABBox, const Vector3& maxAABBox);
bool intersectionTriangleAABBox(const Vector3& t0, const Vector3& t1, const Vector3& t2, const Vector3& boxCenter, const Vector3& halfSizeX, const Vector3& halfSizeY, const Vector3& halfSizeZ);
bool intersectionAABBoxPlane(const Vector3& boxMin, const Vector3& boxMax, const Triangle& triangle);
Vector3 intersectionRayAABBox(const Vector3& orig, const Vector3& dir, const AABBox& box);
Vector3 intersectionRayAABBox(const Vector3& orig, const Vector3& dir, const Vector3& boxMin, const Vector3& boxMax);

Vector3 projectPointOnSegment(const Vector3& point, const Vector3& segmentStart, const Vector3& segmentEnd);
Vector3 projectPointOnSphere(const Vector3& point, const Vector3& sphereCenter, float sphereRadius);

bool pointInPolygon(const Vector3& point, std::vector<Vector3> polygon);

float pointToTriangleDistanceSquared(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& point);
}


#endif // COLLISIONS_H
