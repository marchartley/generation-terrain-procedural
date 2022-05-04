#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "DataStructure/Vector3.h"

namespace Collision {
Vector3 intersectionBetweenTwoSegments(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4);
float shortestDistanceBetweenSegments(Vector3 p11, Vector3 p12, Vector3 p21, Vector3 p22);
float tetrahedronSignedVolume(Vector3 a, Vector3 b, Vector3 c, Vector3 d);
int sign(float a);
Vector3 segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3);
Vector3 intersectionRayPlane(Vector3 rayOrigin, Vector3 rayDir, Vector3 planeCenter, Vector3 planeNormal);
Vector3 intersectionRaySphere(Vector3 rayOrigin, Vector3 rayDir, Vector3 sphereCenter, float sphereRadius, bool returnClosestPoint = true);

}

#endif // COLLISIONS_H
