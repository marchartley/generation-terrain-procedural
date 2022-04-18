#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include "DataStructure/Vector3.h"

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

Vector3 intersectionBetweenTwoSegments(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4);
float shortestDistanceBetweenSegments(Vector3 p11, Vector3 p12, Vector3 p21, Vector3 p22);
float tetrahedronSignedVolume(Vector3 a, Vector3 b, Vector3 c, Vector3 d);
int sign(float a);
int segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3);
Vector3 intersectionRayPlane(Vector3 rayOrigin, Vector3 rayDir, Vector3 planeCenter, Vector3 planeNormal);
Vector3 intersectionRaySphere(Vector3 rayOrigin, Vector3 rayDir, Vector3 sphereCenter, float sphereRadius, bool returnClosestPoint = true);

std::string toUpper(std::string s);
std::string getExtention(std::string file);

namespace interpolation {
    float linear(float x, float _min = 0.0, float _max = 1.0);
    float inv_linear(float x, float _min = 0.0, float _max = 1.0);
    float sigmoid(float _x, float lambda = 10.0, float offset = -0.5, float _min = 0.0, float _max = 1.0);
    float smooth(float _x, float _min = 0.0, float _max = 1.0);
    float quadratic(float _x, float _min = 0.0, float _max = 1.0);
    float cubic(float _x, float _min = 0.0, float _max = 1.0);
    float cosine(float _x, float _min = 0.0, float _max = 1.0);
    float binary(float _x, float _min = 0.0, float _max = 1.0);
    float wyvill(float _x, float _min = 0.0, float _max = 1.0);
}


template <class T>
T lerp(float t, T min, T max) {
    return min + (max - min) * t;
}
template <class T>
float inverseLerp(T val, T min, T max) {
    return (val - min) / (max - min);
}

template <class T, class U>
U remap(T val, T oldMin, T oldMax, U newMin, U newMax)
{
    float oldProgress = inverseLerp(val, oldMin, oldMax);
    return lerp(oldProgress, newMin, newMax);
}

#endif // UTILS_H
