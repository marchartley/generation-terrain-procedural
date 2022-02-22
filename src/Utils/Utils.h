#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include "DataStructure/Vector3.h"

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

Vector3 intersectionPoint(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4);
bool intersection(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4);

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
