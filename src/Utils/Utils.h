#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include "DataStructure/Vector3.h"

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

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
