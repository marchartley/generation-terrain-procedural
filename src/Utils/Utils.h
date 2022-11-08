#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <set>
#include "DataStructure/Vector3.h"

#define PI 3.14159265358979323846

template<class T>
bool isIn(T& elem, std::vector<T> arr) {
    return std::find(arr.begin(), arr.end(), elem) != arr.end();
}

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

std::string toUpper(std::string s);
std::string toLower(std::string s);
std::string getExtention(std::string file);

float rad2deg(float rad);
float deg2rad(float deg);

/// Careful, the order of the vectors are not preserved in these functions
template <class T>
std::vector<T> convertSetToVector(std::set<T> _set) {
    return std::vector<T>(_set.begin(), _set.end());
}
template <class T>
std::set<T> convertVectorToSet(std::vector<T> _vector) {
    return std::set<T>(_vector.begin(), _vector.end());
}
template <class T>
std::vector<T> removeDuplicatesFromVector(std::vector<T> _vector) {
    return convertSetToVector(convertVectorToSet(_vector));
}
template <class T>
std::vector<T> vectorMerge(std::vector<T> v1, std::vector<T> v2) {
    std::vector<T> result = v1;
    result.insert(result.end(), v2.begin(), v2.end());
    return result;
}
template <class T>
std::vector<T> vectorUnion(std::vector<T> v1, std::vector<T> v2) {
    return removeDuplicatesFromVector(vectorMerge(v1, v2));
}
template <class T>
std::vector<T> vectorIntersection(std::vector<T> v1, std::vector<T> v2) {
    std::vector<T> result;
    // Remove duplicates and sort the array by the same time
    v1 = removeDuplicatesFromVector(v1);
    v2 = removeDuplicatesFromVector(v2);
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
    return result;
}

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

    // Found in A Review of Digital Terrain Modeling (Eric Galin, Eric Guérin, Adrien Peytavie, Guillaume Cordonnier, Marie-Paule
    // Cani, Bedrich Benes, James Gain) [2019]
    // Used to generate a terrain from random faults
    float fault_distance(float distance, float impactRadius);
}


template <class T>
T lerp(float t, T min, T max) {
    return min + (max - min) * t;
}
template <class T>
float inverseLerp(T val, T min, T max) {
    return (val - min) / (max - min);
}

template <class T>
T clamp(T val, T min, T max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

template <class T, class U>
U remap(T val, T oldMin, T oldMax, U newMin, U newMax)
{
    float oldProgress = inverseLerp(val, oldMin, oldMax);
    return lerp(oldProgress, newMin, newMax);
}

float gaussian(Vector3 size, Vector3 position, float sigma);
float normalizedGaussian(Vector3 size, Vector3 position, float sigma);

#endif // UTILS_H
