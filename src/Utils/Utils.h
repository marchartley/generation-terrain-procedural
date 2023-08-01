#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <set>
#include "DataStructure/Vector3.h"

#define PI 3.14159265358979323846

template<class T>
bool isIn(T elem, std::vector<T> arr) {
    return std::find(arr.begin(), arr.end(), elem) != arr.end();
}

void displayProgress(float percent, bool displayPercent = true, int consoleWidth = 20, std::string progressSign = "=");

std::string replaceInString(std::string initial, std::string toReplace, std::string replacing);

bool startsWith(std::string text, std::string needle);
bool endsWith(std::string text, std::string needle);
std::vector<std::string> split(std::string str, std::string c);
std::vector<std::string> split(std::string str);
bool makedir(std::string path);
bool checkPathExists(std::string path);
Vector3 HSVtoRGB(float H, float S,float V);

std::string toUpper(std::string s);
std::string toLower(std::string s);
std::string toCapitalize(std::string s);
std::string getExtension(std::string file);
std::string getFilename(std::string path);
std::string simplify(std::string s);

float rad2deg(float rad);
float deg2rad(float deg);

void sleep(int milliseconds);

//template <typename ... FArgs, typename ... Args>
//float timeIt(std::function<void(FArgs...)> func,
//                  Args && ... as) {
//    auto start = std::chrono::system_clock::now();
//    func(std::forward<Args>(as)...);
//    auto end = std::chrono::system_clock::now();
//    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//}
float timeIt(std::function<void()> func);


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

float gaussian(float sigma, float sqrDist);
float gaussian(const Vector3 &size, const Vector3 &position, float sigma);
float normalizedGaussian(float sigma, float sqrDist);
float normalizedGaussian(const Vector3& size, const Vector3& position, float sigma);


template<typename T>
std::vector<T> flattenArray(std::vector<std::vector<T>> arr) {
    std::vector<T> finalArray;
    for (const std::vector<T>& val : arr)
        finalArray.insert(finalArray.end(), val.begin(), val.end());
    return finalArray;
}

// COmpletely stolen from : https://stackoverflow.com/a/54512651
template <typename Iterator>
std::string join(Iterator begin, Iterator end, std::string separator = "")
{
    std::ostringstream o;
    if(begin != end)
    {
        o << *begin++;
        for(;begin != end; ++begin)
            o  << separator << *begin;
    }
    return o.str();
}

template <typename Container>
std::string join(Container const& c, std::string separator = "")
{
    using std::begin;
    using std::end;
    return join(begin(c), end(c), separator);
    // not using std::... directly:
    // there might be a non-std overload that wouldn't be found if we did
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


#endif // UTILS_H
