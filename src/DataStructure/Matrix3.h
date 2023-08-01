#ifndef MATRIX3_H
#define MATRIX3_H

#include <iostream>
#include <vector>
#include <tuple>
#include <memory>
#include <string>
#include <math.h>
#include <iomanip>
#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"
#include "Utils/Collisions.h"
#include "Utils/Utils.h"

#include "Utils/FastNoiseLit.h"

enum RESIZE_MODE {
    LINEAR = 0,
    MAX_VAL = 1,
    MIN_VAL = 2,
    FILL_WITH_DEFAULT = 3,
    NEAREST = 4
};

enum CONVOLUTION_BORDERS {
    ZERO_PAD = 0,
    REPEAT = 1,
    MIRROR = 2,
    WRAPPING = 3,
    IGNORED = 4,
    COPY = 5
};

enum RETURN_VALUE_ON_OUTSIDE {
    DEFAULT_VALUE = 0,
    MIRROR_VALUE = 1,
    WRAPPED_VALUE = 2,
    REPEAT_VALUE = 3
};

// Warning : don't use bool type...
// This class is based on a std::vector, which has specifications on bools
// the [] operator won't work ... Use int or short int instead
template <class T>
class Matrix3
{
public:
    Matrix3();
    Matrix3(size_t sizeX, size_t sizeY, size_t sizeZ = 1, T initValue = T());
    Matrix3(Vector3 size, T initValue = T());
    Matrix3(std::vector<std::vector<std::vector<T>>> data);
    Matrix3(std::vector<std::vector<T>> data);
    Matrix3(std::vector<T> data, size_t sizeX, size_t sizeY, int sizeZ = -1);

    const T& at(int i, int j, int k = 0) const;
    const T& at(Vector3 pos) const;
    const T& at(size_t i) const;
    const T& operator()(size_t x, size_t y, size_t z = 0) const;
    const T& operator[](size_t i) const;
    const T& operator[](Vector3 pos) const;
    T& at(int i, int j, int k = 0);
    T& at(Vector3 pos);
    T& at(size_t i);
    T& operator()(size_t x, size_t y, size_t z = 0);
    T& operator[](size_t i);
    T& operator[](Vector3 pos);

    Vector3 getDimensions() const;
    int width() const;
    int depth() const;
    int height() const;
    int getIndex(size_t x, size_t y, size_t z) const;
    int getIndex(Vector3 coord) const;
    std::tuple<size_t, size_t, size_t> getCoord(size_t index) const;
    Vector3 getCoordAsVector3(size_t index) const;
    bool checkCoord(int x, int y, int z = 0) const;
    bool checkCoord(Vector3 pos) const;
    bool checkIndex(size_t i) const;

    T interpolate(Vector3 coord, RETURN_VALUE_ON_OUTSIDE padding = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE) const;
    T interpolate(float x, float y, float z = 0, RETURN_VALUE_ON_OUTSIDE padding = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE) const;
    Matrix3<T>& addValueAt(T value, Vector3 coord);
    Matrix3<T>& addValueAt(T value, float x, float y, float z = 0.f);

    int getNumberNeighbors(size_t x, size_t y, size_t z, bool using4connect = true) const;
    int getNumberNeighbors(Vector3 pos, bool using4connect = true) const;

    Matrix3<T> resize(float factor, RESIZE_MODE mode = LINEAR) const;
    Matrix3<T> resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode = LINEAR) const;
    Matrix3<T> resize(Vector3 newSize, RESIZE_MODE mode = LINEAR) const;

    Matrix3<T> resizeNearest(float factor) const;
    Matrix3<T> resizeNearest(size_t newX, size_t newY, size_t newZ) const;
    Matrix3<T> resizeNearest(Vector3 newSize) const;

    Matrix3<T> subset(int startX, int endX, int startY, int endY, int startZ = 0, int endZ = -1);
    Matrix3<T> subset(Vector3 start, Vector3 end);
    Matrix3<T>& paste(Matrix3<T>& matrixToPaste, Vector3 upperLeftFrontCorner = Vector3());
    Matrix3<T>& paste(Matrix3<T>& matrixToPaste, int left, int up, int front);
    Matrix3<T>& add(Matrix3<T>& matrixToAdd, Vector3 upperLeftFrontCorner, bool useInterpolation = false);
    Matrix3<T>& add(Matrix3<T> &matrixToAdd, int left, int up, int front, bool useInterpolation = false);
    Matrix3<T> concat(Matrix3<T> matrixToConcat);

    Matrix3<float> toDistanceMap(bool ignoreZlayer = false, bool considerBorders = true);

    Matrix3<T> flip(bool onX, bool onY = false, bool onZ = false);

    template<typename U>
    Matrix3<T> convolution(Matrix3<U>& convMatrix, CONVOLUTION_BORDERS border = ZERO_PAD);

    T min() const;
    T max() const;

    Matrix3<T>& max(Matrix3<T>& otherMatrix, Vector3 upperLeftFrontCorner);
    Matrix3<T>& max(Matrix3<T>& otherMatrix, int left, int up, int front);
    Matrix3<T>& min(Matrix3<T>& otherMatrix, Vector3 upperLeftFrontCorner);
    Matrix3<T>& min(Matrix3<T>& otherMatrix, int left, int up, int front);

    static Matrix3<T> max(Matrix3<T> m1, Matrix3<T> m2);
    static Matrix3<T> min(Matrix3<T> m1, Matrix3<T> m2);

    Matrix3<T> abs() const;
    T sum();

    Matrix3<T> rounded(int precision = 0) const;

    Matrix3<T>& normalize();
    Matrix3<T> normalized() const;

    Matrix3<T> transposeXY();

    Matrix3<Vector3> gradient();
    Matrix3<Vector3> grad();
    Matrix3<float> divergence();
    Matrix3<float> div();
    Matrix3<T> curl();
    Matrix3<T> rot();
    Matrix3<T> laplacian();
    Vector3 gradient(Vector3 position);
    Vector3 gradient(float posX, float posY, float posZ = 0);

    Matrix3<int> skeletonize();
    Matrix3<T> dilate(float t = 1.f);
    Matrix3<T> erode(float t = 1.f);

    T trace();

    Matrix3<T> wrapWith(Matrix3<Vector3> wrapper);
    Matrix3<T> wrapWith(BSpline original, BSpline wrapperCurve);

    Matrix3<T> wrapWithoutInterpolation(Matrix3<Vector3> wrapper);
    Matrix3<T> wrapWithoutInterpolation(BSpline original, BSpline wrapperCurve);

    static Matrix3<float> fbmNoise1D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);
    static Matrix3<Vector3> fbmNoise2D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);
    static Matrix3<Vector3> fbmNoise3D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ = 1);

    static Matrix3<float> gaussian(int sizeOnX, int sizeOnY, int sizeOfZ, float sigma, Vector3 offset = Vector3());
    Matrix3<T> LaplacianOfGaussian(int sizeOnX, int sizeOnY, int sizeOfZ, float sigma);
    Matrix3<T> meanSmooth(int sizeOnX = 3, int sizeOnY = 3, int sizeOfZ = 3, bool ignoreBorders = false);

    Matrix3<T>& insertRow(size_t indexToInsert, int affectedDimension, T newData = T());

    void clear() { this->sizeX = 0; this->sizeY = 0; this->sizeZ = 0; return this->data.clear(); }
    void reset(T newVal = T()) { for (auto& val : data) val = newVal; }

    Matrix3<int> binarize(T limitValue = T(), bool greaterValuesAreSetToOne = true, bool useAlsoTheEqualSign = false);
    Matrix3<int> binarizeBetween(T minValue, T maxValue, bool insideValuesAreSetToOne = true, bool useAlsoTheEqualSign = false);

    template<typename U>
    operator Matrix3<U>() {
        Matrix3<U> returned(this->getDimensions());
        for (size_t i = 0; i < this->size(); i++)
            returned[i] = (U)((*this)[i]);
        return returned;
    }


    static Matrix3<T> random(size_t sizeX, size_t sizeY, size_t sizeZ = 1);
    static Matrix3<T> identity(size_t sizeX, size_t sizeY, size_t sizeZ = 1);

    static Matrix3<Vector3> fromImageRGB(std::string filename);
    static Matrix3<float> fromImageBW(std::string filename);

    Matrix3<T> operator-() const;
    template<typename U>
    Matrix3<T>& operator+=(const Matrix3<U>& o);
    template<typename U>
    Matrix3<T>& operator-=(const Matrix3<U>& o);
    template<typename U>
    Matrix3<T>& operator*=(Matrix3<U>& o);
    template<typename U>
    Matrix3<T>& operator/=(Matrix3<U>& o);
    template<typename U>
    Matrix3<T>& operator*=(U o);
    template<typename U>
    Matrix3<T>& operator/=(U o);
    template<typename U>
    Matrix3<T>& operator+=(U o);
    template<typename U>
    Matrix3<T>& operator-=(U o);
    template<typename U>
    bool operator==(Matrix3<U> o);

    std::string toString() const {return "Matrix3 (" + std::to_string(this->sizeX) + "x" + std::to_string(this->sizeY) + "x" + std::to_string(this->sizeZ) + ")"; }

    std::vector<T> data;
    int sizeX = 0, sizeY = 0, sizeZ = 1;

    bool raiseErrorOnBadCoord = false; // THIS SHOULD CLEARLY BE SET TO TRUE, BUT F**K IT!
    T defaultValueOnBadCoord = T();
    RETURN_VALUE_ON_OUTSIDE returned_value_on_outside = DEFAULT_VALUE;
    bool stillRaiseErrorForX = false;
    bool stillRaiseErrorForY = false;
    bool stillRaiseErrorForZ = false;

    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    std::size_t size() const { return end() - begin(); }
    bool empty() const { return begin() == end(); }

    Vector3 getMirrorPosition(Vector3 pos) const;
    Vector3 getWrappedPosition(Vector3 pos) const;
    Vector3 getRepeatPosition(Vector3 pos) const;

    Matrix3& init(const std::vector<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ);

    std::string displayValues();
    std::string displayAsPlot(T min = 0.f, T max = 0.f, std::vector<std::string> patterns = {}, std::map<T, std::string> specialCharactersAtValue = {}, T specialCharEpsilon = 1e-5, std::string charForError = "X", std::string separator = "");
};
//template<class T>
Matrix3<float> operator-(const float a, Matrix3<float> b);
//template<class T>
Matrix3<float> operator+(const float a, Matrix3<float> b);
#include <sstream>
template<class T>
std::string Matrix3<T>::displayValues()
{
    std::stringstream out;
    for (int z = 0; z < this->sizeZ; z++) {
        out << "[Z-level = " << z << "] : \n";
        for (int y = 0; y < this->sizeY; y++) {
            for (int x = 0; x < this->sizeX; x++) {
                out << std::setw(2) << at(x, y, z) << "\t";
            }
            out << "\n";
        }
    }
    return out.str();
}

template<class T>
std::string Matrix3<T>::displayAsPlot(T min, T max, std::vector<std::string> patterns, std::map<T, std::string> specialCharactersAtValue, T specialCharEpsilon, std::string charForError, std::string separator)
{
    if (patterns.empty())
        patterns = {"#", "=", "-", "."};

    if (min == 0.f && max == 0.f) {
        min = this->min();
        max = this->max();
    }
    std::stringstream out;
    for (int z = 0; z < this->sizeZ; z++) {
        out << "[Z-level = " << z << "] : \n";
        for (int y = 0; y < this->sizeY; y++) {
            for (int x = 0; x < this->sizeX; x++) {
                T val = this->at(x, y, z);
                bool specialCharUsed = false;
                for (const auto& [targetValue, character] : specialCharactersAtValue) {
                    if (targetValue - specialCharEpsilon <= val && val <= targetValue + specialCharEpsilon) {
                        out << character << separator;
                        specialCharUsed = true;
                        break; // Don't add other characters
                    }
                }
                if (!specialCharUsed) {
                    float prop = interpolation::linear(val, min, max);
                    if (prop < 0.f || prop > 1.f || isnan(prop))
                        out << charForError << separator;
                    else
                        out << patterns[int(prop * float(patterns.size() - 1))] << separator;
                }
            }
            out << "\n";
        }
    }
    return out.str();
}



template<class T>
Matrix3<T>::Matrix3()
{
}
template<class T>
Matrix3<T>::Matrix3(Vector3 size, T initValue) : Matrix3<T>(size.x, size.y, size.z, initValue)
{
}
template<class T>
Matrix3<T>::Matrix3(size_t sizeX, size_t sizeY, size_t sizeZ, T initValue)
{
    std::vector<T> data(sizeX * sizeY * sizeZ, initValue);
    init(data, sizeX, sizeY, sizeZ);
}
template<class T>
Matrix3<T>::Matrix3(std::vector<T> data, size_t sizeX, size_t sizeY, int sizeZ)
{
    if (sizeZ == -1) {
        sizeZ = int(data.size()) / (sizeX * sizeY);
    }
    init(data, sizeX, sizeY, sizeZ);
}
template<class T>
Matrix3<T>::Matrix3(std::vector<std::vector<T>> data)
{
    std::vector<T> oneMatrix;
    for (std::vector<T>& row : data)
        oneMatrix.insert(oneMatrix.end(), row.begin(), row.end());
    int sizeX = data[0].size();
    int sizeY = data.size();
    init(oneMatrix, sizeX, sizeY, 1);
}
template<class T>
Matrix3<T>::Matrix3(std::vector<std::vector<std::vector<T>>> data)
{
    std::vector<T> oneMatrix;
    for (std::vector<std::vector<T>>& grid : data)
        for (std::vector<T>& row : grid)
            oneMatrix.insert(oneMatrix.end(), row.begin(), row.end());
    int sizeX = data[0].size();
    int sizeY = data.size();
    int sizeZ = (sizeY > 0 ? data[0][0].size() : 0);
    init(oneMatrix, sizeX, sizeY, sizeZ);
}

template<class T>
bool Matrix3<T>::checkCoord(int x, int y, int z) const
{
    return ((0 <= x && x < sizeX) && (0 <= y && y < sizeY) && (0 <= z && z < sizeZ));
}

template<class T>
bool Matrix3<T>::checkCoord(Vector3 pos) const
{
    return checkCoord(pos.x, pos.y, pos.z);
}

template<class T>
bool Matrix3<T>::checkIndex(size_t i) const
{
    return (0 <= i && i < sizeX * sizeY * sizeZ);
}

template<class T>
T Matrix3<T>::interpolate(Vector3 coord, RETURN_VALUE_ON_OUTSIDE padding) const
{
    Vector3 round = coord.floor();
    Vector3 cellOffset = coord - round;

    bool previousErrorConfig = this->raiseErrorOnBadCoord;
    RETURN_VALUE_ON_OUTSIDE previousOutsideConfig = this->returned_value_on_outside;
//    this->raiseErrorOnBadCoord = false;
//    this->returned_value_on_outside = padding;
    T f000 = this->at(round + Vector3(0, 0, 0));
    T f100 = this->at(round + Vector3(1, 0, 0));
    T f010 = this->at(round + Vector3(0, 1, 0));
    T f110 = this->at(round + Vector3(1, 1, 0));
    T f001 = this->at(round + Vector3(0, 0, 1));
    T f101 = this->at(round + Vector3(1, 0, 1));
    T f011 = this->at(round + Vector3(0, 1, 1));
    T f111 = this->at(round + Vector3(1, 1, 1));
    // Interpolation
    T interpol = ((
                              f000 * (1-cellOffset.x) + f100 * cellOffset.x) * (1-cellOffset.y) + (
                              f010 * (1-cellOffset.x) + f110 * cellOffset.x) * cellOffset.y) * (1 - cellOffset.z) +
                        ((
                             f001 * (1-cellOffset.x) + f101 * cellOffset.x) * (1-cellOffset.y) + (
                             f011 * (1-cellOffset.x) + f111 * cellOffset.x) * cellOffset.y) * cellOffset.z;
//    this->raiseErrorOnBadCoord = previousErrorConfig;
//    this->returned_value_on_outside = previousOutsideConfig;
    return interpol;
}

template<class T>
T Matrix3<T>::interpolate(float x, float y, float z, RETURN_VALUE_ON_OUTSIDE padding) const
{
    return interpolate(Vector3(x, y, z), padding);
}

template<class T>
const T &Matrix3<T>::at(Vector3 pos) const
{
    return this->at(pos.x, pos.y, pos.z);
}
template<class T>
const T &Matrix3<T>::at(int i, int j, int k) const
{
    if (checkCoord(i, j, k)) {
        int index = getIndex(i, j, k);
        return this->data[index];
    }
    bool raiseError = raiseErrorOnBadCoord;
    Vector3 newPos(i, j, k);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (newPos.x < 0 || sizeX <= int(newPos.x)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (newPos.y < 0 || sizeY <= int(newPos.y)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (newPos.z < 0 || sizeZ <= int(newPos.z)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);

    }
    if (!raiseError)
        return this->at(newPos);
    else
        throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template<class T>
const T &Matrix3<T>::at(size_t i) const
{
    if (i >= 0 && i < sizeX * sizeY * sizeZ) {
        return this->data[i];
    }
    int x, y, z;
    std::tie(x, y, z) = this->getCoord(i);

    bool raiseError = raiseErrorOnBadCoord;
    Vector3 newPos(x, y, z);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (newPos.x < 0 || sizeX <= int(newPos.x)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (newPos.y < 0 || sizeY <= int(newPos.y)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (newPos.z < 0 || sizeZ <= int(newPos.z)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);
    }
    if (!raiseError)
        return this->at(newPos);
    else
        throw std::out_of_range("Trying to access index " + std::to_string(i) + " (coord " + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}
template<typename T>
const T& Matrix3<T>::operator()(size_t x, size_t y, size_t z) const {
    return this->at(x, y, z);
}
template<typename T>
const T& Matrix3<T>::operator[](size_t i) const {
    return this->at(i);
}
template<typename T>
const T& Matrix3<T>::operator[](Vector3 pos) const {
    return this->at(pos);
}

template<class T>
T &Matrix3<T>::at(Vector3 pos)
{
    return this->at(pos.x, pos.y, pos.z);
}
template<class T>
T &Matrix3<T>::at(int i, int j, int k)
{
    if (checkCoord(i, j, k)) {
        int index = getIndex(i, j, k);
        return this->data[index];
    }
    bool raiseError = raiseErrorOnBadCoord;
    Vector3 newPos(i, j, k);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (newPos.x < 0 || sizeX <= int(newPos.x)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (newPos.y < 0 || sizeY <= int(newPos.y)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (newPos.z < 0 || sizeZ <= int(newPos.z)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);

    }
    if (!raiseError)
        return this->at(newPos);
    else
        throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template<class T>
T &Matrix3<T>::at(size_t i)
{
    if (i >= 0 && i < sizeX * sizeY * sizeZ) {
        return this->data[i];
    }
    int x, y, z;
    std::tie(x, y, z) = this->getCoord(i);

    bool raiseError = raiseErrorOnBadCoord;
    Vector3 newPos(x, y, z);
    if (!raiseErrorOnBadCoord) {
        if (returned_value_on_outside == DEFAULT_VALUE)
            return defaultValueOnBadCoord;

        if (stillRaiseErrorForX && (newPos.x < 0 || sizeX <= int(newPos.x)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForY && (newPos.y < 0 || sizeY <= int(newPos.y)))
            return defaultValueOnBadCoord;
        if (stillRaiseErrorForZ && (newPos.z < 0 || sizeZ <= int(newPos.z)))
            return defaultValueOnBadCoord;

        if (returned_value_on_outside == MIRROR_VALUE)
            newPos = getMirrorPosition(newPos);
        else if (returned_value_on_outside == WRAPPED_VALUE)
            newPos = getWrappedPosition(newPos);
        else if (returned_value_on_outside == REPEAT_VALUE)
            newPos = getRepeatPosition(newPos);
    }
    if (!raiseError)
        return this->at(newPos);
    else
        throw std::out_of_range("Trying to access index " + std::to_string(i) + " (coord " + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ") on matrix of size "
            + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}
template<typename T>
T& Matrix3<T>::operator()(size_t x, size_t y, size_t z) {
    return this->at(x, y, z);
}
template<typename T>
T& Matrix3<T>::operator[](size_t i) {
    return this->at(i);
}
template<typename T>
T& Matrix3<T>::operator[](Vector3 pos) {
    return this->at(pos);
}

template<class T>
int Matrix3<T>::getIndex(size_t x, size_t y, size_t z) const
{
    return z * (this->sizeX * this->sizeY) + y * (this->sizeX) + x;
}
template<class T>
int Matrix3<T>::getIndex(Vector3 coord) const
{
    return this->getIndex(int(coord.x), int(coord.y), int(coord.z));
}

template<class T>
std::tuple<size_t, size_t, size_t> Matrix3<T>::getCoord(size_t index) const
{
    int z = index / (this->sizeX * this->sizeY);
    int y = (index % (this->sizeX * this->sizeY)) / this->sizeX;
    int x = index % this->sizeX;
    return std::make_tuple(x, y, z);
}
template<class T>
Vector3 Matrix3<T>::getCoordAsVector3(size_t index) const
{
    int z = index / (this->sizeX * this->sizeY);
    int y = (index % (this->sizeX * this->sizeY)) / this->sizeX;
    int x = index % this->sizeX;
    return Vector3(x, y, z);
}

template <class T>
Matrix3<T>& Matrix3<T>::addValueAt(T value, Vector3 coord) {
    bool previousError = this->raiseErrorOnBadCoord;

    this->raiseErrorOnBadCoord = false;

    Vector3 floorPos = coord.floor();
    Vector3 offset = coord - floorPos;

    this->at(floorPos + Vector3(0, 0, 0)) += value * (1 - offset.x) * (1 - offset.y) * (1 - offset.z);
    this->at(floorPos + Vector3(0, 0, 1)) += value * (1 - offset.x) * (1 - offset.y) * (    offset.z);
    this->at(floorPos + Vector3(0, 1, 0)) += value * (1 - offset.x) * (    offset.y) * (1 - offset.z);
    this->at(floorPos + Vector3(0, 1, 1)) += value * (1 - offset.x) * (    offset.y) * (    offset.z);
    this->at(floorPos + Vector3(1, 0, 0)) += value * (    offset.x) * (1 - offset.y) * (1 - offset.z);
    this->at(floorPos + Vector3(1, 0, 1)) += value * (    offset.x) * (1 - offset.y) * (    offset.z);
    this->at(floorPos + Vector3(1, 1, 0)) += value * (    offset.x) * (    offset.y) * (1 - offset.z);
    this->at(floorPos + Vector3(1, 1, 1)) += value * (    offset.x) * (    offset.y) * (    offset.z);

    this->raiseErrorOnBadCoord = previousError;
    return *this;
}
template <class T>
Matrix3<T>& Matrix3<T>::addValueAt(T value, float x, float y, float z) {
    return this->addValueAt(value, Vector3(x, y, z));
}

template<class T>
Vector3 Matrix3<T>::gradient(Vector3 position)
{
    this->raiseErrorOnBadCoord = false;
    Vector3 flooredPos = position.floor();
//    Vector3 offset = position - flooredPos;

    return Vector3(
                (at(flooredPos + Vector3(1, 0, 0)) - at(flooredPos)), // * offset.x,
                (at(flooredPos + Vector3(0, 1, 0)) - at(flooredPos)), // * offset.y,
                (at(flooredPos + Vector3(0, 0, 1)) - at(flooredPos)) // * offset.z
                );
    /*
    float NWB = at(flooredPos + Vector3(0, 0, 0));
    float NEB = at(flooredPos + Vector3(1, 0, 0));
    float SWB = at(flooredPos + Vector3(0, 1, 0));
    float SEB = at(flooredPos + Vector3(1, 1, 0));
    float NWT = at(flooredPos + Vector3(0, 0, 1));
    float NET = at(flooredPos + Vector3(1, 0, 1));
    float SWT = at(flooredPos + Vector3(0, 1, 1));
    float SET = at(flooredPos + Vector3(1, 1, 1));

    return Vector3(
                ()
                );
    return Vector3(
                at(flooredPos + Vector3(1, 0, 0)) * (1 - offset.x) + at(flooredPos) * offset.x,
                at(flooredPos + Vector3(0, 1, 0)) * (1 - offset.y) + at(flooredPos) * offset.y,
                at(flooredPos + Vector3(0, 0, 1)) * (1 - offset.z) + at(flooredPos) * offset.z
                );*/
}

template<class T>
Vector3 Matrix3<T>::gradient(float posX, float posY, float posZ)
{
    return gradient(Vector3(posX, posY, posZ));
}

template<class T>
Matrix3<T> Matrix3<T>::dilate(float t)
{
    Matrix3<T> res = *this;
    while (t > 0.f) {
        Matrix3<T> copy = res;
        copy.raiseErrorOnBadCoord = false;
        copy.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        float dt = (t < 1.f ? t : 1.f);
        for (int x = 0; x < res.sizeX; x++) {
            for (int y = 0; y < res.sizeY; y++) {
                for (int z = 0; z < res.sizeZ; z++) {
                    T maxVal = res.at(x, y, z);
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                maxVal = std::max(maxVal, copy.at(x + dx, y + dy, z + dz));
                            }
                        }
                    }
                    res.at(x, y, z) = (1 - dt) * res.at(x, y, z) + t * maxVal;
                }
            }
        }
        t -= 1.f;
    }
    return res;
}

template<class T>
Matrix3<T> Matrix3<T>::erode(float t)
{
    Matrix3<T> res = *this;
    while (t > 0.f) {
        Matrix3<T> copy = res;
        copy.raiseErrorOnBadCoord = false;
        copy.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        float dt = (t < 1.f ? t : 1.f);
        for (int x = 0; x < res.sizeX; x++) {
            for (int y = 0; y < res.sizeY; y++) {
                for (int z = 0; z < res.sizeZ; z++) {
                    T minVal = res.at(x, y, z);
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                minVal = std::min(minVal, copy.at(x + dx, y + dy, z + dz));
                            }
                        }
                    }
                    res.at(x, y, z) = (1 - dt) * res.at(x, y, z) + t * minVal;;
                }
            }
        }
        t -= 1.f;
    }
    return res;
}

template<class T>
T Matrix3<T>::trace()
{
    if (sizeZ != 1)
        throw std::domain_error("Cannot compute the trace of the matrix : Matrix should be 2D (sizeZ = 1) but here sizeZ = " + std::to_string(this->sizeZ));
    if (sizeX != sizeY)
        throw std::domain_error("Cannot compute the trace of the matrix : Matrix should be squared (sizeX = sizeY) but here sizeX = " + std::to_string(this->sizeX) + " and sizeY = " + std::to_string(this->sizeY));
    T sum = T();
    for (int i = 0; i < sizeX; i++)
        sum += this->at(i, i);
    return sum;
}

template<class T>
Matrix3<T>& Matrix3<T>::init(const std::vector<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ)
{
    this->data = data;
    this->sizeX = sizeX;
    this->sizeY = sizeY;
    this->sizeZ = sizeZ;
    return *this;
}

template<class T>
std::ostream& operator<<(std::ostream& io, const Matrix3<T>& v) {
    io << v.toString();
    return io;
}

template<class T>
std::ostream& operator<<(std::ostream& io, std::shared_ptr<Matrix3<T>> v) {
    io << v->toString();
    return io;
}

template<class T>
T Matrix3<T>::min() const
{
    T min = std::numeric_limits<T>::max();
    for(const T& val : this->data)
        min = std::min(min, val);
    return min;
}
template<class T>
T Matrix3<T>::max() const
{
    T max = std::numeric_limits<T>::min();
    for(const T& val : this->data)
        max = std::max(max, val);
    return max;
}

template <class T>
Matrix3<T> Matrix3<T>::abs() const
{
    Matrix3<T> m = *this;
    for (T& val : m)
        val = std::abs(val);
    return m;
}

template <class T>
T Matrix3<T>::sum()
{
    T sum = T();
    for (T& val : this->data)
        sum += val;
    return sum;
}

template<class T>
Matrix3<T> Matrix3<T>::rounded(int precision) const
{
    for(T& val : this->data)
        val = (val * std::pow(10, precision)) / (int)std::pow(10, precision);
}

template<class T>
Matrix3<T>& Matrix3<T>::normalize() {
    if (this->data.empty()) return *this;
    T min = data[0], max = data[0];
    for (const T& val : data)
    {
        if (min > val) min = val;
        if (max < val) max = val;
    }
    if (min != max) {
        for (T& val : data)
        {
            val = (val - min)/ (max - min);
        }
    } else {
        for (T& val : data)
        {
            val = T();
        }
    }
    return *this;
}
template<class T>
Matrix3<T> Matrix3<T>::normalized() const {
    Matrix3 mat = *this;
    return mat.normalize();
}

template<class T>
Matrix3<T> Matrix3<T>::transposeXY()
{
    Matrix3<T> res(this->getDimensions().yxz());
    for (int x = 0; x < res.sizeX; x++) {
        for (int y = 0; y < res.sizeY; y++) {
            for (int z = 0; z < res.sizeZ; z++) {
                T prevVal = this->at(y, x, z);
                res.at(x, y, z) = prevVal;
            }
        }
    }
    return res;
}

template<class T>
Matrix3<float> Matrix3<T>::gaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma, Vector3 offset) {
    Matrix3<float> gaussian(sizeOnX, sizeOnY, sizeOnZ);
    Vector3 center = Vector3(sizeOnX/2.f, sizeOnY/2.f, sizeOnZ/2.f) + offset;
    center -= Vector3((sizeOnX > 1 ? .5 : 0), (sizeOnY > 1 ? .5 : 0), (sizeOnZ > 1 ? .5 : 0));
    float oneOverSqrt2Pi = 1.f/std::sqrt(2 * 3.141592);
    float sqrSigma = sigma * sigma;
    for (int x = 0; x < gaussian.sizeX; x++) {
        for (int y = 0; y < gaussian.sizeY; y++) {
            for (int z = 0; z < gaussian.sizeZ; z++) {
                Vector3 pos = Vector3(x, y, z) - center;
                gaussian.at(x, y, z) = std::exp(-pos.norm2()/(2*sqrSigma)) * oneOverSqrt2Pi;
            }
        }
    }
    return gaussian;
}

template<class T>
Matrix3<T> Matrix3<T>::LaplacianOfGaussian(int sizeOnX, int sizeOnY, int sizeOnZ, float sigma) {
    Matrix3<float> laplacian = Matrix3<float>(3 + sizeOnX/2 +1, 3 + sizeOnY/2 +1, 3 + sizeOnZ/2 +1, 0.f);
    laplacian.at(sizeOnX/2 + 1, sizeOnY/2 + 1, sizeOnZ/2 + 1) = 1.f;
    laplacian = laplacian.laplacian();
    Matrix3<float> gaussian = Matrix3<T>::gaussian(sizeOnX, sizeOnY, sizeOnZ, sigma);
    return this->convolution(laplacian.convolution(gaussian));
}
template<class T>
Matrix3<T> Matrix3<T>::meanSmooth(int sizeOnX, int sizeOnY, int sizeOfZ, bool ignoreBorders) {
    Matrix3<float> meanMatrix(sizeOnX, sizeOnY, sizeOfZ, 1.f);
    return this->convolution(meanMatrix, (ignoreBorders ? CONVOLUTION_BORDERS::IGNORED : CONVOLUTION_BORDERS::ZERO_PAD));
}

/*
template<class T>
Matrix3<Vector3> Matrix3<T>::gradient() {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3 vec;
                Vector3 allDirections;
                if (x == 0 || x == this->sizeX - 1 || y == 0 || y == this->sizeY - 1
                        || z == 0 || z == this->sizeZ - 1) {
                    returningGrid.at(x, y, z) = vec;
                } else {
                    for (int dx = std::max(x - 1, 0); dx <= std::min(x + 1, this->sizeX - 1); dx++) {
                        for (int dy = std::max(y - 1, 0); dy <= std::min(y + 1, this->sizeY - 1); dy++) {
                            for (int dz = std::max(z - 1, 0); dz <= std::min(z + 1, this->sizeZ - 1); dz++) {
                                if(dx != int(x) || dy != int(y) || dz != int(z)) {
                                    T f = this->at(dx, dy, dz);
                                    vec += Vector3(dx - x, dy - y, dz - z) * f; // * this->at(dx,dy,dz];
                                }
                            }
                        }
                    }
                }

                returningGrid.at(x, y, z) = vec;

            }
        }
    }
    return returningGrid;
}*/
template<class T>
Matrix3<Vector3> Matrix3<T>::gradient() {
    Matrix3<Vector3> returningGrid(this->sizeX, this->sizeY, this->sizeZ);
    bool oldError = this->raiseErrorOnBadCoord;
    this->raiseErrorOnBadCoord = false;
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                returningGrid.at(x, y, z) = Vector3((at(x + 1, y, z) - at(x - 1, y, z)) * .5f,
                                                    (at(x, y + 1, z) - at(x, y - 1, z)) * .5f,
                                                    (at(x, y, z + 1) - at(x, y, z - 1)) * .5f);
            }
        }
    }
    this->raiseErrorOnBadCoord = oldError;
    return returningGrid;
}
template<class T>
Matrix3<Vector3> Matrix3<T>::grad() {return this->gradient(); }

template<class T>
Matrix3<T> Matrix3<T>::laplacian()
{
    Matrix3 returningGrid = *this;
    this->raiseErrorOnBadCoord = false;
    this->defaultValueOnBadCoord = 0;
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                T val = T();
                /*if (x == 0 || x == this->sizeX - 1 || y == 0 || y == this->sizeY - 1
                        || z == 0 || z == this->sizeZ - 1) {
                    returningGrid.at(x, y, z) = val;
                } else {*/
                    val += this->at(x    , y    , z + 1);
                    val += this->at(x    , y    , z - 1);
                    val += this->at(x    , y + 1, z    );
                    val += this->at(x    , y - 1, z    );
                    val += this->at(x + 1, y    , z    );
                    val += this->at(x - 1, y    , z    );
                    val -= this->at(x    , y    , z    ) * 6;
//                }
                returningGrid.at(x, y, z) = val;
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}

template<typename T>
Matrix3<int> Matrix3<T>::binarize(T limitValue, bool greaterValuesAreSetToOne, bool useAlsoTheEqualSign)
{
    Matrix3<int> bin(this->sizeX, this->sizeY, this->sizeZ);
    for (size_t i = 0; i < this->size(); i++) {
        if (greaterValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= limitValue) || (!useAlsoTheEqualSign && this->data[i] > limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= limitValue) || (!useAlsoTheEqualSign && this->data[i] < limitValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    }
    return bin;
}

template<class T>
Matrix3<int> Matrix3<T>::binarizeBetween(T minValue, T maxValue, bool insideValuesAreSetToOne, bool useAlsoTheEqualSign)
{
    Matrix3<int> bin(this->sizeX, this->sizeY, this->sizeZ);
    for (size_t i = 0; i < this->size(); i++) {
        if (insideValuesAreSetToOne) {
            if ((useAlsoTheEqualSign && this->data[i] >= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] > minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        } else {
            if ((useAlsoTheEqualSign && this->data[i] <= minValue && this->data[i] <= maxValue) || (!useAlsoTheEqualSign && this->data[i] < minValue && this->data[i] < maxValue)) {
                bin[i] = 1;
            } else {
                bin[i] = 0;
            }
        }
    }
    return bin;
}

template<typename T>
Matrix3<T>& Matrix3<T>::insertRow(size_t indexToInsert, int affectedDimension, T newData)
{
    auto it = this->data.begin();
    int jumps = 0;
    switch(affectedDimension) {
    case 0:
        jumps = sizeX;
        if (indexToInsert < 0) indexToInsert = sizeX; // If default value, it's last X-index
        it += indexToInsert;
        for (; it <= this->data.end() - (indexToInsert == 0 ? 1 : 0); it += jumps) {
            it = this->data.insert(it, newData) + 1; // Set "it" to the value next to the inserted value
        }
        this->sizeX++;
        break;
    case 1:
        jumps = sizeX * sizeY;
        if (indexToInsert < 0) indexToInsert = sizeY; // If default value, it's last Y-index
        it += (indexToInsert * sizeX);
        for (; it <= this->data.end() - (indexToInsert == 0 ? 1 : 0); it += jumps) {
            for (int i = 0; i < sizeX; i++) {
                it = this->data.insert(it, newData) + 1;
            }
        }
        this->sizeY++;
        break;
    case 2:
        if (indexToInsert < 0) indexToInsert = sizeZ; // If default value, it's last Z-index
        it += (sizeX * sizeY) * indexToInsert;
        for (int i = 0; i < sizeX * sizeY; i++)
            it = this->data.insert(it, newData) + 1;
        this->sizeZ++;
        break;
    default:
        throw std::out_of_range("insertRow can only be processed on dim 0, 1 or 2 (resp. X, Y, Z)");
    }
    return *this;
}

template<typename T>
Matrix3<T> Matrix3<T>::identity(size_t sizeX, size_t sizeY, size_t sizeZ)
{
    static_assert(std::is_arithmetic<T>::value, "");
    Matrix3<T> mat(sizeX, sizeY, sizeZ);
    if (sizeZ == 1) {
        if (sizeX != sizeY)
            throw std::invalid_argument("Identity matrix must be square (dim X == dim Y)");
        for (int i = 0; i < sizeX; i++)
            mat.at(i, i, 0) = 1;
    } else {
        if (sizeX != sizeY || sizeX != sizeZ)
            throw std::invalid_argument("All dimensions must be the same length to do a 3-d identity matrix");
        for (int i = 0; i < sizeX; i++)
            mat.at(i, i, i) = 1;
    }
    return mat;
}

template<typename T>
Matrix3<T> operator+(Matrix3<T> a, Matrix3<T> b) {
    a += b;
    return a;
}
template<typename T>
Matrix3<T> Matrix3<T>::operator-() const {
    return *this * -1.f;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator+=(const Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] += o.data[i];
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator-(Matrix3<T> a, const Matrix3<U>& b) {
    a -= b;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator-=(const Matrix3<U>& o)  {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be substracted (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o.data[i];
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator*(Matrix3<T> a, Matrix3<U> o) {
    a *= o;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator*=(Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be multiplied (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o.data[i];
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator/(Matrix3<T>& a, const Matrix3<U>& b) {
    a /= b;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator/=(Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices must have same sizes to be divided (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o.data[i];
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator*(Matrix3<T> a, U o) {
    a *= o;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator*=(U o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o;
    }
    return *this;
}

template<typename T, typename U>
Matrix3<T> operator/(Matrix3<T> a, U o) {
    a /= o;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator/=(U o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o;
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator+(Matrix3<T> a, U o) {
    a += o;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator+=(U o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] += o;
    }
    return *this;
}
template<typename T, typename U>
Matrix3<T> operator-(Matrix3<T> a, U o) {
    a -= o;
    return a;
}
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator-=(U o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o;
    }
    return *this;
}
template<typename T> template<typename U>
bool Matrix3<T>::operator==(Matrix3<U> o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        return false;
    for (size_t i = 0; i < this->data.size(); i++)
        if (this->data[i] != o.data[i])
            return false;
    return true;
}

template<class T>
int Matrix3<T>::getNumberNeighbors(size_t x, size_t y, size_t z, bool using4connect) const
{
    int neighbors = 0;
    if (using4connect) {
        if (x > 0) neighbors++;
        if (x < sizeX-1) neighbors++;
        if (y > 0) neighbors++;
        if (y < sizeY - 1) neighbors++;
        if (z > 0) neighbors++;
        if (z < sizeZ - 1) neighbors++;
    }
    else {
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                    if (dx != 0 || dy != 0 || dz != 0)
                        neighbors += (checkCoord(x + dx, y + dy, z + dz) ? 1 : 0);
    }
    return neighbors;
}
template<class T>
int Matrix3<T>::getNumberNeighbors(Vector3 pos, bool using4connect) const
{
    return getNumberNeighbors(pos.x, pos.y, pos.z, using4connect);
}

template<typename T>
Matrix3<T> Matrix3<T>::resize(float factor, RESIZE_MODE mode) const
{
    Vector3 newSize = this->getDimensions() * factor;
    if (newSize.x < 1) newSize.x = 1;
    if (newSize.y < 1) newSize.y = 1;
    if (newSize.z < 1) newSize.z = 1;
    return this->resize(newSize, mode);
}
template<typename T>
Matrix3<T> Matrix3<T>::resize(Vector3 newSize, RESIZE_MODE mode) const
{
    return this->resize(newSize.x, newSize.y, newSize.z, mode);
}

template<typename T>
Matrix3<T> Matrix3<T>::resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode) const
{
    Matrix3<T> newMat(newX, newY, newZ);
    newMat.raiseErrorOnBadCoord = false;
    float rx = (this->sizeX - 1) / std::max(1.f, (float)(newX - 1)), ry = (this->sizeY - 1) / std::max(1.f, (float)(newY - 1)), rz = (this->sizeZ - 1) / std::max(1.f, (float)(newZ - 1));

    if (mode == LINEAR) {
        // Apply interpolations
        for (size_t x = 0; x < newX; x++) {
            int x_original = int(x * rx);
            int x_plus_1 = (x_original >= this->sizeX - 1 ? x_original : x_original + 1);
            float d_x = (x * rx) - x_original;
            for (size_t y = 0; y < newY; y++) {
                int y_original = int(y * ry);
                int y_plus_1 = (y_original >= this->sizeY - 1 ? y_original : y_original + 1);
                float d_y = (y * ry) - y_original;
                for (size_t z = 0; z < newZ; z++) {
                    int z_original = int(z * rz);
                    int z_plus_1 = (z_original >= this->sizeZ - 1 ? z_original : z_original + 1);
                    float d_z = (z * rz) - z_original;

                    T f000 = this->at(x_original    , y_original    , z_original    );
                    T f100 = this->at(x_plus_1, y_original    , z_original    );
                    T f010 = this->at(x_original    , y_plus_1, z_original    );
                    T f110 = this->at(x_plus_1, y_plus_1, z_original    );
                    T f001 = this->at(x_original    , y_original    , z_plus_1);
                    T f101 = this->at(x_plus_1, y_original    , z_plus_1);
                    T f011 = this->at(x_original    , y_plus_1, z_plus_1);
                    T f111 = this->at(x_plus_1, y_plus_1, z_plus_1);
                    // Interpolation
                    T res = ((
                                              f000 * (1-d_x) + f100 * d_x) * (1-d_y) + (
                                              f010 * (1-d_x) + f110 * d_x) * d_y) * (1 - d_z) +
                                        ((
                                             f001 * (1-d_x) + f101 * d_x) * (1-d_y) + (
                                             f011 * (1-d_x) + f111 * d_x) * d_y) * d_z;

                    newMat.at(x, y, z) = res;
                }
            }
        }
    } else if (mode == NEAREST) {
        newMat = this->resizeNearest(newX, newY, newZ);

    } else if (mode == MAX_VAL || mode == MIN_VAL) {
//        for (auto& val : newMat)
//            val = (mode == MAX_VAL ? std::numeric_limits<T>::min() : std::numeric_limits<T>::max());
        Matrix3<short int> modifiedMatrix(newX, newY, newZ, 0);
        for (int x = 0; x < this->sizeX; x++) {
            int startX = x / rx;
            int endX = (x + 1) / rx;
            for (int y = 0; y < this->sizeY; y++) {
                int startY = y / ry;
                int endY = (y + 1) / ry;
                for (int z = 0; z < this->sizeZ; z++) {
                    int startZ = z / rz;
                    int endZ = (z + 1) / rz;

                    for (int dx = startX; dx <= endX; dx++) {
                        for (int dy = startY; dy <= endY; dy++) {
                            for (int dz = startZ; dz <= endZ; dz++) {
                                if (modifiedMatrix.checkCoord(dx, dy, dz)) {
                                    // If this cell hasn't been modified yet, we cannot apply the min/max operator
                                    if (modifiedMatrix.at(dx, dy, dz) == 0) {
                                        newMat.at(dx, dy, dz) = this->at(x, y, z);
                                        modifiedMatrix.at(dx, dy, dz) = 1;
                                    } else {
                                        newMat.at(dx, dy, dz) = (mode == MAX_VAL ? std::max(newMat.at(dx, dy, dz), this->at(x, y, z)) : std::min(newMat.at(dx, dy, dz), this->at(x, y, z)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else if (mode == FILL_WITH_DEFAULT) {
        if (rz == 0) rz = 1;
        for (int x = 0; x < this->sizeX; x++) {
            for (int y = 0; y < this->sizeY; y++) {
                for (int z = 0; z < this->sizeZ; z++) {
                    T old = this->at(x, y, z);
                    Vector3 pos(x / rx, y / ry, z / rz);
                    newMat.at(pos) = old;
                }
            }
        }
    }
    newMat.raiseErrorOnBadCoord = this->raiseErrorOnBadCoord;
    return newMat;
}

template<class T>
Matrix3<T> Matrix3<T>::resizeNearest(float factor) const
{
    Vector3 newSize = this->getDimensions() * factor;
    if (newSize.x < 1) newSize.x = 1;
    if (newSize.y < 1) newSize.y = 1;
    if (newSize.z < 1) newSize.z = 1;
    return this->resizeNearest(newSize);
}

template<class T>
Matrix3<T> Matrix3<T>::resizeNearest(size_t newX, size_t newY, size_t newZ) const
{
    Matrix3<T> newMat(newX, newY, newZ);
    newMat.raiseErrorOnBadCoord = false;
    float rx = (this->sizeX - 1) / std::max(1.f, (float)(newX - 1)), ry = (this->sizeY - 1) / std::max(1.f, (float)(newY - 1)), rz = (this->sizeZ - 1) / std::max(1.f, (float)(newZ - 1));

    // Apply interpolations
    for (int x = 0; x < newX; x++) {
        int x_rounded = std::round(x * rx);
        for (int y = 0; y < newY; y++) {
            int y_rounded = std::round(y * ry);
            for (int z = 0; z < newZ; z++) {
                int z_rounded = std::round(z * rz);
                T res = this->at(x_rounded, y_rounded, z_rounded);
                newMat.at(x, y, z) = res;
            }
        }
    }
    return newMat;
}

template<class T>
Matrix3<T> Matrix3<T>::resizeNearest(Vector3 newSize) const
{
    return this->resizeNearest(newSize.x, newSize.y, newSize.z);
}


template<typename T>
Matrix3<T> Matrix3<T>::subset(Vector3 start, Vector3 end)
{
    if (start.z == 0 && end.z == 0)
        end.z = -1; // Give it the default value so it will be managed by the main function
    return this->subset(start.x, end.x, start.y, end.y, start.z, end.z);
}

template<typename T>
Matrix3<T> Matrix3<T>::subset(int startX, int endX, int startY, int endY, int startZ, int endZ)
{
    if (endZ == -1) endZ = this->sizeZ;
    Matrix3<T> croppedMatrix(std::max(endX - startX, 0), std::max(endY - startY, 0), std::max(endZ - startZ, 0));
    for (int x = startX; x < endX; x++) {
        if (x < 0 || this->sizeX <= x) continue;
        for (int y = startY; y < endY; y++) {
            if (y < 0 || this->sizeY <= y) continue;
            for (int z = startZ; z < endZ; z++) {
                if (z < 0 || this->sizeZ <= z) continue;
                croppedMatrix.at(x - startX, y - startY, z - startZ) = this->at(x, y, z);
            }
        }
    }
    return croppedMatrix;
}


template<typename T>
Matrix3<T>& Matrix3<T>::paste(Matrix3<T>& matrixToPaste, Vector3 upperLeftFrontCorner)
{
    return this->paste(matrixToPaste, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z);
}
template<typename T>
Matrix3<T>& Matrix3<T>::paste(Matrix3<T>& matrixToPaste, int left, int up, int front)
{
    for (int x = std::max(left, 0); x < std::min(matrixToPaste.sizeX + left, this->sizeX); x++) {
        for (int y = std::max(up, 0); y < std::min(matrixToPaste.sizeY + up, this->sizeY); y++) {
            for (int z = std::max(front, 0); z < std::min(matrixToPaste.sizeZ + front, this->sizeZ); z++) {
                this->at(x, y, z) = matrixToPaste.at(x - left, y - up, z - front);
            }
        }
    }
    return *this;
}

template<typename T>
Matrix3<T>& Matrix3<T>::add(Matrix3<T>& matrixToAdd, Vector3 upperLeftFrontCorner, bool useInterpolation)
{
    if (useInterpolation) {
        for (int x = 0; x < matrixToAdd.sizeX; x++) {
            for (int y = 0; y < matrixToAdd.sizeY; y++) {
                for (int z = 0; z < matrixToAdd.sizeZ; z++) {
    //                T& val = matrixToAdd.at(x - left, y - up, z - front);
    //                this->at(x, y, z) += val;
                    T& val = matrixToAdd.at(x, y, z);
                    this->addValueAt(val, upperLeftFrontCorner + Vector3(x, y, z));
                }
            }
        }
        return *this;
    } else {
        return this->add(matrixToAdd, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z, useInterpolation);
    }
}
template<typename T>
Matrix3<T>& Matrix3<T>::add(Matrix3<T>& matrixToAdd, int left, int up, int front, bool useInterpolation)
{
    for (int x = std::max(left, 0); x < std::min(matrixToAdd.sizeX + left, this->sizeX); x++) {
        for (int y = std::max(up, 0); y < std::min(matrixToAdd.sizeY + up, this->sizeY); y++) {
            for (int z = std::max(front, 0); z < std::min(matrixToAdd.sizeZ + front, this->sizeZ); z++) {
                T& val = matrixToAdd.at(x - left, y - up, z - front);
                this->at(x, y, z) += val;
            }
        }
    }
    return *this;
}

template<class T>
Matrix3<T> Matrix3<T>::concat(Matrix3<T> matrixToConcat)
{
    Matrix3<T> newMatrix(this->getDimensions() + matrixToConcat.getDimensions() * Vector3(1, 0, 0));
    newMatrix.paste(*this, Vector3());
    newMatrix.paste(matrixToConcat, (newMatrix.getDimensions() - matrixToConcat.getDimensions()) * Vector3(1, 0, 0));
    return newMatrix;
}

template<typename T>
Matrix3<T>& Matrix3<T>::max(Matrix3<T>& otherMatrix, Vector3 upperLeftFrontCorner)
{
    return this->max(otherMatrix, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z);
}
template<typename T>
Matrix3<T>& Matrix3<T>::max(Matrix3<T>& otherMatrix, int left, int up, int front)
{
    for (int x = std::max(left, 0); x < std::min(otherMatrix.sizeX + left, this->sizeX); x++) {
        for (int y = std::max(up, 0); y < std::min(otherMatrix.sizeY + up, this->sizeY); y++) {
            for (int z = std::max(front, 0); z < std::min(otherMatrix.sizeZ + front, this->sizeZ); z++) {
                this->at(x, y, z) = std::max(this->at(x, y, z), otherMatrix.at(x - left, y - up, z - front));
            }
        }
    }
    return *this;
}

template<typename T>
Matrix3<T>& Matrix3<T>::min(Matrix3<T>& otherMatrix, Vector3 upperLeftFrontCorner)
{
    return this->min(otherMatrix, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z);
}
template<typename T>
Matrix3<T>& Matrix3<T>::min(Matrix3<T>& otherMatrix, int left, int up, int front)
{
    for (int x = std::max(left, 0); x < std::min(otherMatrix.sizeX + left, this->sizeX); x++) {
        for (int y = std::max(up, 0); y < std::min(otherMatrix.sizeY + up, this->sizeY); y++) {
            for (int z = std::max(front, 0); z < std::min(otherMatrix.sizeZ + front, this->sizeZ); z++) {
                this->at(x, y, z) = std::min(this->at(x, y, z), otherMatrix.at(x - left, y - up, z - front));
            }
        }
    }
    return *this;
}

template<class T>
Matrix3<T> Matrix3<T>::max(Matrix3<T> m1, Matrix3<T> m2)
{
    if (m1.getDimensions() != m2.getDimensions())
        throw std::domain_error("Matrices must have same sizes to be maxed (M1 = " + m1.toString() + " and M2 = " + m2.toString());
    Matrix3<T> res(m1.getDimensions());
    for (size_t i = 0; i < m1.size(); i++) {
        res[i] = std::max(m1[i], m2[i]);
    }
    return res;
}

template<class T>
Matrix3<T> Matrix3<T>::min(Matrix3<T> m1, Matrix3<T> m2)
{
    if (m1.getDimensions() != m2.getDimensions())
        throw std::domain_error("Matrices must have same sizes to be mined (M1 = " + m1.toString() + " and M2 = " + m2.toString());
    Matrix3<T> res(m1.getDimensions());
    for (size_t i = 0; i < m1.size(); i++) {
        res[i] = std::min(m1[i], m2[i]);
    }
    return res;
}

template<class T>
Matrix3<float> Matrix3<T>::toDistanceMap(bool ignoreZlayer, bool considerBorders)
{
    Matrix3<float> distances(this->sizeX, this->sizeY, this->sizeZ, std::numeric_limits<float>::max() - 10000);
    distances.raiseErrorOnBadCoord = false;
    distances.defaultValueOnBadCoord = (considerBorders ? 0.f : std::numeric_limits<float>::max() - 10000);

    // Using the Chamfer distance -> direct neighbor               => distance = 3
    //                               diagonal on 2 axis neighbor   => distance = 4
    //                               diagonal on all axis neighbor => distance = 5
//    float predefinedDistances[4] = {0, 3, 4, 5};
    // First pass
    for (int x = 0; x < distances.sizeX; x++) {
        for (int y = 0; y < distances.sizeY; y++) {
            for (int z = 0; z < distances.sizeZ; z++) {
                float currentVal = distances.at(x, y, z);
                if (!this->at(x, y, z)) {
                    distances.at(x, y, z) = 0;
                    continue;
                }
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            if (ignoreZlayer && dz != 0) continue;
                            // Weighted distance transform
//                            currentVal = std::min(currentVal, distances.at(dx, dy, dz) + predefinedDistances[std::abs(dx) + std::abs(dy) + std::abs(dz)]);
                            currentVal = std::min(currentVal, distances.at(x+dx, y+dy, z+dz) + (float)std::sqrt(dx*dx + dy*dy + dz*dz));
                        }
                    }
                }
                distances.at(x, y, z) = currentVal;
            }
        }
    }
    // Second pass
    for (int x = distances.sizeX-1; x >= 0; x--) {
        for (int y = distances.sizeY-1; y >= 0; y--) {
            for (int z = distances.sizeZ-1; z >= 0; z--) {
                if (!this->at(x, y, z)) {
                    distances.at(x, y, z) = 0;
                    continue;
                }
                float currentVal = distances.at(x, y, z);
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            if (ignoreZlayer && dz != 0) continue;
                            currentVal = std::min(currentVal, distances.at(x+dx, y+dy, z+dz) + (float)std::sqrt(dx*dx + dy*dy + dz*dz));
                        }
                    }
                }
                distances.at(x, y, z) = currentVal;
            }
        }
    }
//    distances /= 3.f; // We used a weighted distance, go back to normal
    /*for (auto& d : distances) {
        d = std::sqrt(d); // We used an Exact Euclidean Distance, go back to normal
    }*/
    distances.raiseErrorOnBadCoord = true;
    distances.defaultValueOnBadCoord = 0.f;
    return distances; //.normalize();
}

template<class T>
Matrix3<T> Matrix3<T>::flip(bool onX, bool onY, bool onZ)
{
    Matrix3<T> result = *this;
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                int targetX = (onX ? this->sizeX - (x +1) : x);
                int targetY = (onY ? this->sizeY - (y +1) : y);
                int targetZ = (onZ ? this->sizeZ - (z +1) : z);
                result.at(x, y, z) = this->at(targetX, targetY, targetZ);
            }
        }
    }
    return result;
}


template<class T> template<class U>
Matrix3<T> Matrix3<T>::convolution(Matrix3<U>& convMatrix, CONVOLUTION_BORDERS borders)
{
    Matrix3<T> result(this->sizeX, this->sizeY, this->sizeZ);
    this->raiseErrorOnBadCoord = false;
//    this->defaultValueOnBadCoord = 0;

    for (int x = 0; x < result.sizeX; x++) {
        for (int y = 0; y < result.sizeY; y++) {
            for (int z = 0; z < result.sizeZ; z++) {
                Vector3 pos(x, y, z);
                T neighboringSum = T();
                T normalisationValue = T();
                for (int dx = 0; dx < convMatrix.sizeX; dx++) {
                    for (int dy = 0; dy < convMatrix.sizeY; dy++) {
                        for (int dz = 0; dz < convMatrix.sizeZ; dz++) {
                            normalisationValue += convMatrix.at(dx, dy, dz);
                            Vector3 dt(dx - (convMatrix.sizeX/2), dy - (convMatrix.sizeY/2), dz - (convMatrix.sizeZ/2));
                            Vector3 cellValuePosition = pos + dt;
                            if (borders == CONVOLUTION_BORDERS::IGNORED && !result.checkCoord(cellValuePosition))
                                continue;
                            else if (borders == CONVOLUTION_BORDERS::MIRROR && !result.checkCoord(cellValuePosition)) {
                                cellValuePosition = getMirrorPosition(cellValuePosition);
                            }
                            else if (borders == CONVOLUTION_BORDERS::REPEAT && !result.checkCoord(cellValuePosition))
                                cellValuePosition = getRepeatPosition(cellValuePosition);
                            else if (borders == CONVOLUTION_BORDERS::WRAPPING && !result.checkCoord(cellValuePosition))
                                cellValuePosition = this->getWrappedPosition(cellValuePosition);
                            else if (borders == CONVOLUTION_BORDERS::COPY && !result.checkCoord(cellValuePosition))
                                cellValuePosition = pos;

                            T cellValue = this->at(cellValuePosition);
                            neighboringSum += cellValue * convMatrix.at(dx, dy, dz);
                        }
                    }

                }
                result.at(pos) = neighboringSum;
                if (normalisationValue != T())
                    result.at(pos) /= normalisationValue;
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return result;
}

template<class T>
Vector3 Matrix3<T>::getMirrorPosition(Vector3 pos)  const
{
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    x = int(x < 0 ? std::abs(x) : (x >= sizeX ? sizeX - (x - sizeX) -1 : x));
    y = int(y < 0 ? std::abs(y) : (y >= sizeY ? sizeY - (y - sizeY) -1 : y));
    z = int(x < 0 ? std::abs(z) : (z >= sizeZ ? sizeZ - (z - sizeZ) -1 : z));
    return Vector3(x, y, z);
}

template<class T>
Vector3 Matrix3<T>::getWrappedPosition(Vector3 pos) const
{
    Vector3 rounded = pos.roundedDown();
    Vector3 decimals = pos - rounded;
    Vector3 wrap = Vector3(int(rounded.x + sizeX) % sizeX,
                           int(rounded.y + sizeY) % sizeY,
                           int(rounded.z + sizeZ) % sizeZ
                           ) + decimals;
    return wrap;
//    return     Vector3(int(pos.x) % sizeX,
//                       int(pos.y) % sizeY,
//                       int(pos.z) % sizeZ);
}

template<class T>
Vector3 Matrix3<T>::getRepeatPosition(Vector3 pos) const
{
    Vector3 returned;
    returned.x = std::min(std::max(0.f, pos.x), (float)sizeX - 1);
    returned.y = std::min(std::max(0.f, pos.y), (float)sizeY - 1);
    returned.z = std::min(std::max(0.f, pos.z), (float)sizeZ - 1);
    return returned;
}

template <class T>
Matrix3<T> Matrix3<T>::wrapWith(Matrix3<Vector3> wrapper)
{
    // Wrap definition : f(wrap(p)) = f(p + wrap vec)
    // But f(p) != f(p - wrap vec), in the definition. I think it should create
    // better results because we can fetch outside values (with mirror for ex)
    // But that's not the current definition.
    Matrix3<T> result(getDimensions());
    this->raiseErrorOnBadCoord = false;
    this->returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    result.raiseErrorOnBadCoord = false;
    Matrix3<float> unit(result.getDimensions(), 1.f);
    Matrix3<float> unit_ctrl = unit;
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                Vector3 pos(x, y, z);
                Vector3& wrap = wrapper.at(pos);
                result.addValueAt(this->at(pos), pos + wrap);
//                result.addValueAt(this->interpolate(pos - wrap), pos);

                unit_ctrl.addValueAt(1.f, pos + wrap);
//                unit_ctrl.addValueAt(1.f, pos);
            }
        }
    }
    for (auto& val : unit_ctrl) val = (val == 0 ? .00001f : val);
    return result /= unit_ctrl;
}

template<class T>
Matrix3<T> Matrix3<T>::wrapWith(BSpline original, BSpline wrapperCurve)
{
    // For now, start from a straight line on the X-axis
//    BSpline original = BSpline({this->getDimensions() * Vector3(0, .5, .5) + Vector3(1, 0, 0), this->getDimensions() * Vector3(1, .5, .5) - Vector3(1, 0, 0)});
    float pathsResolution = 1000.f;
    std::vector<Vector3> originalCurvePoints = original.getPath(pathsResolution);
    std::vector<Vector3> wrapperCurvePoints = wrapperCurve.getPath(pathsResolution);

    Matrix3<Vector3> wrapper(this->getDimensions());
    wrapper.raiseErrorOnBadCoord = false;
    Matrix3<float> modifications(this->getDimensions(), 0.f);
    modifications.raiseErrorOnBadCoord = false;

    // Vectors along the curve
    for (size_t i = 0; i < originalCurvePoints.size(); i++) {
        Vector3 pos = originalCurvePoints[i];
        Vector3 dir = wrapperCurvePoints[i] - pos;
        float curveWrapLength = dir.norm();
        dir.normalize();

        // In direction of the curve wrapping
        Vector3 endingPropagationPoint = Collision::intersectionRayAABBox(pos + Vector3(.5, .5, .5), dir, Vector3(), getDimensions());
        float distanceToBorder = (dir.norm2() > 0 ? (endingPropagationPoint - pos).norm() : 1.f);

        for (int j = 0; j < distanceToBorder; j++) {
            float wrapLength = (1 - interpolation::linear(j / distanceToBorder)) * curveWrapLength;
//            wrapper.at(pos + dir * (float)j) += dir * wrapLength;
//            modifications.at(pos + dir * (float)j) += 1.f;
            wrapper.addValueAt(dir * wrapLength, pos + dir * (float)j);
            modifications.addValueAt(1.f, pos + dir * (float)j);
        }

        // In opposite direction of the curve wrapping
        endingPropagationPoint = Collision::intersectionRayAABBox(pos + Vector3(.5, .5, .5), dir * -1.f, Vector3(), getDimensions());
        distanceToBorder = (dir.norm2() > 0 ? (endingPropagationPoint - pos).norm() : 1.f);

        for (int j = 1; j < distanceToBorder; j++) {
            float wrapLength = (1 - interpolation::linear(j / distanceToBorder)) * curveWrapLength;
//            wrapper.at(pos - dir * (float)j) += dir * wrapLength;
//            modifications.at(pos - dir * (float)j) += 1.f;
            wrapper.addValueAt(dir * wrapLength, pos - dir * (float)j);
            modifications.addValueAt(1.f, pos - dir * (float)j);
        }
    }

    for (size_t i = 0; i < wrapper.size(); i++) {
        wrapper[i] /= (modifications[i] > 0 ? modifications[i] : 1.f);
        modifications[i] = (modifications[i] > 0 ? 1.f : 0.f);
    }

    // Fill the empty gaps left
    while (modifications.min() <= 0) {
        for (size_t i = 0; i < wrapper.size(); i++) {
            if (modifications[i] == 0) {
                Vector3 pos = wrapper.getCoordAsVector3(i);
                Vector3 replaceValue;
                float divisor = 0.f;
                for (int x = -1; x <= 1; x++) {
                    for (int y = -1; y <= 1; y++) {
                        for (int z = -1; z <= 1; z++) {
                            Vector3 offset(x, y, z);
                            if (modifications.at(pos + offset) > 0) {
                                replaceValue += wrapper.at(pos + offset);
                                divisor++;
                            }
                        }
                    }
                }
                if (divisor > 0.f) {
                    wrapper.at(pos) = replaceValue / divisor;
                    modifications.at(pos) = 1.f;
                }

            }
        }
    }
//    Matrix3<float> tempReturn(wrapper.getDimensions());
//    for (size_t i = 0; i < tempReturn.size(); i++) {
//        tempReturn[i] = wrapper[i].x;
//    }
//    return tempReturn;
    return this->wrapWith(wrapper);

    //    Matrix3<float> unit_ctrl(this->getDimensions(), 1.f);
}

template<class T>
Matrix3<T> Matrix3<T>::wrapWithoutInterpolation(Matrix3<Vector3> wrapper)
{
    Matrix3<T> result = *this; //(getDimensions());
    this->raiseErrorOnBadCoord = false;
    this->returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    result.raiseErrorOnBadCoord = false;
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                Vector3 pos(x, y, z);
                Vector3& wrap = wrapper.at(pos);
                result.at(pos + wrap) = this->at(pos);
            }
        }
    }
    return result;
}

template<class T>
Matrix3<T> Matrix3<T>::wrapWithoutInterpolation(BSpline original, BSpline wrapperCurve)
{
    bool previousRaise = this->raiseErrorOnBadCoord;
    T previousDefault  = this->returned_value_on_outside;

    this->raiseErrorOnBadCoord = false;
    this->returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    Matrix3<Vector3> indices(this->getDimensions());
    for (size_t i = 0; i < this->size(); i++)
        indices[i] = indices.getCoordAsVector3(i);

    indices = indices.wrapWith(original, wrapperCurve);

    Matrix3<T> values(this->getDimensions());
    for (size_t i = 0; i < this->size(); i++) {
        values[i] = this->at(indices[i]);
    }

    this->raiseErrorOnBadCoord = previousRaise;
    this->returned_value_on_outside = previousDefault;

    return values;
}

template<class T>
Matrix3<float> Matrix3<T>::fbmNoise1D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Matrix3<float> values(sizeX, sizeY, sizeZ);
    for (size_t i = 0; i < values.size(); i++) {
        float x, y, z;
        std::tie(x, y, z) = values.getCoord(i);

        values.at(i) = noise.GetNoise(x, y, z);
    }
    return values;
}

template<class T>
Matrix3<Vector3> Matrix3<T>::fbmNoise2D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Matrix3<Vector3> values = Matrix3<Vector3>::fbmNoise3D(noise, sizeX, sizeY, sizeZ);
    for (auto& vec : values)
        vec = vec.xy();
    return values;
}

template<class T>
Matrix3<Vector3> Matrix3<T>::fbmNoise3D(FastNoiseLite noise, int sizeX, int sizeY, int sizeZ)
{
    Vector3 offsetDim1 = Vector3(   0,    0,    0);
    Vector3 offsetDim2 = Vector3(  42,  103, 2048);
    Vector3 offsetDim3 = Vector3(  15,  128, 1000);
    Matrix3<Vector3> values(sizeX, sizeY, sizeZ);
    for (size_t i = 0; i < values.size(); i++) {
        float x, y, z;
        std::tie(x, y, z) = values.getCoord(i);
        // Add some offset (chosen randomly)
        values.at(i) = Vector3(noise.GetNoise(x + offsetDim1.x, y + offsetDim1.y, z + offsetDim1.z),
                               noise.GetNoise(x + offsetDim2.x, y + offsetDim2.y, z + offsetDim2.z),
                               noise.GetNoise(x + offsetDim3.x, y + offsetDim3.y, z + offsetDim3.z));
    }
    return values;
}

template<class T>
int Matrix3<T>::width() const {
    return this->sizeX;
}
template<class T>
int Matrix3<T>::depth() const {
    return this->sizeY;
}
template<class T>
int Matrix3<T>::height() const {
    return this->sizeZ;
}

template<class T>
Vector3 Matrix3<T>::getDimensions() const
{
    return Vector3(this->sizeX, this->sizeY, this->sizeZ);
}



#endif // MATRIX3_H
