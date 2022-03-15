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

enum RESIZE_MODE {
    LINEAR = 0,
    MAX_VAL = 1,
    MIN_VAL = 2
};

// Warning : don't use bool type...
// This class is based on a std::vector, which has specifications on bools
// the [] operator won't work ... Use int or short int instead
template <class T>
class Matrix3
{
public:
    Matrix3();
    Matrix3(size_t sizeX, size_t szeY, size_t sizeZ = 1, T initValue = T());
    Matrix3(Vector3 size, T initValue = T());
    Matrix3(std::vector<std::vector<std::vector<T>>> data);
    Matrix3(std::vector<std::vector<T>> data);
    Matrix3(std::vector<T> data, size_t sizeX, size_t sizeY, int sizeZ = -1);

    T& at(size_t i, size_t j, size_t k = 0);
    T& at(Vector3 pos);
    T& at(size_t i);
    T& operator()(size_t x, size_t y, size_t z = 0);
    T& operator[](size_t i);
    T& operator[](Vector3 pos);

    int getIndex(size_t x, size_t y, size_t z) const;
    std::tuple<size_t, size_t, size_t> getCoord(size_t index) const;
    bool checkCoord(int x, int y, int z = 0) const;
    bool checkCoord(Vector3 pos) const;
    bool checkIndex(size_t i) const;

    int getNumberNeighbors(size_t x, size_t y, size_t z, bool using4connect = true) const;
    int getNumberNeighbors(Vector3 pos, bool using4connect = true) const;

    Matrix3<T> resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode = LINEAR);
    Matrix3 resize(Vector3 newSize, RESIZE_MODE mode = LINEAR);

    Matrix3<T> subset(int startX, int endX, int startY, int endY, int startZ = 0, int endZ = -1);
    Matrix3<T>& paste(Matrix3<T> matrixToPaste, Vector3 upperLeftFrontCorner);
    Matrix3<T>& paste(Matrix3<T> matrixToPaste, int left, int up, int front);
    Matrix3<T>& add(Matrix3<T> matrixToAdd, Vector3 upperLeftFrontCorner);
    Matrix3<T>& add(Matrix3<T> matrixToAdd, int left, int up, int front);

    Matrix3<float> toDistanceMap();

    T min() const;
    T max() const;

    Matrix3<T> abs() const;

    Matrix3<T> rounded(int precision = 0) const;

    Matrix3<T>& normalize();
    Matrix3<T> normalized() const;

    Matrix3<Vector3> gradient();
    Matrix3<Vector3> grad();
    Matrix3<float> divergence();
    Matrix3<float> div();
    Matrix3<T> curl();
    Matrix3<T> rot();
    Matrix3<T> laplacian() const;

    Matrix3<T>& insertRow(size_t indexToInsert, int affectedDimension, T newData = T());

    void clear() { this->sizeX = 0; this->sizeY = 0; this->sizeZ = 0; return this->data.clear(); }


    static Matrix3<T> random(size_t sizeX, size_t sizeY, size_t sizeZ = 1);
    static Matrix3<T> identity(size_t sizeX, size_t sizeY, size_t sizeZ = 1);

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
    int sizeX, sizeY, sizeZ;

    bool raiseErrorOnBadCoord = true;
    T defaultValueOnBadCoord = T();

    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    std::size_t size() const { return end() - begin(); }
    bool empty() const { return begin() == end(); }

    Matrix3& init(std::vector<T> data, size_t sizeX, size_t sizeY, size_t sizeZ);

    std::string displayValues();
};

#include <sstream>
template<class T>
std::string Matrix3<T>::displayValues()
{
    std::stringstream out;
    for (int z = 0; z < this->sizeZ; z++) {
        out << "[Z-level = " << z << "] : \n";
        for (int y = 0; y < this->sizeY; y++) {
            for (int x = 0; x < this->sizeX; x++) {
                out << std::setw(5) << at(x, y, z) << "\t";
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
    int sizeY = (data.size() > 0 ? data[0].size() : 0);
    init(oneMatrix, data.size(), sizeY, 1);
}
template<class T>
Matrix3<T>::Matrix3(std::vector<std::vector<std::vector<T>>> data)
{
    std::vector<T> oneMatrix;
    for (std::vector<std::vector<T>>& grid : data)
        for (std::vector<T>& row : grid)
            oneMatrix.insert(oneMatrix.end(), row.begin(), row.end());
    int sizeY = (data.size() > 0 ? data[0].size() : 0);
    int sizeZ = (sizeY > 0 ? data[0][0].size() : 0);
    init(oneMatrix, data.size(), sizeY, sizeZ);
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
T &Matrix3<T>::at(Vector3 pos)
{
    return this->at(pos.x, pos.y, pos.z);
}
template<class T>
T &Matrix3<T>::at(size_t i, size_t j, size_t k)
{
    if (checkCoord(i, j, k)) {
        int index = getIndex(i, j, k);
        return this->data[index];
    }
    if (!raiseErrorOnBadCoord)
        return defaultValueOnBadCoord;
    throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
        + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template<class T>
T &Matrix3<T>::at(size_t i)
{
    if (i >= 0 && i < sizeX * sizeY * sizeZ) {
        return this->data[i];
    }
    if (!raiseErrorOnBadCoord)
        return defaultValueOnBadCoord;
    int x, y, z;
    std::tie(x, y, z) = this->getCoord(i);
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
std::tuple<size_t, size_t, size_t> Matrix3<T>::getCoord(size_t index) const
{
    int z = index / (this->sizeX * this->sizeY);
    int y = (index % (this->sizeX * this->sizeY)) / this->sizeX;
    int x = index % this->sizeX;
    return std::make_tuple(x, y, z);
}

template<class T>
Matrix3<T>& Matrix3<T>::init(std::vector<T> data, size_t sizeX, size_t sizeY, size_t sizeZ)
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
    for (T& val : data)
    {
        val = (val - min)/ (max - min);
    }
    return *this;
}
template<class T>
Matrix3<T> Matrix3<T>::normalized() const {
    Matrix3 mat = *this;
    return mat.normalize();
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
    this->raiseErrorOnBadCoord = false;
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                returningGrid.at(x, y, z) = Vector3((at(x + 1, y, z) - at(x - 1, y, z)) * .5f,
                                                    (at(x, y + 1, z) - at(x, y - 1, z)) * .5f,
                                                    (at(x, y, z + 1) - at(x, y, z - 1)) * .5f);
            }
        }
    }
    this->raiseErrorOnBadCoord = true;
    return returningGrid;
}
template<class T>
Matrix3<Vector3> Matrix3<T>::grad() {return this->gradient(); }

template<class T>
Matrix3<T> Matrix3<T>::laplacian() const
{
    Matrix3 returningGrid = *this;
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                Vector3 vec = Vector3();
                if (x == 0 || x == this->sizeX - 1 || y == 0 || y == this->sizeY - 1
                        || z == 0 || z == this->sizeZ - 1) {
                    returningGrid.at(x, y, z) = vec;
                } else {
                    vec += this->at(x    , y    , z + 1);
                    vec += this->at(x    , y    , z - 1);
                    vec += this->at(x    , y + 1, z    );
                    vec += this->at(x    , y - 1, z    );
                    vec += this->at(x + 1, y    , z    );
                    vec += this->at(x - 1, y    , z    );
                    vec -= this->at(x    , y    , z    ) * 6;
                }
                returningGrid.at(x, y, z) = vec;
            }
        }
    }
    return returningGrid;
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
template<typename T> template<typename U>
Matrix3<T>& Matrix3<T>::operator+=(const Matrix3<U>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
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
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
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
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
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
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
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
Matrix3<T> Matrix3<T>::resize(Vector3 newSize, RESIZE_MODE mode)
{
    return this->resize(newSize.x, newSize.y, newSize.z, mode);
}
template<typename T>
Matrix3<T> Matrix3<T>::resize(size_t newX, size_t newY, size_t newZ, RESIZE_MODE mode)
{
    Matrix3<T> newMat(newX, newY, newZ);
    float rx = this->sizeX / (float)newX, ry = this->sizeY / (float)newY, rz = this->sizeZ / (float)newZ;

    if (mode == LINEAR) {
        // Apply interpolations
        for (int x = 0; x < newX; x++) {
            int x_original = int(x * rx);
            int x_plus_1 = (x_original >= this->sizeX - 2 ? x_original : x_original + 1);
            float d_x = (x * rx) - x_original;
            for (int y = 0; y < newY; y++) {
                int y_original = int(y * ry);
                int y_plus_1 = (y_original >= this->sizeY - 2 ? y_original : y_original + 1);
                float d_y = (y * ry) - y_original;
                for (int z = 0; z < newZ; z++) {
                    int z_original = int(z * rz);
                    int z_plus_1 = (z_original >= this->sizeZ - 2 ? z_original : z_original + 1);
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
                    newMat.at(x, y, z) = ((
                                              f000 * (1-d_x) + f100 * d_x) * (1-d_y) + (
                                              f010 * (1-d_x) + f110 * d_x) * d_y) * (1 - d_z) +
                                        ((
                                             f001 * (1-d_x) + f101 * d_x) * (1-d_y) + (
                                             f011 * (1-d_x) + f111 * d_x) * d_y) * d_z;
                }
            }
        }
    }
    if (mode == MAX_VAL || mode == MIN_VAL) {
//        for (auto& val : newMat)
//            val = (mode == MAX_VAL ? std::numeric_limits<T>::min() : std::numeric_limits<T>::max());
        newMat.raiseErrorOnBadCoord = false;
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
        newMat.raiseErrorOnBadCoord = this->raiseErrorOnBadCoord;
    }
    return newMat;
}


template<typename T>
Matrix3<T> Matrix3<T>::subset(int startX, int endX, int startY, int endY, int startZ, int endZ)
{
    if (endZ == -1) endZ = this->sizeZ;
    Matrix3<T> croppedMatrix(endX - startX, endY - startY, endZ - startZ);
    for (int x = startX; x < endX; x++) {
        for (int y = startY; y < endY; y++) {
            for (int z = startZ; z < endZ; z++) {
                croppedMatrix.at(x - startX, y - startY, z - startZ) = this->at(x, y, z);
            }
        }
    }
    return croppedMatrix;
}


template<typename T>
Matrix3<T>& Matrix3<T>::paste(Matrix3<T> matrixToPaste, Vector3 upperLeftFrontCorner)
{
    return this->paste(matrixToPaste, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z);
}
template<typename T>
Matrix3<T>& Matrix3<T>::paste(Matrix3<T> matrixToPaste, int left, int up, int front)
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
Matrix3<T>& Matrix3<T>::add(Matrix3<T> matrixToAdd, Vector3 upperLeftFrontCorner)
{
    return this->add(matrixToAdd, upperLeftFrontCorner.x, upperLeftFrontCorner.y, upperLeftFrontCorner.z);
}
template<typename T>
Matrix3<T>& Matrix3<T>::add(Matrix3<T> matrixToAdd, int left, int up, int front)
{
    for (int x = std::max(left, 0); x < std::min(matrixToAdd.sizeX + left, this->sizeX); x++) {
        for (int y = std::max(up, 0); y < std::min(matrixToAdd.sizeY + up, this->sizeY); y++) {
            for (int z = std::max(front, 0); z < std::min(matrixToAdd.sizeZ + front, this->sizeZ); z++) {
                this->at(x, y, z) += matrixToAdd.at(x - left, y - up, z - front);
            }
        }
    }
    return *this;
}

template<class T>
Matrix3<float> Matrix3<T>::toDistanceMap()
{
    Matrix3<float> distances(this->sizeX, this->sizeY, this->sizeZ, std::numeric_limits<float>::max() - 10000);
    distances.raiseErrorOnBadCoord = false;
    distances.defaultValueOnBadCoord = 0.f;

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

#endif // MATRIX3_H
