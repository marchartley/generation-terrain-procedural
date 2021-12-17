#ifndef MATRIX3_H
#define MATRIX3_H

#include <vector>
#include <tuple>
#include <memory>
#include <string>
#include <math.h>
#include "Vector3.h"

template <class T>
class Matrix3
{
public:
    Matrix3();
    Matrix3(int sizeX, int szeY, int sizeZ = 1);
    Matrix3(std::vector<std::vector<std::vector<T>>> data);
    Matrix3(std::vector<std::vector<T>> data);
    Matrix3(std::vector<T> data, int sizeX, int sizeY, int sizeZ = -1);

    T& at(int i, int j, int k = 0);
    int getIndex(int x, int y, int z);
    std::tuple<int, int, int> getCoord(int index);


    Matrix3<T> rounded(int precision = 0);

    Matrix3<T>& normalize();
    Matrix3<T> normalized();

    Matrix3<Vector3> gradient();
    Matrix3<Vector3> grad();
    Matrix3<float> divergence();
    Matrix3<float> div();
    Matrix3<T> curl();
    Matrix3<T> rot();
    Matrix3<T> laplacian();


    static Matrix3<T> random(int sizeX, int sizeY, int sizeZ = 1);
    static Matrix3<T> identity(int sizeX, int sizeY, int sizeZ = 1);

//    operator float*() const { return new float[3]{this->x, this->y, this->z}; }
    // friend Matrix3<T> operator+(Matrix3<T> a, Matrix3<T> b);
    Matrix3<T>& operator+=(const Matrix3<T>& o);
    // friend Matrix3<T> operator-(Matrix3<T> a, const Matrix3<T>& b);
    Matrix3<T>& operator-=(const Matrix3<T>& o);
    // friend Matrix3<T> operator*(Matrix3<T> a, Matrix3<T> o);
    Matrix3<T>& operator*=(Matrix3<T>& o);
    // friend Matrix3<T> operator/(Matrix3<T> a, Matrix3<T> o);
    Matrix3<T>& operator/=(Matrix3<T>& o);
    // friend Matrix3<T> operator*(Matrix3<T> a, float o);
    Matrix3<T>& operator*=(float o);
    // friend Matrix3<T> operator/(Matrix3<T> a, float o);
    Matrix3<T>& operator/=(float o);
    // friend Matrix3<T> operator+(Matrix3<T> a, float o);
    Matrix3<T>& operator+=(float o);
    // friend Matrix3<T> operator-(Matrix3<T> a, float o);
    Matrix3<T>& operator-=(float o);
//    Matrix3<T>& operator=(const Matrix3<T>& o);
    bool operator==(Matrix3<T> o);

    std::string toString() const {return "Matrix3 (" + std::to_string(this->sizeX) + "x" + std::to_string(this->sizeY) + "x" + std::to_string(this->sizeZ) + ")"; }

    std::vector<T> data;
    int sizeX, sizeY, sizeZ;

    Matrix3& init(std::vector<T> data, int sizeX, int sizeY, int sizeZ);
};




template<class T>
Matrix3<T>::Matrix3()
{
}
template<class T>
Matrix3<T>::Matrix3(int sizeX, int sizeY, int sizeZ)
{
    std::vector<T> data(sizeX * sizeY * sizeZ);
    init(data, sizeX, sizeY, sizeZ);
}
template<class T>
Matrix3<T>::Matrix3(std::vector<T> data, int sizeX, int sizeY, int sizeZ)
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
T &Matrix3<T>::at(int i, int j, int k)
{
    if (i >= 0 && i < this->sizeX && j >= 0 && j < this->sizeY && k >= 0 && k < this->sizeZ)
        return this->data[getIndex(i, j, k)];
    throw std::out_of_range("Trying to access coord (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") on matrix of size "
        + std::to_string(sizeX) + "x" + std::to_string(sizeY) + "x" + std::to_string(sizeZ) + ". Max index is " + std::to_string(sizeX * sizeY * sizeZ - 1));
}

template<class T>
int Matrix3<T>::getIndex(int x, int y, int z)
{
    return z * (this->sizeX * this->sizeY) + y * (this->sizeX) + x;
}

template<class T>
std::tuple<int, int, int> Matrix3<T>::getCoord(int index)
{
    int z = index / (this->sizeX * this->sizeY);
    int y = (index % (this->sizeX * this->sizeY)) / this->sizeX;
    int x = index % this->sizeX;
    return std::make_tuple(x, y, z);
}

template<class T>
Matrix3<T>& Matrix3<T>::init(std::vector<T> data, int sizeX, int sizeY, int sizeZ)
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
Matrix3<T> Matrix3<T>::rounded(int precision)
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
Matrix3<T> Matrix3<T>::normalized() {
    Matrix3 mat = *this;
    return mat.normalize();
}

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
}
template<class T>
Matrix3<Vector3> Matrix3<T>::grad(){return this->gradient(); }

template<class T>
Matrix3<T> Matrix3<T>::laplacian()
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
Matrix3<T> Matrix3<T>::identity(int sizeX, int sizeY, int sizeZ)
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
Matrix3<T>& Matrix3<T>::operator+=(const Matrix3<T>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] += o.data[i];
    }
    return *this;
}
template<typename T>
Matrix3<T> operator-(Matrix3<T> a, const Matrix3<T>& b) {
    a -= b;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator-=(const Matrix3<T>& o)  {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o.data[i];
    }
    return *this;
}
template<typename T>
Matrix3<T> operator*(Matrix3<T> a, Matrix3<T> o) {
    a *= o;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator*=(Matrix3<T>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o.data[i];
    }
    return *this;
}
template<typename T>
Matrix3<T> operator/(Matrix3<T>& a, const Matrix3<T>& b) {
    a /= b;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator/=(Matrix3<T>& o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        throw std::domain_error("Matrices maust have same sizes to be added (M1 = " + this->toString() + " and M2 = " + o.toString());
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o.data[i];
    }
    return *this;
}
template<typename T>
Matrix3<T> operator*(Matrix3<T> a, float o) {
    a *= o;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator*=(float o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= o;
    }
    return *this;
}

template<typename T>
Matrix3<T> operator/(Matrix3<T> a, float o) {
    a /= o;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator/=(float o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] /= o;
    }
    return *this;
}
template<typename T>
Matrix3<T> operator+(Matrix3<T> a, float o) {
    a += o;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator+=(float o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] += o;
    }
    return *this;
}
template<typename T>
Matrix3<T> operator-(Matrix3<T> a, float o) {
    a -= o;
    return a;
}
template<typename T>
Matrix3<T>& Matrix3<T>::operator-=(float o) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= o;
    }
    return *this;
}
template<typename T>
bool Matrix3<T>::operator==(Matrix3<T> o) {
    if (this->sizeX != o.sizeX || this->sizeY != o.sizeY || this->sizeZ != o.sizeZ)
        return false;
    for (size_t i = 0; i < this->data.size(); i++)
        if (this->data[i] != o.data[i])
            return false;
    return true;
}
#endif // MATRIX3_H
