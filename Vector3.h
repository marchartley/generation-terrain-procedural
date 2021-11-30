#ifndef VECTOR3_H
#define VECTOR3_H

#include "Globals.h"
#include <iostream>
#include <vector>

class Vector3 {
public:
    float x, y, z;
    Vector3();
    Vector3(float x, float y, float z = 0.f);
    Vector3(const Vector3& copy);
    Vector3(Vector3* copy);

    static std::vector<float> toArray(Vector3 v);
    static std::vector<float> toArray(std::vector<Vector3> vs);
    std::tuple<int, int, int> toIntTuple() {return std::make_tuple<int, int, int>(int(this->x), int(this->y), int(this->z)); }


    friend std::ostream& operator<<(std::ostream& io, const Vector3& v);
    friend std::ostream& operator<<(std::ostream& io, Vector3* v);

    float dot(Vector3& o);
    Vector3 cross(Vector3 o);
    Vector3 rounded(int precision = 0);

    float norm();
    Vector3& normalize();

    static Vector3 random();

    operator float*() const { return new float[3]{this->x, this->y, this->z}; }
//    friend Vector3 operator+(Vector3 a, Vector3& b);
    friend Vector3 operator+(Vector3 a, Vector3 b);
    Vector3& operator+=(const Vector3& o);
    friend Vector3 operator-(Vector3 a, const Vector3& b);
    Vector3& operator-=(const Vector3& o);
    friend Vector3 operator*(Vector3 a, const Vector3& o);
    Vector3& operator*=(Vector3& o);
    Vector3 operator/(Vector3& o);
    Vector3& operator/=(Vector3& o);
    Vector3 operator*(float o);
    Vector3& operator*=(float o);
    Vector3 operator/(float o);
    Vector3& operator/=(float o);
    Vector3 operator+(float o);
    Vector3& operator+=(float o);
    Vector3 operator-(float o);
    Vector3& operator-=(float o);
    Vector3& operator=(const Vector3& o);

};

#endif // VECTOR3_H
