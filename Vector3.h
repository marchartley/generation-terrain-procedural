#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>

class Vector3 {
public:
    float x, y, z;
    Vector3();
    Vector3(float x, float y, float z = 0.f);
    Vector3(const Vector3& copy);
    Vector3(Vector3* copy);

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

    friend std::ostream& operator<<(std::ostream& io, const Vector3& v);
    friend std::ostream& operator<<(std::ostream& io, Vector3* v);

    float dot(Vector3& o);
    Vector3 cross(Vector3 o);

    float norm();
    Vector3& normalize();

};

#endif // VECTOR3_H
