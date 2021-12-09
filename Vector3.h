#ifndef VECTOR3_H
#define VECTOR3_H

#include "Globals.h"
#include <iostream>
#include <vector>

#define GRID1D_float std::vector<float>
#define GRID1D_vec3 std::vector<Vector3>
#define GRID2D_float std::vector<std::vector<float>>
#define GRID2D_vec3 std::vector<std::vector<Vector3>>
#define GRID3D_float std::vector<std::vector<std::vector<float>>>
#define GRID3D_vec3 std::vector<std::vector<std::vector<Vector3>>>

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

    float divergence() { return x + y + z; }

    static GRID3D_vec3 gradient(GRID3D_float field);
    static GRID3D_vec3 grad(GRID3D_float field);
    static GRID3D_float divergence(GRID3D_vec3 field);
    static GRID3D_float div(GRID3D_vec3 field);
    static GRID3D_vec3 curl(GRID3D_vec3 field);
    static GRID3D_vec3 rot(GRID3D_vec3 field);
    static GRID3D_float laplacian(GRID3D_float field);
    static GRID3D_vec3 laplacian(GRID3D_vec3 field);


    static Vector3 random();

    static Vector3 nabla;

    operator float*() const { return new float[3]{this->x, this->y, this->z}; }
//    friend Vector3 operator+(Vector3 a, Vector3& b);
    friend Vector3 operator+(Vector3 a, Vector3 b);
    Vector3& operator+=(const Vector3& o);
    friend Vector3 operator-(Vector3 a, const Vector3& b);
    Vector3& operator-=(const Vector3& o);
    friend Vector3 operator*(Vector3 a, Vector3 o);
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

    std::string toString() const {return "Vector3 (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")"; }

};

#endif // VECTOR3_H
