#ifndef VECTOR3_H
#define VECTOR3_H

#include "Utils/Globals.h"
#include <iostream>
#include <vector>
#include <QGLViewer/vec.h>
#include "DataStructure/Matrix.h"

class Vector3 {
public:
    float x, y, z;
    bool valid = true;
    Vector3();
    Vector3(float x, float y, float z = 0.f, bool valid = true);
    Vector3(const Vector3& copy);
    Vector3(Vector3* copy);
    Vector3(qglviewer::Vec other);
    Vector3(bool valid);
    Vector3(const float* coords, bool valid = true);

    static std::vector<float> toArray(Vector3 v);
    static std::vector<float> toArray(std::vector<Vector3> vs);
    std::tuple<int, int, int> toIntTuple() {return std::make_tuple<int, int, int>(int(this->x), int(this->y), int(this->z)); }


    friend std::ostream& operator<<(std::ostream& io, const Vector3& v);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Vector3> v);

    float dot(Vector3 o);
    Vector3 cross(Vector3 o);
    Vector3 rounded(int precision = 0) const;

    float norm();
    float norm2();
    Vector3& normalize();
    Vector3 normalized() const;
    Vector3 abs();

    bool isAlmostVertical();

    Vector3& rotate(float angle_x, float angle_y, float angle_z);
    Vector3& rotate(Vector3 eulerAngles);
    Vector3 rotated(float angle_x, float angle_y, float angle_z);
    Vector3 rotated(Vector3 eulerAngles);
    Vector3& rotate(float angle, float dir_x, float dir_y, float dir_z);
    Vector3& rotate(float angle, Vector3 direction);
    Vector3 rotated(float angle, float dir_x, float dir_y, float dir_z);
    Vector3 rotated(float angle, Vector3 direction);
    Vector3& translate(float move_x, float move_y, float move_z);
    Vector3& translate(Vector3 move);
    Vector3 translated(float move_x, float move_y, float move_z);
    Vector3 translated(Vector3 move);
//    Vector3& setDirection(Vector3 dir, Vector3 upVector = Vector3(0, 0, 1));

    Vector3& changeBasis(Vector3 newX, Vector3 newY, Vector3 newZ);
    Vector3 changedBasis(Vector3 newX, Vector3 newY, Vector3 newZ);

    float divergence() { return x + y + z; }

    Matrix toMatrix();

    static Vector3 lerp(float t, Vector3 min, Vector3 max) {
        return min + (max - min) * t;
    }
    static float inverseLerp(Vector3 val, Vector3 min, Vector3 max) {
        return (val - min).dot(max - min) / ((max - min).norm2() > 0 ? (max - min).norm2() : 0.0001f);
    }
    float inverseLerp(Vector3 min, Vector3 max) {
        return Vector3::inverseLerp(*this, min, max);
    }
    static float remap(Vector3 val, Vector3 oldMin, Vector3 oldMax, float newMin, float newMax)
    {
        float oldProgress = Vector3::inverseLerp(val, oldMin, oldMax);
        return newMin + (newMax - newMin) * oldProgress;
    }
    float remap(Vector3 oldMin, Vector3 oldMax, float newMin, float newMax)
    {
        return Vector3::remap(*this, oldMin, oldMax, newMin, newMax);
    }

    static Vector3 random();

    static Vector3 nabla;

    bool isValid() const { return this->valid; }
    operator qglviewer::Vec() const { return qglviewer::Vec(this->x, this->y, this->z); }
    operator float*() const { return new float[3]{this->x, this->y, this->z}; }
//    friend Vector3 operator+(Vector3 a, Vector3& b);
    friend Vector3 operator+(Vector3 a, Vector3 b);
    Vector3& operator+=(const Vector3& o);
    friend Vector3 operator-(Vector3 a, const Vector3& b);
    Vector3& operator-=(const Vector3& o);
    friend Vector3 operator*(Vector3 a, Vector3 o);
    Vector3& operator*=(Vector3 o);
    Vector3 operator/(Vector3 o);
    Vector3& operator/=(Vector3 o);
    Vector3 operator*(float o);
    Vector3& operator*=(float o);
    Vector3 operator/(float o);
    Vector3& operator/=(float o);
    Vector3 operator+(float o);
    Vector3& operator+=(float o);
    Vector3 operator-(float o);
    Vector3& operator-=(float o);
    Vector3& operator=(const Vector3& o);
    friend bool operator==(Vector3 a, Vector3 b);
    friend bool operator!=(Vector3 a, Vector3 b);
    friend bool operator<(Vector3 a, Vector3 b);
    friend bool operator<=(Vector3 a, Vector3 b);
    friend bool operator>(Vector3 a, Vector3 b);
    friend bool operator>=(Vector3 a, Vector3 b);
//    using std::abs;
//    friend Vector3 abs(Vector3 o) { return o.abs(); }

    std::string toString() const {return "Vector3 (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")"; }
//    const char* toHashString() const {return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z); }

};

/*class Vector3Hash
{

};*/

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
namespace std {
    Vector3 abs(Vector3 o);
    template <> struct hash<Vector3>
    {
        size_t operator()(const Vector3 & x) const
        {
            size_t seed = 0;
            hash_combine(seed, int(x.x * 100));
            hash_combine(seed, int(x.y * 100));
            hash_combine(seed, int(x.z * 100));
            return seed;
        }
    };
}

#endif // VECTOR3_H
