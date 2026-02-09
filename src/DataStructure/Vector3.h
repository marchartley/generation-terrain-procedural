#ifndef VECTOR3_H
#define VECTOR3_H

#include "Utils/Globals.h"
#include <iostream>
#include <vector>
#include <QGLViewer/vec.h>
#include <QGLViewer/quaternion.h>
#include "DataStructure/Matrix.h"
//#include <glm/glm.hpp>
#include <cmath>

class AABBox;
template <class T>
class Vec3;

template <class T>
class Vec3 {
public:
    Vec3();
    constexpr Vec3(T x, T y, T z = 0.f, bool valid = true);
//    Vec3(const Vec3& copy);
//    Vec3(Vec3* copy);
    Vec3(qglviewer::Vec other);
    explicit Vec3(bool valid);
    // explicit Vec3(const T* coords, bool valid = true);

    ~Vec3() = default;

    static std::vector<T> toArray(const Vec3& v);
    static std::vector<T> toArray(const std::vector<Vec3> &vs);
    std::tuple<int, int, int> toIntTuple() const {return {int(x()), int(y()), int(z())}; }

    static Vec3 min();
    static Vec3 max();
    static Vec3 min(const Vec3<T>& a, const Vec3<T>& b);
    static Vec3 max(const Vec3<T>& a, const Vec3<T>& b);
    static Vec3 min(const std::vector<Vec3>& allVectors);
    static Vec3 max(const std::vector<Vec3>& allVectors);

    static std::vector<Vec3> getAABBoxVertices(const Vec3& mini, const Vec3& maxi);

    // friend std::ostream& operator<<(std::ostream& io, const Vec3<T>& bbox);
    // friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Vec3<T>> bbox);

    T dot(const Vec3& o) const;
    Vec3 cross(const Vec3& o) const;
    Vec3 rounded(int precision = 0) const;
    Vec3 roundedUp(int precision = 0) const;
    Vec3 roundedDown(int precision = 0) const;
    Vec3 floor() const;
    Vec3 ceil() const;
    static Vec3 wrap(const Vec3& p, const Vec3& mini, const Vec3& maxi);

    float magnitude() const;
    float length() const;
    float lengthSquared() const;
    float norm() const;
    float norm2() const;
    Vec3& normalize();
    Vec3 normalized() const;
    Vec3 abs() const;
    Vec3& setMag(float newMag);
    Vec3& maxMagnitude(float maxMag);
    Vec3& minMagnitude(float minMag);
    Vec3& clamp(float minMag, float maxMag);
    Vec3 clamped(float minMag, float maxMag) const;
    Vec3& clamp(const Vec3 &minBound, const Vec3 &maxBound);
    Vec3 clamped(const Vec3& minBound, const Vec3& maxBound) const;

    bool isAlmostVertical();

    Vec3& rotate(float angle_x, float angle_y, float angle_z);
    Vec3& rotate(const Vec3& eulerAngles);
    Vec3 rotated(float angle_x, float angle_y, float angle_z) const;
    Vec3 rotated(const Vec3& eulerAngles) const;
    Vec3& rotate(float angle, float dir_x, float dir_y, float dir_z);
    Vec3& rotate(float angle, const Vec3& direction);
    Vec3 rotated(float angle, float dir_x, float dir_y, float dir_z) const;
    Vec3 rotated(float angle, const Vec3& direction) const;
    Vec3& translate(T move_x, T move_y, T move_z);
    Vec3& translate(const Vec3& move);
    Vec3 translated(T move_x, T move_y, T move_z);
    Vec3 translated(const Vec3& move);
    Vec3& applyTransform(Matrix transformMatrix);
//    Vec3& setDirection(const Vec3& dir, const Vec3& upVector = Vec3(0, 0, 1));
    Vec3 rotated90XY() const;

    Vec3& changeBasis(const Vec3& newX, const Vec3& newY, const Vec3& newZ);
    Vec3 changedBasis(const Vec3& newX, const Vec3& newY, const Vec3& newZ);

    Vec3 reflexion(const Vec3& normal);

    Vec3 toPolar();
    Vec3 fromPolar();

    T divergence() { return x() + y() + z(); }

    static Vec3 fromMatrix(Matrix mat);
    Matrix toMatrix() const;
    Matrix toRotationMatrix() const;


    Vec3 toEulerAngles();
    Vec3 eulerAnglesWith(const Vec3& other);

    Vec3 getAllAnglesWith(const Vec3& otherVector) const;
    float getAngleWith(const Vec3& otherVector) const;
    float getSignedAngleWith(const Vec3& otherVector) const;
    float getSignedAngleAroundAxisWith(const Vec3& otherVector, const Vec3& axis) const;

    static Vec3 quaternionToEuler(qglviewer::Quaternion quaternion);
    static Vec3 quaternionToEuler(float x, float y, float z, float w);

    static Vec3 slerp(float t, const Vec3<T>& a, const Vec3<T>& b);

    static Vec3 lerp(float t, const Vec3& min, const Vec3& max) {
        return min + (max - min) * t;
    }
    static float inverseLerp(const Vec3& val, const Vec3& min, const Vec3& max) {
        Vec3 d = max - min;
        float d2 = d.norm2();
        return (val - min).dot(d) / (d2 > 0 ? d2 : 1e-4f);

    }
    float inverseLerp(const Vec3& min, const Vec3& max) {
        return Vec3::inverseLerp(*this, min, max);
    }
    static float remap(const Vec3& val, const Vec3& oldMin, const Vec3& oldMax, float newMin, float newMax)
    {
        float oldProgress = Vec3::inverseLerp(val, oldMin, oldMax);
        return newMin + (newMax - newMin) * oldProgress;
    }
    float remap(const Vec3& oldMin, const Vec3& oldMax, float newMin, float newMax)
    {
        return Vec3::remap(*this, oldMin, oldMax, newMin, newMax);
    }


    template <class U>
    operator Vec3<U>() const {
        return Vec3<U>(this->x(), this->y(), this->z());
    }

    static Vec3 random();
    static Vec3 random(float norm);
    static Vec3 random(float minNorm, float maxNorm);
    static Vec3 random(const Vec3& maxValues);
    static Vec3 random(const Vec3& minValues, const Vec3& maxValues);
    static Vec3 random(const AABBox& bounds);

    T maxComp() const { return std::max(x(), std::max(y(), z())); }
    T minComp() const { return std::min(x(), std::min(y(), z())); }


    template <class U, typename V>
    static bool isInBox(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<V>& maxPos);
    template <class U>
    static float signedManhattanDistanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension = false);
    template <class U>
    static float manhattanDistanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension = false);
    template <class U>
    static float signedDistanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension = false);
    template <class U>
    static float distanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension = false);

    bool isValid() const { return this->valid && std::isfinite(x()) && std::isfinite(y()) && std::isfinite(z()); }
    void setValid(bool newValidValue) { this->valid = newValidValue; }
    static Vec3 invalid() { return Vec3(false); }
    static Vec3 origin() { return Vec3(); }
    operator qglviewer::Vec() const { return qglviewer::Vec(this->x(), this->y(), this->z()); }
    T* data() { return v; }
    const T* data() const { return v; }
//    operator glm::vec3() const { return glm::vec3(this->x, this->y, this->z); }
//    friend Vec3 operator+(const Vec3<T>& a, Vec3& b);
    template <class U>
    Vec3& operator+=(const Vec3<U>& o);
    template <class U>
    Vec3& operator-=(const Vec3<U>& o);
    template <class U>
    Vec3& operator*=(const Vec3<U>& o);
    // template <class U>
//    Vec3 operator/(const Vec3<U>& o);
    template <class U>
    Vec3& operator/=(const Vec3<U>& o);
    // template <class U>
//    Vec3 operator*(U o) const;
    template <class U>
    Vec3& operator*=(U o);
    // template <class U>
//    Vec3 operator/(U o);
    template <class U>
    Vec3& operator/=(U o);
    // template <class U>
//    Vec3 operator+(U o);
    template <class U>
    Vec3& operator+=(U o);
    // template <class U>
//    Vec3 operator-(U o);
    template <class U>
    Vec3& operator-=(U o);
    // Friend functions
    /*friend Vec3 operator+(const Vec3<T>& a, const Vec3<T>& b);
    friend Vec3 operator-(const Vec3<T>& a, const Vec3<T>& b);
    friend Vec3 operator*(const Vec3<T>& a, const Vec3& o);
    friend Vec3 operator/(const Vec3<T>& a, const Vec3& o);

//    friend Vec3 operator+(T a, const Vec3<T>& b);
//    friend Vec3 operator-(T a, const Vec3<T>& b);
    friend Vec3 operator*(T a, const Vec3<T>& b);
    friend Vec3 operator/(T a, const Vec3<T>& b);

//    friend Vec3 operator+(const Vec3& b, T a);
//    friend Vec3 operator-(const Vec3& b, T a);
    friend Vec3 operator*(const Vec3& b, T a);
    friend Vec3 operator/(const Vec3& b, T a);

    friend Vec3 operator-(const Vec3& v);

//    Vec3& operator=(const Vec3& o);
    friend bool operator==(const Vec3<T>& a, const Vec3<T>& b);
    friend bool operator!=(const Vec3<T>& a, const Vec3<T>& b);
    friend bool operator<(const Vec3<T>& a, const Vec3<T>& b);
    friend bool operator<=(const Vec3<T>& a, const Vec3<T>& b);
    friend bool operator>(const Vec3<T>& a, const Vec3<T>& b);
    friend bool operator>=(const Vec3<T>& a, const Vec3<T>& b);
//    using std::abs;
//    friend Vec3 abs(const Vec3& o) { return o.abs(); }
    */
    // End friend functions

    T& operator[](size_t i);
    const T& operator[](size_t i) const;

    std::string toString() const {return "Vec3 (" + std::to_string(x()) + ", " + std::to_string(y()) + ", " + std::to_string(z()) + ")" + (this->isValid() ? "" : "[INVALID]"); }
//    const char* toHashString() const {return std::to_string(x()) + "," + std::to_string(y) + "," + std::to_string(z); }

    Vec3 xxx() const { return Vec3(x(), x(), x()); }
    Vec3 xx()  const { return Vec3(x(), x(), 0); }
    Vec3 yyy() const { return Vec3(y(), y(), y()); }
    Vec3 yy()  const { return Vec3(y(), y(), 0); }
    Vec3 xy()  const { return Vec3(x(), y(), 0); }
    Vec3 yz()  const { return Vec3(y(), z(), 0); }
    Vec3 xz()  const { return Vec3(x(), z(), 0); }
    Vec3 xyz() const { return Vec3(x(), y(), z()); }
    Vec3 yxz() const { return Vec3(y(), x(), z()); }
    Vec3 zyx() const { return Vec3(z(), y(), x()); }
    Vec3 xzy() const { return Vec3(x(), z(), y()); }
    T& x() { return v[0]; }
    T& y() { return v[1]; }
    T& z() { return v[2]; }
    const T& x() const {return v[0]; }
    const T& y() const {return v[1]; }
    const T& z() const {return v[2]; }

    T& r() { return v[0]; }
    T& g() { return v[1]; }
    T& b() { return v[2]; }
    const T& r() const {return v[0]; }
    const T& g() const {return v[1]; }
    const T& b() const {return v[2]; }



    inline static constexpr Vec3 white { T(1), T(1), T(1) };
    inline static constexpr Vec3 black { T(0), T(0), T(0) };
    inline static constexpr Vec3 red   { T(1), T(0), T(0) };
    inline static constexpr Vec3 green { T(0), T(1), T(0) };
    inline static constexpr Vec3 blue  { T(0), T(0), T(1) };

// protected:
    T v[3];
    bool valid = true;

};
template <class T, class U>
Vec3<T> operator+(const Vec3<T>& a, const Vec3<U>& b);
template <class T, class U>
Vec3<T> operator-(const Vec3<T>& a, const Vec3<U>& b);
template <class T, class U>
Vec3<T> operator*(const Vec3<T>& a, const Vec3<U>& b);
template <class T, class U>
Vec3<T> operator/(const Vec3<T>& a, const Vec3<U>& b);


template <class T, class U>
Vec3<T> operator+(U a, const Vec3<T>& b);
template <class T, class U>
Vec3<T> operator-(U a, const Vec3<T>& b);
template <class T, class U>
Vec3<T> operator*(U a, const Vec3<T>& b);
template <class T, class U>
Vec3<T> operator/(U a, const Vec3<T>& b);

template <class T, class U>
Vec3<T> operator+(const Vec3<T>& b, U a);
template <class T, class U>
Vec3<T> operator-(const Vec3<T>& b, U a);
template <class T, class U>
Vec3<T> operator*(const Vec3<T>& b, U a);
template <class T, class U>
Vec3<T> operator/(const Vec3<T>& b, U a);

using Vector3 = Vec3<float>;
using Vector3i = Vec3<int>;


class AABBox { //: public std::pair<Vec3, Vec3> {
public:
    AABBox();
    AABBox(const Vector3& mini, const Vector3& maxi);
    AABBox(const std::tuple<Vector3, Vector3>& minMax);
    AABBox(const std::vector<Vector3>& allPointsToContain);
    Vector3 min() const { return this->mini; }
    Vector3 max() const { return this->maxi; }
    Vector3 dimensions() const { return max() - min(); }
    bool contains(const Vector3& position) const { return Vector3::isInBox(position, min(), max()); }
    bool containsXY(const Vector3& position) const { return Vector3::isInBox(position.xy(), min().xy(), max().xy()); }
    Vector3 center() const { return (this->min() + this->max()) * .5f; }
    Vector3 normalize(const Vector3& p);

    AABBox& expand(const Vector3& newPoint);
    AABBox& expand(const std::vector<Vector3>& newPoints);

    float distanceTo(const Vector3& p);

    Vector3 mini;
    Vector3 maxi;

    static Vector3 random(const Vector3& mini, const Vector3& maxi);

    Vector3 intersects(const Vector3& rayStart, const Vector3& rayEnd);

    std::string toString() const {return "AABBox from " + min().toString() + " to " + max().toString(); }
    friend std::ostream& operator<<(std::ostream& io, const AABBox& v);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<AABBox> v);
};

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
namespace std {
template <class T>
Vec3<T> abs(const Vec3<T>& o);
    /*template <> struct hash<Vec3>
    {
        size_t operator()(const Vec3&  x) const
        {
            size_t seed = 0;
            hash_combine(seed, int(x.x() * 100));
            hash_combine(seed, int(x.y() * 100));
            hash_combine(seed, int(x.z() * 100));
            return seed;
        }
    };*/
}



#include "Utils/Utils.h"
#include "Utils/json.h"




template <class T>
constexpr Vec3<T>::Vec3(T x, T y, T z, bool valid) : v{x, y, z}, valid(valid) {

}
template <class T>
Vec3<T>::Vec3() : Vec3(0.f, 0.f, 0.f) {

}/*
template <class T>
Vec3<T>::Vec3(const Vec3& copy) : Vec3(copy.x, copy.y, copy.z, copy.valid) {

}
template <class T>
Vec3<T>::Vec3(Vec3* copy) : Vec3(copy->x, copy->y, copy->z, copy->valid) {

}*/

template <class T>
Vec3<T>::Vec3(qglviewer::Vec other)
    : Vec3(other.x, other.y, other.z)
{

}

template <class T>
Vec3<T>::Vec3(bool valid) : Vec3()
{
    this->valid = valid;
}
/*
template <class T>
Vec3<T>::Vec3(const T *coords, bool valid)
    : Vec3(coords[0], coords[1], coords[2], valid)
{

}
*/
template <class T>
float Vec3<T>::norm() const {
    // if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return sqrt(this->x() * this->x() + this->y() * this->y() + this->z() * this->z());
}
template <class T>
float Vec3<T>::norm2() const {
    // if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return this->x() * this->x() + this->y() * this->y() + this->z() * this->z();
}


template <class T>
Vec3<T>& Vec3<T>::normalize() {
    float n2 = x()*x() + y()*y() + z()*z();
    if (n2 < 1e-10f) return *this;   // (1e-5)^2
    float inv = 1.0f / std::sqrt(n2);
    x() *= inv; y() *= inv; z() *= inv;
    return *this;
}
template <class T>
Vec3<T> Vec3<T>::normalized() const {
    Vec3 a = *this;
    return a.normalize();
}

template <class T>
Vec3<T> Vec3<T>::abs() const
{
    Vec3 a = *this;
    a.x() = std::abs(x());
    a.y() = std::abs(y());
    a.z() = std::abs(z());
    return a;
}

template <class T>
Vec3<T>& Vec3<T>::setMag(float newMag)
{
    this->normalize() *= newMag;
    return *this;
}

template <class T>
Vec3<T>& Vec3<T>::maxMagnitude(float maxMag)
{
    if (this->norm2() > maxMag*maxMag)
        this->setMag(maxMag);
    return *this;
}

template <class T>
Vec3<T>& Vec3<T>::minMagnitude(float minMag)
{
    if (this->norm2() < minMag*minMag)
        this->setMag(minMag);
    return *this;
}

template <class T>
Vec3<T>& Vec3<T>::clamp(float minMag, float maxMag)
{
    float mag2 = this->norm2();
    if (mag2 > maxMag*maxMag)
        this->setMag(maxMag);
    else if (mag2 < minMag*minMag)
        this->setMag(minMag);
    return *this;
}

template <class T>
Vec3<T> Vec3<T>::clamped(float minMag, float maxMag) const
{
    Vec3 v = *this;
    return v.clamp(minMag, maxMag);
}

template <class T>
Vec3<T> &Vec3<T>::clamp(const Vec3 &minBound, const Vec3 &maxBound)
{
    this->x() = (this->x() < minBound.x() ? minBound.x() : this->x() > maxBound.x() ? maxBound.x() : this->x());
    this->y() = (this->y() < minBound.y() ? minBound.y() : this->y() > maxBound.y() ? maxBound.y() : this->y());
    this->z() = (this->z() < minBound.z() ? minBound.z() : this->z() > maxBound.z() ? maxBound.z() : this->z());
    return *this;
}

template <class T>
Vec3<T> Vec3<T>::clamped(const Vec3 &minBound, const Vec3 &maxBound) const
{
    Vec3 v = *this;
    return v.clamp(minBound, maxBound);
}

template <class T>
bool Vec3<T>::isAlmostVertical()
{
    static const Vec3 up(0, 0, 1);
    return std::abs(this->dot(up)) > 0.999;
}

template <class T>
Matrix Vec3<T>::toMatrix() const
{
    return Matrix(3, 1, this->data());
}

template <class T>
Matrix Vec3<T>::toRotationMatrix() const
{
    Matrix Rx (3, 3, std::vector<float>({
                                           1, 0, 0,
                                           0, cos(this->x()), -sin(this->x()),
                                           0, sin(this->x()), cos(this->x())
                                       }).data());
    Matrix Rz (3, 3, std::vector<float>({
                                           cos(this->y()), 0, -sin(this->y()),
                                           0, 1, 0,
                                           sin(this->y()), 0, cos(this->y())
                                       }).data());
    Matrix Ry (3, 3, std::vector<float>({
                                           cos(this->z()), -sin(this->z()), 0,
                                           sin(this->z()), cos(this->z()), 0,
                                           0, 0, 1
                                       }).data());
    return Rx.product(Ry).product(Rz);
}

template <class T>
Vec3<T> Vec3<T>::toEulerAngles()
{
    Vec3 self = this->normalized();
    return Vec3(
        std::acos(self.y()),
        std::acos(self.z()),
        std::acos(self.x())
        );
}

template <class T>
Vec3<T> Vec3<T>::eulerAnglesWith(const Vec3& other)
{
    return Vec3(false);
    /*Vec3 v1 = this->normalized();
    Vec3 v2 = other.normalized();
    Vec3 v3 = v1.cross(v2);
    v2 = v3.cross(v1); // make it perpendicular

    float R11 = v1.x, R12 = v2.x, R13 = v3.x,
            R21 = v1.y, R22 = v2.y, R23 = v3.y,
            R31 = v1.z, R32 = v2.z, R33 = v3.z;
    Vec3 result;
    result.x() = std::atan2(R32, R33) * (180.0 / M_PI);
    result.y() = std::atan2(-1 * R31, std::sqrt(R32 * R32 + R33 * R33)) * (180.0 / M_PI);
    result.z() = std::atan2(R21, R11) * (180.0 / M_PI);
    return result;*/
}

template <class T>
Vec3<T> Vec3<T>::getAllAnglesWith(const Vec3 &otherVector) const
{
    if (*this == otherVector)
        return Vec3(0.f, 0.f, 0.f);
    float onX = this->getSignedAngleAroundAxisWith(otherVector, Vec3(1, 0, 0));
    float onY = this->getSignedAngleAroundAxisWith(otherVector, Vec3(0, 1, 0));
    float onZ = this->getSignedAngleAroundAxisWith(otherVector, Vec3(0, 0, 1));
    return Vec3(onX, onY, onZ);
}

template <class T>
float Vec3<T>::getAngleWith(const Vec3& otherVector) const
{
    float denom = std::sqrt(this->norm2() * otherVector.norm2());
    if (denom < 1e-10) return 0;
    float c = dot(otherVector) / denom;
    c = std::clamp(c, -1.0f, 1.0f);
    return std::acos(c);

}

template <class T>
float Vec3<T>::getSignedAngleWith(const Vec3 &otherVector) const
{
    return this->getAngleWith(otherVector) * sign(this->cross(otherVector).z());
}

template <class T>
float Vec3<T>::getSignedAngleAroundAxisWith(const Vec3 &otherVector, const Vec3 &axis) const
{
    Vec3 normalizedAxis = axis.normalized();
    Vec3 vA = *this - normalizedAxis * this->dot(normalizedAxis); //this->normalized();
    Vec3 vB = otherVector - normalizedAxis * otherVector.dot(normalizedAxis); //.normalized();
    if (vA == vB)
        return 0.f;
    //    float dot = vA.dot(vB);
    Vec3 cross = vA.cross(vB);
    if (cross.norm2() < 1e-5)
        return 0.f;
    float angle = vA.getAngleWith(vB);
    if (cross.dot(axis) < 0)
        angle = -angle;
    return angle;
}

template <class T>
Vec3<T> Vec3<T>::quaternionToEuler(qglviewer::Quaternion quaternion)
{
    return Vec3::quaternionToEuler(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
}

template <class T>
Vec3<T> Vec3<T>::quaternionToEuler(float x, float y, float z, float w)
{
    Vec3 angles;

    // roll (x-axis rotation)
    float sinr_cosp = 2 * (w * x + y * z);
    float cosr_cosp = 1 - 2 * (x * x + y * y);
    angles.x() = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    float sinp = std::sqrt(1 + 2 * (w * y - x * z));
    float cosp = std::sqrt(1 - 2 * (w * y - x * z));
    angles.y() = 2 * std::atan2(sinp, cosp) - M_PI / 2;

    // yaw (z-axis rotation)
    float siny_cosp = 2 * (w * z + x * y);
    float cosy_cosp = 1 - 2 * (y * y + z * z);
    angles.z() = std::atan2(siny_cosp, cosy_cosp);

    return angles;
}

template <class T>
Vec3<T> Vec3<T>::slerp(float t, const Vec3 &A, const Vec3 &B) {
    float lengthA = A.norm();
    float lengthB = B.norm();
    auto start = (lengthA > 0 ? A / lengthA : A);
    auto end = (lengthB > 0 ? B / lengthB : B);
    float dot = start.dot(end);
    dot = std::clamp(dot, -1.0f, 1.0f);

    float theta = std::acos(dot) * t;
    Vec3 relativeVec = end - start * dot;
    relativeVec.normalize(); // Ensure it's a unit vector

    return ((start * std::cos(theta)) + (relativeVec * std::sin(theta))) * interpolation::inv_linear(t, lengthA, lengthB);
}

template <class T>
T Vec3<T>::dot(const Vec3& o) const {
    return (this->x() * o.x()) + (this->y() * o.y()) + (this->z() * o.z());
}
template <class T>
Vec3<T> Vec3<T>::cross(const Vec3& o) const {
    Vec3 v(this->y() * o.z() - this->z() * o.y(),
           this->z() * o.x() - this->x() * o.z(),
           this->x() * o.y() - this->y() * o.x());
    return v;
}
template <class T>
Vec3<T> Vec3<T>::rounded(int precision) const
{
    return this->roundedDown(precision);
    /*Vec3 v = *this;
    v.x() = (int)(v.x() * pow(10, precision)) / (float)(pow(10, precision));
    v.y() = (int)(v.y() * pow(10, precision)) / (float)(pow(10, precision));
    v.z() = (int)(v.z() * pow(10, precision)) / (float)(pow(10, precision));
    return v;*/
}

template <class T>
Vec3<T> Vec3<T>::roundedUp(int precision) const
{
    float power = std::pow(10, precision);
    return ((*this) * power).ceil() / power;
}

template <class T>
Vec3<T> Vec3<T>::roundedDown(int precision) const
{
    float power = std::pow(10, precision);
    return ((*this) * power).floor() / power;
}
template <class T>
Vec3<T> Vec3<T>::floor() const
{
    Vec3 v = *this;
    v.x() = std::floor(v.x());
    v.y() = std::floor(v.y());
    v.z() = std::floor(v.z());
    return v;
}
template <class T>
Vec3<T> Vec3<T>::ceil() const
{
    Vec3 v = *this;
    v.x() = std::ceil(v.x());
    v.y() = std::ceil(v.y());
    v.z() = std::ceil(v.z());
    return v;
}

template <class T>
Vec3<T> Vec3<T>::wrap(const Vec3& p, const Vec3& mini, const Vec3& maxi)
{
    Vec3 newMaxi = maxi - mini;
    Vec3 newP = p - mini;
    /*
//    p -= mini;
//    maxi -= mini;
    Vec3 rounded = newP.roundedDown();
    Vec3 decimals = newP - rounded;
    Vec3 wrap = Vec3(int(rounded.x() + newMaxi.x) % int(newMaxi.x),
                           int(rounded.y() + newMaxi.y) % int(newMaxi.y),
                           int(rounded.z() + newMaxi.z) % int(newMaxi.z)
                           ) + decimals + mini;
    return wrap;*/
    return Vec3(fmod(newP.x(), newMaxi.x()) + mini.x(), fmod(newP.y(), newMaxi.y()) + mini.y(), fmod(newP.z(), newMaxi.z()) + mini.z());
}

template <class T>
float Vec3<T>::magnitude() const
{
    return this->norm();
}

template <class T>
float Vec3<T>::length() const
{
    return this->norm();
}

template <class T>
float Vec3<T>::lengthSquared() const
{
    return this->norm2();
}

template <class T>
Vec3<T> Vec3<T>::random() {
    Vec3 v(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0));
    v.normalize();
    return v;
}

template <class T>
Vec3<T> Vec3<T>::random(float norm)
{
    return Vec3::random() * norm;
}

template <class T>
Vec3<T> Vec3<T>::random(float minNorm, float maxNorm)
{
    return Vec3::random() * random_gen::generate(minNorm, maxNorm);
}

template <class T>
Vec3<T> Vec3<T>::random(const Vec3& maxValues)
{
    return Vec3::random(Vec3(), maxValues);
}

template <class T>
Vec3<T> Vec3<T>::random(const Vec3& minValues, const Vec3& maxValues)
{
    return Vec3(random_gen::generate(minValues.x(), maxValues.x()),
                random_gen::generate(minValues.y(), maxValues.y()),
                random_gen::generate(minValues.z(), maxValues.z()));
}

template <class T>
Vec3<T> Vec3<T>::random(const AABBox &bounds)
{
    return Vec3::random(bounds.min(), bounds.max());
}

template <class T>
std::vector<T> Vec3<T>::toArray(const Vec3& v)
{
    return {v.x(), v.y(), v.z()};
}
template <class T>
std::vector<T> Vec3<T>::toArray(const std::vector<Vec3>& vs)
{
    std::vector<T> arr;
    arr.reserve(vs.size() * 3);
    for (const auto& v : vs) { arr.push_back(v.x()); arr.push_back(v.y()); arr.push_back(v.z()); }
    return arr;
}

template <class T>
Vec3<T> Vec3<T>::min()
{
    return Vec3(std::numeric_limits<T>::lowest(),
                std::numeric_limits<T>::lowest(),
                std::numeric_limits<T>::lowest());
}

template <class T>
Vec3<T> Vec3<T>::max()
{
    return Vec3(std::numeric_limits<T>::max(),
                std::numeric_limits<T>::max(),
                std::numeric_limits<T>::max());
}

template <class T>
Vec3<T> Vec3<T>::min(const Vec3<T>& a, const Vec3<T>& b)
{
    if (!a.isValid()) return b;
    if (!b.isValid()) return a;
    return Vec3(std::min(a.x(), b.x()), std::min(a.y(), b.y()), std::min(a.z(), b.z()));
}

template <class T>
Vec3<T> Vec3<T>::max(const Vec3<T>& a, const Vec3<T>& b)
{
    if (!a.isValid()) return b;
    if (!b.isValid()) return a;
    return Vec3(std::max(a.x(), b.x()), std::max(a.y(), b.y()), std::max(a.z(), b.z()));
}

template <class T>
Vec3<T> Vec3<T>::min(const std::vector<Vec3>& allVectors)
{
    if (allVectors.empty())
        return Vec3(false);
    Vec3 res = allVectors[0];
    for (size_t i = 1; i < allVectors.size(); i++)
        res = Vec3::min(res, allVectors[i]);
    return res;
}

template <class T>
Vec3<T> Vec3<T>::max(const std::vector<Vec3>& allVectors)
{
    if (allVectors.empty())
        return Vec3(false);
    Vec3 res = allVectors[0];
    for (size_t i = 1; i < allVectors.size(); i++)
        res = Vec3::max(res, allVectors[i]);
    return res;
}

template <class T>
std::vector<Vec3<T>> Vec3<T>::getAABBoxVertices(const Vec3& mini, const Vec3& maxi)
{
    T minX = mini.x(), maxX = maxi.x(),
        minY = mini.y(), maxY = maxi.y(),
        minZ = mini.z(), maxZ = maxi.z();
    return {
        Vec3(minX, minY, minZ),
        Vec3(minX, minY, maxZ),
        Vec3(minX, maxY, minZ),
        Vec3(minX, maxY, maxZ),
        Vec3(maxX, minY, minZ),
        Vec3(maxX, minY, maxZ),
        Vec3(maxX, maxY, minZ),
        Vec3(maxX, maxY, maxZ)
    };
}


template <class T>
Vec3<T>& Vec3<T>::rotate(float angle_x, float angle_y, float angle_z) {
    return this->rotate(Vec3(angle_x, angle_y, angle_z));
}
template <class T>
Vec3<T>& Vec3<T>::rotate(const Vec3& eulerAngles) {
    return this->applyTransform(eulerAngles.toRotationMatrix());
    /*Matrix newCoords = R.product(this->toMatrix());
    this->x = newCoords[0][0];
    this->y = newCoords[1][0];
    this->z = newCoords[2][0];
    return *this;*/
}
template <class T>
Vec3<T> Vec3<T>::rotated(float angle_x, float angle_y, float angle_z) const {
    return this->rotated(Vec3(angle_x, angle_y, angle_z));
}
template <class T>
Vec3<T> Vec3<T>::rotated(const Vec3& eulerAngles) const {
    Vec3 v = *this;
    return v.rotate(eulerAngles);
}
template <class T>
Vec3<T>& Vec3<T>::rotate(float angle, float dir_x, float dir_y, float dir_z) {
    return this->rotate(angle, Vec3(dir_x, dir_y, dir_z));
}
template <class T>
Vec3<T>& Vec3<T>::rotate(float angle, const Vec3& direction) {
    float c = cos(angle), s = sin(angle);
    Vec3 v = direction.normalized(); // alias
    Matrix R (3, 3, std::vector<float>({
                                          v.x()*v.x()*(1-c)+c, v.x()*v.y()*(1-c)-v.z()*s, v.x()*v.z()*(1-c)+v.y()*s,
                                          v.x()*v.y()*(1-c)+v.z()*s, v.y()*v.y()*(1-c)+c, v.y()*v.z()*(1-c)-v.x()*s,
                                          v.x()*v.z()*(1-c)-v.y()*s, v.y()*v.z()*(1-c)+v.x()*s, v.z()*v.z()*(1-c)+c
                                      }).data());
    return this->applyTransform(R);
    /*Matrix newCoords = R.product(this->toMatrix());
    this->x = newCoords[0][0];
    this->y = newCoords[1][0];
    this->z = newCoords[2][0];
    return *this;*/
}
template <class T>
Vec3<T> Vec3<T>::rotated(float angle, float dir_x, float dir_y, float dir_z) const {
    return this->rotated(angle, Vec3(dir_x, dir_y, dir_z));
}
template <class T>
Vec3<T> Vec3<T>::rotated(float angle, const Vec3& direction) const {
    Vec3 v = *this;
    return v.rotate(angle, direction);
}
template <class T>
Vec3<T>& Vec3<T>::translate(T move_x, T move_y, T move_z) {
    return this->translate(Vec3(move_x, move_y, move_z));
}
template <class T>
Vec3<T>& Vec3<T>::translate(const Vec3& move) {
    return (*this += move);
}
template <class T>
Vec3<T> Vec3<T>::translated(T move_x, T move_y, T move_z) {
    return this->translated(Vec3(move_x, move_y, move_z));
}
template <class T>
Vec3<T> Vec3<T>::translated(const Vec3& move) {
    Vec3 v = *this;
    return v.translate(move);
}

template <class T>
Vec3<T>& Vec3<T>::applyTransform(Matrix transformMatrix)
{
    // Matrix newCoords = transformMatrix.product(this->toMatrix());
    // this->x() = newCoords[0][0];
    // this->y() = newCoords[1][0];
    // this->z() = newCoords[2][0];
    // return *this;
    const Vec3 newP(
        x() * transformMatrix[0][0] + y() * transformMatrix[0][1] + z() * transformMatrix[0][2],
        x() * transformMatrix[1][0] + y() * transformMatrix[1][1] + z() * transformMatrix[1][2],
        x() * transformMatrix[2][0] + y() * transformMatrix[2][1] + z() * transformMatrix[2][2]
        );
    this->x() = newP.x(); this->y() = newP.y(); this->z() = newP.z();
    return *this;
}

template <class T>
Vec3<T> Vec3<T>::rotated90XY() const
{
    return Vec3(y(), -x());
}

template <class T>
Vec3<T>& Vec3<T>::changeBasis(const Vec3& newX, const Vec3& newY, const Vec3& newZ)
{
    Vec3 newVec = this->changedBasis(newX, newY, newZ);
    this->x() = newVec.x();
    this->y() = newVec.y();
    this->z() = newVec.z();
    return *this;
}

template <class T>
Vec3<T> Vec3<T>::changedBasis(const Vec3& newX, const Vec3& newY, const Vec3& newZ)
{
    Vec3 newVec;
    newVec += newX * this->x() + newY * this->y() + newZ * this->z();
    return newVec;
}

template <class T>
Vec3<T> Vec3<T>::reflexion(const Vec3& normal)
{
    T dot = normal.dot(*this);
    Vec3 n2d = normal * 2.f * dot;
    return *this - n2d;
}

template <class T>
Vec3<T> Vec3<T>::toPolar()
{
    static const Vec3 right(1, 0);
    return Vec3(this->getAngleWith(right) / (2.f * M_PI), this->norm());
}

template <class T>
Vec3<T> Vec3<T>::fromPolar()
{
    return Vec3(1, 0).rotate(0, 0, this->x() * (2.f * M_PI)) * this->y();
}

template <class T>
Vec3<T> Vec3<T>::fromMatrix(Matrix mat)
{
    if (mat.size() == 1) {
        return Vec3(mat[0][0], mat[0][1], mat[0][2]);
    } else if (mat[0].size() == 1) {
        return Vec3(mat[0][0], mat[1][0], mat[2][0]);
    } else {
        throw std::domain_error("Cannot transform Matrix " + mat.toString() + " to Vec3");
    }
}

template <class T, class U>
Vec3<T> operator+(const Vec3<T>& a, const Vec3<U>& b) {
    Vec3 res = a;
    res += b;
    return res;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator+=(const Vec3<U>& o) {
    this->x() += o.x();
    this->y() += o.y();
    this->z() += o.z();
    this->valid &= o.valid;
    return *this;
}
template <class T, class U>
Vec3<T> operator-(const Vec3<T>& a, const Vec3<U>& b) {
    Vec3 res = a;
    res -= b;
    return res;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator-=(const Vec3<U>& o) {
    this->x() -= o.x();
    this->y() -= o.y();
    this->z() -= o.z();
    this->valid &= o.valid;
    return *this;
}

template <class T>
T &Vec3<T>::operator[](size_t i)
{
    return ((T*)(this))[i];
}
template <class T>
const T &Vec3<T>::operator[](size_t i) const
{
    return ((T*)(this))[i];
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator*=(const Vec3<U>& o) {
    this->x() *= o.x();
    this->y() *= o.y();
    this->z() *= o.z();
    this->valid &= o.valid;
    return *this;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator/=(const Vec3<U>& o) {
    this->x() /= o.x();
    this->y() /= o.y();
    this->z() /= o.z();
    this->valid &= o.valid;
    return *this;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator*=(U o) {
    this->x() *= o;
    this->y() *= o;
    this->z() *= o;
    return *this;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator/=(U o) {
    this->x() /= o;
    this->y() /= o;
    this->z() /= o;
    return *this;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator+=(U o) {
    this->x() += o;
    this->y() += o;
    this->z() += o;
    return *this;
}
template <class T> template <class U>
Vec3<T>& Vec3<T>::operator-=(U o) {
    this->x() -= o;
    this->y() -= o;
    this->z() -= o;
    return *this;
}
template <class T, class U>
Vec3<T> operator/(const Vec3<T>& a, const Vec3<U>& b) {
    Vec3 res = a;
    res /= b;
    return res;
}
template <class T, class U>
Vec3<T> operator*(const Vec3<T>& a, const Vec3<U>& b) {
    Vec3 res = a;
    res *= b;
    return res;
}
template <class T, class U>
Vec3<T> operator*(U a, const Vec3<T>& b) {
    Vec3 res = b;
    res *= a;
    return res;
}

template <class T, class U>
Vec3<T> operator/(U a, const Vec3<T>& b) {
    return Vec3(a / b.x(), a / b.y(), a / b.z());
}

template <class T, class U>
Vec3<T> operator*(const Vec3<T>& b, U a) {
    Vec3 res = b;
    res *= a;
    return res;
}
template <class T, class U>
Vec3<T> operator/(const Vec3<T>& b, U a) {
    Vec3 res = b;
    res /= a;
    return res;
}

template <class T, class U>
bool operator==(const Vec3<T>& a, const Vec3<U>& b)
{
    float epsilon = 1e-8;
    return std::abs(a.x() - b.x()) < epsilon && std::abs(a.y() - b.y()) < epsilon && std::abs(a.z() - b.z()) < epsilon; //int(a.x() * 1e8) == int(b.x() * 1e8) && int(a.y() * 1e8) == int(b.y() * 1e8) && int(a.z() * 1e8) == int(b.z() * 1e8);
}

template <class T>
Vec3<T> operator-(const Vec3<T>& v) {
    return Vec3<T>() - v;
}

template <class T, class U>
bool operator!=(const Vec3<T>& a, const Vec3<U>& b)
{
    return !(a == b);
}

template <class T, class U>
bool operator<(const Vec3<T>& a, const Vec3<U>& b)
{
    return a.norm2() < b.norm2();
}

template <class T, class U>
bool operator<=(const Vec3<T>& a, const Vec3<U>& b)
{
    return !(a > b);
}

template <class T, class U>
bool operator>(const Vec3<T>& a, const Vec3<U>& b)
{
    return a.norm2() > b.norm2();
}

template <class T, class U>
bool operator>=(const Vec3<T>& a, const Vec3<U>& b)
{
    return !(a < b);
}

template <class T>
std::ostream& operator<<(std::ostream& io, const Vec3<T>& v) {
    io << v.toString();
    return io;
}

template <class T>
std::ostream& operator<<(std::ostream& io, std::shared_ptr<Vec3<T>> v) {
    io << v->toString();
    return io;
}

template <class T>
Vec3<T> std::abs(const Vec3<T>& o) { return o.abs(); }

template <class T> template <class U, class V>
bool Vec3<T>::isInBox(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<V> &maxPos) {
    return (minPos.x() <= pos.x() && pos.x() <= maxPos.x()) && (minPos.y() <= pos.y() && pos.y() <= maxPos.y()) && (minPos.z() <= pos.z() && pos.z() <= maxPos.z());
    //    return (pos - minPos).minComp() >= 0.f && (pos - (minPos + maxPos)).maxComp() <= 0.f;
}

// Inside : positive
// Outside: negative
template <class T> template <class U>
float Vec3<T>::signedManhattanDistanceToBoundaries(const Vec3& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension)
{
    Vec3<T> newPos = pos - minPos;
    Vec3<T> newMax = maxPos - minPos;

    if (ignoreZdimension) {
        newPos = newPos.xy();
        newMax = newMax.xy();
    }

    static const Vec3 origin;
    if (Vec3::isInBox(newPos, origin, newMax)) {
        if (ignoreZdimension)
            return -std::min({newPos.x(), newPos.y(), newMax.x() - newPos.x(), newMax.y() - newPos.y()});
        else
            return -std::min({newPos.x(), newPos.y(), newPos.z(), newMax.x() - newPos.x(), newMax.y() - newPos.y(), newMax.z() - newPos.z()});
    } else {
        if (ignoreZdimension)
            return -std::min({newPos.x(), newPos.y(), newMax.x() - newPos.x(), newMax.y() - newPos.y()});
        else
            return -std::min({newPos.x(), newPos.y(), newPos.z(), newMax.x() - newPos.x(), newMax.y() - newPos.y(), newMax.z() - newPos.z()});
    }
}

template <class T> template <class U>
float Vec3<T>::manhattanDistanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension)
{
    return std::abs(Vec3::signedManhattanDistanceToBoundaries(pos, minPos, maxPos, ignoreZdimension));
}

template <class T> template <class U>
float Vec3<T>::signedDistanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension)
{
    Vec3 boxDim = maxPos - minPos;
    Vec3 halfDim = boxDim * .5f;
    //    pos = pos - (minPos + halfDim); // ((maxPos - minPos) * .5f);
    //    Vec3 q = pos.abs() - halfDim; // - ((maxPos - minPos) * .5f);
    Vec3 q = (pos - (minPos + halfDim)).abs() - halfDim; // - ((maxPos - minPos) * .5f);
    if (ignoreZdimension)
        q.z() = std::numeric_limits<T>::lowest();

    static const Vec3 origin(0, 0, 0);
    float d = Vec3::max(q, origin).norm() + std::min(q.maxComp(), 0.f);
    return d;

}

template <class T> template <class U>
float Vec3<T>::distanceToBoundaries(const Vec3<T>& pos, const Vec3<U>& minPos, const Vec3<U>& maxPos, bool ignoreZdimension)
{
    return std::abs(Vec3::signedDistanceToBoundaries(pos, minPos, maxPos, ignoreZdimension));
}



template <class T>
nlohmann::json vec3_to_json(const Vec3<T>& vec);
template <class T>
Vec3<T> json_to_vec3(nlohmann::json json);

std::vector<float> json_to_color(nlohmann::json json);
nlohmann::json color_to_json(const std::vector<float>& color);
nlohmann::json color_to_json(const Vector3& color);






template <class T>
nlohmann::json vec3_to_json(const Vec3<T>& vec) {
    return nlohmann::json({{"x", vec.x()}, {"y", vec.y()}, {"z", vec.z()}});
}
template <class T>
Vec3<T> json_to_vec3(nlohmann::json json)
{
    return Vec3<T>(json.at("x").get<T>(), json.at("y").get<T>(), json.at("z").get<T>());
}


#endif // VECTOR3_H
