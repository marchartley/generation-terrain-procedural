#ifndef VECTOR3_H
#define VECTOR3_H

#include "Utils/Globals.h"
#include <iostream>
#include <vector>
#include <QGLViewer/vec.h>
#include <QGLViewer/quaternion.h>
#include "DataStructure/Matrix.h"
//#include <glm/glm.hpp>

class Vector3 {
public:
    Vector3();
    Vector3(float x, float y, float z = 0.f, bool valid = true);
//    Vector3(const Vector3& copy);
//    Vector3(Vector3* copy);
    Vector3(qglviewer::Vec other);
    explicit Vector3(bool valid);
    explicit Vector3(const float* coords, bool valid = true);

    static std::vector<float> toArray(const Vector3& v);
    static std::vector<float> toArray(std::vector<Vector3> vs);
    std::tuple<int, int, int> toIntTuple() {return std::make_tuple<int, int, int>(int(this->x), int(this->y), int(this->z)); }

    static Vector3 min();
    static Vector3 max();
    static Vector3 min(const Vector3& a, const Vector3& b);
    static Vector3 max(const Vector3& a, const Vector3& b);
    static Vector3 min(std::vector<Vector3> allVectors);
    static Vector3 max(std::vector<Vector3> allVectors);

    static std::vector<Vector3> getAABBoxVertices(const Vector3& mini, const Vector3& maxi);

    friend std::ostream& operator<<(std::ostream& io, const Vector3& bbox);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<Vector3> bbox);

    float dot(const Vector3& o) const;
    Vector3 cross(const Vector3& o) const;
    Vector3 rounded(int precision = 0) const;
    Vector3 roundedUp(int precision = 0) const;
    Vector3 roundedDown(int precision = 0) const;
    Vector3 floor() const;
    Vector3 ceil() const;
    static Vector3 wrap(const Vector3& p, const Vector3& mini, const Vector3& maxi);

    float magnitude() const;
    float length() const;
    float lengthSquared() const;
    float norm() const;
    float norm2() const;
    Vector3& normalize();
    Vector3 normalized() const;
    Vector3 abs() const;
    Vector3& setMag(float newMag);
    Vector3& maxMagnitude(float maxMag);
    Vector3& minMagnitude(float minMag);
    Vector3& clamp(float minMag, float maxMag);
    Vector3 clamped(float minMag, float maxMag) const;

    bool isAlmostVertical();

    Vector3& rotate(float angle_x, float angle_y, float angle_z);
    Vector3& rotate(const Vector3& eulerAngles);
    Vector3 rotated(float angle_x, float angle_y, float angle_z) const;
    Vector3 rotated(const Vector3& eulerAngles) const;
    Vector3& rotate(float angle, float dir_x, float dir_y, float dir_z);
    Vector3& rotate(float angle, const Vector3& direction);
    Vector3 rotated(float angle, float dir_x, float dir_y, float dir_z) const;
    Vector3 rotated(float angle, const Vector3& direction) const;
    Vector3& translate(float move_x, float move_y, float move_z);
    Vector3& translate(const Vector3& move);
    Vector3 translated(float move_x, float move_y, float move_z);
    Vector3 translated(const Vector3& move);
    Vector3& applyTransform(Matrix transformMatrix);
//    Vector3& setDirection(const Vector3& dir, const Vector3& upVector = Vector3(0, 0, 1));
    Vector3 rotated90XY() const;

    Vector3& changeBasis(const Vector3& newX, const Vector3& newY, const Vector3& newZ);
    Vector3 changedBasis(const Vector3& newX, const Vector3& newY, const Vector3& newZ);

    Vector3 reflexion(const Vector3& normal);

    Vector3 toPolar();
    Vector3 fromPolar();

    float divergence() { return x + y + z; }

    static Vector3 fromMatrix(Matrix mat);
    Matrix toMatrix() const;
    Matrix toRotationMatrix() const;


    Vector3 toEulerAngles();
    Vector3 eulerAnglesWith(const Vector3& other);

    Vector3 getAllAnglesWith(const Vector3& otherVector) const;
    float getAngleWith(const Vector3 &otherVector) const;
    float getSignedAngleWith(const Vector3 &otherVector) const;
    float getSignedAngleAroundAxisWith(const Vector3 &otherVector, const Vector3& axis) const;

    static Vector3 quaternionToEuler(qglviewer::Quaternion quaternion);
    static Vector3 quaternionToEuler(float x, float y, float z, float w);

    static Vector3 slerp(float t, const Vector3& a, const Vector3& b);

    static Vector3 lerp(float t, const Vector3& min, const Vector3& max) {
        return min + (max - min) * t;
    }
    static float inverseLerp(const Vector3& val, const Vector3& min, const Vector3& max) {
        return (val - min).dot(max - min) / ((max - min).norm2() > 0 ? (max - min).norm2() : 0.0001f);
    }
    float inverseLerp(const Vector3& min, const Vector3& max) {
        return Vector3::inverseLerp(*this, min, max);
    }
    static float remap(const Vector3& val, const Vector3& oldMin, const Vector3& oldMax, float newMin, float newMax)
    {
        float oldProgress = Vector3::inverseLerp(val, oldMin, oldMax);
        return newMin + (newMax - newMin) * oldProgress;
    }
    float remap(const Vector3& oldMin, const Vector3& oldMax, float newMin, float newMax)
    {
        return Vector3::remap(*this, oldMin, oldMax, newMin, newMax);
    }

    static Vector3 random();
    static Vector3 random(float norm);
    static Vector3 random(float minNorm, float maxNorm);
    static Vector3 random(const Vector3& maxValues);
    static Vector3 random(const Vector3& minValues, const Vector3& maxValues);

    static Vector3 nabla;

    float maxComp() const { return std::max({x, y, z}); };
    float minComp() const { return std::min({x, y, z}); };


    static bool isInBox(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos);
    static float signedManhattanDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension = false);
    static float manhattanDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension = false);
    static float signedDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension = false);
    static float distanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension = false);

    bool isValid() const { return this->valid && (this->x == this->x && this->y == this->y && this->z == this->z); }
    void setValid(bool newValidValue) { this->valid = newValidValue; }
    static Vector3 invalid() { return Vector3(false); }
    operator qglviewer::Vec() const { return qglviewer::Vec(this->x, this->y, this->z); }
    explicit operator float*() const { return new float[3]{this->x, this->y, this->z}; }
//    operator glm::vec3() const { return glm::vec3(this->x, this->y, this->z); }
//    friend Vector3 operator+(const Vector3& a, Vector3& b);
    Vector3& operator+=(const Vector3& o);
    Vector3& operator-=(const Vector3& o);
    Vector3& operator*=(const Vector3& o);
//    Vector3 operator/(const Vector3& o);
    Vector3& operator/=(const Vector3& o);
//    Vector3 operator*(float o) const;
    Vector3& operator*=(float o);
//    Vector3 operator/(float o);
    Vector3& operator/=(float o);
//    Vector3 operator+(float o);
    Vector3& operator+=(float o);
//    Vector3 operator-(float o);
    Vector3& operator-=(float o);
    friend Vector3 operator+(const Vector3& a, const Vector3& b);
    friend Vector3 operator-(const Vector3& a, const Vector3& b);
    friend Vector3 operator*(const Vector3& a, const Vector3& o);
    friend Vector3 operator/(const Vector3& a, const Vector3& o);

//    friend Vector3 operator+(float a, const Vector3& b);
//    friend Vector3 operator-(float a, const Vector3& b);
    friend Vector3 operator*(float a, const Vector3& b);
    friend Vector3 operator/(float a, const Vector3& b);

//    friend Vector3 operator+(const Vector3& b, float a);
//    friend Vector3 operator-(const Vector3& b, float a);
    friend Vector3 operator*(const Vector3& b, float a);
    friend Vector3 operator/(const Vector3& b, float a);

    friend Vector3 operator-(const Vector3& v);

//    Vector3& operator=(const Vector3& o);
    friend bool operator==(const Vector3& a, const Vector3& b);
    friend bool operator!=(const Vector3& a, const Vector3& b);
    friend bool operator<(const Vector3& a, const Vector3& b);
    friend bool operator<=(const Vector3& a, const Vector3& b);
    friend bool operator>(const Vector3& a, const Vector3& b);
    friend bool operator>=(const Vector3& a, const Vector3& b);
//    using std::abs;
//    friend Vector3 abs(const Vector3& o) { return o.abs(); }

    float& operator[](size_t i);
    const float& operator[](size_t i) const;

    std::string toString() const {return "Vector3 (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")" + (this->isValid() ? "" : "[INVALID]"); }
//    const char* toHashString() const {return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z); }

    Vector3 xxx() const { return Vector3(x, x, x); }
    Vector3 xx()  const { return Vector3(x, x, 0); }
    Vector3 yyy() const { return Vector3(y, y, y); }
    Vector3 yy()  const { return Vector3(y, y, 0); }
    Vector3 xy()  const { return Vector3(x, y, 0); }
    Vector3 yz()  const { return Vector3(y, z, 0); }
    Vector3 xz()  const { return Vector3(x, z, 0); }
    Vector3 xyz() const { return Vector3(x, y, z); }
    Vector3 yxz() const { return Vector3(y, x, z); }
    Vector3 zyx() const { return Vector3(z, y, x); }
    Vector3 xzy() const { return Vector3(x, z, y); }
    float x, y, z;
protected:
    bool valid = true;

};

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
namespace std {
    Vector3 abs(const Vector3& o);
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

class AABBox { //: public std::pair<Vector3, Vector3> {
public:
    AABBox();
    AABBox(const Vector3& mini, const Vector3& maxi);
    AABBox(std::tuple<Vector3, Vector3> minMax);
    AABBox(std::vector<Vector3> allPointsToContain);
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

#include "Utils/json.h"
nlohmann::json vec3_to_json(const Vector3& vec);
Vector3 json_to_vec3(nlohmann::json json);

std::vector<float> json_to_color(nlohmann::json json);
nlohmann::json color_to_json(const std::vector<float>& color);
nlohmann::json color_to_json(const Vector3& color);

#endif // VECTOR3_H
