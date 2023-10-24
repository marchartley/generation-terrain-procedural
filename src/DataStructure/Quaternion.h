#ifndef QUATERNION_H
#define QUATERNION_H

#include "DataStructure/Vector3.h"

class Quaternion {
public:
    float w, x, y, z;

    Quaternion(float w, float x, float y, float z);

    // ...

    static Quaternion AxisAngle(const Vector3& axis, float angle);

    Quaternion operator*(const Quaternion& b) const;

    Vector3 rotate(const Vector3& v) const;

    static Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t);

    float dot(const Quaternion& b) const;

    Quaternion operator*(float s) const;
    Quaternion operator/(float s) const;
    Quaternion operator+(const Quaternion& b) const;

    Quaternion operator-(const Quaternion& b) const;
    Quaternion operator-() const;

    Vector3 toVector3() const;
    float norm() const;

    Quaternion normalized() const;

    Quaternion conjugate() const;

    Quaternion inverse() const;

    static Quaternion lerp(const Quaternion& q1, const Quaternion& q2, float t);
    static Quaternion squad(const Quaternion& q1, const Quaternion& q2, const Quaternion& a, const Quaternion& b, float t);
    Quaternion log() const;
    Quaternion exp() const;
    Vector3 toAxisAngle(float& angle) const;
    static Quaternion fromAxisAngle(const Vector3& axis, float angle);

    Matrix toRotationMatrix() const;
    Quaternion fromRotationMatrix(const Matrix& m);



};

Quaternion operator*(float s, const Quaternion& q);




#endif // QUATERNION_H
