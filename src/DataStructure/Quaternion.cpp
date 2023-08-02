#include "Quaternion.h"

//Quaternion::Quaternion()
//{

//}

Quaternion operator*(float s, const Quaternion &q) {
    return q * s;
}

Quaternion::Quaternion(float w, float x, float y, float z)
    : w(w), x(x), y(y), z(z) {}

Quaternion Quaternion::AxisAngle(const Vector3 &axis, float angle) {
    float halfAngle = angle * 0.5f;
    float s = std::sin(halfAngle);
    return Quaternion(std::cos(halfAngle), axis.x * s, axis.y * s, axis.z * s);
}

Vector3 Quaternion::rotate(const Vector3 &v) const {
    Quaternion qv(0, v.x, v.y, v.z);
    Quaternion rotated = (*this) * qv * conjugate();
    return rotated.toVector3();
}

Quaternion Quaternion::slerp(const Quaternion &q1, const Quaternion &q2, float t) {
    float dot = q1.dot(q2);
    if (dot < 0) {
        return q1 * (-1) + t * ((-q2) - q1 * (-1));
    } else {
        return q1 + t * (q2 - q1);
    }
}

float Quaternion::dot(const Quaternion &b) const {
    return w*b.w + x*b.x + y*b.y + z*b.z;
}

Quaternion Quaternion::operator*(float s) const {
    return Quaternion(w*s, x*s, y*s, z*s);
}

Quaternion Quaternion::operator/(float s) const
{
    return *this * (1.f / s);
}

Quaternion Quaternion::operator+(const Quaternion &b) const {
    return Quaternion(w + b.w, x + b.x, y + b.y, z + b.z);
}

Quaternion Quaternion::operator-(const Quaternion &b) const {
    return Quaternion(w - b.w, x - b.x, y - b.y, z - b.z);
}

Quaternion Quaternion::operator-() const {
    return *this * -1.f;
}

Vector3 Quaternion::toVector3() const {
    return Vector3(x, y, z);
}

float Quaternion::norm() const {
    return std::sqrt(w*w + x*x + y*y + z*z);
}

Quaternion Quaternion::normalized() const {
    float n = norm();
    return Quaternion(w/n, x/n, y/n, z/n);
}

Quaternion Quaternion::conjugate() const {
    return Quaternion(w, -x, -y, -z);
}

Quaternion Quaternion::inverse() const {
    float n2 = norm() * norm();
    return conjugate() / n2;
}

Quaternion Quaternion::lerp(const Quaternion &q1, const Quaternion &q2, float t) {
    return (1 - t) * q1 + t * q2;
}

Quaternion Quaternion::squad(const Quaternion &q1, const Quaternion &q2, const Quaternion &a, const Quaternion &b, float t) {
    Quaternion c = slerp(q1, q2, t);
    Quaternion d = slerp(a, b, t);
    return slerp(c, d, 2.0f * t * (1.0f - t));
}

Quaternion Quaternion::log() const {
    float a = std::acos(w);
    float s = std::sin(a);

    if (std::abs(s) < std::numeric_limits<float>::epsilon())
        return Quaternion(0, x, y, z);
    else
        return Quaternion(0, a * x / s, a * y / s, a * z / s);
}

Quaternion Quaternion::exp() const {
    float a = std::sqrt(x * x + y * y + z * z);
    float s = std::sin(a);

    if (std::abs(s) < std::numeric_limits<float>::epsilon())
        return Quaternion(std::cos(a), x, y, z);
    else
        return Quaternion(std::cos(a), s * x / a, s * y / a, s * z / a);
}

void Quaternion::toAxisAngle(Vector3 &axis, float &angle) const {
    angle = 2 * std::acos(w);
    float s = std::sqrt(1 - w*w);

    if (s < 0.001) {
        axis.x = x;
        axis.y = y;
        axis.z = z;
    } else {
        axis.x = x / s;
        axis.y = y / s;
        axis.z = z / s;
    }
}

Quaternion Quaternion::fromAxisAngle(const Vector3 &axis, float angle) {
    float s = std::sin(angle / 2);
    return Quaternion(std::cos(angle / 2), axis.x * s, axis.y * s, axis.z * s);
}

Quaternion Quaternion::operator*(const Quaternion &b) const {
    return Quaternion(
                w * b.w - x * b.x - y * b.y - z * b.z,  // new w
                w * b.x + x * b.w + y * b.z - z * b.y,  // new x
                w * b.y - x * b.z + y * b.w + z * b.x,  // new y
                w * b.z + x * b.y - y * b.x + z * b.w   // new z
                );
}
