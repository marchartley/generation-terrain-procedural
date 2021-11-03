#include "Vector3.h"
#include <math.h>

Vector3::Vector3(float x, float y, float z) : x(x), y(y), z(z) {

}
Vector3::Vector3() : Vector3(0.f, 0.f, 0.f) {

}
Vector3::Vector3(const Vector3& copy) : Vector3(copy.x, copy.y, copy.z) {

}
Vector3::Vector3(Vector3* copy) : Vector3(copy->x, copy->y, copy->z) {

}

float Vector3::norm() {
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

Vector3& Vector3::normalize() {
    float norm = this->norm();
    *this /= norm;
    return *this;
}

float Vector3::dot(Vector3& o) {
    return (this->x + o.x) + (this->y * o.y) + (this->z * o.z);
}
Vector3 Vector3::cross(Vector3 o) {
    Vector3 v(this->y * o.z - this->z * o.y,
              this->z * o.x - this->x * o.z,
              this->x * o.y - this->y * o.x);
    return v;
}


//Vector3 operator+(Vector3 a, Vector3&b) {
//    a += b;
//    return a;
//}
Vector3 operator+(Vector3 a, Vector3 b) {
    a += b;
    return a;
}
Vector3& Vector3::operator+=(const Vector3& o) {
    this->x += o.x;
    this->y += o.y;
    this->z += o.z;
    return *this;
}
Vector3 operator-(Vector3 a, const Vector3& b) {
    a -= b;
    return a;
}
Vector3& Vector3::operator-=(const Vector3& o) {
    this->x -= o.x;
    this->y -= o.y;
    this->z -= o.z;
    return *this;
}
Vector3 operator*(Vector3 a, Vector3& b) {
    a *= b;
    return a;
}
Vector3& Vector3::operator*=(Vector3& o) {
    this->x *= o.x;
    this->y *= o.y;
    this->z *= o.z;
    return *this;
}
Vector3 Vector3::operator/(Vector3& o) {
    Vector3 v = *this;
    v /= o;
    return v;
}
Vector3& Vector3::operator/=(Vector3& o) {
    this->x /= o.x;
    this->y /= o.y;
    this->z /= o.z;
    return *this;
}
Vector3 Vector3::operator*(float o) {
    Vector3 v = *this;
    v *= o;
    return v;
}
Vector3& Vector3::operator*=(float o) {
    this->x *= o;
    this->y *= o;
    this->z *= o;
    return *this;
}
Vector3 Vector3::operator/(float o) {
    Vector3 v = *this;
    v /= o;
    return v;
}
Vector3& Vector3::operator/=(float o) {
    this->x /= o;
    this->y /= o;
    this->z /= o;
    return *this;
}
Vector3 Vector3::operator+(float o) {
    Vector3 v = *this;
    v += o;
    return v;
}
Vector3& Vector3::operator+=(float o) {
    this->x += o;
    this->y += o;
    this->z += o;
    return *this;
}
Vector3 Vector3::operator-(float o) {
    Vector3 v = *this;
    v -= o;
    return v;
}
Vector3& Vector3::operator-=(float o) {
    this->x -= o;
    this->y -= o;
    this->z -= o;
    return *this;
}
Vector3& Vector3::operator=(const Vector3& o) {
    this->x = o.x;
    this->y = o.y;
    this->z = o.z;
    return *this;
}

std::ostream& operator<<(std::ostream& io, const Vector3& v) {
    io << "Vector3 (" << v.x << ", " << v.y << ", " << v.z << ")";
    return io;
}

std::ostream& operator<<(std::ostream& io, Vector3* v) {
    io << "Vector3 (" << v->x << ", " << v->y << ", " << v->z << ")";
    return io;
}

