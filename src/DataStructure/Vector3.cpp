#include "DataStructure/Vector3.h"
#include <math.h>
#include "Utils/Globals.h"

Vector3 Vector3::nabla = Vector3(1.f, 1.f, 1.f).normalize();

Vector3::Vector3(float x, float y, float z) : x(x), y(y), z(z) {

}
Vector3::Vector3() : Vector3(0.f, 0.f, 0.f) {

}
Vector3::Vector3(const Vector3& copy) : Vector3(copy.x, copy.y, copy.z) {

}
Vector3::Vector3(Vector3* copy) : Vector3(copy->x, copy->y, copy->z) {

}

float Vector3::norm() {
    if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}
float Vector3::norm2() {
    if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return this->x * this->x + this->y * this->y + this->z * this->z;
}

Vector3& Vector3::normalize() {
    if(this->norm() == 0)
        return *this;
    float norm = this->norm();
    *this /= norm;
    return *this;
}
Vector3 Vector3::normalized() const {
    Vector3 a = *this;
    return a.normalize();
}

Vector3 Vector3::abs()
{
    Vector3 a = *this;
    a.x = std::abs(x);
    a.y = std::abs(y);
    a.z = std::abs(z);
    return a;
}

float Vector3::dot(Vector3 o) {
    return (this->x * o.x) + (this->y * o.y) + (this->z * o.z);
}
Vector3 Vector3::cross(Vector3 o) {
    Vector3 v(this->y * o.z - this->z * o.y,
              this->z * o.x - this->x * o.z,
              this->x * o.y - this->y * o.x);
    return v;
}
Vector3 Vector3::rounded(int precision)
{
    Vector3 v = *this;
    v.x = (int)(v.x * pow(10, precision)) / (float)(pow(10, precision));
    v.y = (int)(v.y * pow(10, precision)) / (float)(pow(10, precision));
    v.z = (int)(v.z * pow(10, precision)) / (float)(pow(10, precision));
    return v;
}

Vector3 Vector3::random() {
    Vector3 v(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0));
    v.normalize();
    return v;
}

std::vector<float> Vector3::toArray(Vector3 v)
{
    std::vector<float> arr;
    arr.insert(arr.end(), {v.x, v.y, v.z});
    return arr;
}
std::vector<float> Vector3::toArray(std::vector<Vector3> vs)
{
    std::vector<float> arr;
    for (Vector3& v : vs)
        arr.insert(arr.end(), {v.x, v.y, v.z});
    return arr;
}
/*
Matrix3<Vector3> Vector3::gradient(Matrix3<float> field)
{
    Matrix3<Vector3> returningGrid(field.sizeX, field.sizeY, field.sizeZ);
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                Vector3 vec = Vector3();
                Vector3 allDirections;
                if (x == 0 || x == field.sizeX - 1 || y == 0 || y == field.sizeY - 1
                        || z == 0 || z == field.sizeZ - 1) {
                    returningGrid(x, y, z) = vec;
                } else {
                    for (int dx = std::max(x-1, 0); dx <= std::min(x+1, field.sizeX - 1); dx++) {
                        for (int dy = std::max(y-1, 0); dy <= std::min(y+1, field.sizeY - 1); dy++) {
                            for (int dz = std::max(z-1, 0); dz <= std::min(z+1, field.sizeZ - 1); dz++) {
                                if(dx != x || dy != y || dz != z) {
                                    Vector3 a = Vector3(dx - x, dy - y, dz - z);
//                                    allDirections += a;
                                    float f = field(dx, dy, dz);
                                    vec += a * f; //Vector3(dx - x, dy - y, dz - z) * field(dx, dy, dz);
                                }
                            }
                        }
                    }
                }

                returningGrid(x, y, z) = vec;

            }
        }
    }
    return returningGrid;
}

Matrix3<Vector3> Vector3::grad(Matrix3<float> field) {
    return Vector3::gradient(field);
}

Matrix3<float> Vector3::divergence(Matrix3<Vector3> field)
{
    Matrix3<float> returningGrid(field.sizeX, field.sizeY, field.sizeZ);
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                returningGrid(x, y, z) = field(x, y, z).divergence();
            }
        }
    }
    return returningGrid;
}

Matrix3<float> Vector3::div(Matrix3<Vector3> field)
{
    return Vector3::divergence(field);
}

Matrix3<Vector3> Vector3::curl(Matrix3<Vector3> field)
{
    Matrix3<Vector3> returningGrid(field.sizeX, field.sizeY, field.sizeZ);
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                Vector3& vec = field(x, y, z);
                returningGrid(x, y, z) = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
            }
        }
    }
    return returningGrid;
}

Matrix3<Vector3> Vector3::rot(Matrix3<Vector3> field)
{
    return Vector3::curl(field);
}

Matrix3<Vector3> Vector3::laplacian(Matrix3<Vector3> field)
{
    Matrix3<Vector3> returningGrid(field.sizeX, field.sizeY, field.sizeZ);
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                Vector3 vec = Vector3();
                if (x == 0 || x == field.sizeX - 1 || y == 0 || y == field.sizeY - 1
                        || z == 0 || z == field.sizeZ - 1) {
                    returningGrid(x, y, z) = vec;
                } else {
                    vec += field(x    , y    , z + 1);
                    vec += field(x    , y    , z - 1);
                    vec += field(x    , y + 1, z    );
                    vec += field(x    , y - 1, z    );
                    vec += field(x + 1, y    , z    );
                    vec += field(x - 1, y    , z    );
                    vec -= field(x    , y    , z    ) * 6;
                }
                returningGrid(x, y, z) = vec;
            }
        }
    }
    return returningGrid;
}

Matrix3<float> Vector3::laplacian(Matrix3<float> field)
{
    Matrix3<float> returningGrid(field.sizeX, field.sizeY, field.sizeZ);
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                float laplace = 0;
                if (x == 0 || x == field.sizeX - 1 || y == 0 || y == field.sizeY - 1
                        || z == 0 || z == field.sizeZ - 1) {
                    // Keep laplacian value to 0
                } else {
                    laplace += field(x    , y    , z + 1);
                    laplace += field(x    , y    , z - 1);
                    laplace += field(x    , y + 1, z    );
                    laplace += field(x    , y - 1, z    );
                    laplace += field(x + 1, y    , z    );
                    laplace += field(x - 1, y    , z    );
                    laplace -= field(x    , y    , z    ) * 6;
                }
                returningGrid(x, y, z) = laplace;
            }
        }
    }
    return returningGrid;
}
*/








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
Vector3 operator*(Vector3 a, Vector3 b) {
    a *= b;
    return a;
}
Vector3& Vector3::operator*=(Vector3 o) {
    this->x *= o.x;
    this->y *= o.y;
    this->z *= o.z;
    return *this;
}
Vector3 Vector3::operator/(Vector3 o) {
    Vector3 v = *this;
    v /= o;
    return v;
}
Vector3& Vector3::operator/=(Vector3 o) {
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
bool operator==(Vector3 a, Vector3 b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool operator!=(Vector3 a, Vector3 b)
{
    return !(a == b);
}

bool operator<(Vector3 a, Vector3 b)
{
    return a.norm2() < b.norm2();
}

bool operator<=(Vector3 a, Vector3 b)
{
    return !(a > b);
}

bool operator>(Vector3 a, Vector3 b)
{
    return a.norm2() > b.norm2();
}

bool operator>=(Vector3 a, Vector3 b)
{
    return !(a < b);
}

std::ostream& operator<<(std::ostream& io, const Vector3& v) {
    io << v.toString();
    return io;
}

std::ostream& operator<<(std::ostream& io, std::shared_ptr<Vector3> v) {
    io << v->toString();
    return io;
}

