#include "DataStructure/Vector3.h"
#include <math.h>
#include "Utils/Globals.h"
#include "Utils/Utils.h"

Vector3 Vector3::nabla = Vector3(1.f, 1.f, 1.f).normalize();

Vector3::Vector3(float x, float y, float z, bool valid) : x(x), y(y), z(z), valid(valid) {

}
Vector3::Vector3() : Vector3(0.f, 0.f, 0.f) {

}/*
Vector3::Vector3(const Vector3& copy) : Vector3(copy.x, copy.y, copy.z, copy.valid) {

}
Vector3::Vector3(Vector3* copy) : Vector3(copy->x, copy->y, copy->z, copy->valid) {

}*/

Vector3::Vector3(qglviewer::Vec other)
    : Vector3(other.x, other.y, other.z)
{

}

Vector3::Vector3(bool valid) : Vector3()
{
    this->valid = valid;
}

Vector3::Vector3(const float *coords, bool valid)
    : Vector3(coords[0], coords[1], coords[2], valid)
{

}

float Vector3::norm() const {
    if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}
float Vector3::norm2() const {
    if(this->x == 0 && this->y == 0 && this->z == 0) return 0;
    return this->x * this->x + this->y * this->y + this->z * this->z;
}


Vector3& Vector3::normalize() {
    if(this->norm() < 1e-5)
        return *this;
    float norm = this->norm();
    *this /= norm;
    return *this;
}
Vector3 Vector3::normalized() const {
    Vector3 a = *this;
    return a.normalize();
}

Vector3 Vector3::abs() const
{
    Vector3 a = *this;
    a.x = std::abs(x);
    a.y = std::abs(y);
    a.z = std::abs(z);
    return a;
}

Vector3& Vector3::setMag(float newMag)
{
    this->normalize() *= newMag;
    return *this;
}

Vector3& Vector3::maxMagnitude(float maxMag)
{
    if (this->norm2() > maxMag*maxMag)
        this->setMag(maxMag);
    return *this;
}

Vector3 &Vector3::minMagnitude(float minMag)
{
    if (this->norm2() < minMag*minMag)
        this->setMag(minMag);
    return *this;
}

Vector3 &Vector3::clamp(float minMag, float maxMag)
{
    float mag2 = this->norm2();
    if (mag2 > maxMag*maxMag)
        this->setMag(maxMag);
    else if (mag2 < minMag*minMag)
        this->setMag(minMag);
    return *this;
}

Vector3 Vector3::clamped(float minMag, float maxMag) const
{
    Vector3 v = *this;
    return v.clamp(minMag, maxMag);
}

bool Vector3::isAlmostVertical()
{
    return std::abs(this->dot(Vector3(0, 0, 1))) > 0.999;
}

Matrix Vector3::toMatrix() const
{
    return Matrix(3, 1, (float*)(*this));
}

Matrix Vector3::toRotationMatrix() const
{
    Matrix Rx (3, 3, std::vector<float>({
                   1, 0, 0,
                   0, cos(this->x), -sin(this->x),
                   0, sin(this->x), cos(this->x)
               }).data());
    Matrix Rz (3, 3, std::vector<float>({
                    cos(this->y), 0, -sin(this->y),
                    0, 1, 0,
                    sin(this->y), 0, cos(this->y)
                }).data());
    Matrix Ry (3, 3, std::vector<float>({
                    cos(this->z), -sin(this->z), 0,
                    sin(this->z), cos(this->z), 0,
                    0, 0, 1
                }).data());
    return Rx.product(Ry).product(Rz);
}

Vector3 Vector3::toEulerAngles()
{
    Vector3 self = this->normalized();
    return Vector3(
                    std::acos(self.y),
                    std::acos(self.z),
                    std::acos(self.x)
                );
}

Vector3 Vector3::eulerAnglesWith(const Vector3& other)
{
    return Vector3(false);
    /*Vector3 v1 = this->normalized();
    Vector3 v2 = other.normalized();
    Vector3 v3 = v1.cross(v2);
    v2 = v3.cross(v1); // make it perpendicular

    float R11 = v1.x, R12 = v2.x, R13 = v3.x,
            R21 = v1.y, R22 = v2.y, R23 = v3.y,
            R31 = v1.z, R32 = v2.z, R33 = v3.z;
    Vector3 result;
    result.x = std::atan2(R32, R33) * (180.0 / M_PI);
    result.y = std::atan2(-1 * R31, std::sqrt(R32 * R32 + R33 * R33)) * (180.0 / M_PI);
    result.z = std::atan2(R21, R11) * (180.0 / M_PI);
    return result;*/
}

Vector3 Vector3::getAllAnglesWith(const Vector3 &otherVector) const
{
    if (*this == otherVector)
        return Vector3(0.f, 0.f, 0.f);
    float onX = this->getSignedAngleAroundAxisWith(otherVector, Vector3(1, 0, 0));
    float onY = this->getSignedAngleAroundAxisWith(otherVector, Vector3(0, 1, 0));
    float onZ = this->getSignedAngleAroundAxisWith(otherVector, Vector3(0, 0, 1));
    return Vector3(onX, onY, onZ);
}

float Vector3::getAngleWith(const Vector3& otherVector) const
{
    Vector3 vA = this->normalized();
    Vector3 vB = otherVector.normalized();
    return std::acos(vA.dot(vB));
}

float Vector3::getSignedAngleWith(const Vector3 &otherVector) const
{
    return this->getAngleWith(otherVector) * sign(this->cross(otherVector).z);
}

float Vector3::getSignedAngleAroundAxisWith(const Vector3 &otherVector, const Vector3 &axis) const
{
    Vector3 normalizedAxis = axis.normalized();
    Vector3 vA = *this - normalizedAxis * this->dot(normalizedAxis); //this->normalized();
    Vector3 vB = otherVector - normalizedAxis * otherVector.dot(normalizedAxis); //.normalized();
    if (vA == vB)
        return 0.f;
//    float dot = vA.dot(vB);
    Vector3 cross = vA.cross(vB);
    if (cross.norm2() < 1e-5)
        return 0.f;
    float angle = vA.getAngleWith(vB);
    if (cross.dot(axis) < 0)
        angle = -angle;
    return angle;
}

Vector3 Vector3::quaternionToEuler(qglviewer::Quaternion quaternion)
{
    return Vector3::quaternionToEuler(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
}

Vector3 Vector3::quaternionToEuler(float x, float y, float z, float w)
{
    Vector3 angles;

    // roll (x-axis rotation)
    float sinr_cosp = 2 * (w * x + y * z);
    float cosr_cosp = 1 - 2 * (x * x + y * y);
    angles.x = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    float sinp = std::sqrt(1 + 2 * (w * y - x * z));
    float cosp = std::sqrt(1 - 2 * (w * y - x * z));
    angles.y = 2 * std::atan2(sinp, cosp) - M_PI / 2;

    // yaw (z-axis rotation)
    float siny_cosp = 2 * (w * z + x * y);
    float cosy_cosp = 1 - 2 * (y * y + z * z);
    angles.z = std::atan2(siny_cosp, cosy_cosp);

    return angles;
}

Vector3 Vector3::slerp(float t, const Vector3 &A, const Vector3 &B) {
    float lengthA = A.norm();
    float lengthB = B.norm();
    auto start = A / lengthA, end = B / lengthB;
    float dot = start.dot(end);
    dot = std::clamp(dot, -1.0f, 1.0f);

    float theta = std::acos(dot) * t;
    Vector3 relativeVec = end - start * dot;
    relativeVec.normalize(); // Ensure it's a unit vector

    return ((start * std::cos(theta)) + (relativeVec * std::sin(theta))) * interpolation::inv_linear(t, lengthA, lengthB);
}

float Vector3::dot(const Vector3& o) const {
    return (this->x * o.x) + (this->y * o.y) + (this->z * o.z);
}
Vector3 Vector3::cross(const Vector3& o) const {
    Vector3 v(this->y * o.z - this->z * o.y,
              this->z * o.x - this->x * o.z,
              this->x * o.y - this->y * o.x);
    return v;
}
Vector3 Vector3::rounded(int precision) const
{
    return this->roundedDown(precision);
    /*Vector3 v = *this;
    v.x = (int)(v.x * pow(10, precision)) / (float)(pow(10, precision));
    v.y = (int)(v.y * pow(10, precision)) / (float)(pow(10, precision));
    v.z = (int)(v.z * pow(10, precision)) / (float)(pow(10, precision));
    return v;*/
}

Vector3 Vector3::roundedUp(int precision) const
{
    float power = std::pow(10, precision);
    return ((*this) * power).ceil() / power;
}

Vector3 Vector3::roundedDown(int precision) const
{
    float power = std::pow(10, precision);
    return ((*this) * power).floor() / power;
}
Vector3 Vector3::floor() const
{
    Vector3 v = *this;
    v.x = std::floor(v.x);
    v.y = std::floor(v.y);
    v.z = std::floor(v.z);
    return v;
}
Vector3 Vector3::ceil() const
{
    Vector3 v = this->floor();
    if (v.x != x) v.x += 1;
    if (v.y != y) v.y += 1;
    if (v.z != z) v.z += 1;
    return v;
}

Vector3 Vector3::wrap(const Vector3& p, const Vector3& mini, const Vector3& maxi)
{
    Vector3 newP = p - mini;
    Vector3 newMaxi = maxi - mini;
//    p -= mini;
//    maxi -= mini;
    Vector3 rounded = newP.roundedDown();
    Vector3 decimals = newP - rounded;
    Vector3 wrap = Vector3(int(rounded.x + newMaxi.x) % int(newMaxi.x),
                           int(rounded.y + newMaxi.y) % int(newMaxi.y),
                           int(rounded.z + newMaxi.z) % int(newMaxi.z)
                           ) + decimals + mini;
    return wrap;
}

float Vector3::magnitude() const
{
    return this->norm();
}

float Vector3::length() const
{
    return this->norm();
}

Vector3 Vector3::random() {
    Vector3 v(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0));
    v.normalize();
    return v;
}

Vector3 Vector3::random(float norm)
{
    return Vector3::random() * norm;
}

Vector3 Vector3::random(float minNorm, float maxNorm)
{
    return Vector3::random() * random_gen::generate(minNorm, maxNorm);
}

Vector3 Vector3::random(const Vector3& maxValues)
{
    return Vector3::random(Vector3(), maxValues);
}

Vector3 Vector3::random(const Vector3& minValues, const Vector3& maxValues)
{
    return Vector3(random_gen::generate(minValues.x, maxValues.x),
                   random_gen::generate(minValues.y, maxValues.y),
                   random_gen::generate(minValues.z, maxValues.z));
}

std::vector<float> Vector3::toArray(const Vector3& v)
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

Vector3 Vector3::min()
{
    return Vector3(std::numeric_limits<float>::lowest(),
                   std::numeric_limits<float>::lowest(),
                   std::numeric_limits<float>::lowest());
}

Vector3 Vector3::max()
{
    return Vector3(std::numeric_limits<float>::max(),
                   std::numeric_limits<float>::max(),
                   std::numeric_limits<float>::max());
}

Vector3 Vector3::min(const Vector3& a, const Vector3& b)
{
    if (!a.isValid()) return b;
    if (!b.isValid()) return a;
    return Vector3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

Vector3 Vector3::max(const Vector3& a, const Vector3& b)
{
    if (!a.isValid()) return b;
    if (!b.isValid()) return a;
    return Vector3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

Vector3 Vector3::min(std::vector<Vector3> allVectors)
{
    if (allVectors.empty())
        return Vector3(false);
    Vector3 res = allVectors[0];
    for (size_t i = 1; i < allVectors.size(); i++)
        res = Vector3::min(res, allVectors[i]);
    return res;
}

Vector3 Vector3::max(std::vector<Vector3> allVectors)
{
    if (allVectors.empty())
        return Vector3(false);
    Vector3 res = allVectors[0];
    for (size_t i = 1; i < allVectors.size(); i++)
        res = Vector3::max(res, allVectors[i]);
    return res;
}

std::vector<Vector3> Vector3::getAABBoxVertices(const Vector3& mini, const Vector3& maxi)
{
    float minX = mini.x, maxX = maxi.x,
            minY = mini.y, maxY = maxi.y,
            minZ = mini.z, maxZ = maxi.z;
    return {
        Vector3(minX, minY, minZ),
        Vector3(minX, minY, maxZ),
        Vector3(minX, maxY, minZ),
        Vector3(minX, maxY, maxZ),
        Vector3(maxX, minY, minZ),
        Vector3(maxX, minY, maxZ),
        Vector3(maxX, maxY, minZ),
        Vector3(maxX, maxY, maxZ)
    };
}


Vector3& Vector3::rotate(float angle_x, float angle_y, float angle_z) {
    return this->rotate(Vector3(angle_x, angle_y, angle_z));
}
Vector3& Vector3::rotate(const Vector3& eulerAngles) {
    return this->applyTransform(eulerAngles.toRotationMatrix());
    /*Matrix newCoords = R.product(this->toMatrix());
    this->x = newCoords[0][0];
    this->y = newCoords[1][0];
    this->z = newCoords[2][0];
    return *this;*/
}
Vector3 Vector3::rotated(float angle_x, float angle_y, float angle_z) const {
    return this->rotated(Vector3(angle_x, angle_y, angle_z));
}
Vector3 Vector3::rotated(const Vector3& eulerAngles) const {
    Vector3 v = *this;
    return v.rotate(eulerAngles);
}
Vector3& Vector3::rotate(float angle, float dir_x, float dir_y, float dir_z) {
    return this->rotate(angle, Vector3(dir_x, dir_y, dir_z));
}
Vector3& Vector3::rotate(float angle, const Vector3& direction) {
    float c = cos(angle), s = sin(angle);
    Vector3 v = direction.normalized(); // alias
    Matrix R (3, 3, std::vector<float>({
                   v.x*v.x*(1-c)+c, v.x*v.y*(1-c)-v.z*s, v.x*v.z*(1-c)+v.y*s,
                   v.x*v.y*(1-c)+v.z*s, v.y*v.y*(1-c)+c, v.y*v.z*(1-c)-v.x*s,
                   v.x*v.z*(1-c)-v.y*s, v.y*v.z*(1-c)+v.x*s, v.z*v.z*(1-c)+c
               }).data());
    return this->applyTransform(R);
    /*Matrix newCoords = R.product(this->toMatrix());
    this->x = newCoords[0][0];
    this->y = newCoords[1][0];
    this->z = newCoords[2][0];
    return *this;*/
}
Vector3 Vector3::rotated(float angle, float dir_x, float dir_y, float dir_z) const {
    return this->rotated(angle, Vector3(dir_x, dir_y, dir_z));
}
Vector3 Vector3::rotated(float angle, const Vector3& direction) const {
    Vector3 v = *this;
    return v.rotate(angle, direction);
}
Vector3& Vector3::translate(float move_x, float move_y, float move_z) {
    return this->translate(Vector3(move_x, move_y, move_z));
}
Vector3& Vector3::translate(const Vector3& move) {
    return (*this += move);
}
Vector3 Vector3::translated(float move_x, float move_y, float move_z) {
    return this->translated(Vector3(move_x, move_y, move_z));
}
Vector3 Vector3::translated(const Vector3& move) {
    Vector3 v = *this;
    return v.translate(move);
}

Vector3& Vector3::applyTransform(Matrix transformMatrix)
{
    Matrix newCoords = transformMatrix.product(this->toMatrix());
    this->x = newCoords[0][0];
    this->y = newCoords[1][0];
    this->z = newCoords[2][0];
    return *this;
}

Vector3 &Vector3::changeBasis(const Vector3& newX, const Vector3& newY, const Vector3& newZ)
{
    Vector3 newVec = this->changedBasis(newX, newY, newZ);
    this->x = newVec.x;
    this->y = newVec.y;
    this->z = newVec.z;
    return *this;
}

Vector3 Vector3::changedBasis(const Vector3& newX, const Vector3& newY, const Vector3& newZ)
{
    Vector3 newVec;
    newVec += newX * this->x + newY * this->y + newZ * this->z;
    return newVec;
}

Vector3 Vector3::reflexion(const Vector3& normal)
{
    float dot = normal.dot(*this);
    Vector3 n2d = normal * 2.f * dot;
    return *this - n2d;
}

Vector3 Vector3::fromMatrix(Matrix mat)
{
    if (mat.size() == 1) {
        return Vector3(mat[0][0], mat[0][1], mat[0][2]);
    } else if (mat[0].size() == 1) {
        return Vector3(mat[0][0], mat[1][0], mat[2][0]);
    } else {
        throw std::domain_error("Cannot transform Matrix " + mat.toString() + " to Vector3");
    }
}


Vector3 operator+(const Vector3& a, const Vector3& b) {
    Vector3 res = a;
    res += b;
    return res;
}
Vector3& Vector3::operator+=(const Vector3& o) {
    this->x += o.x;
    this->y += o.y;
    this->z += o.z;
    return *this;
}
Vector3 operator-(const Vector3& a, const Vector3& b) {
    Vector3 res = a;
    res -= b;
    return res;
}
Vector3& Vector3::operator-=(const Vector3& o) {
    this->x -= o.x;
    this->y -= o.y;
    this->z -= o.z;
    return *this;
}

float &Vector3::operator[](size_t i)
{
    return ((float*)(this))[i];
}
const float &Vector3::operator[](size_t i) const
{
    return ((float*)(this))[i];
}
//Vector3 operator*(const Vector3& a, const Vector3& b) {
//    a *= b;
//    return a;
//}
//Vector3 operator/(const Vector3& a, const Vector3& b) {
//    a /= b;
//    return a;
//}
Vector3& Vector3::operator*=(const Vector3& o) {
    this->x *= o.x;
    this->y *= o.y;
    this->z *= o.z;
    return *this;
}
//Vector3 Vector3::operator/(const Vector3& o) {
//    Vector3 v = *this;
//    v /= o;
//    return v;
//}
Vector3& Vector3::operator/=(const Vector3& o) {
    this->x /= o.x;
    this->y /= o.y;
    this->z /= o.z;
    return *this;
}/*
Vector3 operator/(const Vector3& a, const Vector3& b) {
    a /= b;
    return a;
}*/
//Vector3 Vector3::operator*(float o) const {
//    Vector3 v = *this;
//    v *= o;
//    return v;
//}
Vector3& Vector3::operator*=(float o) {
    this->x *= o;
    this->y *= o;
    this->z *= o;
    return *this;
}
//Vector3 Vector3::operator/(float o) {
//    Vector3 v = *this;
//    v /= o;
//    return v;
//}
Vector3& Vector3::operator/=(float o) {
    this->x /= o;
    this->y /= o;
    this->z /= o;
    return *this;
}
//Vector3 Vector3::operator+(float o) {
//    Vector3 v = *this;
//    v += o;
//    return v;
//}
Vector3& Vector3::operator+=(float o) {
    this->x += o;
    this->y += o;
    this->z += o;
    return *this;
}
//Vector3 Vector3::operator-(float o) {
//    Vector3 v = *this;
//    v -= o;
//    return v;
//}
Vector3& Vector3::operator-=(float o) {
    this->x -= o;
    this->y -= o;
    this->z -= o;
    return *this;
}
/*Vector3& Vector3::operator=(const Vector3& o) {
    this->x = o.x;
    this->y = o.y;
    this->z = o.z;
    this->valid = o.valid;
    return *this;
}*/

Vector3 operator/(const Vector3& a, const Vector3& b) {
    Vector3 res = a;
    res /= b;
    return res;
}
Vector3 operator*(const Vector3& a, const Vector3& b) {
    Vector3 res = a;
    res *= b;
    return res;
}
//Vector3 operator+(float a, const Vector3& b) {
//    return b + a;
//}
//Vector3 operator-(float a, const Vector3& b) {
//    return b * -1.f + a;
//}
Vector3 operator*(float a, const Vector3& b) {
    Vector3 res = b;
    res *= a;
    return res;
}

Vector3 operator/(float a, const Vector3& b) {
    return Vector3(a / b.x, a / b.y, a / b.z);
}

//Vector3 operator+(const Vector3& b, float a) {
//    return b + a;
//}
//Vector3 operator-(const Vector3& b, float a) {
//    return b * -1.f + a;
//}
Vector3 operator*(const Vector3& b, float a) {
    Vector3 res = b;
    res *= a;
    return res;
}
Vector3 operator/(const Vector3& b, float a) {
    Vector3 res = b;
    res /= a;
    return res;
}

bool operator==(const Vector3& a, const Vector3& b)
{
    float epsilon = 1e-8;
    return std::abs(a.x - b.x) < epsilon && std::abs(a.y - b.y) < epsilon && std::abs(a.z - b.z) < epsilon; //int(a.x * 1e8) == int(b.x * 1e8) && int(a.y * 1e8) == int(b.y * 1e8) && int(a.z * 1e8) == int(b.z * 1e8);
}

Vector3 operator-(const Vector3& v) {
    return v * -1.f;
}

bool operator!=(const Vector3& a, const Vector3& b)
{
    return !(a == b);
}

bool operator<(const Vector3& a, const Vector3& b)
{
    return a.norm2() < b.norm2();
}

bool operator<=(const Vector3& a, const Vector3& b)
{
    return !(a > b);
}

bool operator>(const Vector3& a, const Vector3& b)
{
    return a.norm2() > b.norm2();
}

bool operator>=(const Vector3& a, const Vector3& b)
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

Vector3 std::abs(const Vector3& o) { return o.abs(); }

bool Vector3::isInBox(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos) {
    return (minPos.x <= pos.x && pos.x <= maxPos.x) && (minPos.y <= pos.y && pos.y <= maxPos.y) && (minPos.z <= pos.z && pos.z <= maxPos.z);
//    return (pos - minPos).minComp() >= 0.f && (pos - (minPos + maxPos)).maxComp() <= 0.f;
}

// Inside : positive
// Outside: negative
float Vector3::signedManhattanDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension)
{
    Vector3 newPos = pos - minPos;
    Vector3 newMax = maxPos - minPos;

    if (ignoreZdimension) {
        newPos = newPos.xy();
        newMax = newMax.xy();
    }
    if (Vector3::isInBox(newPos, Vector3(), newMax)) {
        if (ignoreZdimension)
            return -std::min({newPos.x, newPos.y, newMax.x - newPos.x, newMax.y - newPos.y});
        else
            return -std::min({newPos.x, newPos.y, newPos.z, newMax.x - newPos.x, newMax.y - newPos.y, newMax.z - newPos.z});
    } else {
        if (ignoreZdimension)
            return -std::min({newPos.x, newPos.y, newMax.x - newPos.x, newMax.y - newPos.y});
        else
            return -std::min({newPos.x, newPos.y, newPos.z, newMax.x - newPos.x, newMax.y - newPos.y, newMax.z - newPos.z});
    }
}

float Vector3::manhattanDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension)
{
    return std::abs(Vector3::signedManhattanDistanceToBoundaries(pos, minPos, maxPos, ignoreZdimension));
}

float Vector3::signedDistanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension)
{
    Vector3 boxDim = maxPos - minPos;
    Vector3 halfDim = boxDim * .5f;
//    pos = pos - (minPos + halfDim); // ((maxPos - minPos) * .5f);
//    Vector3 q = pos.abs() - halfDim; // - ((maxPos - minPos) * .5f);
    Vector3 q = (pos - (minPos + halfDim)).abs() - halfDim; // - ((maxPos - minPos) * .5f);
    if (ignoreZdimension)
        q.z = 0.f;
    float d = Vector3::max(q, Vector3(0, 0, 0)).norm() + std::min(q.maxComp(), 0.f);
    return d;

}

float Vector3::distanceToBoundaries(const Vector3& pos, const Vector3& minPos, const Vector3& maxPos, bool ignoreZdimension)
{
    return std::abs(Vector3::signedDistanceToBoundaries(pos, minPos, maxPos, ignoreZdimension));
}



AABBox::AABBox()
    : AABBox(Vector3(0, 0, 0), Vector3(0, 0, 0))
{

}

AABBox::AABBox(const Vector3& mini, const Vector3& maxi) : mini(mini), maxi(maxi)
{

}

AABBox::AABBox(std::tuple<Vector3, Vector3> minMax)
    : AABBox(std::get<0>(minMax), std::get<1>(minMax))
{

}

AABBox::AABBox(std::vector<Vector3> allPointsToContain)
    : AABBox(Vector3::min(allPointsToContain), Vector3::max(allPointsToContain))
{

}

Vector3 AABBox::normalize(const Vector3& p)
{
    auto mini = this->min();
    auto maxi = this->max();
    Vector3 ret = (p - mini) / (maxi - mini);
    return ret;
}

AABBox &AABBox::expand(const Vector3 &newPoint)
{
    this->mini = Vector3::min(this->mini, newPoint);
    this->maxi = Vector3::max(this->maxi, newPoint);
    return *this;
}

AABBox& AABBox::expand(const std::vector<Vector3>& newPoints)
{
    for (const auto& p : newPoints)
        this->expand(p);
    return *this;
}

float AABBox::distanceTo(const Vector3 &p)
{
    return Vector3::distanceToBoundaries(p, this->min(), this->max());
}

Vector3 AABBox::random(const Vector3& mini, const Vector3& maxi)
{
    return Vector3::random(mini, maxi);
}

Vector3 AABBox::intersects(const Vector3& rayStart, const Vector3& rayEnd) {
    if (this->contains(rayStart)) return rayStart;

    Vector3 direction = (rayEnd - rayStart).normalized();
    Vector3 tMin = (min() - rayStart) / direction;
    Vector3 tMax = (max() - rayStart) / direction;

    if (tMin.x > tMax.x) std::swap(tMin.x, tMax.x);
    if (tMin.y > tMax.y) std::swap(tMin.y, tMax.y);
    if (tMin.z > tMax.z) std::swap(tMin.z, tMax.z);

    float tNear = std::max({tMin.x, tMin.y, tMin.z});
    float tFar = std::min({tMax.x, tMax.y, tMax.z});

    if (tNear > tFar || tFar < 0) {
        return Vector3::invalid();  // return invalid Vector3 if no intersection
    }

    Vector3 intersectionPoint = rayStart + direction * tNear; // calculate intersection point
    return intersectionPoint;
}

std::ostream& operator<<(std::ostream& io, const AABBox& bbox) {
    io << bbox.toString();
    return io;
}

std::ostream& operator<<(std::ostream& io, std::shared_ptr<AABBox> bbox) {
    io << bbox->toString();
    return io;
}





nlohmann::json vec3_to_json(const Vector3& vec) {
    return nlohmann::json({{"x", vec.x}, {"y", vec.y}, {"z", vec.z}});
}

Vector3 json_to_vec3(nlohmann::json json)
{
    return Vector3(json.at("x").get<float>(), json.at("y").get<float>(), json.at("z").get<float>());
}

std::vector<float> json_to_color(nlohmann::json json)
{
    return std::vector<float>({json.at("r").get<float>(), json.at("g").get<float>(), json.at("b").get<float>(), (json.contains("a") ? json.at("a").get<float>() : 1.f)});
}

nlohmann::json color_to_json(const std::vector<float> &color)
{
    return nlohmann::json({{"r", color[0]}, {"g", color[1]}, {"b", color[2]}, {"a", (color.size() > 3 ? color[3] : 1.f)}});
}

nlohmann::json color_to_json(const Vector3 &color)
{
    return color_to_json(std::vector<float>({color.x, color.y, color.z, 1.f}));
}
