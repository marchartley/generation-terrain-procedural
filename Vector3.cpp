#include "Vector3.h"
#include <math.h>
#include "Globals.h"

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

Vector3& Vector3::normalize() {
    if(this->norm() == 0)
        return *this;
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

GRID3D_vec3 Vector3::gradient(GRID3D_float field)
{
    GRID3D_vec3 returningGrid = GRID3D_vec3(field.size(), GRID2D_vec3(field[0].size(), GRID1D_vec3(field[0][0].size())));
    for (size_t x = 0; x < field.size(); x++) {
        for (size_t y = 0; y < field[x].size(); y++) {
            for (size_t z = 0; z < field[x][y].size(); z++) {
                Vector3 vec = Vector3();
                if (x == 0 || x == field.size() - 1 || y == 0 || y == field[x].size() - 1
                        || z == 0 || z == field[x][y].size() - 1) {
                    returningGrid[x][y][z] = vec;
                } else {
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if(dx != 0 && dy != 0 && dz != 0)
                                    vec += Vector3(dx, dy, dz) * field[x + dx][y + dy][z + dz];
                            }
                        }
                    }
                }
                returningGrid[x][y][z] = vec;

            }
        }
    }
    return returningGrid;
}

GRID3D_vec3 Vector3::grad(GRID3D_float field) {
    return Vector3::gradient(field);
}

GRID3D_float Vector3::divergence(GRID3D_vec3 field)
{
    GRID3D_float returningGrid = GRID3D_float(field.size(), GRID2D_float(field[0].size(), GRID1D_float(field[0][0].size())));
    for (size_t x = 0; x < field.size(); x++) {
        for (size_t y = 0; y < field[x].size(); y++) {
            for (size_t z = 0; z < field[x][y].size(); z++) {
                returningGrid[x][y][z] = field[x][y][z].divergence();
            }
        }
    }
    return returningGrid;
}

GRID3D_float Vector3::div(GRID3D_vec3 field)
{
    return Vector3::divergence(field);
}

GRID3D_vec3 Vector3::curl(GRID3D_vec3 field)
{
    GRID3D_vec3 returningGrid = GRID3D_vec3(field.size(), GRID2D_vec3(field[0].size(), GRID1D_vec3(field[0][0].size())));
    for (size_t x = 0; x < field.size(); x++) {
        for (size_t y = 0; y < field[x].size(); y++) {
            for (size_t z = 0; z < field[x][y].size(); z++) {
                Vector3& vec = field[x][y][z];
                returningGrid[x][y][z] = Vector3(vec.z - vec.y, vec.x - vec.z, vec.y - vec.x);
            }
        }
    }
    return returningGrid;
}

GRID3D_vec3 Vector3::rot(GRID3D_vec3 field)
{
    return Vector3::curl(field);
}

GRID3D_vec3 Vector3::laplacian(GRID3D_vec3 field)
{
    GRID3D_vec3 returningGrid = GRID3D_vec3(field.size(), GRID2D_vec3(field[0].size(), GRID1D_vec3(field[0][0].size())));
    for (size_t x = 0; x < field.size(); x++) {
        for (size_t y = 0; y < field[x].size(); y++) {
            for (size_t z = 0; z < field[x][y].size(); z++) {
                Vector3 vec = Vector3();
                if (x == 0 || x == field.size() - 1 || y == 0 || y == field[x].size() - 1
                        || z == 0 || z == field[x][y].size() - 1) {
                    returningGrid[x][y][z] = vec;
                } else {
                    vec += field[x    ][y    ][z + 1];
                    vec += field[x    ][y    ][z - 1];
                    vec += field[x    ][y + 1][z    ];
                    vec += field[x    ][y - 1][z    ];
                    vec += field[x + 1][y    ][z    ];
                    vec += field[x - 1][y    ][z    ];
                    vec -= field[x    ][y    ][z    ] * 6;
                }
                returningGrid[x][y][z] = vec;
            }
        }
    }
    return returningGrid;
}

GRID3D_float Vector3::laplacian(GRID3D_float field)
{
    GRID3D_float returningGrid = GRID3D_float(field.size(), GRID2D_float(field[0].size(), GRID1D_float(field[0][0].size())));
    for (size_t x = 0; x < field.size(); x++) {
        for (size_t y = 0; y < field[x].size(); y++) {
            for (size_t z = 0; z < field[x][y].size(); z++) {
                float laplace = 0;
                if (x == 0 || x == field.size() - 1 || y == 0 || y == field[x].size() - 1
                        || z == 0 || z == field[x][y].size() - 1) {
                    // Keep laplacian value to 0
                } else {
                    laplace += field[x    ][y    ][z + 1];
                    laplace += field[x    ][y    ][z - 1];
                    laplace += field[x    ][y + 1][z    ];
                    laplace += field[x    ][y - 1][z    ];
                    laplace += field[x + 1][y    ][z    ];
                    laplace += field[x - 1][y    ][z    ];
                    laplace -= field[x    ][y    ][z    ] * 6;
                }
                returningGrid[x][y][z] = laplace;
            }
        }
    }
    return returningGrid;
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
Vector3 operator*(Vector3 a, Vector3 b) {
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
    io << v.toString();
    return io;
}

std::ostream& operator<<(std::ostream& io, Vector3* v) {
    io << v->toString();
    return io;
}

