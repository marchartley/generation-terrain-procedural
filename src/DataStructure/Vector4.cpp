#include "Vector4.h"

Vector4::Vector4() : Vector4(0, 0, 0, 0)
{

}

Vector4::Vector4(float x, float y, float z, float w)
    : Vector3(x, y, z), w(w)
{

}

Vector4::Vector4(const Vector3 &xyz, float w)
    : Vector4(xyz.x, xyz.y, xyz.z, w)
{

}
