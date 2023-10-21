#ifndef VECTOR4_H
#define VECTOR4_H

#include "Vector3.h"
#include <glm/vec4.hpp>

class Vector4 : public Vector3
{
public:
    Vector4();
    Vector4(float x, float y, float z, float w);
    Vector4(const Vector3& xyz, float w);

    operator glm::vec4() const { return glm::vec4(this->x, this->y, this->z, this->w); }
    explicit operator float*() const { return new float[4]{this->x, this->y, this->z, this->w}; }

    float w;
};

#endif // VECTOR4_H
