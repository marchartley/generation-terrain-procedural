#ifndef ENVMATERIAL_H
#define ENVMATERIAL_H

#include "DataStructure/Matrix3.h"

class EnvMaterial
{
public:
    EnvMaterial();
    EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, float mass, float decay, const Vector3& gridSize);

    std::string name;
    float diffusionSpeed;
    float waterTransport;
    float mass = 0.f;
    float decay;

    GridF currentState;
};

#endif // ENVMATERIAL_H
