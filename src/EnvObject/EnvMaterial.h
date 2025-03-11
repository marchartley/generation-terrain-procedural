#ifndef ENVMATERIAL_H
#define ENVMATERIAL_H

#include "DataStructure/Matrix3.h"

class EnvMaterial
{
public:
    EnvMaterial();
    EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, float mass, float decay, float virtualHeight, const Vector3& gridSize);

    std::string name;
    float diffusionSpeed;
    float waterTransport;
    float mass = 0.f;
    float decay;
    float virtualHeight = 1.f;

    bool isStable = false;

    GridF currentState;

    nlohmann::json toJSON() const;
    bool fromJSON(nlohmann::json json);

    void update(const GridV3& waterCurrents, const GridV3& heightsGradients, float dt);
};

#endif // ENVMATERIAL_H
