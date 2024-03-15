#include "EnvMaterial.h"

EnvMaterial::EnvMaterial()
{

}

EnvMaterial::EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, const Vector3 &gridSize)
    : name(name), diffusionSpeed(diffusionSpeed), waterTransport(waterTransport), currentState(GridF(gridSize))
{

}
