#include "EnvMaterial.h"

EnvMaterial::EnvMaterial()
{

}

EnvMaterial::EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, float mass, const Vector3 &gridSize)
    : name(name), diffusionSpeed(diffusionSpeed), waterTransport(waterTransport), mass(mass), currentState(GridF(gridSize))
{

}
