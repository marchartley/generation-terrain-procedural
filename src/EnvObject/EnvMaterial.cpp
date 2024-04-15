#include "EnvMaterial.h"

EnvMaterial::EnvMaterial()
{

}

EnvMaterial::EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, float mass, float decay, const Vector3 &gridSize)
    : name(name), diffusionSpeed(diffusionSpeed), waterTransport(waterTransport), mass(mass), decay(decay), currentState(GridF(gridSize))
{

}

nlohmann::json EnvMaterial::toJSON() const
{
    nlohmann::json json;
    json["name"] = this->name;
    json["data"] = stringifyGridF(this->currentState, false);
    return json;
}

bool EnvMaterial::fromJSON(nlohmann::json json)
{
    if (json["name"] == this->name) {
        this->currentState = loadGridF(json["data"], false);
        return true;
    }
    return false;
}
