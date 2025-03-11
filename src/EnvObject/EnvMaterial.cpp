#include "EnvMaterial.h"

EnvMaterial::EnvMaterial()
{

}

EnvMaterial::EnvMaterial(std::string name, float diffusionSpeed, float waterTransport, float mass, float decay, float virtualHeight, const Vector3 &gridSize)
    : name(name), diffusionSpeed(diffusionSpeed), waterTransport(waterTransport), mass(mass), decay(decay), virtualHeight(virtualHeight), currentState(GridF(gridSize))
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

void EnvMaterial::update(const GridV3 &waterCurrents, const GridV3 &heightsGradients, float dt)
{
    if (this->diffusionSpeed < 1.f) {
        if (random_gen::generate() < this->diffusionSpeed) {
            this->currentState = this->currentState.meanSmooth(3, 3, 1, false); // Diffuse a little bit
        }
    } else {
        this->currentState = this->currentState.meanSmooth(this->diffusionSpeed * 2 + 1, this->diffusionSpeed * 2 + 1, 1, false); // Diffuse
    }
    this->currentState = (this->currentState.warpWith((waterCurrents * this->waterTransport) - (heightsGradients * this->mass))) * std::pow(this->decay, dt);
}
