#include "EnvMaterial.h"

EnvMaterial::EnvMaterial()
{

}

EnvMaterial::EnvMaterial(const std::string& name, float diffusionSpeed, float waterTransport, float mass, float decay, float virtualHeight, const Vector3 &gridSize)
    : name(name), diffusionSpeed(diffusionSpeed), waterTransport(waterTransport), mass(mass), decay(decay), virtualHeight(virtualHeight), currentState(GridF(gridSize))
{

}

void EnvMaterial::update(const GridV3 &waterCurrents, const GridV3 &heightsGradients, float dt)
{
    // if (!this->isStable) {
        if (this->diffusionSpeed < 1.f) {
            if (random_gen::generate() < this->diffusionSpeed) {
                this->currentState = this->currentState.meanSmooth(3, 3, 1, false); // Diffuse a little bit
            }
        } else {
            this->currentState = this->currentState.meanSmooth(this->diffusionSpeed * 2 + 1, this->diffusionSpeed * 2 + 1, 1, false); // Diffuse
        }
        this->currentState = this->currentState.warpWith((waterCurrents * this->waterTransport) - (heightsGradients * this->mass)).max(0.f) * std::pow(this->decay, dt);
    // }
}
