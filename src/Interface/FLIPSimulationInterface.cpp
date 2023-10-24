#include "FLIPSimulationInterface.h"

FLIPSimulationInterface::FLIPSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("flipsimulation", "FLIP simulation", "physics", "FLIP Simulation", "flip_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[FLIP]; // = dynamic_cast<FLIPSimulation*>(_simulation);
//    *_simulation = FLIPSimulation(1000, 20, 20, 1, 1, .2f, 20000, 0.1f);
}

void FLIPSimulationInterface::updateParticlesMesh()
{
    auto simulation = dynamic_cast<FLIPSimulation*>(_simulation);
    auto& particles = simulation->particles;
    std::vector<Vector3> particlePositions;
    std::vector<Vector3> colors;

    GridF pressure = GridF(simulation->dimensions);
    pressure.data = simulation->p;
    pressure.normalize();

    std::vector<Vector3> positions(simulation->particles.size());
    std::vector<Vector3> velocities(simulation->particles.size());
    Vector3 minVel = Vector3::max(), maxVel = Vector3::min();
    for (size_t i = 0; i < simulation->particles.size(); i++) {
        velocities[i] = simulation->particles[i].velocity;
        minVel = std::min(minVel, velocities[i]);
        maxVel = std::max(maxVel, velocities[i]);
        positions[i] = simulation->particles[i].position / simulation->dimensions;
    }
    for (size_t i = 0; i < simulation->particles.size(); i++)
        velocities[i] = (velocities[i] - minVel) / (maxVel - minVel);


    for (size_t i = 0; i < simulation->particles.size(); i++) {
            particlePositions.push_back(positions[i] * voxelGrid->getDimensions());
            colors.push_back(Vector3(particles[i].velocity.norm(), 0, 1));
    }
    particlesMesh.colorsArray = colors;
    particlesMesh.fromArray(particlePositions);
}


void FLIPSimulationInterface::updateParticles()
{
}
