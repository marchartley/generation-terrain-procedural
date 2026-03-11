#include "StableFluidSimulationInterface.h"
#include "FluidSimulation/OpenFoamParser.h"

StableFluidSimulationInterface::StableFluidSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("stablefluidssimulation", "Stable Fluids Simulation", "physics", "Stable fluids simulation", "stable_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[STABLE]; // = dynamic_cast<StableFluidsSimulation*>(_simulation);
//    *_simulation = StableFluidsSimulation();
    this->nbComputationsPerFrame = 5;
}

void StableFluidSimulationInterface::computeSimulation(int nbSteps)
{
    displayProcessTime("OpenFoam computation... ", [&]() {
        OpenFoamParser::createSimulationFile("OpenFoam/OF_Sim_Marcos/", this->voxelGrid->getVoxelValues());
    });
    this->_simulation->_cachedStep = -1;
    this->updateVectorsMesh();
//    AbstractFluidSimulationInterface::computeSimulation(nbSteps);
}
