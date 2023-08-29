#include "LBMFluidSimulationInterface.h"

LBMFluidSimulationInterface::LBMFluidSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("lbmluidsimulation", "LBM Fluid Simulation", "physics", "LBM fluid Simulation", "lbm_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[LBM]; // = dynamic_cast<LBMFluidSimulation*>(_simulation);
//    _simulation = new LBMFluidSimulation();
}
