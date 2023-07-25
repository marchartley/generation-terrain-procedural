#include "LBMFluidSimulationInterface.h"

LBMFluidSimulationInterface::LBMFluidSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("LBM Fluids", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[LBM]; // = dynamic_cast<LBMFluidSimulation*>(_simulation);
//    _simulation = new LBMFluidSimulation();
}
