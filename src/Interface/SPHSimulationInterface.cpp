#include "SPHSimulationInterface.h"

SPHSimulationInterface::SPHSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("sphsimulation", "SPH Fluid simulation", "physics", "SPH Simulation", "sph_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[SPH];// = dynamic_cast<SPHSimulation*>(_simulation);
//    _simulation = new SPHSimulation();
}
