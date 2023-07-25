#ifndef LBMFLUIDSIMULATIONINTERFACE_H
#define LBMFLUIDSIMULATIONINTERFACE_H

#include "AbstractFluidSimulationInterface.h"
#include "FluidSimulation/LBMFluidSimulation.h"

class LBMFluidSimulationInterface : public AbstractFluidSimulationInterface
{
public:
    LBMFluidSimulationInterface(QWidget* parent);
};

#endif // LBMFLUIDSIMULATIONINTERFACE_H
