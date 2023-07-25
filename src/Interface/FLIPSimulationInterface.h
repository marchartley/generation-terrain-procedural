#ifndef FLIPSIMULATIONINTERFACE_H
#define FLIPSIMULATIONINTERFACE_H

class FLIPSimulationInterface;
#include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "FluidSimulation/FLIPSimulation.h"

class FLIPSimulationInterface : public AbstractFluidSimulationInterface
{
    Q_OBJECT
public:
    FLIPSimulationInterface(QWidget *parent = nullptr);

    void updateParticlesMesh();

public Q_SLOTS:
    void updateParticles();

protected:
};

#endif // FLIPSIMULATIONINTERFACE_H
