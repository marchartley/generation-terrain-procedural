#ifndef FLOWFIELDINTERFACE_H
#define FLOWFIELDINTERFACE_H


class StableFluidSimulationInterface;
// #include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "FluidSimulation/StableFluidsFluidSimulation.h"

class StableFluidSimulationInterface : public AbstractFluidSimulationInterface
{
    Q_OBJECT
public:
    StableFluidSimulationInterface(QWidget *parent = nullptr);

//    void updateParticlesMesh();

public Q_SLOTS:
    virtual void computeSimulation(int nbSteps = 1);

protected:
};

#endif // FLOWFIELDINTERFACE_H
