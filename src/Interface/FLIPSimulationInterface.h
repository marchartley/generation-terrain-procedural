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
    QLayout* createGUI();

    void resetParticles();

    virtual void display(const Vector3 &camPos = Vector3(false));

public Q_SLOTS:
    void updateParticles();
    void computeSimulation(int nbSteps);

    void keyPressEvent(QKeyEvent* e);

protected:
    FLIPSimulation* simulation;
    std::vector<float> changesHistory;

    GridV3 _currentVelocities;
};

#endif // FLIPSIMULATIONINTERFACE_H
