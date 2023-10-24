#ifndef SPHSIMULATIONINTERFACE_H
#define SPHSIMULATIONINTERFACE_H

class SPHSimulationInterface;
#include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "FluidSimulation/SPHSimulation.h"

class SPHSimulationInterface : public AbstractFluidSimulationInterface
{
    Q_OBJECT
public:
    SPHSimulationInterface(QWidget *parent = nullptr);
/*
    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

    void updateSimulationMesh();

public Q_SLOTS:
    void show();
    void hide();

    virtual void afterTerrainUpdated();

    void updateParticles();

protected:
    Mesh particlesMesh;
    Mesh ghostsMesh;
    Mesh connectionsMesh;

    SPH::SPHSimulation simulation;

    std::vector<BSpline> particlePaths;
    */
};

#endif // SPHSIMULATIONINTERFACE_H
