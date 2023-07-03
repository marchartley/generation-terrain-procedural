#ifndef FLIPSIMULATIONINTERFACE_H
#define FLIPSIMULATIONINTERFACE_H

class FLIPSimulationInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "FluidSimulation/FLIPSimulation.h"

class FLIPSimulationInterface : public ActionInterface
{
    Q_OBJECT
public:
    FLIPSimulationInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(Vector3 camPos = Vector3(false));

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

    FLIP::FLIPSimulation simulation;

    std::vector<BSpline> particlePaths;


    Mesh debugObstacleMesh;
};

#endif // FLIPSIMULATIONINTERFACE_H
