#ifndef FLOWFIELDINTERFACE_H
#define FLOWFIELDINTERFACE_H


class FlowFieldInterface;
#include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "FluidSimulation/StableFluidsFluidSimulation.h"

class FlowFieldInterface : public AbstractFluidSimulationInterface
{
    Q_OBJECT
public:
    FlowFieldInterface(QWidget *parent = nullptr);

//    void updateParticlesMesh();

public Q_SLOTS:
//    void updateParticles();

protected:
};

/*
class FlowFieldInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"

class FlowFieldInterface : public ActionInterface
{
    Q_OBJECT
public:
    FlowFieldInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch);

    void display(Vector3 camPos = Vector3(false));
    void displayPressureDensities(Vector3 camPos = Vector3(false));
    void displayFlows(Vector3 camPos = Vector3(false));
    void displaySumOfFlows(Vector3 camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    void recomputeFlowfield(int steps = 30);
    void updateFlowfieldDebugMesh();


protected:
//    Mesh flowMesh;

    std::vector<Mesh> flowMeshes;
    Mesh sumOfFlowsMesh;

    GridF pressureDensityVoxels;
    GridV3 totalFlow;
    Mesh pressureDensityMesh;

    std::vector<bool> displayedFlowfields;
    bool displayingPressure = false;
    bool displayingSumOfFlows = false;
    bool computeAtEachFrame = true;
};
*/
#endif // FLOWFIELDINTERFACE_H
