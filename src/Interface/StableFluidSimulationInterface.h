#ifndef FLOWFIELDINTERFACE_H
#define FLOWFIELDINTERFACE_H


class StableFluidSimulationInterface;
#include <QWidget>
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

    void display(const Vector3& camPos = Vector3(false));
    void displayPressureDensities(const Vector3& camPos = Vector3(false));
    void displayFlows(const Vector3& camPos = Vector3(false));
    void displaySumOfFlows(const Vector3& camPos = Vector3(false));

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
