#ifndef ENVOBJECTFLUIDSIMULATION_H
#define ENVOBJECTFLUIDSIMULATION_H


class WarpFluidSimulationInterface;

#include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "EnvObject/EnvObject.h"

class EnvObjectFluidSimulation : public ActionInterface
{
    Q_OBJECT
public:
    EnvObjectFluidSimulation(QWidget* parent = nullptr);

    virtual void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    virtual void display(const Vector3& camPos = Vector3(false));

    virtual void replay(nlohmann::json action);

    virtual QLayout* createGUI();

    virtual void updateSimulationMeshes();

    virtual void updateVectorsMesh();
    virtual void updateParticlesMesh();
    virtual void updateBoundariesMesh();


public Q_SLOTS:
    virtual void computeSimulation(int nbSteps = 1);
    virtual void show();
    virtual void hide();

    virtual void afterTerrainUpdated();

protected:
    Mesh vectorsMesh;
    Mesh particlesMesh;
    Mesh boundariesMesh;
    Mesh gridBoundaryMesh;
    Mesh otherMeshToDisplay;

    FluidSimulation* _simulation;

    bool displayBoundaries = false;
    bool displayGridBoundaries = false;
    bool displayVectors = true;
    bool displayParticles = false;
    bool computeAtEachFrame = false;

    int nbComputationsPerFrame = 1;
};


#endif // ENVOBJECTFLUIDSIMULATION_H
