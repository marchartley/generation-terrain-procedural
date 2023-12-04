#ifndef WARPFLUIDSIMULATIONINTERFACE_H
#define WARPFLUIDSIMULATIONINTERFACE_H


class WarpFluidSimulationInterface;

#include <QWidget>
#include "Interface/AbstractFluidSimulationInterface.h"
#include "FluidSimulation/WarpedFluidSimulation.h"

class WarpFluidSimulationInterface : public AbstractFluidSimulationInterface
{
    Q_OBJECT
public:
    WarpFluidSimulationInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void computeFromTerrain(TerrainModel* terrain);

    void display(const Vector3 &camPos);
public Q_SLOTS:

    virtual void afterTerrainUpdated();

protected:
};

#endif // WARPFLUIDSIMULATIONINTERFACE_H
