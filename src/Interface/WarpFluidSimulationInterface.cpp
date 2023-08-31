#include "WarpFluidSimulationInterface.h"

WarpFluidSimulationInterface::WarpFluidSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("warpfluidsimulation", "Warp Fluids simulation", "physics", "Warp fluid Simulation", "warp_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[WARP];
}

void WarpFluidSimulationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    AbstractFluidSimulationInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    this->computeFromTerrain(voxelGrid.get());
}

void WarpFluidSimulationInterface::display(const Vector3& camPos)
{
    AbstractFluidSimulationInterface::display(camPos);
}


void WarpFluidSimulationInterface::computeFromTerrain(TerrainModel *terrain)
{
    auto simulation = dynamic_cast<WarpedFluidSimulation*>(_simulation);
    simulation->recomputeVelocities();
}

void WarpFluidSimulationInterface::afterTerrainUpdated()
{
    this->computeFromTerrain(voxelGrid.get());
}
