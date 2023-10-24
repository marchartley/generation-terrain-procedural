#ifndef GLOBALTERRAINPROPERTIES_H
#define GLOBALTERRAINPROPERTIES_H

class GlobalTerrainProperties;

//#include "TerrainGen/Heightmap.h"
//#include "TerrainGen/LayerBasedGrid.h"
//#include "TerrainGen/VoxelGrid.h"
//#include "TerrainGen/ImplicitPatch.h"

#include "FluidSimulation/FLIPSimulation.h"
#include "FluidSimulation/LBMFluidSimulation.h"
#include "FluidSimulation/SPHSimulation.h"
#include "FluidSimulation/ShallowWaterSimulation.h"
#include "FluidSimulation/StableFluidsFluidSimulation.h"
#include "FluidSimulation/WarpedFluidSimulation.h"


enum FluidSimType {
    FLIP,
    LBM,
    SHALLOW,
    SPH,
    STABLE,
    WARP
};

class GlobalTerrainProperties
{
public:
    GlobalTerrainProperties();
    ~GlobalTerrainProperties();

    static GlobalTerrainProperties* get() {
        if (instance == nullptr) {
            instance = new GlobalTerrainProperties;
        }
        return instance;
    }

    static GlobalTerrainProperties* instance;

    // Terrains models
//    std::shared_ptr<Heightmap> heightmap;
//    std::shared_ptr<VoxelGrid> voxelGrid;
//    std::shared_ptr<LayerBasedGrid> layerGrid;
//    std::shared_ptr<ImplicitNaryOperator*> implicitTerrain;

    // Fluid sim models
    std::map<FluidSimType, FluidSimulation*> simulations;
//    FLIPSimulation* flipSimulation;
//    LBMFluidSimulation* lbmSimulation;
//    ShallowWaterSimulation* shallowSimulation;
//    SPHSimulation* sphSimulation;
//    StableFluidsSimulation* stableSimulation;
//    WarpedFluidSimulation* warpSimulation;


    GridF environmentalDensities;
};

std::string stringFromFluidSimType(FluidSimType type);

FluidSimType FluidSimTypeFromString(std::string type);

#endif // GLOBALTERRAINPROPERTIES_H
