#include "GlobalTerrainProperties.h"

GlobalTerrainProperties* GlobalTerrainProperties::instance = nullptr;

GlobalTerrainProperties::GlobalTerrainProperties()
{

    float dt = 0.05f;
    float diffusion = 0.8f;
    float viscosity = 0.01f;
    int fluidSolverIterations = 5;
    float density = 1000.f;
    Vector3 fluidSimRescale = Vector3(2, 2, 2);

    Vector3 simulationSize = Vector3(100, 100, 40) / fluidSimRescale;
    this->simulations[FluidSimType::FLIP] = new FLIPSimulation(density, simulationSize.x, simulationSize.y, simulationSize.z, .5f, 1.0f, 750, dt);
    this->simulations[FluidSimType::SPH] = new SPHSimulation();
    this->simulations[FluidSimType::LBM] = new LBMFluidSimulation(false);
//    this->simulations[FluidSimType::SHALLOW] = new ShallowWaterSimulation();
    this->simulations[FluidSimType::STABLE] = new StableFluidsSimulation(simulationSize.x, simulationSize.y, simulationSize.z, dt, diffusion, viscosity, fluidSolverIterations);
    this->simulations[FluidSimType::WARP] = new WarpedFluidSimulation(100, 100, 100); //simulationSize.x, simulationSize.y, simulationSize.z);
    /*
    Vector3 simulationSize(1, 1, 1);
    this->simulations[FluidSimType::FLIP] = new FLIPSimulation(1000, simulationSize.x, simulationSize.y, simulationSize.z, 1, .2f, 20000, 1.f); // dt);
    this->simulations[FluidSimType::SPH] = new SPHSimulation();
    this->simulations[FluidSimType::LBM] = new LBMFluidSimulation(false);
//    *properties->simulations[FluidSimType::SHALLOW] = ShallowWaterSimulation();
    this->simulations[FluidSimType::STABLE] = new StableFluidsSimulation(simulationSize.x, simulationSize.y, simulationSize.z, 1.f, 1.f, 1.f, 1); //dt, diffusion, viscosity, fluidSolverIterations);
    this->simulations[FluidSimType::WARP] = new WarpedFluidSimulation(simulationSize.x, simulationSize.y, simulationSize.z);
    */
}

GlobalTerrainProperties::~GlobalTerrainProperties()
{
    for (auto& simu : simulations) {
        if (simu.second)
            delete simu.second;
    }
//    if (flipSimulation)
//        delete flipSimulation;
//    if (lbmSimulation)
//        delete lbmSimulation;
//    if (shallowSimulation)
//        delete shallowSimulation;
//    if (sphSimulation)
//        delete sphSimulation;
//    if (stableSimulation)
//        delete stableSimulation;
//    if (warpSimulation)
//        delete warpSimulation;
}

std::string stringFromFluidSimType(FluidSimType type)
{
   if (type == FLIP) return "FLIP";
   if (type == LBM) return "LBM";
   if (type == SHALLOW) return "Shallow";
   if (type == SPH) return "SPH";
   if (type == STABLE) return "Stable";
   if (type == WARP) return "Warp";
   return "";
}

FluidSimType FluidSimTypeFromString(std::string type)
{
    type = toUpper(type);
    if (type == "FLIP") return FluidSimType::FLIP;
    if (type == "LBM") return FluidSimType::LBM;
    if (type == "SHALLOW") return FluidSimType::SHALLOW;
    if (type == "SPH") return FluidSimType::SPH;
    if (type == "STABLE") return FluidSimType::STABLE;
    if (type == "WARP") return FluidSimType::WARP;
    return FluidSimType::LBM;
}
