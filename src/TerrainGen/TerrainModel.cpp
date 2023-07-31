#include "TerrainModel.h"
#include "TerrainGen/ImplicitPatch.h"
#include "Utils/ShapeCurve.h"

TerrainModel::TerrainModel()
{
    this->properties = GlobalTerrainProperties::get();
}

TerrainModel::~TerrainModel()
{

}

float TerrainModel::getHeight(Vector3 pos)
{
    return this->getHeight(pos.x, pos.y);
}

bool TerrainModel::contains(Vector3 v)
{
    return Vector3::isInBox(v, Vector3(), this->getDimensions());
}

bool TerrainModel::contains(float x, float y, float z)
{
    return this->contains(Vector3(x, y, z));
}

size_t TerrainModel::getCurrentHistoryIndex() const
{
    return this->_cachedHistoryIndex; // Not good
}

void TerrainModel::initFluidSim()
{
    float dt = 0.1f;
    float diffusion = 0.8f;
    float viscosity = 0.01f;
    int fluidSolverIterations = 5;
    float density = 1000.f;
    this->fluidSimRescale = Vector3(4, 4, 4);

    /*
    Vector3 simulationSize = this->getDimensions() / fluidSimRescale;
    *properties->simulations[FluidSimType::FLIP] = FLIPSimulation(density, simulationSize.x, simulationSize.y, simulationSize.z, 1, .2f, 20000, dt);
    *properties->simulations[FluidSimType::SPH] = SPHSimulation();
    *properties->simulations[FluidSimType::LBM] = LBMFluidSimulation(true);
//    *properties->simulations[FluidSimType::SHALLOW] = ShallowWaterSimulation();
    *properties->simulations[FluidSimType::STABLE] = StableFluidsSimulation(simulationSize.x, simulationSize.y, simulationSize.z, dt, diffusion, viscosity, fluidSolverIterations);
    *properties->simulations[FluidSimType::WARP] = WarpedFluidSimulation(simulationSize.x, simulationSize.y, simulationSize.z);
    */
    /*
    this->fluidSimulation = StableFluidsSimulation(this->getSizeX() / this->fluidSimRescale.x, this->getSizeY() / this->fluidSimRescale.y, this->getSizeZ() / this->fluidSimRescale.z, dt, diffusion, viscosity, fluidSolverIterations);

    this->multipleFluidSimulations.resize(4);
    this->multipleFlowFields.resize(4);
    this->multipleSeaCurrents.resize(4);
    for (size_t i = 0; i < this->multipleFluidSimulations.size(); i++) {
        this->multipleFluidSimulations[i] = StableFluidsSimulation(this->getSizeX() / this->fluidSimRescale.x, this->getSizeY() / this->fluidSimRescale.y, this->getSizeZ() / this->fluidSimRescale.z, dt, diffusion, viscosity, fluidSolverIterations);
        this->multipleFlowFields[i] = Matrix3<Vector3>(this->getDimensions());
    }
    float waterStrength = 1.f;
    this->multipleSeaCurrents = {
        Vector3(1.f, 0.f, 0.f).normalized() * waterStrength,
        Vector3(0.f, 1.f, 0.f).normalized() * waterStrength,
        Vector3(-1.f, 0.f, 0.f).normalized() * waterStrength,
        Vector3(0.f, -1.f, 0.f).normalized() * waterStrength
    };
    */
}

void TerrainModel::initEnvironmentalDensities()
{
    this->properties->environmentalDensities = Matrix3<float>(this->getDimensions(), 1); // Fill with air density for now
}


Matrix3<Vector3> TerrainModel::getFlowfield(FluidSimType simu)
{
//    return this->getFlowfield(0);
    return properties->simulations[simu]->getVelocities(this->getDimensions());
}
/*
Matrix3<Vector3> TerrainModel::getFlowfield(size_t flowIndex)
{
    Vector3 dimensions = this->getDimensions();
    return properties->simulations[LBM]->getVelocities(dimensions.x, dimensions.y, dimensions.z);
}
*/
void TerrainModel::computeFlowfield(FluidSimType simu)
{
//    return properties->simulations[simu]->step();
    this->computeMultipleFlowfields(simu);
}

void TerrainModel::computeMultipleFlowfields(FluidSimType simu, int steps, TerrainModel *implicit)
{
    FluidSimulation* simulation = this->properties->simulations[simu];
    auto primitives = dynamic_cast<ImplicitNaryOperator*>(implicit);
    auto smallerVoxelGrid = this->getVoxelized().resize(simulation->dimensions);
//    auto geom = Mesh().applyMarchingCubes(smallerVoxelGrid).getTriangles(); //this->getGeometry().getTriangles();
    Matrix3<float> obstacleMap = smallerVoxelGrid.binarize();
    auto densities = this->getEnvironmentalDensities().resize(simulation->dimensions);
    float maxDensity = densities.max();
    for (size_t i = 0; i < obstacleMap.size(); i++)
        obstacleMap[i] = ((maxDensity > 10 && densities[i] < 1000) || obstacleMap.getCoordAsVector3(i).z < 1 ? 1 : obstacleMap[i]);
    size_t nbCurrentsToCompute = 1; //this->multipleFluidSimulations.size();
    fluidSimRescale = this->getDimensions() / simulation->dimensions;


    /***
     * Implicit affect flowfield part :
     ***/
    std::vector<BSpline> allTunnelsCurves;
    std::vector<Matrix3<float>> rasterizedCurves;
    auto tunnelsPatches = primitives->findAll(ImplicitPatch::ParametricTunnel);
    for (auto& tunnelPatch : tunnelsPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(tunnelPatch);
        if (asPrimitive && asPrimitive->material == WATER) {
            BSpline curve = asPrimitive->optionalCurve;
            for (auto& p : curve.points) {
                p = asPrimitive->getGlobalPositionOf(p);
                p /= this->fluidSimRescale; // Rescale the curves to fit the simulation process
            }
            allTunnelsCurves.push_back(curve);
            Matrix3<float> rasterizedCurve(smallerVoxelGrid.getDimensions(), -1.f);
            for (int x = 0; x < rasterizedCurve.sizeX; x++) {
                for (int y = 0; y < rasterizedCurve.sizeY; y++) {
                    for (int z = 0; z < rasterizedCurve.sizeZ; z++) {
                        Vector3 pos(x, y, z);
                        float closestTime = curve.estimateClosestTime(pos);
                        Vector3 closest = curve.getPoint(closestTime);
                        if ((pos - closest).norm2() < (asPrimitive->parametersProvided[0] / fluidSimRescale.x) * (asPrimitive->parametersProvided[0] / fluidSimRescale.x))
                            rasterizedCurve.at(pos) = closestTime;
                    }
                }
            }
            rasterizedCurves.push_back(rasterizedCurve);
        }
    }

    std::vector<Matrix3<float>> rasterizedLagoonAreas;
    auto lagoonPatches = primitives->findAll(ImplicitPatch::Polygon);
    for (auto& lagoonPatch : lagoonPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(lagoonPatch);
        if (asPrimitive /* && asPrimitive->material == WATER*/) {
            ShapeCurve curve = asPrimitive->optionalCurve;
            for (auto& p : curve.points) {
                p = asPrimitive->getGlobalPositionOf(p);
                p /= this->fluidSimRescale; // Rescale the curves to fit the simulation process
            }
            Matrix3<float> rasterizedArea(smallerVoxelGrid.getDimensions(), -1.f);
            for (int x = 0; x < rasterizedArea.sizeX; x++) {
                for (int y = 0; y < rasterizedArea.sizeY; y++) {
                    for (int z = 0; z < rasterizedArea.sizeZ; z++) {
                        Vector3 pos(x, y, z);
                        rasterizedArea.at(pos) = (curve.contains(pos, false) ? 1 : 0);
                    }
                }
            }
            rasterizedArea = 1.f - rasterizedArea.toDistanceMap().normalize();
            rasterizedLagoonAreas.push_back(rasterizedArea);
        }
    }


    /***
     * End Implicit affect flowfield
     ***/

#pragma omp parallel for
    for (size_t iCurrent = 0; iCurrent < nbCurrentsToCompute; iCurrent++) {
//        if (iCurrent > 0) continue;
        Vector3 simulationDimensions = simulation->dimensions;
        simulation->setObstacles(obstacleMap);
//        this->multipleFluidSimulations[iCurrent].setObstacles(geom);
        for (int x = 0; x < simulationDimensions.x; x++) {
            for (int y = 0; y < simulationDimensions.y; y++) {
                for (int z = 0; z < simulationDimensions.z; z++) {
                    Vector3 pos(x, y, z);
                    if (!Vector3::isInBox(pos - (sea_current.normalized() * 2.f), Vector3(), simulationDimensions))
                        simulation->setVelocity(x, y, z, sea_current);
//                        this->multipleFluidSimulations[iCurrent].velocity(x, y, z) = this->multipleSeaCurrents[iCurrent];// / this->fluidSimRescale;
                }
            }
        }
        for (int i = 0; i < steps; i++) {

            /***
             * Implicit affect flowfield part :
             ***/
            for (size_t iTunnel = 0; iTunnel < allTunnelsCurves.size(); iTunnel++) {
                auto& curve = allTunnelsCurves[iTunnel];
                Vector3 inputFlow = simulation->getVelocity(curve.points.front());
                float inputStrength = inputFlow.norm();

                Matrix3<float>& rasterizedCurve = rasterizedCurves[iTunnel];
                for (int x = 0; x < rasterizedCurve.sizeX; x++) {
                    for (int y = 0; y < rasterizedCurve.sizeY; y++) {
                        for (int z = 0; z < rasterizedCurve.sizeZ; z++) {
                            Vector3 direction = curve.getDirection(rasterizedCurve.at(x, y, z));
                            simulation->addVelocity(x, y, z, direction * inputStrength);
                        }
                    }
                }
            }

            for (size_t iLagoon = 0; iLagoon < rasterizedLagoonAreas.size(); iLagoon++) {

                Matrix3<float>& rasterizedArea = rasterizedLagoonAreas[iLagoon];
                for (int x = 0; x < rasterizedArea.sizeX; x++) {
                    for (int y = 0; y < rasterizedArea.sizeY; y++) {
                        for (int z = 0; z < rasterizedArea.sizeZ; z++) {
                            simulation->setVelocity(x, y, z, simulation->getVelocity(x, y, z) * rasterizedArea.at(x, y, z));
//                            this->multipleFluidSimulations[iCurrent].velocity.at(x, y, z) *=
                        }
                    }
                }
            }

            /***
             * End Implicit affect flowfield
             ***/
//            this->multipleFluidSimulations[iCurrent].step();
            simulation->step();
        }
//        std::cout << "Max " << iCurrent << " " << this->multipleFluidSimulations[iCurrent].velocity.max() << " sum " << this->multipleFluidSimulations[iCurrent].velocity.sum() << "\n" << std::flush;
//        this->multipleFlowFields[iCurrent] = this->multipleFluidSimulations[iCurrent].getVelocities(this->getSizeX(), this->getSizeY(), this->getSizeZ());
//        std::cout << "Max " << iCurrent << " " << this->multipleFlowFields[iCurrent].max() << " " << this->multipleFlowFields[iCurrent].max().norm() << " sum " << this->multipleFlowFields[iCurrent].sum() << "\n" << std::flush;
    }

}


Matrix3<float> &TerrainModel::getEnvironmentalDensities()
{
    if (this->properties->environmentalDensities.size() < 2) {
        this->properties->environmentalDensities = Matrix3<float>(this->getDimensions());
        this->updateEnvironmentalDensities(0.f);
    }
    return this->properties->environmentalDensities;
}

void TerrainModel::updateEnvironmentalDensities(float waterLevel)
{
    this->storedWaterLevel = waterLevel;
    float waterLevelOnVoxels = waterLevel * this->getSizeZ();
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            for (int z = 0; z < this->getSizeZ(); z++) {
//                this->environmentalDensities.at(x, y, z) = (z < waterLevel ? (1000 + (1 - float(z)/waterLevel) * 1000)
//                                                                           : (1000 - (z - waterLevel)/float(this->getSizeZ()) * 999));
                this->properties->environmentalDensities.at(x, y, z) = (z < waterLevelOnVoxels ? 1000 : 1);
            }
        }
    }
}
