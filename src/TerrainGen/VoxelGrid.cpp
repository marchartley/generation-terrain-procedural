#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Utils/ShapeCurve.h"

VoxelGrid::VoxelGrid(int nx, int ny, int nz, float noise_shifting)
//    : /*blockSize(blockSize),*/ noise_shifting(noise_shifting)
{
    this->_cachedVoxelValues = Matrix3<float>(nx /* * chunkSize*/, ny /* * chunkSize*/, nz /* * chunkSize + 1*/, 1.f);
    this->initMap();
    this->fromCachedData();
}

VoxelGrid::VoxelGrid(Heightmap& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), 1.f) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(1, 1, 1, 1.f) {

}
VoxelGrid::~VoxelGrid()
{
//    this->chunks.clear();
}
void VoxelGrid::from2DGrid(Heightmap grid, Vector3 subsectionStart, Vector3 subsectionEnd, float scaleFactor) {
    float zScale = 2.f;
    this->_cachedVoxelValues = Matrix3<float>(grid.getSizeX(), grid.getSizeY(), grid.getSizeZ() * zScale);
    this->initMap();
    if (subsectionEnd == subsectionStart) {
        // If they are not set, we want the entire map
        subsectionStart = Vector3();
        subsectionEnd = Vector3(grid.getSizeX(), grid.getSizeY());
//        this->getSizeX() = grid.getSizeX();
//        this->getSizeY() = grid.getSizeY();
//        this->getSizeZ() = grid.getMaxHeight() /* * scaleFactor*/ * 2;
    } else {
        // Otherwise, we want a subset of the map, we just need to clamp the dimensions
        subsectionStart.x = std::max(subsectionStart.x, 0.f);
        subsectionStart.y = std::max(subsectionStart.y, 0.f);
        subsectionEnd.x = std::min(subsectionEnd.x, (float)grid.getSizeX());
        subsectionEnd.y = std::min(subsectionEnd.y, (float)grid.getSizeY());
    }
    this->_cachedVoxelValues = Matrix3<float>((subsectionEnd.x - subsectionStart.x) * scaleFactor, (subsectionEnd.y - subsectionStart.y) * scaleFactor, grid.getMaxHeight(), 0.f);
    this->initMap();

    Matrix3<float> gridHeights = grid.getHeights().subset(subsectionStart.xy(), subsectionEnd.xy()).resize(this->getDimensions());
    gridHeights.raiseErrorOnBadCoord = false;
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float grid_height = gridHeights.at(x, y) * (this->getSizeZ() / (grid.getMaxHeight()));
            int z = int(std::max(grid_height, 2.f));
            // Positive values
            for (int i = 0; i < int(z); i++) {
                _cachedVoxelValues.at(x, y, i) = .5f;
            }
            if (z < this->getSizeZ()) {
                _cachedVoxelValues.at(x, y, z) = interpolation::inv_linear(z - int(z), -.5f, .5f);
                for (int i = z+1; i < this->getSizeZ(); i++) {
                    _cachedVoxelValues.at(x, y, i) = -.5f;
                }
            }
        }
    }
    this->fromCachedData();
    this->smoothVoxels();
    this->smoothVoxels();
}

void VoxelGrid::fromLayerBased(LayerBasedGrid layerBased, int fixedHeight)
{
    this->setVoxelValues(layerBased.voxelize(fixedHeight == -1 ? layerBased.getSizeZ() * 2.f : fixedHeight));
}

void VoxelGrid::fromImplicit(ImplicitPatch *implicitTerrain, int fixedHeight)
{
    this->setVoxelValues(implicitTerrain->getVoxelized(implicitTerrain->getDimensions().xy() + Vector3(0, 0, fixedHeight == -1 ? 40.f : fixedHeight)));
    this->smoothVoxels();
}

VoxelGrid* VoxelGrid::fromCachedData()
{
    this->voxelsValuesStack.clear();
    this->voxelsValuesAnchorStack.clear();
    this->currentHistoryIndex = 0;
    this->applyModification(this->_cachedVoxelValues);

    /*
    for (auto& vc : this->chunks) {
        vc->currentHistoryIndex = 0;
        vc->voxelsValuesStack.clear();
        vc->voxelsValuesAnchorStack.clear();
        vc->applyModification(this->_cachedVoxelValues.subset(vc->x, vc->x + vc->sizeX, vc->y, vc->y + vc->sizeY, vc->z, vc->z + vc->sizeZ));
    }
    */
    if (this->_smoothingNeeded)  {
        this->smoothVoxels();
        this->_smoothingNeeded = false;
    }
    this->_cachedHistoryIndex = -1; // Force refresh of "getVoxelValues"
    return this;
}

void VoxelGrid::setVoxelValues(const Matrix3<float> &values)
{
    this->_cachedVoxelValues = values;
    this->initMap();
    this->fromCachedData();
}
/*
void VoxelGrid::computeFlowfield(FluidSimType type)
{
    this->computeMultipleFlowfields(type);
//    this->flowField = this->getFlowfield(type);
}
*/
/*
void VoxelGrid::computeMultipleFlowfields(FluidSimType type, int steps, ImplicitNaryOperator *primitives)
{
    auto smallerVoxelGrid = this->getVoxelValues().resize(this->multipleFluidSimulations[0].dimensions);
//    auto geom = Mesh().applyMarchingCubes(smallerVoxelGrid).getTriangles(); //this->getGeometry().getTriangles();
    Matrix3<float> obstacleMap = smallerVoxelGrid.binarize();
    auto densities = this->getEnvironmentalDensities().resize(this->multipleFluidSimulations[0].dimensions);
    float maxDensity = densities.max();
    for (size_t i = 0; i < obstacleMap.size(); i++)
        obstacleMap[i] = ((maxDensity > 10 && densities[i] < 1000) || obstacleMap.getCoordAsVector3(i).z < 1 ? 1 : obstacleMap[i]);
    size_t nbCurrentsToCompute = this->multipleFluidSimulations.size();


    ///
     /// Implicit affect flowfield part :
     ///
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
        if (asPrimitive) {
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


    ///
     /// End Implicit affect flowfield
     ///

#pragma omp parallel for
    for (size_t iCurrent = 0; iCurrent < nbCurrentsToCompute; iCurrent++) {
//        if (iCurrent > 0) continue;
        Vector3 simulationDimensions = this->multipleFluidSimulations[iCurrent].dimensions;
        this->multipleFluidSimulations[iCurrent].setObstacles(obstacleMap);
//        this->multipleFluidSimulations[iCurrent].setObstacles(geom);
        for (int x = 0; x < simulationDimensions.x; x++) {
            for (int y = 0; y < simulationDimensions.y; y++) {
                for (int z = 0; z < simulationDimensions.z; z++) {
                    Vector3 pos(x, y, z);
                    if (!Vector3::isInBox(pos - (this->multipleSeaCurrents[iCurrent].normalized() * 2.f), Vector3(), simulationDimensions))
                        this->multipleFluidSimulations[iCurrent].velocity(x, y, z) = this->multipleSeaCurrents[iCurrent];// / this->fluidSimRescale;
                }
            }
        }
        for (int i = 0; i < steps; i++) {

            ///
             /// Implicit affect flowfield part :
             ///
            for (size_t iTunnel = 0; iTunnel < allTunnelsCurves.size(); iTunnel++) {
                auto& curve = allTunnelsCurves[iTunnel];
                Vector3 inputFlow = this->multipleFluidSimulations[iCurrent].velocity.at(curve.points.front());
                float inputStrength = inputFlow.norm();

                Matrix3<float>& rasterizedCurve = rasterizedCurves[iTunnel];
                for (int x = 0; x < rasterizedCurve.sizeX; x++) {
                    for (int y = 0; y < rasterizedCurve.sizeY; y++) {
                        for (int z = 0; z < rasterizedCurve.sizeZ; z++) {
                            Vector3 direction = curve.getDirection(rasterizedCurve.at(x, y, z));
                            this->multipleFluidSimulations[iCurrent].velocity.at(x, y, z) += direction * inputStrength;
                        }
                    }
                }
            }

            for (size_t iLagoon = 0; iLagoon < rasterizedLagoonAreas.size(); iLagoon++) {

                Matrix3<float>& rasterizedArea = rasterizedLagoonAreas[iLagoon];
                for (int x = 0; x < rasterizedArea.sizeX; x++) {
                    for (int y = 0; y < rasterizedArea.sizeY; y++) {
                        for (int z = 0; z < rasterizedArea.sizeZ; z++) {
                            this->multipleFluidSimulations[iCurrent].velocity.at(x, y, z) *= rasterizedArea.at(x, y, z);
                        }
                    }
                }
            }

            ///
             /// End Implicit affect flowfield
             ///
            this->multipleFluidSimulations[iCurrent].step();
        }
        std::cout << "Max " << iCurrent << " " << this->multipleFluidSimulations[iCurrent].velocity.max() << " sum " << this->multipleFluidSimulations[iCurrent].velocity.sum() << "\n" << std::flush;
        this->multipleFlowFields[iCurrent] = this->multipleFluidSimulations[iCurrent].getVelocities(this->getSizeX(), this->getSizeY(), this->getSizeZ());
        std::cout << "Max " << iCurrent << " " << this->multipleFlowFields[iCurrent].max() << " " << this->multipleFlowFields[iCurrent].max().norm() << " sum " << this->multipleFlowFields[iCurrent].sum() << "\n" << std::flush;
    }
}*/
/*
void VoxelGrid::affectFlowfieldAround(Vector3 pos, Vector3 newVal, int kernelSize)
{
    this->affectFlowfieldAround(pos.x, pos.y, pos.z, newVal, kernelSize);
}
void VoxelGrid::affectFlowfieldAround(float x, float y, float z, Vector3 newVal, int kernelSize)
{
    Vector3 pos(x, y, z);
    for (int dx = -(kernelSize/2); dx <= (kernelSize/2); dx++) {
        for (int dy = -(kernelSize/2); dy <= (kernelSize/2); dy++) {
            for (int dz = -(kernelSize/2); dz <= (kernelSize/2); dz++) {
                Vector3 offset(dx, dy, dz);
                Vector3 currFlow = this->getFlowfield(pos + offset);

//                currFlow += ((offset * -1.0) + newVal) * (1.0 - (offset.norm() / (kernelSize/2.0)));
                currFlow += (newVal * (1.0 - (offset.norm() / (kernelSize/2.0)))); // Change proportionnal to distance
                this->setFlowfield(pos + offset, currFlow); // / 2.0); // Store the mean of the previous and new value
            }
        }
    }
}
void VoxelGrid::affectFlowfieldAround(Vector3 pos, float alphaEffect, int kernelSize)
{
    this->affectFlowfieldAround(pos.x, pos.y, pos.z, alphaEffect, kernelSize);
}
void VoxelGrid::affectFlowfieldAround(float x, float y, float z, float alphaEffect, int kernelSize)
{
    Vector3 pos(x, y, z);
    for (int dx = -(kernelSize/2); dx <= (kernelSize/2); dx++) {
        for (int dy = -(kernelSize/2); dy <= (kernelSize/2); dy++) {
            for (int dz = -(kernelSize/2); dz <= (kernelSize/2); dz++) {
                Vector3 offset(dx, dy, dz);
                Vector3 currFlow = this->getFlowfield(pos + offset);

                currFlow += ((offset * -1.0) * alphaEffect);
//                currFlow += (newVal * (1.0 - (offset.norm() / (kernelSize/2.0)))); // Change proportionnal to distance
                this->setFlowfield(pos + offset, currFlow / 2.0); // Store the mean of the previous and new value
            }
        }
    }
}
*/
void VoxelGrid::initMap()
{
//    this->chunkSize = std::min(int(this->getSizeX()), this->chunkSize);

//    this->chunks.clear();
//    this->flowField.clear();
    this->distanceField.clear();
    this->pressureField.clear();
    this->_cachedHistoryIndex = -1;

    initFluidSim();
    initEnvironmentalDensities();
/*
    // Create and configure FastNoise object
    this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    this->noise.SetFrequency(1.f / (float) this->getSizeX());
    this->noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    this->noise.SetFractalLacunarity(2.0);
    this->noise.SetFractalGain(0.7);
    this->noise.SetFractalWeightedStrength(0.5);
    this->noise.SetFractalOctaves(10);
    for (int x = 0; x < this->getSizeX(); x++)
        for (int y = 0; y < this->getSizeY(); y++)
            for (int h = 0; h < this->getSizeZ(); h++)
                noiseMinMax.update(this->noise.GetNoise((float)x, (float)y, (float)h));

    this->chunks = std::vector<std::shared_ptr<VoxelChunk>>(this->numberOfChunksX() * this->numberOfChunksY());

    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            this->chunks[iChunk] = std::make_shared<VoxelChunk>(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), Matrix3<float>(chunkSize, chunkSize, this->getSizeZ(), 0.f), this);
            this->chunks[iChunk]->lastChunkOnX = (xChunk == this->numberOfChunksX() - 1);
            this->chunks[iChunk]->lastChunkOnY = (yChunk == this->numberOfChunksY() - 1);
            iChunk++;
        }
    }
    for(size_t i = 0; i < this->chunks.size(); i++) {
        if (int(i) > this->numberOfChunksY() - 1) {
            this->chunks[i]->neighboring_chunks[LEFT] = this->chunks[i - int(this->getSizeY() / chunkSize)];
            this->chunks[i - int(this->getSizeY()/chunkSize)]->neighboring_chunks[RIGHT] = this->chunks[i];
        }
        if (i % this->numberOfChunksY() >= 1) {
            this->chunks[i]->neighboring_chunks[FRONT] = this->chunks[i - 1];
            this->chunks[i - 1]->neighboring_chunks[BACK] = this->chunks[i];
        }
        this->chunks[i]->updateLoDsAvailable();
        this->chunks[i]->LoDIndex = 1; // std::min(i % this->numberOfChunksY() + i / this->numberOfChunksX(), this->chunks[i]->LoDs.size() - 1);
    }
*//*
    float dt = 0.1f;
    float diffusion = 0.8f;
    float viscosity = 0.01f;
    int fluidSolverIterations = 5;
    this->fluidSimRescale = Vector3(4, 4, 4);
    this->fluidSimulation = StableFluidsSimulation(this->getSizeX() / this->fluidSimRescale.x, this->getSizeY() / this->fluidSimRescale.y, this->getSizeZ() / this->fluidSimRescale.z, dt, diffusion, viscosity, fluidSolverIterations);
    this->environmentalDensities = Matrix3<float>(this->getDimensions(), 1); // Fill with air density for now

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
    };*/
}

void VoxelGrid::makeItFall(float erosionStrength)
{
    /*
    computeVoxelGroups();
    for(int i = 0; i < 1; i++) {
        for(std::shared_ptr<VoxelChunk>& vc : this->chunks) {
            vc->makeItFall();
        }
    }
//    remeshAll();
    if (erosionStrength > 0.0) {
//        UnderwaterErosion erod(this, 10, erosionStrength, 100); // (this->shared_from_this(), 10, erosionStrength, 100);
//        float _a, _b;
//        erod.Apply(0, _a, _b);
    }
    */
}
void VoxelGrid::letGravityMakeSandFall(bool remesh)
{
    auto start = std::chrono::system_clock::now();/*
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks) {
        vc->letGravityMakeSandFall();
    }*/
    this->letGravityMakeSandFallWithFlow();
    this->applyModification(this->shareSandWithNeighbors());
//    if (remesh)
//        remeshAll();
    auto duration = std::chrono::duration<float>(std::chrono::system_clock::now() - start);
    std::cout << duration.count() << " s for making sand fall once" << std::endl;
}
void VoxelGrid::letGravityMakeSandFallWithFlow(bool remesh)
{
    Matrix3<float> transportMatrix(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    Matrix3<float> currentTransportMatrix(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    Matrix3<float> voxelValues = this->getVoxelValues();
    Matrix3<Vector3> flow(this->getSizeX(), this->getSizeY(), this->getSizeZ(), this->sea_current); // To change with real flow
    float gravityForce = 1.f;
//    Matrix3<Vector3> flow = this->flowField;
    float sandLowerLimit = 0.0, sandUpperLimit = 1.0;

    int iter = 0;
    do {
        currentTransportMatrix = Matrix3<float>(this->getSizeX(), this->getSizeY(), this->getSizeZ());
        std::vector<Vector3> sandyPositions;
        for (size_t i = 0; i < voxelValues.size(); i++) {
            if (sandLowerLimit <= voxelValues[i] + transportMatrix[i] && voxelValues[i] + transportMatrix[i] < sandUpperLimit) {
                sandyPositions.push_back(voxelValues.getCoordAsVector3(i));
            }
        }
        std::sort(sandyPositions.begin(), sandyPositions.end(),
                  [](const Vector3& a, const Vector3& b) -> bool { return a.z < b.z; });

        int a = 0;
        for (const Vector3& sandPos : sandyPositions) {
            // If some sand has fall there, don't update... Maybe.. Maybe not... Maybe remove it
            if (voxelValues.at(sandPos) + transportMatrix.at(sandPos) + currentTransportMatrix.at(sandPos) >= sandUpperLimit)
                continue;
            Vector3 currentPos = sandPos;
            Vector3 currentDir(0, 0, -1.0);
            while ((currentPos + currentDir).rounded() != sandPos.rounded() && this->contains(currentPos + currentDir) &&
                   voxelValues.at(currentPos + currentDir) + transportMatrix.at(currentPos + currentDir) + currentTransportMatrix.at(currentPos + currentDir) <= 0) {
                currentPos += currentDir;
                currentDir += flow.at(currentPos) + Vector3(0, 0, -1 * gravityForce);
                currentDir.normalize();
            }
            a++;
            // If it goes away
            if (!this->contains(currentPos)) {
                currentTransportMatrix.at(sandPos) -= voxelValues.at(sandPos) + transportMatrix.at(sandPos) + currentTransportMatrix.at(sandPos);
                continue;
            } else {
                /*if (voxelValues.at(currentPos) + transportMatrix.at(currentPos) + currentTransportMatrix.at(currentPos) > sandUpperLimit)
                    currentPos -= currentDir; // Go back to previous pos if it saturates the landing voxel
                float transportedValue = voxelValues.at(sandPos) + transportMatrix.at(sandPos) + currentTransportMatrix.at(sandPos);
                currentTransportMatrix.at(currentPos) += transportedValue;
                currentTransportMatrix.at(sandPos) -= transportedValue;*/
                if (voxelValues.at(currentPos) + transportMatrix.interpolate(currentPos) + currentTransportMatrix.interpolate(currentPos) > sandUpperLimit)
                    currentPos -= currentDir; // Go back to previous pos if it saturates the landing voxel
                float transportedValue = voxelValues.at(sandPos) + transportMatrix.at(sandPos) + currentTransportMatrix.at(sandPos);
                currentTransportMatrix.addValueAt(transportedValue, currentPos);
                currentTransportMatrix.addValueAt(transportedValue, sandPos);
            }
        }
        transportMatrix += currentTransportMatrix;
        std::cout << currentTransportMatrix.abs().sum() << " quantity moved at iteration " << ++iter << std::endl;
    } while(currentTransportMatrix.abs().sum() > 1e5);
    this->applyModification(transportMatrix);
}

Matrix3<float> VoxelGrid::shareSandWithNeighbors()
{
    Matrix3<float> allVoxels = this->getVoxelValues();
    Matrix3<float> transport(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    allVoxels.raiseErrorOnBadCoord = false;
    allVoxels.defaultValueOnBadCoord = std::numeric_limits<float>::max();
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            for (int z = 0; z < this->getSizeZ(); z++) {
                // Maybe check if val is movable (0 < val < 1.0)
                float voxelValue = allVoxels.at(x, y, z);
                if (voxelValue <= 0.0 || 1.0 <= voxelValue) continue;
                float coefficient_to_share = 0.f;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        if (dx == 0 && dy == 0) continue;
                        float neighbor = allVoxels.at(x+dx, y+dy, z);
                        if (0.0 <= neighbor && neighbor < 1.0) {
                            coefficient_to_share += 1/sqrt((float)(dx*dx)+(dy*dy));
                        }
                    }
                }
                float value_to_share = voxelValue / (coefficient_to_share + 1.0);
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        if (dx == 0 && dy == 0) continue;
                        float neighbor = allVoxels.at(x+dx, y+dy, z);
                        if (0.0 <= neighbor && neighbor < 1.0) {
                            transport.at(x+dx, y+dy, z) += value_to_share / sqrt((float)(dx*dx)+(dy*dy)) * value_to_share;
                        }
                    }
                }
                transport.at(x, y, z) -= value_to_share * coefficient_to_share;
            }
        }
    }
    return transport;
}

void VoxelGrid::applyModification(Matrix3<float> modifications, Vector3 anchor)
{
    if (currentHistoryIndex < this->voxelsValuesStack.size()) {
        this->voxelsValuesStack.erase(this->voxelsValuesStack.begin() + currentHistoryIndex, this->voxelsValuesStack.end());
        this->voxelsValuesAnchorStack.erase(this->voxelsValuesAnchorStack.begin() + currentHistoryIndex, this->voxelsValuesAnchorStack.end());
    }
    this->currentHistoryIndex = std::max(1, this->currentHistoryIndex + 1);

    this->voxelsValuesStack.push_back(modifications);
    this->voxelsValuesAnchorStack.push_back(anchor);
}

void VoxelGrid::add2DHeightModification(Matrix3<float> heightmapModifier, float factor, Vector3 anchor)
{
    /// Two possibilities : either inverse the voxel values when needed, or just add random values based on Perlin noise
    Matrix3<float> modification(this->getDimensions(), 0.f);
    Matrix3<float> previousValues = this->getVoxelValues();

    /// Third possibility : Apply deformation on the Z-axis
    Matrix3<Vector3> deformation(this->getDimensions());
    float maxDepthAllowed = this->getSizeZ();
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float currentHeight = this->getHeight(x, y);
            float heightmapValue = heightmapModifier.at(x, y);
            for (int z = 0; z < this->getSizeZ(); z++) {
                float coef = z <= currentHeight ?
                            1 - (currentHeight - z) / currentHeight :
                            1 - (z - currentHeight) / (maxDepthAllowed - currentHeight);
                deformation.at(x, y, z) = Vector3(0, 0, heightmapValue * coef);
            }
        }
    }
    Matrix3<float> newVoxels = previousValues.wrapWith(deformation) - previousValues;

    this->applyModification(newVoxels, anchor);

}

bool VoxelGrid::undo()
{
    if (this->voxelsValuesStack.size() > 1 && this->currentHistoryIndex > 0) {
        this->currentHistoryIndex --;
        return true;
    }
    return false;
    /*
    bool effective = false;
    for (auto& vc : this->chunks)
        effective = vc->undo();
    return effective;
    */
}

bool VoxelGrid::redo()
{
    if (this->currentHistoryIndex < this->voxelsValuesStack.size()) {
        this->currentHistoryIndex ++;
        return true;
    }
    return false;
    /*
    bool effective = false;
    for (auto& vc : this->chunks)
        effective = vc->redo();
    return effective;
    */
}

size_t VoxelGrid::getCurrentHistoryIndex() const
{
    return this->currentHistoryIndex;
//    return this->chunks.front()->currentHistoryIndex;
}

float VoxelGrid::getHeight(float x, float y) {
    for (int z = this->getSizeZ() - 1; z >= 0; z--)
        if (this->getVoxelValue(x, y, z) > 0.f)
            return z;
    return 0;
}

bool VoxelGrid::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelGrid::contains(float x, float y, float z) {
    return (0 <= x && x < this->getSizeX() && 0 <= y && y < this->getSizeY() && 0 <= z && z < this->getSizeZ());
}

bool VoxelGrid::checkIsInGround(Vector3 position)
{
    return this->_cachedVoxelValues.at(position) > 0.f;
//    return this->getVoxelValues().at(position) > 0.f;
//    return this->getVoxelValue(position) > 0.f;
}

void VoxelGrid::limitVoxelValues(float limitedValue)
{
    auto values = this->getVoxelValues();
    for (int z = 0; z < values.sizeZ - 1; z++) {
        #pragma omp parallel for collapse(2)
        for (int x = 0; x < values.sizeX; x++) {
            for (int y = 0; y < values.sizeY; y++) {
                if (values.at(x, y, z) > limitedValue) {
                    values.at(x, y, z + 1) += values.at(x, y, z) - limitedValue;
                    values.at(x, y, z) = limitedValue;
                }
            }
        }
    }
    this->applyModification(values - this->getVoxelValues()); // Add the difference with initial values
}

void VoxelGrid::computeVoxelGroups()
{
    int currentMarker = 0;
    std::set<int> neighbors_coords;
    Matrix3<int> connected(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    Matrix3<int> labels(this->getSizeX(), this->getSizeY(), this->getSizeZ(), -1);
    Matrix3<float> voxelValues = this->getVoxelValues();

    for (size_t i = 0; i < voxelValues.size(); i++) {
        if (voxelValues[i] > 0.f) connected[i] = 1;
    }
    connected.raiseErrorOnBadCoord = false; // Allow to search on out-of-bounds
    connected.defaultValueOnBadCoord = 0; // But mark it as background

    for (size_t i = 0; i < connected.size(); i++) {
        if (connected.at(i) == 1) {
            currentMarker++;
            neighbors_coords.insert(i);
            while (!neighbors_coords.empty()) {
                int nx, ny, nz;
                size_t i_neighbor = *neighbors_coords.begin();
                neighbors_coords.erase(neighbors_coords.begin());
                labels.at(i_neighbor) = currentMarker;
                std::tie(nx, ny, nz) = connected.getCoord(i_neighbor);
                if (connected.at(i_neighbor) == 1) {
                    connected.at(i_neighbor) = 0;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (connected(nx + dx, ny + dy, nz + dz) == 1) {
                                    neighbors_coords.insert(connected.getIndex(nx + dx, ny + dy, nz + dz));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

//    for (auto& vc : this->chunks)
//        vc->voxelGroups = labels.subset(vc->x, vc->x + vc->sizeX, vc->y, vc->y + vc->sizeY, 0, 0 + vc->sizeZ);
}

Matrix3<float> VoxelGrid::getVoxelValues()
{
    if (this->_cachedHistoryIndex != int(this->getCurrentHistoryIndex())) {
        this->_cachedHistoryIndex = this->getCurrentHistoryIndex();

        this->_cachedVoxelValues = Matrix3<float>(this->getSizeX(), this->getSizeY(), this->getSizeZ());
        this->_cachedVoxelValues.raiseErrorOnBadCoord = false;
        for (int i = 0; i < this->currentHistoryIndex; i++)
            _cachedVoxelValues.add(this->voxelsValuesStack[i], this->voxelsValuesAnchorStack[i]);
        /*
         for (auto& vc : this->chunks) {
            this->_cachedVoxelValues.paste(vc->getVoxelValues(), vc->x, vc->y, 0);
        }
        */
//        this->updateLayersRepresentation();
    }
    return this->_cachedVoxelValues;
}

Matrix3<float> VoxelGrid::getVoxelized(Vector3 dimensions, Vector3 scale)
{
    return this->getVoxelValues();
}

Mesh VoxelGrid::getGeometry(Vector3 dimensions)
{
    Vector3 originalDimensions = this->getDimensions();
    if (!dimensions.isValid())
        dimensions = originalDimensions;
    Mesh m = Mesh::applyMarchingCubes(this->getVoxelValues().resize(dimensions).meanSmooth(5, 5, 5));
    m.scale(originalDimensions / dimensions);
    return m;
    /*
    auto triTable = MarchingCubes::triangleTable;
//    auto edges = MarchingCubes::cubeEdges;
    auto values = this->getVoxelValues().meanSmooth(5, 5, 5);
    values.defaultValueOnBadCoord = -1;

    float offsetX = 0.f;
    float offsetY = 0.f;
    float offsetZ = 0.f;
    Vector3 scale(1.f, 1.f, 1.f);
    float isolevel = 0.f;

    Vector3 vertDecals[8] = {
        Vector3(0.0, 0.0, 0.0),
        Vector3(1.0, 0.0, 0.0),
        Vector3(1.0, 1.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, 1.0),
        Vector3(1.0, 0.0, 1.0),
        Vector3(1.0, 1.0, 1.0),
        Vector3(0.0, 1.0, 1.0)
    };
    std::function cubePos = [&](Vector3 voxelPos, int i) -> Vector3 {
        return voxelPos + vertDecals[i];
    };

    //Get vertex i value within current marching cube
    std::function cubeVal = [&](Vector3 pos) -> float {
        if (!values.checkCoord(pos)) return -1.f;
        return values.at(pos);
    };
    //Get vertex i value within current marching cube
    std::function cubeVali = [&](Vector3 voxelPos, int i) -> float {
        return cubeVal(cubePos(voxelPos, i));
    };

    //Get triangle table value
    std::function triTableValue = [&](int i, int j) -> int{
        return triTable[i][j];
    };

    //Compute interpolated vertex along an edge
    std::function vertexInterp = [&](float isolevel, Vector3 v0, float l0, Vector3 v1, float l1) -> Vector3 {
        float iso = std::clamp((isolevel-l0)/(l1-l0), 0.f, 1.f);
        return v0 * (1.f - iso) + v1 * iso;
//        return mix(v0, v1, clamp() * scale + vec3(offsetX, offsetY, offsetZ);
    };

    std::function getPosition = [&](Vector3 position, Vector3 _offset) -> Vector3 {
    //    return position + vec4(_offset, 0.0);
        _offset += Vector3(offsetX, offsetY, offsetZ);
        position *= scale;

//        float distToLimits = (voxels_displayed_on_borders > 1 ? min(mincomp(abs(position.xyz - min_vertice_positions)), mincomp(abs(position.xyz + vec3(1.0) - max_vertice_positions))) : 1.0);
        Vector3 off = _offset * 1.f; // (clamp(distToLimits / float(voxels_displayed_on_borders), 0.0, 1.0));
        return position + off; //clamp(position + vec4(off, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
    //    return clamp (position + vec4(_offset, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
    };

    std::function getCubeIndex = [&](Vector3 voxPos, Vector3 normal) -> int {
        int cubeindex = 0;
        float cubeVal0 = cubeVali(voxPos, 0);
        float cubeVal1 = cubeVali(voxPos, 1);
        float cubeVal2 = cubeVali(voxPos, 2);
        float cubeVal3 = cubeVali(voxPos, 3);
        float cubeVal4 = cubeVali(voxPos, 4);
        float cubeVal5 = cubeVali(voxPos, 5);
        float cubeVal6 = cubeVali(voxPos, 6);
        float cubeVal7 = cubeVali(voxPos, 7);
        float refined_isolevel = isolevel + 0.0001;
        //Determine the index into the edge table which
        //tells us which vertices are inside of the surface
        cubeindex  = int(cubeVal0 < refined_isolevel);
        cubeindex += int(cubeVal1 < refined_isolevel)*2;
        cubeindex += int(cubeVal2 < refined_isolevel)*4;
        cubeindex += int(cubeVal3 < refined_isolevel)*8;
        cubeindex += int(cubeVal4 < refined_isolevel)*16;
        cubeindex += int(cubeVal5 < refined_isolevel)*32;
        cubeindex += int(cubeVal6 < refined_isolevel)*64;
        cubeindex += int(cubeVal7 < refined_isolevel)*128;

        normal = Vector3(0, 0, 0);

        if (cubeindex != 0 && cubeindex != 255) {
            Vector3 vertlist[12];

            //Find the vertices where the surface intersects the cube
            vertlist[0] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 1), cubeVal1);
            vertlist[1] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 2), cubeVal2);
            vertlist[2] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 3), cubeVal3);
            vertlist[3] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 0), cubeVal0);
            vertlist[4] = vertexInterp(refined_isolevel, cubePos(voxPos, 4), cubeVal4, cubePos(voxPos, 5), cubeVal5);
            vertlist[5] = vertexInterp(refined_isolevel, cubePos(voxPos, 5), cubeVal5, cubePos(voxPos, 6), cubeVal6);
            vertlist[6] = vertexInterp(refined_isolevel, cubePos(voxPos, 6), cubeVal6, cubePos(voxPos, 7), cubeVal7);
            vertlist[7] = vertexInterp(refined_isolevel, cubePos(voxPos, 7), cubeVal7, cubePos(voxPos, 4), cubeVal4);
            vertlist[8] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 4), cubeVal4);
            vertlist[9] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 5), cubeVal5);
            vertlist[10] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 6), cubeVal6);
            vertlist[11] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 7), cubeVal7);


//            vec3 edge1 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 1)];
//            vec3 edge2 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 2)];
//            normal = normalize(cross(edge1, edge2));
        }
        return cubeindex;
    };

    Mesh marched;
    marched.useIndices = false;
    float refined_isolevel = isolevel + 0.0001;
    Vector3 vertlist[12];
    for (int x = -1; x < values.sizeX; x++) {
        for (int y = -1; y < values.sizeY; y++) {
            for (int z = -1; z < values.sizeZ; z++) {
                Vector3 position = Vector3(x, y, z);
                Vector3 voxPos = position;

                float cubeVal0 = cubeVali(voxPos, 0);
                float cubeVal1 = cubeVali(voxPos, 1);
                float cubeVal2 = cubeVali(voxPos, 2);
                float cubeVal3 = cubeVali(voxPos, 3);
                float cubeVal4 = cubeVali(voxPos, 4);
                float cubeVal5 = cubeVali(voxPos, 5);
                float cubeVal6 = cubeVali(voxPos, 6);
                float cubeVal7 = cubeVali(voxPos, 7);

                Vector3 normal;
                int cubeindex = getCubeIndex(voxPos, normal);

                //Cube is entirely in/out of the surface
                if (cubeindex == 0 || cubeindex == 255)
                    continue;

                //Find the vertices where the surface intersects the cube
                vertlist[0] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 1), cubeVal1);
                vertlist[1] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 2), cubeVal2);
                vertlist[2] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 3), cubeVal3);
                vertlist[3] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 0), cubeVal0);
                vertlist[4] = vertexInterp(refined_isolevel, cubePos(voxPos, 4), cubeVal4, cubePos(voxPos, 5), cubeVal5);
                vertlist[5] = vertexInterp(refined_isolevel, cubePos(voxPos, 5), cubeVal5, cubePos(voxPos, 6), cubeVal6);
                vertlist[6] = vertexInterp(refined_isolevel, cubePos(voxPos, 6), cubeVal6, cubePos(voxPos, 7), cubeVal7);
                vertlist[7] = vertexInterp(refined_isolevel, cubePos(voxPos, 7), cubeVal7, cubePos(voxPos, 4), cubeVal4);
                vertlist[8] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 4), cubeVal4);
                vertlist[9] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 5), cubeVal5);
                vertlist[10] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 6), cubeVal6);
                vertlist[11] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 7), cubeVal7);

                int i = 0;
                while(true){
                    if(triTableValue(cubeindex, i) != -1){
                        position = vertlist[triTableValue(cubeindex, i+0)];
                        marched.vertexArray.push_back(position);

                        position = vertlist[triTableValue(cubeindex, i+1)];
                        marched.vertexArray.push_back(position);

                        position = vertlist[triTableValue(cubeindex, i+2)];
                        marched.vertexArray.push_back(position);
                    }else{
                        break;
                    }

                    i = i + 3;
                }
            }
        }
    }
    return marched;*/
}
/*
std::tuple<int, int, int, int> VoxelGrid::getChunksAndVoxelIndices(Vector3 pos) {
    return getChunksAndVoxelIndices(pos.x, pos.y, pos.z);
}
std::tuple<int, int, int, int> VoxelGrid::getChunksAndVoxelIndices(float x, float y, float z) {

    if(z < 0 || x < 0 || y < 0)
        return std::make_tuple(-1, -1, -1, -1); // -1 = not good coord
    int _x = x, _y = y, _z = z;
    int xChunk = int(_x / chunkSize);
    int voxPosX = _x % chunkSize;
    int yChunk = _y / chunkSize;
    int voxPosY = _y % chunkSize;

    if (xChunk < 0 || yChunk < 0 || _x >= getSizeX() || _y >= getSizeY() || _z >= this->getSizeZ())
        return std::make_tuple(-1, -1, -1, -1); // -1 = not good coord
    return std::make_tuple(xChunk * (this->getSizeY() / chunkSize) + yChunk, voxPosX, voxPosY, _z);
}*/
float VoxelGrid::getVoxelValue(Vector3 pos) {
    return this->getVoxelValue(pos.x, pos.y, pos.z);
}

float VoxelGrid::getVoxelValue(float x, float y, float z) {
    /*int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->getVoxelValue(voxPosX, voxPosY, _z);
    return -1;*/
    return this->getVoxelValues().at(x, y, z);
}
void VoxelGrid::setVoxelValue(Vector3 pos, float newVal) {
    this->setVoxelValue(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setVoxelValue(float x, float y, float z, float newVal)
{
    this->applyModification(Matrix3<float>(1, 1, 1, -this->getVoxelValue(x, y, z) + newVal), Vector3(x, y, z));
    /*
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1) {
        Matrix3<float> setterMatrix(this->chunks[iChunk]->sizeX, this->chunks[iChunk]->sizeY, this->chunks[iChunk]->sizeZ);
        setterMatrix.at(voxPosX, voxPosY, _z) = -this->chunks[iChunk]->getVoxelValue(voxPosX, voxPosY, _z) + newVal;
        this->chunks[iChunk]->applyModification(setterMatrix);
//        this->chunks[iChunk]->voxelValues(voxPosX, voxPosY, _z) = newVal;
        this->chunks[iChunk]->needRemeshing = true;
    }
    */
}
/*
float VoxelGrid::getOriginalVoxelValue(Vector3 pos) {
    return this->getOriginalVoxelValue(pos.x, pos.y, pos.z);
}

float VoxelGrid::getOriginalVoxelValue(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1 && !this->chunks[iChunk]->voxelsValuesStack.empty())
        return this->chunks[iChunk]->voxelsValuesStack[0].at(voxPosX, voxPosY, _z);
    return -1;
}*//*
Matrix3<Vector3> VoxelGrid::getFlowfield()
{
    return this->getFlowfield(0);
}

Matrix3<Vector3> VoxelGrid::getFlowfield(size_t flowIndex)
{
    Matrix3<float> binary = this->getVoxelValues().binarize(0.f, false);
//    return Matrix3<Vector3>(this->getDimensions(), Vector3(1, 0, 0)) * binary;
    return this->multipleFlowFields[flowIndex] * binary;
}*/
/*
Vector3 VoxelGrid::getFlowfield(Vector3 pos) {
    return this->getFlowfield(pos.x, pos.y, pos.z);
}

Vector3 VoxelGrid::getFlowfield(float x, float y, float z) {
    return this->getFlowfield().at(x, y, z);
}
void VoxelGrid::setFlowfield(Vector3 pos, Vector3 newVal) {
    this->setFlowfield(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setFlowfield(float x, float y, float z, Vector3 newVal)
{
    this->multipleFlowFields[0].at(x, y, z) = newVal;
}
*/
/*
int VoxelGrid::getVoxelGroup(Vector3 pos) {
    return getVoxelGroup(pos.x, pos.y, pos.z);
}
int VoxelGrid::getVoxelGroup(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z);
    return -1;
}*//*
void VoxelGrid::setVoxelGroup(Vector3 pos, int newVal) {
    this->setVoxelGroup(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setVoxelGroup(float x, float y, float z, int newVal)
{
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z) = newVal;
}*//*
bool VoxelGrid::getVoxelIsOnGround(Vector3 pos) {
    return getVoxelIsOnGround(pos.x, pos.y, pos.z);
}
bool VoxelGrid::getVoxelIsOnGround(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z) == 0;
    return -1;
}
void VoxelGrid::setVoxelIsOnGround(Vector3 pos, bool newVal) {
    this->setVoxelIsOnGround(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setVoxelIsOnGround(float x, float y, float z, bool newVal)
{
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z) = (newVal ? 0 : -1);
}
*/


/*
int VoxelGrid::getMaxLoD()
{
    if (chunks.size() > 0 && chunks[0]) {
        int maxLoD = this->getSizeX();
        for (std::shared_ptr<VoxelChunk>& vc : this->chunks)
            if (int(vc->LoDs.size() - 1) < maxLoD)
                maxLoD = vc->LoDs.size() - 1;
        return maxLoD;
    } else {
        return 8;
    }
}
*/
void VoxelGrid::saveState()
{
    Matrix3<float> voxelValues = this->getVoxelValues();
    this->voxelsValuesStack.clear();
    this->voxelsValuesAnchorStack.clear();
    this->currentHistoryIndex = 0;
    /*for (auto& vc : chunks) {
        vc->voxelsValuesStack.clear();
        vc->voxelsValuesAnchorStack.clear();
        vc->currentHistoryIndex = 0;
    }*/
    this->applyModification(voxelValues);

}

void VoxelGrid::saveMap(std::string filename)
{
    std::ofstream out;
    out.open(filename);
    out << this->getSizeX() << " " << this->getSizeY() << " " << this->getSizeZ() << "\n";
    auto values = this->getVoxelValues();
    for (size_t i = 0; i < values.size(); i++)
        out << values[i] << " ";
    /*
    out << this->getSizeX() << " " << this->getSizeY() << " " << this->getSizeZ() << " " << this->chunkSize << "\n";
    for (auto& vc : chunks) {
        for (int x = 0; x < vc->sizeX; x++)
            for (int y = 0; y < vc->sizeY; y++)
                for (int z = 0; z < vc->sizeZ; z++)
                    out << vc->getVoxelValue(x, y, z) << " ";
    }*/
    out.close();
}
/*
void VoxelGrid::retrieveMap(std::string filename)
{
    std::ifstream in;
    in.open(filename);
    if (in.fail()) {
        std::cerr << "Unable to open file " << filename << "..." << std::endl;
        return;
    }
    int _x, _y, _z;
    int _chunkSize;
    in >> _x >> _y >> _z >> _chunkSize;
    this->chunkSize = _chunkSize;
    Vector3 finalSize = this->getDimensions();
    this->_cachedVoxelValues = Matrix3<float>(_x, _y, _z, 0.f);
    initMap();

    if (_chunkSize > 0) {
        float map_val;
        int iChunk = 0;
        for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
            for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
                Vector3 offset(xChunk * chunkSize, yChunk * chunkSize, 0.f);
    //            this->chunks[iChunk] = std::make_shared<VoxelChunk>(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), Matrix3<float>(chunkSize, chunkSize, _z), this); //->shared_from_this());
    //            this->chunks[iChunk]->lastChunkOnX = (xChunk == this->numberOfChunksX() - 1);
    //            this->chunks[iChunk]->lastChunkOnY = (yChunk == this->numberOfChunksY() - 1);
    //            Matrix3<float> values = Matrix3<float>(this->chunkSize, this->chunkSize, _z);
                for (int x = 0; x < chunkSize; x++) {
                    for (int y = 0; y < chunkSize; y++) {
                        for (int z = 0; z < this->getSizeZ(); z++) {
                            in >> map_val;
                            this->_cachedVoxelValues.at(Vector3(x, y, z) + offset) = map_val;
                        }
                    }
                }
    //            this->applyModification(values);
            }
        }
    } else {
        for (size_t i = 0; i < _cachedVoxelValues.size(); i++)
            in >> _cachedVoxelValues[i];
    }
    finalSize.z *= 2.f;
    _cachedVoxelValues.resize(finalSize);
    this->fromCachedData();
}*/

void VoxelGrid::retrieveMap(std::string filename)
{
    std::ifstream in;
    in.open(filename);
    if (in.fail()) {
        std::cerr << "Unable to open file " << filename << "..." << std::endl;
        return;
    }
    int _x, _y, _z;
    int _chunkSize;
    in >> _x >> _y >> _z >> _chunkSize;
    int chunkSize = _chunkSize;
    Vector3 finalSize = this->getDimensions();
    this->_cachedVoxelValues = Matrix3<float>(_x, _y, _z, 0.f);
    initMap();

    if (_chunkSize > 0) { // Compatibility with previous terrain saving system
        float map_val;
        int iChunk = 0;
        for (int xChunk = 0; xChunk < std::ceil(this->getSizeX() / (float)chunkSize); xChunk++) {
            for (int yChunk = 0; yChunk < std::ceil(this->getSizeY() / (float)chunkSize); yChunk++) {
                Vector3 offset(xChunk * chunkSize, yChunk * chunkSize, 0.f);
                for (int x = 0; x < chunkSize; x++) {
                    for (int y = 0; y < chunkSize; y++) {
                        for (int z = 0; z < this->getSizeZ(); z++) {
                            in >> map_val;
                            this->_cachedVoxelValues.at(Vector3(x, y, z) + offset) = map_val;
                        }
                    }
                }
            }
        }
    } else {
        for (size_t i = 0; i < _cachedVoxelValues.size(); i++)
            in >> _cachedVoxelValues[i];
    }
    finalSize.z *= 2.f;
    _cachedVoxelValues.resize(finalSize);
    this->fromCachedData();
}

Vector3 VoxelGrid::getFirstIntersectingVoxel(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    if (!minPos.isValid()) minPos = Vector3();
    if (!maxPos.isValid()) maxPos = this->getDimensions();

    Vector3 currPos = origin;
    auto values = this->getVoxelValues();
    values.raiseErrorOnBadCoord = false;
    float distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
    float distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    // Continue while we are in the grid or we are heading towards the grid
    bool beenInBox = false;
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {

        if (Vector3::isInBox(currPos, minPos, maxPos)) {
            float isoval = values.at(currPos);
            if (isoval > 0.0) {
                return currPos;
            }
            beenInBox = true;
        }
        currPos += dir;
        distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
        distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    }
    if (beenInBox && Vector3::isInBox(currPos.xy(), minPos.xy(), maxPos.xy())) {
        currPos.z = 0;
        return currPos;
    }
    return Vector3(false);
}

Vector3 VoxelGrid::getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    return this->getFirstIntersectingVoxel(origin, dir, minPos, maxPos);
}

float VoxelGrid::getNoiseValue(int x, int y, int z, float noise_shift)
{
    return noiseMinMax.remap(this->noise.GetNoise((float)x, (float)y, (float)z), -1.0 + noise_shift, 1.0 + noise_shift);
}

#include <sstream>
std::string VoxelGrid::toString()
{
    /*
    std::ostringstream ret;
    ret << "{";
    ret << "\n\t\"sizeX\": \"" << sizeX << "\",\n\t\"sizeY\": \"" << sizeY << "\",\n\t\"sizeZ\": \"" << sizeZ;
    ret << "\",\n\t\"blockSize\": \"" << blockSize << "\",\n\t\"noise_shifting\": \"" << noise_shifting;
    ret << "\",\n\t\"chunkSize\": \"" << chunkSize << "\",\n\t\"numberOfVoxels\": \"" << sizeX*sizeY*sizeZ;
    ret << "\",\n\t\"numberOfChunks\": \"" << chunks.size() << "\"\n}";
    return ret.str();*/
    return "";
}
std::string VoxelGrid::toShortString()
{
    /*
    std::ostringstream ret;
    ret << "X" << sizeX << "Y" << sizeY << "Z" << sizeZ;
    ret << "BS" << blockSize << "NS" << noise_shifting;
    ret << "CS" << chunkSize;
    return ret.str();
    */
    return "";
}

void VoxelGrid::smoothVoxels()
{
    Matrix3<float> voxelValues = this->getVoxelValues();
    this->applyModification(voxelValues.meanSmooth(3, 3, 3) - voxelValues);
}

void VoxelGrid::saveHeightmap(std::string heightmap_filename)
{
    std::string ext = toUpper(getExtension(heightmap_filename));
    int width = this->getSizeX();
    int height = this->getSizeY();
    // To heightmap
    std::vector<float> toFloatData(width*height);
    std::vector<uint8_t> toIntData(width*height);

//    float newHeight = std::max(this->maxHeight, this->heights.max());

    Matrix3<float> heights(width, height);
    for (int x = 0; x < heights.sizeX; x++) {
        for (int y = 0; y < heights.sizeY; y++) {
            for (int z = this->getSizeZ() - 1; z >= 0; z--) {
                if (this->getVoxelValues().at(x, y, z) > 0.f) {
                    heights.at(x, y) = float(z) / float(this->getSizeZ());
                    break;
                }
            }
        }
    }

    toFloatData = heights.data;
    for (size_t i = 0; i < heights.size(); i++) {
        toIntData[i] = toFloatData[i] * 255;
    }
    if (ext == "PNG")
        stbi_write_png(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data(), this->getSizeX() * 1);
    else if (ext == "JPG")
        stbi_write_jpg(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data(), 95);
    else if (ext == "BMP")
        stbi_write_bmp(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data());
    else if (ext == "TGA")
        stbi_write_tga(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toIntData.data());
    else if (ext == "HDR")
        stbi_write_hdr(heightmap_filename.c_str(), this->getSizeX(), this->getSizeY(), 1, toFloatData.data());
    else {
        std::cerr << "Trying to save map without valid extension. Possible extensions :\n\t- png\n\t- jpg\n\t- tga\n\t- bmp\n\t- hdr" << std::endl;
    }
}
