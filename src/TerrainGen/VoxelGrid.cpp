#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/UnderwaterErosion.h"

VoxelGrid::VoxelGrid(int nx, int ny, int nz/*, float blockSize*/, float noise_shifting)
    : /*blockSize(blockSize),*/ noise_shifting(noise_shifting) {
    this->_cachedVoxelValues = Matrix3<float>(nx * chunkSize, ny * chunkSize, nz * chunkSize + 1, 1.f);
    this->initMap();
    /*
    this->getSizeZ() = nz * (chunkSize) + 1; // The Z axis is not "increased" by "nz" as it is not splitted in chunks
//    chunkSize += 1; // More useful for the LoDs // TODO : CORRECT THIS TO BE ABLE TO FILL IN THE LAYERED TERRAIN
    this->getSizeX() = nx * (chunkSize);
    this->getSizeY() = ny * (chunkSize);
    this->initMap();



    std::vector<Matrix3<float>> data(this->chunks.size(), Matrix3<float>(chunkSize, chunkSize, this->getSizeZ()));
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            for (int x = 0; x < chunkSize; x++) {
                for (int y = 0; y < chunkSize; y++) {
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float noise_val = getNoiseValue(xChunk * chunkSize + x, yChunk * chunkSize + y, h, noise_shifting);
                        noise_val = 1.1; // To remove
//                        noise_val = (xChunk == 1 && yChunk == 1 && h / chunkSize  == 1? 1.0 : -1.0); // To remove
//                        noise_val -= (h > 10 ? 1.0 : 0.0); // To remove
//                        noise_val = 0.5; // To remove
//                        noise_val = (h < getSizeZ()/4 ? 1.0 : -1.0); // To remove
                        data[iChunk].at(x, y, h) = noise_val;
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;*/
}

VoxelGrid::VoxelGrid(Heightmap& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), 1.f/*, grid.getTileSize()*/) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(1, 1, 1, 1.0) {

}
VoxelGrid::~VoxelGrid()
{
    this->chunks.clear();
}
void VoxelGrid::from2DGrid(Heightmap grid, Vector3 subsectionStart, Vector3 subsectionEnd, float scaleFactor) {
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
//    this->getSizeX() = (subsectionEnd.x - subsectionStart.x) * scaleFactor; //grid.getSizeX();
//    this->getSizeY() = (subsectionEnd.y - subsectionStart.y) * scaleFactor; //grid.getSizeY();
//    this->getSizeZ() = grid.getMaxHeight() /* * scaleFactor*/ * 2; // Give space for arches or things
    this->initMap();
//    std::cout << "New map size : " << this->getDimensions() << std::endl;

    Matrix3<float> gridHeights = grid.getHeights().subset(subsectionStart.xy(), subsectionEnd.xy()).resize(this->getDimensions());
    gridHeights.raiseErrorOnBadCoord = false;
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            float grid_height = gridHeights.at(x, y) * (this->getSizeZ() / (grid.getMaxHeight() * 2));
            int z = int(grid_height);
            // Add some positive noise for h < height
            int i = 0;
            for (; i < z-2; i++) {
                float noise_val = getNoiseValue(x, y, i);
                _cachedVoxelValues.at(x, y, i) = abs(noise_val);
            }
            for (; i < z; i++) {
                float noise_val = .3f; //= getNoiseValue((xChunk * chunkSize + x), (yChunk * chunkSize + y), i);
                _cachedVoxelValues.at(x, y, i) = abs(noise_val);
            }
            // Add negative noise for h > height
            for (; i < this->getSizeZ(); i++) {
                float noise_val = getNoiseValue(x, y, i);
                _cachedVoxelValues.at(x, y, i) = -abs(noise_val);
            }
        }
    }
    /*
    std::vector<Matrix3<float>> data(this->chunks.size(), Matrix3<float>(this->chunkSize, this->chunkSize, this->getSizeZ()));
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            for (int x = 0; x < chunkSize; x++) {
                for (int y = 0; y < chunkSize; y++) {
                    float grid_height = gridHeights.at((xChunk * chunkSize + x), (yChunk * chunkSize + y)) * (this->getSizeZ() / (grid.getMaxHeight() * 2));
                    int z = int(grid_height);
                    // Add some positive noise for h < height
                    int i = 0;
                    for (; i < z-2; i++) {
                        float noise_val = getNoiseValue((xChunk * chunkSize + x), (yChunk * chunkSize + y), i);
                        data[iChunk].at(x, y, i) = abs(noise_val);
                    }
                    for (; i < z; i++) {
                        float noise_val = .3f; //= getNoiseValue((xChunk * chunkSize + x), (yChunk * chunkSize + y), i);
                        data[iChunk].at(x, y, i) = abs(noise_val);
                    }
                    // Add negative noise for h > height
                    for (; i < this->getSizeZ(); i++) {
                        float noise_val = getNoiseValue((xChunk * chunkSize + x), (yChunk * chunkSize + y), i);
                        data[iChunk].at(x, y, i) = -abs(noise_val);
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;
    this->_smoothingNeeded = true;
    */
}

void VoxelGrid::fromLayerBased(LayerBasedGrid layerBased, int fixedHeight)
{
    this->_cachedVoxelValues = layerBased.voxelize(fixedHeight == -1 ? layerBased.getSizeZ() * 2.f : fixedHeight);
//    this->_cachedVoxelValues = Matrix3<float>(layerBased.getSizeX(), layerBased.getSizeY(), (fixedHeight == -1 ? layerBased.getSizeZ() * 2.f : fixedHeight));
//    this->getSizeX() = layerBased.getSizeX();
//    this->getSizeY() = layerBased.getSizeY();
//    this->getSizeZ() = (fixedHeight == -1 ? layerBased.getSizeZ() * 2.f : fixedHeight);
    this->initMap();
/*
    std::vector<Matrix3<float>> data(this->chunks.size(), Matrix3<float>(this->chunkSize, this->chunkSize, this->getSizeZ()));
    Matrix3<float> precomputedVoxelGrid = layerBased.voxelize(fixedHeight);
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            data[iChunk] = precomputedVoxelGrid.subset((xChunk * chunkSize), (xChunk + 1) * chunkSize, (yChunk * chunkSize), (yChunk + 1) * chunkSize, 0, this->getSizeZ());
            /*for (int x = 0; x < chunkSize; x++) {
                for (int y = 0; y < chunkSize; y++) {
                    for (int z = 0; z < this->getSizeZ(); z++) {
                        TerrainTypes type = layerBased.getValue((xChunk * chunkSize + x), (yChunk * chunkSize + y), z);
                        float density = LayerBasedGrid::densityFromMaterial(type);
                        data[iChunk].at(x, y, z) = density;
                    }
                }
            }*//*
            iChunk ++;
        }
    }
    this->tempData = data;
    this->_smoothingNeeded = false;
    */
}

//std::shared_ptr<VoxelGrid> VoxelGrid::fromIsoData()
VoxelGrid* VoxelGrid::fromIsoData()
{
//    std::vector<Matrix3<float>> isoData = this->tempData;
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
//            this->chunks[iChunk] = std::make_shared<VoxelChunk>(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), isoData[iChunk], this); //->shared_from_this());

            this->chunks[iChunk] = std::make_shared<VoxelChunk>(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), this->_cachedVoxelValues.subset((xChunk * chunkSize), (xChunk + 1) * chunkSize, (yChunk * chunkSize), (yChunk + 1) * chunkSize, 0, this->getSizeZ()), this); //->shared_from_this());
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
    if (this->_smoothingNeeded)  {
        this->smoothVoxels();
        this->_smoothingNeeded = false;
    }
    this->_cachedHistoryIndex = -1; // Force refresh of "getVoxelValues"
//    this->computeFlowfield();
    return this; // ->shared_from_this();
}

void VoxelGrid::computeFlowfield()
{
//    DyeField constantDyeSource(fluidSystem->fullDim);
    Matrix3<float> obstacleMap(this->getSizeX(), this->getSizeY(), this->getSizeZ(), false);
    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            for (int z = 0; z < this->getSizeZ(); z++) {
                obstacleMap(x, y, z) = (this->getVoxelValue(x, y, z) > 0.0 ? 1.f : 0.f);
            }
        }
    }

    this->fluidSimulation.setObstacles(obstacleMap);
    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < this->fluidSimulation.sizeY; y++) {
            for (int z = 0; z < this->fluidSimulation.sizeZ; z++) {
                this->fluidSimulation.velocity(x, y, z) = this->sea_current / (float)this->fluidSimRescale;
            }
        }
    }
    for (int i = 0; i < 30; i++)
        this->fluidSimulation.step();
    std::cout << "Fluid sim step : " << this->fluidSimulation.currentStep << ". Max fluid velocity : " << this->fluidSimulation.velocity.max() << std::endl;
    this->flowField = this->fluidSimulation.getVelocities(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    this->flowField.raiseErrorOnBadCoord = false;

    for (auto& vc : this->chunks) {
        for (int x = vc->x; x < vc->x + vc->sizeX; x++) {
            for (int y = vc->y; y < vc->y + vc->sizeY; y++) {
                for (int z = 0; z < vc->sizeZ; z++) {
                    vc->flowField(x - vc->x, y - vc->y, z - 0) = flowField(x, y, z);
                }
            }
        }
    }
    this->flowField.raiseErrorOnBadCoord = true;
}
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

void VoxelGrid::initMap()
{
    this->chunkSize = std::min(int(this->getSizeX()), this->chunkSize);

//    sizeX -= (sizeX % chunkSize);
//    sizeY -= (sizeY % chunkSize);

    this->chunks.clear();
    this->flowField.clear();
    this->distanceField.clear();
    this->pressureField.clear();

    // Create and configure FastNoise object
    this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    this->noise.SetFrequency(/*blockSize*/ 1.f / (float) this->getSizeX());
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

    float dt = 0.1f;
    float diffusion = 0.1f;
    float viscosity = 0.01f;
    this->fluidSimulation = FluidSimulation(this->getSizeX() / this->fluidSimRescale, this->getSizeY() / this->fluidSimRescale, this->getSizeZ() / this->fluidSimRescale, dt, diffusion, viscosity, 10);

    this->environmentalDensities = Matrix3<float>(this->getDimensions(), 1000); // Fill with water density
}
/*
void VoxelGrid::createMesh()
{
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        vc->needRemeshing = true;
    }
    remeshAll();
}
*/
void VoxelGrid::makeItFall(float erosionStrength)
{
    computeVoxelGroups();
    for(int i = 0; i < 1; i++) {
        for(std::shared_ptr<VoxelChunk>& vc : this->chunks) {
            vc->makeItFall();
        }
    }
//    remeshAll();
    if (erosionStrength > 0.0) {
        UnderwaterErosion erod(this, 10, erosionStrength, 100); // (this->shared_from_this(), 10, erosionStrength, 100);
        erod.Apply();
    }
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
                if (a == 0){
                    std::cout << currentPos << "\n";
                }
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
//    std::cout << transportMatrix.abs().sum() << " quantity of matter moved." << std::endl;
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
//    std::cout << "Full matrix weight : " << modifications.sum() << std::endl;
    for (auto& vc : this->chunks) {
        // Check if the modification is affecting the chunk
        if ((vc->x + vc->sizeX < anchor.x || vc->y + vc->sizeY < anchor.y || vc->z + vc->sizeZ < anchor.z) ||
                (vc->x > anchor.x + modifications.sizeX || vc->y > anchor.y + modifications.sizeY || vc->z > anchor.z + modifications.sizeZ)) {
            vc->applyModification(Matrix3<float>(0, 0, 0));
        } else {
            // We are sure that a part of the matrix is in this chunk
            Vector3 O1 = Vector3(vc->x, vc->y, vc->z);
            Vector3 chunkAnchor = anchor - O1;
            Vector3 O2 = chunkAnchor;
            Vector3 A = Vector3(std::max(chunkAnchor.x, 0.f), std::max(chunkAnchor.y, 0.f), std::max(chunkAnchor.z, 0.f));
            Vector3 B = Vector3(std::min(chunkAnchor.x + modifications.sizeX, (float)vc->sizeX), std::min(chunkAnchor.y + modifications.sizeY, (float)vc->sizeY), std::min(chunkAnchor.z + modifications.sizeZ, (float)vc->sizeZ));
            Vector3 relativeA = A - O2;
            Vector3 relativeB = B - O2;
            Matrix3<float> subset = modifications.subset(relativeA.x, relativeB.x, relativeA.y, relativeB.y, relativeA.z, relativeB.z);
//            std::cout << "Submatrix total weight : " << subset.sum() << std::endl;
            vc->applyModification(subset, A);
//            vc->applyModification(modifications.subset(vc->x, vc->x + this->chunkSize, vc->y, vc->y + this->chunkSize, 0, 0 + this->getSizeZ()));
        }
    }
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
//            float desiredDepth = heightmapValue / maxDepthAllowed;
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
    /*
    return;



    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
//            int currentHeight = this->getHeight(x, y);
            float startingHeight = this->getHeight(x, y);
            float height = heightmapModifier.at(x, y);
//            if (height != 10)
//                std::cout << height << " at (" << x << ", " << y << ")" << std::endl;
//            float start = std::max(0.f, std::min((float)currentHeight, currentHeight + height));
//            float end   = std::min((float)getSizeZ(), std::max((float)currentHeight, currentHeight + height));
            if (height > 1) {
//                std::cout << "Adding matter? " << height << "\n";
                for (int z = std::max((int)(startingHeight-2), 0); z < std::min((int)(startingHeight + height), getSizeZ()); z++) {
                    // Method #1 : inverse the values of the voxels
                    modification.at(x, y, z) = previousValues.at(x, y, z) * -2.f * (z < startingHeight ? -1.f : 1.f);
                    // Method #2 : use the noise values
    //                modification.at(x, y, z) = getNoiseValue(x, y, z);
                }
                if (0 < startingHeight + height && startingHeight + height < getSizeZ()) {
                    // Special case on last voxel : take into account the decimal value of the height
                    modification.at(x, y, startingHeight + height) = previousValues.at(x, y, startingHeight + height) * -2.f * (height - (int)height);
                    // Method #2 : use the noise values
        //            modification.at(x, y, startingHeight + height) = getNoiseValue(x, y, startingHeight + height) * (height - (int)height);
                }
            } else if (height < 0) {
//                std::cout << "With start at " << startingHeight << " and height at " << height << " we go from " << std::max((int)std::ceil(startingHeight + height), 0) << " to " << std::min((int)(startingHeight + 1), getSizeZ()) << "\n";
                for (int z = std::max((int)std::ceil(startingHeight + height), 0); z < std::min((int)(startingHeight + 1), getSizeZ()); z++) {
                    // Method #1 : inverse the values of the voxels
//                    modification.at(x, y, z) = previousValues.at(x, y, z) * -2.f * (z == startingHeight-1 ? -1.f : 1.f);
                    // Method #2 : use the noise values
                    modification.at(x, y, z) = -getNoiseValue(x, y, z);
//                    std::cout << "Remove from " << Vector3(x, y, z) << "\n";
                }
                if (0 < startingHeight + height && startingHeight + height < getSizeZ()) {
                    // Special case on last voxel : take into account the decimal value of the height
//                    modification.at(x, y, startingHeight + height) = previousValues.at(x, y, startingHeight + height) * -2.f * (height - (int)height);
                    // Method #2 : use the noise values
                    modification.at(x, y, startingHeight + height) = -getNoiseValue(x, y, startingHeight + height) * ((int)(height) - height);
                }
            }
        }
    }
    this->applyModification(modification * factor, anchor);
    */
}

bool VoxelGrid::undo()
{
    bool effective = false;
    for (auto& vc : this->chunks)
        effective = vc->undo();
//    this->remeshAll();
    return effective;
}

bool VoxelGrid::redo()
{
    bool effective = false;
    for (auto& vc : this->chunks)
        effective = vc->redo();
//    this->remeshAll();
    return effective;
}

size_t VoxelGrid::getCurrentHistoryIndex() const
{
    return this->chunks.front()->currentHistoryIndex;
}
/*
void VoxelGrid::display() {

    for (std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        vc->display();
    }
}*/

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

//void VoxelGrid::remeshAll()
//{
    /*
    for (std::shared_ptr<VoxelChunk>& vc : this->chunks) {
        vc->createMesh(this->displayWithMarchingCubes);
    }*/
/*
    Matrix3<int> cubeEdges(256, 1);
    for (int i = 0; i < 256; i++)
        cubeEdges[i] = MarchingCubes::cubeEdges[i];
    Matrix3<int> triTable(16, 256);
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < 256; j++)
            triTable.at(i, j) = MarchingCubes::triangleTable[j][i];


    Matrix3<float> isoData = this->getVoxelValues().resize(30, 30, 30);
    std::vector<Vector3> points(isoData.size());
    for (size_t i = 0; i < points.size(); i++) {
        points[i] = isoData.getCoordAsVector3(i);
    }

    mesh.useIndices = false;
    mesh.fromArray(points);
//    std::cout << "Setting voxel grid with " << points.size()*3 << " points because we use " << isoData << std::endl;
    mesh.shader->setTexture3D("dataFieldTex", 0, isoData);
    mesh.shader->setTexture2D("edgeTableTex", 1, cubeEdges);
    mesh.shader->setTexture2D("triTableTex", 2, triTable);
    mesh.shader->setFloat("isolevel", 0.f);
    mesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    mesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    mesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    mesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    mesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    mesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    mesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    mesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
    */
//}

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

    for (auto& vc : this->chunks)
        vc->voxelGroups = labels.subset(vc->x, vc->x + vc->sizeX, vc->y, vc->y + vc->sizeY, 0, 0 + vc->sizeZ);
}

Matrix3<float> VoxelGrid::getVoxelValues()
{
    if (this->_cachedHistoryIndex != int(this->getCurrentHistoryIndex())) {
        this->_cachedHistoryIndex = this->getCurrentHistoryIndex();

        this->_cachedVoxelValues = Matrix3<float>(this->getSizeX(), this->getSizeY(), this->getSizeZ());
        for (auto& vc : this->chunks) {
            this->_cachedVoxelValues.paste(vc->getVoxelValues(), vc->x, vc->y, 0);
        }
//        this->updateLayersRepresentation();
    }
    return this->_cachedVoxelValues;
}

Mesh VoxelGrid::getGeometry()
{
    auto triTable = MarchingCubes::triangleTable;
//    auto edges = MarchingCubes::cubeEdges;
    auto values = this->getVoxelValues();

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

    std::function getCubeIndex = [&](Vector3 voxPos, vec3 normal) -> int {
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

        normal = vec3(0, 0, 0);

        if (cubeindex != 0 && cubeindex != 255) {
            vec3 vertlist[12];

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
    for (int x = 0; x < values.sizeX; x++) {
        for (int y = 0; y < values.sizeY; y++) {
            for (int z = 0; z < values.sizeZ; z++) {
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
    return marched;
}

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
}
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
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1) {
        Matrix3<float> setterMatrix(this->chunks[iChunk]->sizeX, this->chunks[iChunk]->sizeY, this->chunks[iChunk]->sizeZ);
        setterMatrix.at(voxPosX, voxPosY, _z) = -this->chunks[iChunk]->getVoxelValue(voxPosX, voxPosY, _z) + newVal;
        this->chunks[iChunk]->applyModification(setterMatrix);
//        this->chunks[iChunk]->voxelValues(voxPosX, voxPosY, _z) = newVal;
        this->chunks[iChunk]->needRemeshing = true;
    }
}

float VoxelGrid::getOriginalVoxelValue(Vector3 pos) {
    return this->getOriginalVoxelValue(pos.x, pos.y, pos.z);
}

float VoxelGrid::getOriginalVoxelValue(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1 && !this->chunks[iChunk]->voxelsValuesStack.empty())
        return this->chunks[iChunk]->voxelsValuesStack[0].at(voxPosX, voxPosY, _z);
    return -1;
}
Matrix3<Vector3> VoxelGrid::getFlowfield()
{
    Matrix3<Vector3> values(this->getSizeX(), this->getSizeY(), this->getSizeZ());
    for (auto& vc : this->chunks) {
        values.paste(vc->flowField, vc->x, vc->y, 0);
    }
    return values;
}

Vector3 VoxelGrid::getFlowfield(Vector3 pos) {
    return this->getFlowfield(pos.x, pos.y, pos.z);
}

Vector3 VoxelGrid::getFlowfield(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->flowField(voxPosX, voxPosY, _z);
    return Vector3();
}
void VoxelGrid::setFlowfield(Vector3 pos, Vector3 newVal) {
    this->setFlowfield(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setFlowfield(float x, float y, float z, Vector3 newVal)
{
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        this->chunks[iChunk]->flowField(voxPosX, voxPosY, _z) = newVal;
}

int VoxelGrid::getVoxelGroup(Vector3 pos) {
    return getVoxelGroup(pos.x, pos.y, pos.z);
}
int VoxelGrid::getVoxelGroup(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z);
    return -1;
}
void VoxelGrid::setVoxelGroup(Vector3 pos, int newVal) {
    this->setVoxelGroup(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setVoxelGroup(float x, float y, float z, int newVal)
{
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        this->chunks[iChunk]->voxelGroups(voxPosX, voxPosY, _z) = newVal;
}
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

void VoxelGrid::saveState()
{
    Matrix3<float> voxelValues = this->getVoxelValues();
    for (auto& vc : chunks) {
        vc->voxelsValuesStack.clear();
        vc->voxelsValuesAnchorStack.clear();
        vc->currentHistoryIndex = 0;
    }
    this->applyModification(voxelValues);
}

void VoxelGrid::saveMap(std::string filename)
{
    std::ofstream out;
    out.open(filename);
    out << this->getSizeX() << " " << this->getSizeY() << " " << this->getSizeZ() << " " << this->chunkSize << "\n";
    for (auto& vc : chunks) {
        for (int x = 0; x < vc->sizeX; x++)
            for (int y = 0; y < vc->sizeY; y++)
                for (int z = 0; z < vc->sizeZ; z++)
                    out << vc->getVoxelValue(x, y, z) << " ";
    }
    out.close();
}
void VoxelGrid::retrieveMap(std::string filename)
{
    std::ifstream in;
    in.open(filename);
    int _x = this->getSizeX(), _y = this->getSizeY(), _z = this->getSizeZ();
    in >> _x >> _y >> _z >> this->chunkSize;
    this->_cachedVoxelValues = Matrix3<float>(_x, _y, _z, 0.f);
    initMap();

    for (int x = 0; x < _x; x++) {
        for (int y = 0; y < _y; y++) {
            for (int z = 0; z < _z; z++) {
                float map_val;
                in >> map_val;
                this->_cachedVoxelValues.at(x, y, z) = map_val;
            }
        }
    }
    /*
    std::vector<Matrix3<float>> data(this->chunks.size(), Matrix3<float>(this->chunkSize, this->chunkSize, this->getSizeZ()));
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            for (int x = 0; x < chunkSize; x++) {
                for (int y = 0; y < chunkSize; y++) {
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float map_val;
                        in >> map_val;
                        data[iChunk].at(x, y, h) = map_val;
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;*/
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
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {
        float isoval = values.at(currPos);
        if (isoval > 0.0) {
            return currPos;
        }
        currPos += dir;
        distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
        distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
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

/*
std::pair<Matrix3<int>, Matrix3<float> > VoxelGrid::getLayersRepresentations()
{
    return {_materialLayers, _heightLayers};

}

int getMaterialPerVoxelValue(float value, const std::vector<float>& materials) {
    int material = 0;
    while (value > materials[material] && material < int(materials.size() - 1))
        material ++;
    return material;
}
void VoxelGrid::updateLayersRepresentation()
{
    std::vector<float> materialsLimits = {.0f, 1.f}; // .5f, 1.5f, 2.5f};

    Matrix3<std::vector<int> > materials(this->getSizeX(), this->getSizeY(), 1, {0}); // Add a "null" material with height 0
    Matrix3<std::vector<float> > heights(this->getSizeX(), this->getSizeY(), 1, {0});
    int maxStackSize = 0;

    auto voxelValues = getVoxelValues();

    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            auto start = std::chrono::system_clock::now();
            int currentMaterial = getMaterialPerVoxelValue(voxelValues.at(x, y, 0), materialsLimits);
            float previousHeight = 0.f;
            for (int z = 0; z < this->getSizeZ(); z++) {
                int newMaterial = getMaterialPerVoxelValue(voxelValues.at(x, y, z), materialsLimits);
                // If there is a change of material or that we reached the top, register
                if (currentMaterial != newMaterial || z == this->getSizeZ() - 1) {
                    materials.at(x, y).push_back(currentMaterial);
                    heights.at(x, y).push_back(z - previousHeight);
                    previousHeight = z;
                    currentMaterial = newMaterial;
                    maxStackSize = std::max(maxStackSize, int(materials.at(x, y).size()));
                }
            }
            auto end = std::chrono::system_clock::now();
//            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
        }
    }

    // Transform the vectors into matrix

    _materialLayers = Matrix3<int>(this->getSizeX(), this->getSizeY(), maxStackSize);
    _heightLayers = Matrix3<float>(this->getSizeX(), this->getSizeY(), maxStackSize);

    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            for (size_t i = 1; i < materials.at(x, y).size(); i++) { // Start at 1 to remove the "null" layer
                _materialLayers.at(x, y, i) = materials.at(x, y)[i];
                _heightLayers.at(x, y, i) = heights.at(x, y)[i];
            }
        }
    }
}
*/


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
