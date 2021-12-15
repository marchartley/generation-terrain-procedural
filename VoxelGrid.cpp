#include "VoxelGrid.h"
#include "UnderwaterErosion.h"

#include "Shader.h"



typedef std::vector<float> DoubleVec;

int findNearestNeighbourIndex(const float ac_dfValue, DoubleVec x)
{
  float lv_dfDistance = 10000;
  int lv_nIndex = -1;

  for (unsigned int i = 0; i < x.size(); i++) {
    float newDist = ac_dfValue - x[i];
    if (newDist >= 0 && newDist < lv_dfDistance) {
      lv_dfDistance = newDist;
      lv_nIndex = i;
    }
  }

  return lv_nIndex;
}

DoubleVec interpolation(DoubleVec x, DoubleVec y,
  DoubleVec xx)
{
  float dx, dy;
  DoubleVec slope, intercept, result;
  slope.resize(x.size());
  intercept.resize(x.size());
  result.resize(xx.size());
  int indiceEnVector;

  for (unsigned i = 0; i < x.size(); i++) {
    if (i < x.size() - 1) {
      dx = x[i + 1] - x[i];
      dy = y[i + 1] - y[i];
      slope[i] = dy / dx;
      intercept[i] = y[i] - x[i] * slope[i];
    }
    else {
      slope[i] = slope[i - 1];
      intercept[i] = intercept[i - 1];
    }
  }

  for (unsigned i = 0; i < xx.size(); i++) {
    indiceEnVector = findNearestNeighbourIndex(xx[i], x);
    if (indiceEnVector != -1) {
      result[i] = slope[indiceEnVector] *
        xx[i] + intercept[indiceEnVector];
    }
    else
      result[i] = 10000;
  }
  return result;
}

std::vector<std::vector<std::vector<Vector3>>> resizeArray(std::vector<std::vector<std::vector<Vector3>>> arr, int newX, int newY, int newZ)
{
    /*int oldX = int(arr.size()), oldY = int(arr[0].size()), oldZ = int(arr[0][0].size());
    float resX = newX/(float)oldX, resY = newY/(float)oldY, resZ = newZ/(float)oldZ;
    float invX = 1/resX, invY = 1/resY, invZ = 1/resZ;
    for (int x = 0; x < newX; x++)
    {
        Vector3 newVal = arr[x/invX] * () + arr[x/invX + 1] * (1 - (x%invX));
    }*/
}






VoxelGrid::VoxelGrid(int nx, int ny, int nz, float blockSize, float noise_shifting)
    : blockSize(blockSize), noise_shifting(noise_shifting) {
    this->sizeZ = nz * (chunkSize) + 1; // The Z axis is not "increased" by "nz" as it is not splitted in chunks
    chunkSize += 1; // More useful for the LoDs
    this->sizeX = nx * (chunkSize);
    this->sizeY = ny * (chunkSize);
    this->initMap();



    std::vector<std::vector<std::vector<std::vector<float>>>> data(this->chunks.size());
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            data[iChunk] = std::vector<std::vector<std::vector<float>>>(this->chunkSize);
            for (int x = 0; x < chunkSize; x++) {
                data[iChunk][x] = std::vector<std::vector<float>>(this->chunkSize);
                for (int y = 0; y < chunkSize; y++) {
                    data[iChunk][x][y] = std::vector<float>(this->getSizeZ());
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)h),
                                                            -1.0 + noise_shifting, 1.0 + noise_shifting);
//                        noise_val = 1.0; // To remove
//                        noise_val = (xChunk == 1 && yChunk == 1 && h / chunkSize  == 1? 1.0 : -1.0); // To remove
//                        noise_val -= (h > 10 ? 1.0 : 0.0); // To remove
                        noise_val = 0.5;
                        data[iChunk][x][y][h] = noise_val;
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;
}
VoxelGrid::VoxelGrid(Grid& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), grid.getTileSize()) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(10, 10, 10, 1.0) {

}
VoxelGrid::~VoxelGrid()
{
    this->chunks.clear();
}
void VoxelGrid::from2DGrid(Grid grid) {
    this->sizeX = grid.getSizeX();
    this->sizeY = grid.getSizeY();
    this->sizeZ = grid.getMaxHeight();
    this->initMap();

    std::vector<std::vector<std::vector<std::vector<float>>>> data(this->chunks.size());
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            data[iChunk] = std::vector<std::vector<std::vector<float>>>(this->chunkSize);
            for (int x = 0; x < chunkSize; x++) {
                data[iChunk][x] = std::vector<std::vector<float>>(this->chunkSize);
                for (int y = 0; y < chunkSize; y++) {
                    data[iChunk][x][y] = std::vector<float>(this->getSizeZ());
                    float grid_height = grid.getHeight(xChunk * chunkSize + x, yChunk * chunkSize + y) * (this->sizeZ / grid.getMaxHeight());
                    int z = int(grid_height);
                    for (int i = 0; i < z; i++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)i),
                                                            -2.0, 2.0);
                        data[iChunk][x][y][i] = abs(noise_val);
                    }
                    for (int i = z; i < this->getSizeZ()+1; i++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)i),
                                                            -2.0, 2.0);
                        data[iChunk][x][y][i] = -abs(noise_val);
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;
//    this->fromIsoData(data);
}

std::shared_ptr<VoxelGrid> VoxelGrid::fromIsoData()
{
    std::vector<std::vector<std::vector<std::vector<float>>>> isoData = this->tempData;
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            this->chunks[iChunk] = std::make_shared<VoxelChunk>(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), isoData[iChunk], this->shared_from_this());
            this->chunks[iChunk]->lastChunkOnX = (xChunk == this->numberOfChunksX() - 1);
            this->chunks[iChunk]->lastChunkOnY = (yChunk == this->numberOfChunksY() - 1);
            iChunk++;
        }
    }
    for(size_t i = 0; i < this->chunks.size(); i++) {
        if (int(i) > this->numberOfChunksY() - 1) {
            this->chunks[i]->neighboring_chunks[LEFT] = this->chunks[i - int(this->sizeY / chunkSize)];
            this->chunks[i - int(this->sizeY/chunkSize)]->neighboring_chunks[RIGHT] = this->chunks[i];
        }
        if (i % this->numberOfChunksY() >= 1) {
            this->chunks[i]->neighboring_chunks[FRONT] = this->chunks[i - 1];
            this->chunks[i - 1]->neighboring_chunks[BACK] = this->chunks[i];
        }
//        this->chunks[i]->resetVoxelsNeighbors();
        this->chunks[i]->updateLoDsAvailable();
        this->chunks[i]->LoDIndex = 1; // std::min(i % this->numberOfChunksY() + i / this->numberOfChunksX(), this->chunks[i]->LoDs.size() - 1);

//        this->chunks[i]->computeFlowfield();
    }
    this->computeFlowfield();
//    this->createMesh();
    return this->shared_from_this();
}

void VoxelGrid::computeFlowfield(Vector3 sea_current)
{
    /*
    this->distanceField = std::vector<std::vector<std::vector<int>>>(this->sizeX, std::vector<std::vector<int>>(this->sizeY, std::vector<int>(this->sizeZ, 9999999)));
    this->pressureField = std::vector<std::vector<std::vector<float>>>(this->sizeX, std::vector<std::vector<float>>(this->sizeY, std::vector<float>(this->sizeZ, 0)));
    this->flowField = std::vector<std::vector<std::vector<Vector3>>>(this->sizeX, std::vector<std::vector<Vector3>>(this->sizeY, std::vector<Vector3>(this->sizeZ)));

    for (auto& vc : this->chunks) {
        vc->computeDistanceField();
        for (int x = vc->x; x < vc->x + vc->sizeX; x++) {
            for (int y = vc->y; y < vc->y + vc->sizeY; y++) {
                for (int z = 0; z < vc->height; z++) {
                    distanceField[x][y][z] = vc->distanceField[x - vc->x][y - vc->y][z - 0];
                    pressureField[x][y][z] = vc->pressureField[x - vc->x][y - vc->y][z - 0];
//                    pressureField[x][y][z] = 0.0;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<Vector3>>> pressionGradient = Vector3::grad(this->pressureField);
    this->flowField = pressionGradient;
    for (size_t x = 0; x < flowField.size(); x++) {
        for (size_t y = 0; y < flowField[x].size(); y++) {
            for (size_t z = 0; z < flowField[x][y].size(); z++) {
                flowField[x][y][z] *= -1.0; // Inverse the gradient
                flowField[x][y][z] += sea_current;
            }
        }
    }*/
    std::vector<std::vector<std::vector<bool>>> obstacleMap = std::vector<std::vector<std::vector<bool>>>(this->sizeX, std::vector<std::vector<bool>>(this->sizeY, std::vector<bool>(this->sizeZ, false)));
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                obstacleMap[x][y][z] = (this->getVoxelValue(x, y, z) > 0.0);
            }
        }
    }
    Vector3 center (5, 45, 45);
    this->fluidSimulation.setObstacles(obstacleMap);/*
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
    float dx = 0.0, dy = 0.0, dz = 0.0;
                this->fluidSimulation.addDensity((center.x + dx)/this->fluidSimRescale, (center.y + dy)/this->fluidSimRescale, (center.z + dz)/this->fluidSimRescale, 100.0);
            }
        }
    }*/
    for (int y = 0; y < this->sizeY-1; y++) {
        for (int z = 0; z < this->sizeZ-1; z++) {
            this->fluidSimulation.addDensity (1, y/this->fluidSimRescale, z/this->fluidSimRescale, 2.0/(float)(this->sizeY * this->sizeZ));
            this->fluidSimulation.addVelocity(1, y/this->fluidSimRescale, z/this->fluidSimRescale, this->sea_current.normalized()/(float)(this->sizeY * this->sizeZ)); //(this->sea_current + Vector3::random()).normalized() * 1.0);
        }
    }
    /*
    this->fluidSimulation.addDensity(5, 28, 60, 10.0);
    this->fluidSimulation.addVelocity(5, 28, 60, this->sea_current);*/
    /*
    for (int x = 0; x < this->sizeX; x++)
        for (int y = 0; y < this->sizeY; y++)
            for (int z = 0; z < this->sizeZ; z++)
                if (this->getVoxelValue(x, y, z) > 0.0) {
                    this->fluidSimulation.addVelocity(x, y, z, this->sea_current);
                    this->fluidSimulation.addDensity(x, y,z, 100.0);
                }*/
    for (int i = 0; i < 2; i++)
        this->fluidSimulation.step();
    this->flowField = this->fluidSimulation.getVelocities(this->sizeX, this->sizeY, this->sizeZ);

    for (auto& vc : this->chunks) {
        for (int x = vc->x; x < vc->x + vc->sizeX; x++) {
            for (int y = vc->y; y < vc->y + vc->sizeY; y++) {
                for (int z = 0; z < vc->height; z++) {
                    vc->flowField[x - vc->x][y - vc->y][z - 0] = flowField[x][y][z];
                }
            }
        }
    }
}

void VoxelGrid::affectFlowfieldAround(Vector3 pos, Vector3 newVal, int kernelSize)
{
    this->affectFlowfieldAround(pos.x, pos.y, pos.z, newVal, kernelSize);
}
void VoxelGrid::affectFlowfieldAround(float x, float y, float z, Vector3 newVal, int kernelSize)
{
//    this->setFlowfield(x, y, z, this->getFlowfield(x, y, z) + newVal);

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
//    this->setFlowfield(x, y, z, this->getFlowfield(x, y, z) + newVal);

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
    this->chunkSize = std::min(this->sizeX, this->chunkSize);

    sizeX -= (sizeX % chunkSize);
    sizeY -= (sizeY % chunkSize);

    this->voxels.clear();
    this->chunks.clear();

    // Create and configure FastNoise object
    this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    this->noise.SetFrequency(blockSize / (float) sizeX);
    this->noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    this->noise.SetFractalLacunarity(2.0);
    this->noise.SetFractalGain(0.7);
    this->noise.SetFractalWeightedStrength(0.5);
    this->noise.SetFractalOctaves(10);
    for (int x = 0; x < this->sizeX; x++)
        for (int y = 0; y < this->sizeY; y++)
            for (int h = 0; h < this->sizeZ; h++)
                noiseMinMax.update(this->noise.GetNoise((float)x, (float)y, (float)h));

    this->chunks = std::vector<std::shared_ptr<VoxelChunk>>(this->numberOfChunksX() * this->numberOfChunksY());

    this->fluidSimulation = FluidSimulation(this->sizeX / this->fluidSimRescale, this->sizeY / this->fluidSimRescale, this->sizeZ / this->fluidSimRescale, 0.2, 0.01, 0.0001, 10);
}

void VoxelGrid::createMesh()
{
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        vc->needRemeshing = true;
//        vc->LoDIndex = vc->LoDs.size() - 1; // To remove
    }
    remeshAll();
}

void VoxelGrid::makeItFall(float erosionStrength)
{
    computeVoxelGroups();
    for(int i = 0; i < 1; i++) {
        for(std::shared_ptr<VoxelChunk>& vc : this->chunks) {
            vc->makeItFall();
        }
    }
    remeshAll();
    if (erosionStrength > 0.0) {
        UnderwaterErosion erod(this->shared_from_this(), 10, erosionStrength, 100);
        erod.Apply();
    }
}
void VoxelGrid::letGravityMakeSandFall(bool remesh)
{
    computeVoxelGroups();
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks) {
        vc->letGravityMakeSandFall();
    }
    if (remesh)
        remeshAll();
}

void VoxelGrid::display() {
    for (std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        vc->display();
    }
}

int VoxelGrid::getHeight(int x, int y) {
    int maxHeight = -1;
    for (Voxel v : this->voxels) {
        if (v.getX() == x && v.getY() == y && v.getZ() > maxHeight)
            maxHeight = v.getZ();
    }
    return maxHeight;
}

bool VoxelGrid::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool VoxelGrid::contains(float x, float y, float z) {
    return (0 <= x && x < this->sizeX && 0 <= y && y < this->sizeY && 0 <= z && z < this->sizeZ);
}
std::shared_ptr<Voxel> VoxelGrid::getVoxel(Vector3 pos) {
    return this->getVoxel(pos.x, pos.y, pos.z);
}

std::shared_ptr<Voxel> VoxelGrid::getVoxel(float x, float y, float z) {
    if(z < 0 || x < 0 || y < 0)
        return nullptr;
    int _x = x, _y = y, _z = z;
    int xChunk = int(_x / chunkSize);
    int voxPosX = _x % chunkSize;
    int yChunk = _y / chunkSize;
    int voxPosY = _y % chunkSize;

    if (xChunk < 0 || yChunk < 0 || _x >= getSizeX() || _y >= getSizeY() || _z >= this->getSizeZ())
        return nullptr;
    return this->chunks[xChunk * (this->sizeY / chunkSize) + yChunk]->voxels[voxPosX][voxPosY][_z];
}

void VoxelGrid::remeshAll()
{
    this->computeVoxelGroups();
    for (std::shared_ptr<VoxelChunk>& vc : this->chunks) {
        vc->createMesh(this->displayWithMarchingCubes);
    }
}


std::vector<std::vector<std::vector<float>>> VoxelGrid::toFloat()
{
    std::vector<std::vector<std::vector<float>>> arr(this->sizeX, std::vector<std::vector<float>>(this->sizeY, std::vector<float>(this->sizeZ, 0.0)));
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        for(int x = 0; x < vc->sizeX; x++) {
            for(int y = 0; y < vc->sizeY; y++) {
                for(int z = 0; z < vc->height; z++) {
                    arr[x + vc->x][y + vc->y][z] = vc->voxels[x][y][z]->getIsosurface();
                }
            }
        }
    }
    return arr;
}
void VoxelGrid::toVoxels(std::vector<std::vector<std::vector<float>>> arr)
{
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        for(int x = 0; x < vc->sizeX; x++) {
            for(int y = 0; y < vc->sizeY; y++) {
                vc->voxelValues[x][y] = arr[vc->x + x][vc->y + y];
            }
        }
        vc->toVoxels();
    }
}
void VoxelGrid::toVoxels()
{
    for(std::shared_ptr<VoxelChunk>& vc : this->chunks)
    {
        vc->toVoxels();
    }
}

void VoxelGrid::computeVoxelGroups()
{
    std::vector<std::vector<std::vector<Vector3*>>> vecs =
            std::vector<std::vector<std::vector<Vector3*>>>(this->sizeX,
            std::vector<std::vector<Vector3*>>(this->sizeY, std::vector<Vector3*>(this->sizeZ)));
    for (int x = 0; x < this->sizeX; x ++) {
        for(int y = 0; y < this->sizeY; y++) {
            for (int z = 0; z < this->sizeZ; z++) {
                if (z > 0)
                    setVoxelGroup(x, y, z, -1);
                vecs[x][y][z] = new Vector3(x, y, z);
            }
        }
    }
    std::unordered_set<Vector3*> groundNeighbors;
    std::unordered_set<Vector3*> done;
    for (int x = 0; x < this->sizeX; x ++) {
        for(int y = 0; y < this->sizeY; y++) {
            int group = getVoxelValue(x, y, 0) > 0.0 ? 0 : -1;
            if (group == 0) {
                setVoxelGroup(x, y, 0, group); // If it's touching the ground, it's directly in first group
                groundNeighbors.insert(vecs[x][y][0]);
                while (groundNeighbors.size() != 0) {
                    Vector3* n = (*groundNeighbors.begin());
                    setVoxelIsOnGround(n, true);
                    done.insert(n);
                    groundNeighbors.erase(groundNeighbors.begin());
                    for (int n_x = -1; n_x <= 1; n_x++)
                        for (int n_y = -1; n_y <= 1; n_y++)
                            for (int n_z = -1; n_z <= 1; n_z++)
                                if (getVoxelGroup(n->x + n_x, n->y + n_y, n->z + n_z) != 0 && getVoxelValue(n->x + n_x, n->y + n_y, n->z + n_z) > 0.0)
                                    groundNeighbors.insert(vecs[n->x + n_x][n->y + n_y][n->z + n_z]);
                }
            }
        }
    }
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
    return std::make_tuple(xChunk * (this->sizeY / chunkSize) + yChunk, voxPosX, voxPosY, _z);
}
float VoxelGrid::getVoxelValue(Vector3 pos) {
    return this->getVoxelValue(pos.x, pos.y, pos.z);
}

float VoxelGrid::getVoxelValue(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelValues[voxPosX][voxPosY][_z];
    return -1;
}
void VoxelGrid::setVoxelValue(Vector3 pos, float newVal) {
    this->setVoxelValue(pos.x, pos.y, pos.z, newVal);
}
void VoxelGrid::setVoxelValue(float x, float y, float z, float newVal)
{
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1) {
        this->chunks[iChunk]->voxelValues[voxPosX][voxPosY][_z] = newVal;
        this->chunks[iChunk]->needRemeshing = true;
    }
}

float VoxelGrid::getOriginalVoxelValue(Vector3 pos) {
    return this->getOriginalVoxelValue(pos.x, pos.y, pos.z);
}

float VoxelGrid::getOriginalVoxelValue(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->originalVoxelValues[voxPosX][voxPosY][_z];
    return -1;
}
Vector3 VoxelGrid::getFlowfield(Vector3 pos) {
    return this->getFlowfield(pos.x, pos.y, pos.z);
}

Vector3 VoxelGrid::getFlowfield(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->flowField[voxPosX][voxPosY][_z];
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
        this->chunks[iChunk]->flowField[voxPosX][voxPosY][_z] = newVal;
}

int VoxelGrid::getVoxelGroup(Vector3 pos) {
    return getVoxelGroup(pos.x, pos.y, pos.z);
}
int VoxelGrid::getVoxelGroup(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelGroups[voxPosX][voxPosY][_z];
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
        this->chunks[iChunk]->voxelGroups[voxPosX][voxPosY][_z] = newVal;
}
bool VoxelGrid::getVoxelIsOnGround(Vector3 pos) {
    return getVoxelIsOnGround(pos.x, pos.y, pos.z);
}
bool VoxelGrid::getVoxelIsOnGround(float x, float y, float z) {
    int iChunk, voxPosX, voxPosY, _z;
    std::tie(iChunk, voxPosX, voxPosY, _z) = this->getChunksAndVoxelIndices(x, y, z);
    if (iChunk != -1)
        return this->chunks[iChunk]->voxelGroups[voxPosX][voxPosY][_z] == 0;
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
        this->chunks[iChunk]->voxelGroups[voxPosX][voxPosY][_z] = (newVal ? 0 : -1);
}




int VoxelGrid::getMaxLoD()
{
    if (chunks.size() > 0 && chunks[0]) {
        int maxLoD = this->sizeX;
        for (std::shared_ptr<VoxelChunk>& vc : this->chunks)
            if (int(vc->LoDs.size() - 1) < maxLoD)
                maxLoD = vc->LoDs.size() - 1;
        return maxLoD;
    } else {
        return 8;
    }
}

void VoxelGrid::saveMap(std::string filename)
{
    std::ofstream out;
    out.open(filename);
    out << this->sizeX << " " << this->sizeY << " " << this->sizeZ << " " << this->chunkSize << "\n";
    for (auto& vc : chunks) {
        for (int x = 0; x < vc->sizeX; x++)
            for (int y = 0; y < vc->sizeY; y++)
                for (int z = 0; z < vc->height; z++)
                    out << vc->voxelValues[x][y][z] << " ";
    }
    out.close();
}
void VoxelGrid::retrieveMap(std::string filename)
{
    std::ifstream in;
    in.open(filename);
    in >> this->sizeX >> this->sizeY >> this->sizeZ >> this->chunkSize;

    std::vector<std::vector<std::vector<std::vector<float>>>> data(this->chunks.size());
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            data[iChunk] = std::vector<std::vector<std::vector<float>>>(this->chunkSize);
            for (int x = 0; x < chunkSize; x++) {
                data[iChunk][x] = std::vector<std::vector<float>>(this->chunkSize);
                for (int y = 0; y < chunkSize; y++) {
                    data[iChunk][x][y] = std::vector<float>(this->getSizeZ());
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float map_val;
                        in >> map_val;
                        data[iChunk][x][y][h] = map_val;
                    }
                }
            }
            iChunk ++;
        }
    }
    this->tempData = data;
}



#include <sstream>
std::string VoxelGrid::toString()
{
    std::ostringstream ret;
    ret << "{";
    ret << "\n\t\"sizeX\": \"" << sizeX << "\",\n\t\"sizeY\": \"" << sizeY << "\",\n\t\"sizeZ\": \"" << sizeZ;
    ret << "\",\n\t\"blockSize\": \"" << blockSize << "\",\n\t\"noise_shifting\": \"" << noise_shifting;
    ret << "\",\n\t\"chunkSize\": \"" << chunkSize << "\",\n\t\"numberOfVoxels\": \"" << sizeX*sizeY*sizeZ;
    ret << "\",\n\t\"numberOfChunks\": \"" << chunks.size() << "\"\n}";
    return ret.str();
}
std::string VoxelGrid::toShortString()
{
    std::ostringstream ret;
    ret << "X" << sizeX << "Y" << sizeY << "Z" << sizeZ;
    ret << "BS" << blockSize << "NS" << noise_shifting;
    ret << "CS" << chunkSize;
    return ret.str();
}
