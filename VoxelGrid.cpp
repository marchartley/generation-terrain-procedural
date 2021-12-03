#include "VoxelGrid.h"
#include "UnderwaterErosion.h"

#include "Shader.h"

VoxelGrid::VoxelGrid(int nx, int ny, int nz, float blockSize, float noise_shifting)
    : sizeX(nx), sizeY(ny), sizeZ(nz), blockSize(blockSize), noise_shifting(noise_shifting) {
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
                        data[iChunk][x][y][h] = noise_val;
                    }
                }
            }
            iChunk ++;
        }
    }
    this->fromIsoData(data);
    /*

    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {

            std::vector<std::vector<std::vector<float>>> iso_data(chunkSize);
            for (int x = 0; x < chunkSize; x++) {
                iso_data[x] = std::vector<std::vector<float>>(chunkSize);
                for (int y = 0; y < chunkSize; y++) {
                    iso_data[x][y] = std::vector<float>(this->getSizeZ());
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)h),
                                                            -1.0 + noise_shifting, 1.0 + noise_shifting);
                        iso_data[x][y][h] = noise_val;
                    }
                }
            }
            this->chunks[xChunk + yChunk * this->numberOfChunksX()] =
                    new VoxelChunk(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), iso_data, this);
            if (xChunk == this->numberOfChunksX() - 1) {
                this->chunks[this->chunks.size() - 1]->lastChunkOnX = true;
            }
            if (yChunk == this->numberOfChunksY() - 1) {
                this->chunks[this->chunks.size() - 1]->lastChunkOnY = true;
            }
        }
    }
    for(size_t i = 0; i < this->chunks.size(); i++) {
        if (int(i) > this->sizeY/chunkSize - 1) {
            this->chunks[i]->neighboring_chunks[LEFT] = this->chunks[i - int(this->sizeY / chunkSize)];
            this->chunks[i - int(this->sizeY/chunkSize)]->neighboring_chunks[RIGHT] = this->chunks[i];
        }
        if (i % int(this->sizeY/chunkSize) >= 1) {
            this->chunks[i]->neighboring_chunks[FRONT] = this->chunks[i - 1];
            this->chunks[i - 1]->neighboring_chunks[BACK] = this->chunks[i];
        }
        this->chunks[i]->resetVoxelsNeighbors();
    }
    this->createMesh();
*/
}
VoxelGrid::VoxelGrid(Grid& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), grid.getTileSize()) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(10, 10, 10, 1.0) {

}
VoxelGrid::~VoxelGrid()
{
    for(VoxelChunk* vc : this->chunks)
        delete vc;
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
    this->fromIsoData(data);
/*
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {

            std::vector<std::vector<std::vector<float>>> data(this->chunkSize);
            for (int x = 0; x < chunkSize; x++) {
                data[x] = std::vector<std::vector<float>>(this->chunkSize);
                for (int y = 0; y < chunkSize; y++) {
                    data[x][y] = std::vector<float>(this->getSizeZ());
                    float grid_height = grid.getHeight(xChunk * chunkSize + x, yChunk * chunkSize + y) * (this->sizeZ / grid.getMaxHeight());
                    int z = int(grid_height);
                    for (int i = 0; i < z; i++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)i),
                                                            -2.0, 2.0);
                        data[x][y][i] = abs(noise_val);
                    }
                    for (int i = z; i < this->getSizeZ()+1; i++) {
                        float noise_val = noiseMinMax.remap(this->noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)i),
                                                            -2.0, 2.0);
                        data[x][y][i] = -abs(noise_val);
                    }
                }
            }
            this->chunks[xChunk + yChunk * this->numberOfChunksX()] =
                    new VoxelChunk();
            if (xChunk == this->numberOfChunksX() - 1) {
                this->chunks[this->chunks.size() - 1]->lastChunkOnX = true;
            }
            if (yChunk == this->numberOfChunksY() - 1) {
                this->chunks[this->chunks.size() - 1]->lastChunkOnY = true;
            }
        }
    }
    std::function<void(Voxel*)> resetVoxelsInChunk([](Voxel* v) -> void { v->resetNeighbors(); });
    for(size_t i = 0; i < this->chunks.size(); i++) {
        if (i > this->numberOfChunksY() - 1) {
            this->chunks[i]->neighboring_chunks[LEFT] = this->chunks[i - int(this->sizeY / chunkSize)];
            this->chunks[i - int(this->sizeY/chunkSize)]->neighboring_chunks[RIGHT] = this->chunks[i];
        }
        if (i % this->numberOfChunksY() >= 1) {
            this->chunks[i]->neighboring_chunks[FRONT] = this->chunks[i - 1];
            this->chunks[i - 1]->neighboring_chunks[BACK] = this->chunks[i];
        }
        this->chunks[i]->applyToVoxels(resetVoxelsInChunk);
    }
    this->createMesh();*/
}

VoxelGrid* VoxelGrid::fromIsoData(std::vector<std::vector<std::vector<std::vector<float>>>>& isoData)
{
    int iChunk = 0;
    for (int xChunk = 0; xChunk < this->numberOfChunksX(); xChunk++) {
        for (int yChunk = 0; yChunk < this->numberOfChunksY(); yChunk++) {
            this->chunks[iChunk] =
                    new VoxelChunk(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), isoData[iChunk], this);
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
        this->chunks[i]->resetVoxelsNeighbors();
    }
    this->createMesh();
    return this;
}

void VoxelGrid::initMap()
{
    this->chunkSize = std::min(this->sizeX, this->chunkSize);

    sizeX -= (sizeX % chunkSize);
    sizeY -= (sizeY % chunkSize);

    this->voxels.clear();
    for(auto& chunk : this->chunks)
        delete chunk;
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

    this->chunks = std::vector<VoxelChunk*>(this->numberOfChunksX() * this->numberOfChunksY());
}

void VoxelGrid::createMesh()
{
    for(VoxelChunk* vc : this->chunks)
    {
        vc->needRemeshing = true;
    }
    remeshAll();
}

void VoxelGrid::makeItFall(float erosionStrength)
{
    for(int i = 0; i < 1; i++) {
        for(VoxelChunk* vc : this->chunks) {
            vc->makeItFall();
        }
    }
    remeshAll();
    if (erosionStrength > 0.0) {
        UnderwaterErosion erod(this, 10, erosionStrength, 100);
        erod.Apply();
    }
}
void VoxelGrid::letGravityMakeSandFall(bool remesh)
{
    for(VoxelChunk* vc : this->chunks) {
        vc->letGravityMakeSandFall();
    }
    if (remesh)
        remeshAll();
}

void VoxelGrid::display() {
    for (VoxelChunk* vc : this->chunks)
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
Voxel* VoxelGrid::getVoxel(Vector3 pos) {
    return this->getVoxel(pos.x, pos.y, pos.z);
}

Voxel* VoxelGrid::getVoxel(float x, float y, float z) {
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
    for (VoxelChunk* vc : this->chunks) {
        vc->createMesh(this->displayWithMarchingCubes);
    }
}


std::vector<std::vector<std::vector<float>>> VoxelGrid::toFloat()
{
    std::vector<std::vector<std::vector<float>>> arr(this->sizeX, std::vector<std::vector<float>>(this->sizeY, std::vector<float>(this->sizeZ, 0.0)));
    for(VoxelChunk* vc : this->chunks)
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
    for(VoxelChunk* vc : this->chunks)
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
    for(VoxelChunk* vc : this->chunks)
    {
        vc->toVoxels();
    }
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
