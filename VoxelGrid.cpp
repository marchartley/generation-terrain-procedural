#include "VoxelGrid.h"
#include "FastNoiseLit.h"

VoxelGrid::VoxelGrid(int nx, int ny, int nz, float blockSize)
    : sizeX(nx), sizeY(ny), sizeZ(nz), blockSize(blockSize) {

    // Create and configure FastNoise object
    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFrequency((blockSize/sizeX));
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(10);

    int chunkSize = (nx > 40 ? 40 : nx);
    this->voxels.clear();
    this->chunks.clear();
    for (int xChunk = 0; xChunk < sizeX / chunkSize; xChunk++) {
        for (int yChunk = 0; yChunk < sizeY / chunkSize; yChunk++) {
            std::vector<std::vector<std::vector<TerrainTypes>>> data;
            for (int x = 0; x < chunkSize; x++) {
                data.push_back(std::vector<std::vector<TerrainTypes>>());
                for (int y = 0; y < chunkSize; y++) {
                    data[x].push_back(std::vector<TerrainTypes>());
                    for(int h = 0; h < this->getSizeZ(); h++) {
                        float noise_val = noise.GetNoise((float)xChunk * chunkSize + x, (float)yChunk * chunkSize + y, (float)h);
                        if (noise_val < 0.0) {
                            data[x][y].push_back(TerrainTypes::DIRT);
                        } else {
                            data[x][y].push_back(TerrainTypes::AIR);
                        }
                    }
                }
            }
            this->chunks.push_back(VoxelChunk(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), data));
        }
    }

    for(int i = 0; i < this->chunks.size() - 1; i++) {
        if (i > this->sizeY/chunkSize - 1) {
            this->chunks[i].neighboring_chunks[LEFT] = &this->chunks[i - int(this->sizeY / chunkSize)];
            this->chunks[i - int(this->sizeY/chunkSize)].neighboring_chunks[RIGHT] = &this->chunks[i];
        }
        if (i % int(this->sizeY/chunkSize) >= 1) {
            this->chunks[i].neighboring_chunks[FRONT] = &this->chunks[i - 1];
            this->chunks[i - 1].neighboring_chunks[BACK] = &this->chunks[i];
        }
    }
    for (VoxelChunk& c: this->chunks)
        c.createMesh();

}
VoxelGrid::VoxelGrid(Grid& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), grid.getTileSize()) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(10, 10, 10, 1.0) {

}
void VoxelGrid::from2DGrid(Grid grid) {
    int chunkSize = (sizeX > 40 ? 40 : sizeX);
    this->voxels.clear();
    this->chunks.clear();
    this->sizeX = grid.getSizeX();
    this->sizeY = grid.getSizeY();
    this->sizeZ = grid.getMaxHeight();
    for (int xChunk = 0; xChunk < sizeX / chunkSize; xChunk++) {
        for (int yChunk = 0; yChunk < sizeY / chunkSize; yChunk++) {
            std::vector<std::vector<std::vector<TerrainTypes>>> data;
            for (int x = 0; x < chunkSize; x++) {
                data.push_back(std::vector<std::vector<TerrainTypes>>());
                for (int y = 0; y < chunkSize; y++) {
                    data[x].push_back(std::vector<TerrainTypes>());
                    float grid_height = grid.getHeight(xChunk * chunkSize + x, yChunk * chunkSize + y) * (this->sizeZ / grid.getMaxHeight());
                    int z = int(grid_height)+1;
                    for (int i = 0; i < z; i++)
                        data[x][y].push_back(TerrainTypes::DIRT);
                    for (int i = z; i < this->getSizeZ(); i++)
                        data[x][y].push_back(TerrainTypes::AIR);
                }
            }
//            data[1][1][0] = TerrainTypes::AIR;
            this->chunks.push_back(VoxelChunk(xChunk * chunkSize, yChunk * chunkSize, chunkSize, chunkSize, this->getSizeZ(), data));
        }
    }
    for(int i = 0; i < this->chunks.size() - 1; i++) {
        if (i > this->sizeY/chunkSize - 1) {
            this->chunks[i].neighboring_chunks[LEFT] = &this->chunks[i - int(this->sizeY / chunkSize)];
            this->chunks[i - int(this->sizeY/chunkSize)].neighboring_chunks[RIGHT] = &this->chunks[i];
        }
        if (i % int(this->sizeY/chunkSize) >= 1) {
            this->chunks[i].neighboring_chunks[FRONT] = &this->chunks[i - 1];
            this->chunks[i - 1].neighboring_chunks[BACK] = &this->chunks[i];
        }
    }
    for (VoxelChunk& c: this->chunks)
        c.createMesh();
}

void VoxelGrid::display(bool apply_marching_cubes, bool display_vertices) {
    glPushMatrix();
    glScalef(1/this->blockSize, 1/this->blockSize, 1/this->blockSize);
    glTranslatef(-this->getSizeX()/2.0, -this->getSizeY()/2.0, -this->getSizeZ()/2.0);
    glColor3f(1.0, 1.0, 1.0);
    for (VoxelChunk& vc : this->chunks) {
        glPushMatrix();
        glTranslatef(vc.x, vc.y, 0.0);
        vc.display(apply_marching_cubes, display_vertices);
        glPopMatrix();
    }
    glPopMatrix();
}

int VoxelGrid::getHeight(int x, int y) {
    int maxHeight = -1;
    for (Voxel v : this->voxels) {
        if (v.getX() == x && v.getY() == y && v.getZ() > maxHeight)
            maxHeight = v.getZ();
    }
    return maxHeight;
}
