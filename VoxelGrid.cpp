#include "VoxelGrid.h"

VoxelGrid::VoxelGrid(int nx, int ny, int nz, float blockSize)
    : sizeX(nx), sizeY(ny), sizeZ(nz), blockSize(blockSize) {

}
VoxelGrid::VoxelGrid(Grid& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), grid.getMaxHeight(), .5) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(10, 10, 10, 1.0) {

}
void VoxelGrid::from2DGrid(Grid grid) {
    this->voxels.clear();
    this->sizeX = grid.getSizeX();
    this->sizeY = grid.getSizeY();
    this->sizeZ = grid.getMaxHeight();

    for (int x = 0; x < grid.getSizeX(); x++) {
        for (int y = 0; y < grid.getSizeY(); y++) {
            float grid_height = grid.getHeight(x, y) * (this->sizeZ / grid.getMaxHeight());
            int z = int(grid_height)+1;
            std::vector<TerrainTypes> data;
            for (int i = 0; i < z; i++)
                data.push_back(TerrainTypes::DIRT);
            for (int i = z; i < this->getSizeZ(); i++)
                data.push_back(TerrainTypes::AIR);

            this->chunks.push_back(VoxelChunk(x, y, this->getSizeZ(), data));
////            Voxel upperVox(x, y, z, TerrainTypes::DIRT, this->blockSize);
////            Voxel lowerVox(x, y, z - 1, TerrainTypes::DIRT, this->blockSize);
////            this->voxels.push_back(upperVox);
////            this->voxels.push_back(lowerVox);
//            for (int _z = 0; _z < z; _z++) {
//                this->voxels.push_back(Voxel(x, y, _z, TerrainTypes::DIRT, this->blockSize));
//            }
        }
    }
}

void VoxelGrid::display() {
    glPushMatrix();
    glScalef(this->blockSize, this->blockSize, this->blockSize);
    glTranslatef(this->getSizeX()/2, this->getSizeY()/2, -this->getSizeZ()/2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(1.0, 1.0, 1.0);
    for (VoxelChunk& vc : this->chunks) {
        glPushMatrix();
        glTranslatef(-vc.x, -vc.y, 0.0);
        vc.display();
        glPopMatrix();
    }
    glPopMatrix();
//    glPushMatrix();
//    glTranslatef(this->getSizeX() * this->blockSize, this->getSizeY() * this->blockSize, this->getSizeZ() * this->blockSize);
//    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//    glColor3f(1.0, 1.0, 1.0);
//    for (Voxel& vox : this->voxels ) {
//        if (vox.getZ() == this->getHeight(vox.getX(), vox.getY()))
//            vox.display();
//    }
//    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//    glColor3f(0.0, 0.0, 0.0);
//    for (Voxel& vox : this->voxels ) {
//        vox.display();
//    }
//    glPopMatrix();
}

int VoxelGrid::getHeight(int x, int y) {
    int maxHeight = -1;
    for (Voxel v : this->voxels) {
        if (v.getX() == x && v.getY() == y && v.getZ() > maxHeight)
            maxHeight = v.getZ();
    }
    return maxHeight;
}
