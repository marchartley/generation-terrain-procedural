#include "VoxelGrid.h"

Voxel::Voxel(int x, int y, int z, TerrainTypes type, float blockSize)
    : x(x), y(y), z(z), type(type), blockSize(blockSize) {

}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0) {

}

void Voxel::display() {
    glPushMatrix();
    glTranslatef(-this->x, -this->y, -this->z);
    glScalef(1/this->blockSize, 1/this->blockSize, 1/this->blockSize);

    // Front
    glBegin(GL_QUADS);
    glVertex3f(0, 0, 0);
    glVertex3f(1, 0, 0);
    glVertex3f(1, 1, 0);
    glVertex3f(0, 1, 0);
    glEnd();

    // Left
    glBegin(GL_QUADS);
    glVertex3f(1, 0, 0);
    glVertex3f(1, 1, 0);
    glVertex3f(1, 1, 1);
    glVertex3f(1, 0, 1);
    glEnd();

    // Back
    glColor3f(0.3, 1.0, 0.2);
    glBegin(GL_QUADS);
    glVertex3f(1, 0, 1);
    glVertex3f(1, 1, 1);
    glVertex3f(0, 1, 1);
    glVertex3f(0, 0, 1);
    glEnd();
    glColor3f(1.0, 1.0, 1.0);

    // Right
    glBegin(GL_QUADS);
    glVertex3f(0, 0, 1);
    glVertex3f(0, 1, 1);
    glVertex3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glEnd();

    // Top
    glBegin(GL_QUADS);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1);
    glVertex3f(1, 0, 1);
    glVertex3f(1, 0, 0);
    glEnd();

    // Bottom
    glBegin(GL_QUADS);
    glVertex3f(0, 1, 0);
    glVertex3f(1, 1, 0);
    glVertex3f(1, 1, 1);
    glVertex3f(0, 1, 1);
    glEnd();
    glPopMatrix();
}


VoxelGrid::VoxelGrid(int nx, int ny, int nz, float blockSize)
    : sizeX(nx), sizeY(ny), sizeZ(nz), blockSize(blockSize) {

}
VoxelGrid::VoxelGrid(Grid& grid) : VoxelGrid(grid.getSizeX(), grid.getSizeY(), 100, .5) {
    this->from2DGrid(grid);
}
VoxelGrid::VoxelGrid() : VoxelGrid(10, 10, 10, 1.0) {

}
void VoxelGrid::from2DGrid(Grid grid) {
    this->voxels.clear();
    float maxi = 1.0;
    float mini = -1.0;
    for (int x = 0; x < grid.getSizeX(); x++) {
        for (int y = 0; y < grid.getSizeY(); y++) {
            this->voxels.push_back(Voxel(x, y, int((grid.getHeight(x, y) + 1)/30 * this->sizeZ/2), TerrainTypes::DIRT, this->blockSize));
            this->voxels.push_back(Voxel(x, y, int((grid.getHeight(x, y) + 1)/30 * this->sizeZ/2)-1, TerrainTypes::DIRT, this->blockSize));
            /*for (int z = 0; z < (grid.getHeight(x, y) + 1) * this->sizeZ/2; z++) {
                this->voxels.push_back(Voxel(x, y, z, TerrainTypes::DIRT));
            }*/
        }
    }
}

void VoxelGrid::display() {
    glPushMatrix();
    glTranslatef(this->getSizeX() * this->blockSize, this->getSizeY() * this->blockSize, this->getSizeZ() * this->blockSize);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(1.0, 1.0, 1.0);
    for (Voxel& vox : this->voxels ) {
        vox.display();
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(0.0, 0.0, 0.0);
    for (Voxel& vox : this->voxels ) {
        vox.display();
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
