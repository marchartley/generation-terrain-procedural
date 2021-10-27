#include "VoxelChunk.h"

VoxelChunk::VoxelChunk(int x, int y, int sizeX, int sizeY, int height, std::vector<std::vector<std::vector<TerrainTypes>>> data)
    : data(data), height(height), x(x), y(y), sizeX(sizeX), sizeY(sizeY) {
}

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, 0, 0, std::vector<std::vector<std::vector<TerrainTypes>>>())
{

}

void VoxelChunk::draw_part(int start_height, int end_height, TerrainTypes type, bool draw_bottom) {
    if (type == TerrainTypes::AIR)
        return;
    else if (type == TerrainTypes::DIRT)
        glColor3f(1.0, 0.5, 0.5);
    else if (type == TerrainTypes::WATER)
        glColor3f(0.2, 0.6, 1.0);


    // Bottom
    if (draw_bottom) {
        glVertex3f(0.0, 0.0, start_height);
        glVertex3f(0.0, 1.0, start_height);
        glVertex3f(1.0, 1.0, start_height);
        glVertex3f(1.0, 0.0, start_height);
    }

    // Left
    glVertex3f(0.0, 0.0, start_height);
    glVertex3f(0.0, 1.0, start_height);
    glVertex3f(0.0, 1.0, end_height);
    glVertex3f(0.0, 0.0, end_height);

    // Back
    glColor3f(1.0, .0, .0);
    glVertex3f(0.0, 1.0, start_height);
    glVertex3f(1.0, 1.0, start_height);
    glVertex3f(1.0, 1.0, end_height);
    glVertex3f(0.0, 1.0, end_height);
    glColor3f(1.0, 1.0, 1.0);

    // Right
    glVertex3f(1.0, 1.0, start_height);
    glVertex3f(1.0, 0.0, start_height);
    glVertex3f(1.0, 0.0, end_height);
    glVertex3f(1.0, 1.0, end_height);

    // Front
    glVertex3f(1.0, 0.0, start_height);
    glVertex3f(0.0, 0.0, start_height);
    glVertex3f(0.0, 0.0, end_height);
    glVertex3f(1.0, 0.0, end_height);

    // Top
    glVertex3f(0.0, 0.0, end_height);
    glVertex3f(0.0, 1.0, end_height);
    glVertex3f(1.0, 1.0, end_height);
    glVertex3f(1.0, 0.0, end_height);

}
void VoxelChunk::createMesh() {
    this->voxels = std::vector<std::vector<std::vector<Voxel>>>();
    for(int x = 0; x < sizeX; x++) {
        this->voxels.push_back(std::vector<std::vector<Voxel>>());
        for(int y = 0; y < sizeY; y++) {
            this->voxels[x].push_back(std::vector<Voxel>());
            for(int h = 0; h < height; h++) {
//                if (this->data[x][y][h] != TerrainTypes::AIR) {
                    Voxel v(x, y, h, this->data[x][y][h], 1.0);
                    v.has_neighbors[LEFT] = !(x == 0 || this->data[x-1][y][h] == TerrainTypes::AIR);
                    v.has_neighbors[RIGHT] = !(x == sizeX-1 || this->data[x+1][y][h] == TerrainTypes::AIR);
                    v.has_neighbors[FRONT] = !(y == 0 || this->data[x][y-1][h] == TerrainTypes::AIR);
                    v.has_neighbors[BACK] = !(y == sizeY-1 || this->data[x][y+1][h] == TerrainTypes::AIR);
                    v.has_neighbors[BOTTOM] = !(h == 0 || this->data[x][y][h-1] == TerrainTypes::AIR);
                    v.has_neighbors[TOP] = !(h == height - 1 || this->data[x][y][h+1] == TerrainTypes::AIR);

                    if(x == 0 && this->neighboring_chunks.find(LEFT) != this->neighboring_chunks.end())
                    {
                        VoxelChunk* n = this->neighboring_chunks[LEFT];
                        v.has_neighbors[LEFT] = !(n->data[n->sizeX-1][y][h] == TerrainTypes::AIR);
                    }
                    if(y == 0 && this->neighboring_chunks.find(FRONT) != this->neighboring_chunks.end())
                    {
                        VoxelChunk* n = this->neighboring_chunks[FRONT];
                        v.has_neighbors[FRONT] = !(n->data[x][n->sizeY - 1][h] == TerrainTypes::AIR);
                    }
                    if(x == sizeX - 1 && this->neighboring_chunks.find(RIGHT) != this->neighboring_chunks.end())
                    {
                        VoxelChunk* n = this->neighboring_chunks[RIGHT];
                        v.has_neighbors[RIGHT] = !(n->data[0][y][h] == TerrainTypes::AIR);
                    }
                    if(y == sizeY - 1 && this->neighboring_chunks.find(BACK) != this->neighboring_chunks.end())
                    {
                        VoxelChunk* n = this->neighboring_chunks[BACK];
                        v.has_neighbors[BACK] = !(n->data[x][0][h] == TerrainTypes::AIR);
                    }
                    this->voxels[x][y].push_back(v); // v.display();
//                }
            }
        }
    }





    for(int x = 0; x < sizeX; x++) {
        for(int y = 0; y < sizeY; y++) {
            for(int h = 0; h < height; h++) {
                Voxel& v = this->voxels[x][y][h];
                // Compute isosurfaces
                if (x == 0){
                    v.vertices[0].isosurface = 1.0;
                    v.vertices[3].isosurface = 1.0;
                    v.vertices[7].isosurface = 1.0;
                    v.vertices[4].isosurface = 1.0;
                } else if (x == sizeX - 1) {
                    v.vertices[1].isosurface = 1.0;
                    v.vertices[2].isosurface = 1.0;
                    v.vertices[5].isosurface = 1.0;
                    v.vertices[6].isosurface = 1.0;
                }
                if (y == 0){
                    v.vertices[2].isosurface = 1.0;
                    v.vertices[3].isosurface = 1.0;
                    v.vertices[6].isosurface = 1.0;
                    v.vertices[7].isosurface = 1.0;
                } else if (y == sizeY - 1) {
                    v.vertices[0].isosurface = 1.0;
                    v.vertices[1].isosurface = 1.0;
                    v.vertices[4].isosurface = 1.0;
                    v.vertices[5].isosurface = 1.0;
                }
                if (h == 0){
                    v.vertices[0].isosurface = 1.0;
                    v.vertices[1].isosurface = 1.0;
                    v.vertices[2].isosurface = 1.0;
                    v.vertices[3].isosurface = 1.0;
                } else if (h == height - 1) {
                    v.vertices[4].isosurface = 1.0;
                    v.vertices[5].isosurface = 1.0;
                    v.vertices[6].isosurface = 1.0;
                    v.vertices[7].isosurface = 1.0;
                }
                if (x > 0) {
                    if (y > 0) {
                        if (h > 0) {
                            if (!this->voxels[x - 1][y - 1][h - 1])
                                v.vertices[3].isosurface = 1.0;
                        }
                        if (h < height - 1) {
                            if (!this->voxels[x - 1][y - 1][h + 1])
                                v.vertices[7].isosurface = 1.0;
                        }
                    }
                    if (y < sizeY - 1) {
                        if (h > 0) {
                            if (!this->voxels[x - 1][y + 1][h - 1])
                                v.vertices[0].isosurface = 1.0;
                        }
                        if (h < height - 1) {
                            if (!this->voxels[x - 1][y + 1][h + 1])
                                v.vertices[4].isosurface = 1.0;
                        }
                    }
                }
                if (x < sizeX - 1) {
                    if (y > 0) {
                        if (h > 0) {
                            if (!this->voxels[x + 1][y - 1][h - 1])
                                v.vertices[2].isosurface = 1.0;
                        }
                        if (h < height - 1) {
                            if (!this->voxels[x + 1][y - 1][h + 1])
                                v.vertices[6].isosurface = 1.0;
                        }
                    }
                    if (y < sizeY - 1) {
                        if (h > 0) {
                            if (!this->voxels[x + 1][y + 1][h - 1])
                                v.vertices[1].isosurface = 1.0;
                        }
                        if (h < height - 1) {
                            if (!this->voxels[x + 1][y + 1][h + 1])
                                v.vertices[5].isosurface = 1.0;
                        }
                    }
                }
            }
        }
    }
}
void VoxelChunk::display()
{
    for(int x = 0; x < sizeX; x++) {
        for(int y = 0; y < sizeY; y++) {
            for(int h = 0; h < height; h++) {
                this->voxels[x][y][h].display();
//                if (this->data[x][y][h] != TerrainTypes::AIR) {
//                    Voxel v(x, y, h, this->data[x][y][h], 1.0);
//                    v.has_neighbors[LEFT] = !(x == 0 || this->data[x-1][y][h] == TerrainTypes::AIR);
//                    v.has_neighbors[RIGHT] = !(x == sizeX-1 || this->data[x+1][y][h] == TerrainTypes::AIR);
//                    v.has_neighbors[FRONT] = !(y == 0 || this->data[x][y-1][h] == TerrainTypes::AIR);
//                    v.has_neighbors[BACK] = !(y == sizeY-1 || this->data[x][y+1][h] == TerrainTypes::AIR);
//                    v.has_neighbors[BOTTOM] = !(h == 0 || this->data[x][y][h-1] == TerrainTypes::AIR);
//                    v.has_neighbors[TOP] = !(h == height || this->data[x][y][h+1] == TerrainTypes::AIR);

//                    v.display();
//                }
            }
        }
    }
    /*
    glBegin(GL_QUADS);
    int start_height = -1;
    int end_height = 0;
    bool draw = false;
    TerrainTypes prev_type = this->data[0];
    for(int h = 0; h < this->height; h++){
        if (this->data[h] != prev_type && start_height != -1) {
            if (!draw)
                continue;
            end_height = h - 1;
            this->draw_part(start_height, end_height, prev_type, (start_height == 0));
            prev_type = this->data[h];
            start_height = -1;
            draw = false;
        } else {
            draw = true;
            if (start_height == -1)
                start_height = h;

        }
    }
    if (this->data[this->data.size() -1] == prev_type){
        this->draw_part(start_height, end_height, prev_type, (start_height <= 0));


    }
    glEnd();
//    for (int h = 0; h < this->height; h++)
//    {
//        Voxel v(0, 0, h, this->data[h], 1.0);
//        v.has_neighbors[VOXEL_NEIGHBOR::BOTTOM] = (h != this->height && this->data[h+1] != TerrainTypes::AIR);
//        v.has_neighbors[VOXEL_NEIGHBOR::TOP] = (h > 0);
//        v.has_neighbors[VOXEL_NEIGHBOR::LEFT] = (this->x > 0);
//        v.has_neighbors[VOXEL_NEIGHBOR::RIGHT] = (this->x < 99);
//        v.has_neighbors[VOXEL_NEIGHBOR::FRONT] = (this-> y > 0);
//        v.has_neighbors[VOXEL_NEIGHBOR::BACK] = (this-> y < 99);
//        v.display();
//    }
*/
}
