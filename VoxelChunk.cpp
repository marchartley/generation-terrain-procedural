#include "VoxelChunk.h"

VoxelChunk::VoxelChunk(int x, int y, int height, std::vector<TerrainTypes> data, int* surrounding_heights)
    : data(data), height(height), x(x), y(y) {
    for (int i = 0; i < 4; i++)
        this->surrounding_heights[i] = surrounding_heights[i];
}

VoxelChunk::VoxelChunk() : VoxelChunk(0, 0, 0, std::vector<TerrainTypes>(), NULL)
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
void VoxelChunk::display()
{
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
}
