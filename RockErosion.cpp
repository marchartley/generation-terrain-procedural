#include "RockErosion.h"

RockErosion::RockErosion()
{

}
RockErosion::RockErosion(int size, float maxStrength)
    : size(size), maxStrength(maxStrength)
{
    std::cout.precision(2);
    attackMask = new float**[size];
    float maxVal = sqrt(3);
    for (int i = 0; i < size; i++) {
        attackMask[i] = new float*[size];
        for (int j = 0; j < size; j++) {
            attackMask[i][j] = new float[size];
            for (int k = 0; k < size; k++) {
                float i_i = (i - (size-1)/2.0) / ((size-1)/2.0);
                float j_i = (j - (size-1)/2.0) / ((size-1)/2.0);
                float k_i = (k - (size-1)/2.0) / ((size-1)/2.0);
                attackMask[i][j][k] = -(maxVal - sqrt(i_i*i_i + j_i*j_i + k_i*k_i))/(maxVal) * maxStrength;
//                std::cout << attackMask[i][j][k] << "\t";
            }
//            std::cout << std::endl;
        }
//        std::cout << std::endl;
    }
}

void RockErosion::Apply(Voxel* main_v, bool addingMatterMode, bool applyRemeshing) {
    for (int x = -size/2; x < size/2; x++) {
        for (int y = -size/2; y < size/2; y++) {
            for (int z = -size/2; z < size/2; z++){
                int v_x = main_v->x + x, v_y = main_v->y + y, v_z = main_v->z + z;
                VoxelChunk* current_parent = main_v->parent;
                if ((main_v->x + x < 0 && current_parent->x == 0) || (main_v->x + x >= current_parent->sizeX && current_parent->lastChunkOnX))
                    continue;
                if ((main_v->y + y < 0 && current_parent->y == 0) || (main_v->y + y >= current_parent->sizeY && current_parent->lastChunkOnY))
                    continue;
                if (main_v->z + z < 0 || main_v->z + z >= current_parent->height)
                    continue;

                if (main_v->x + x < 0) {
                    current_parent = current_parent->neighboring_chunks[LEFT];
                    v_x += current_parent->sizeX;
                }
                if (main_v->y + y < 0) {
                    current_parent = current_parent->neighboring_chunks[FRONT];
                    v_y += current_parent->sizeY;
                }
                if (main_v->x + x >= current_parent->sizeX) {
                    current_parent = current_parent->neighboring_chunks[RIGHT];
                    v_x -= current_parent->sizeX;
                }
                if (main_v->y + y >= current_parent->sizeY) {
                    current_parent = current_parent->neighboring_chunks[BACK];
                    v_y -= current_parent->sizeY;
                }
                Voxel* v = current_parent->voxels[v_x][v_y][v_z];
                v->manual_isosurface += this->attackMask[x+size/2][y+size/2][z+size/2] * (addingMatterMode ? -1 : 1);
                v->manual_isosurface = std::max(v->manual_isosurface, -2.0f);
                v->manual_isosurface = std::min(v->manual_isosurface, 2.0f);
                if (v->getIsosurface() < 0.0)
                    v->type = TerrainTypes::AIR;
                else
                    v->type = TerrainTypes::DIRT;
                v->parent->data[v->x][v->y][v->z] = v->type;
                v->parent->needRemeshing = true;
            }
        }
    }
    if (applyRemeshing)
        main_v->parent->parent->remeshAll();
}
