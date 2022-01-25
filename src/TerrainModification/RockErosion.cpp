#include "TerrainModification/RockErosion.h"

RockErosion::RockErosion()
{

}
RockErosion::RockErosion(int size, float maxStrength)
    : size(size), maxStrength(maxStrength)
{
    std::cout.precision(2);
    attackMask = Matrix3<float>(size, size, size);
    float radius = size / 2.f;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                float t_i = (i - radius) / radius;
                float t_j = (j - radius) / radius;
                float t_k = (k - radius) / radius;

                float dist = std::min(std::sqrt(t_i*t_i + t_j*t_j + t_k*t_k), 1.f);
                attackMask.at(i, j, k) = -(1 - dist) * maxStrength;
            }
        }
    }
    /*
    float maxVal = sqrt(3);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                float i_i = (i - (size-1)/2.0) / ((size-1)/2.0);
                float j_i = (j - (size-1)/2.0) / ((size-1)/2.0);
                float k_i = (k - (size-1)/2.0) / ((size-1)/2.0);
                attackMask(i, j, k) = -(maxVal - sqrt(i_i*i_i + j_i*j_i + k_i*k_i))/(maxVal) * maxStrength;
            }
        }
    }*/
}

void RockErosion::Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode, bool applyRemeshing) {
    Matrix3<float> erosionMatrix(grid->sizeX, grid->sizeY, grid->sizeZ);
    grid->applyModification(this->computeErosionMatrix(erosionMatrix, pos, addingMatterMode));
    /*
    for (int x = -size/2; x < size/2; x++) {
        for (int y = -size/2; y < size/2; y++) {
            for (int z = -size/2; z < size/2; z++){
                float currentIsovalue = grid->getVoxelValue(pos + Vector3(x, y, z) + Vector3(.5, .5, .5));
//                if (currentIsovalue > 0.0) {
                    grid->setVoxelValue(pos + Vector3(x, y, z) + Vector3(.5, .5, .5), currentIsovalue + (this->attackMask(x+size/2, y+size/2, z+size/2) * (addingMatterMode ? -1 : 1)));

//                }
*/
                /*
                if(v != nullptr) {
                    v->isosurface += this->attackMask[x+size/2][y+size/2][z+size/2] * (addingMatterMode ? -1 : 1);
                    v->isosurface = std::max(v->isosurface, -2.0f);
                    v->isosurface = std::min(v->isosurface, 2.0f);
                    v->parent->voxelValues[v->x][v->y][v->z] += this->attackMask[x+size/2][y+size/2][z+size/2] * (addingMatterMode ? -1 : 1);
                    v->parent->needRemeshing = true;
                }*/
/*
            }
        }
    }*/
    if (applyRemeshing)
        grid->remeshAll();
}

Matrix3<float>& RockErosion::computeErosionMatrix(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode)
{
    blankMatrix.add(this->attackMask * (addingMatterMode ? -1.f : 1.f), pos - Vector3(size/2.f, size/2.f, size/2.f));
    return blankMatrix;
}
