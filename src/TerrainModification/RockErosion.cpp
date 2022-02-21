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
}

void RockErosion::Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode, bool applyRemeshing) {
    Matrix3<float> erosionMatrix(grid->sizeX, grid->sizeY, grid->sizeZ);
    grid->applyModification(this->computeErosionMatrix(erosionMatrix, pos, addingMatterMode));
    if (applyRemeshing)
        grid->remeshAll();
}

Matrix3<float>& RockErosion::computeErosionMatrix(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode)
{
    blankMatrix.add(this->attackMask * (addingMatterMode ? -1.f : 1.f), pos - Vector3(size/2.f, size/2.f, size/2.f));
    return blankMatrix;
}
