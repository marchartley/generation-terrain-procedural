#include "TerrainModification/RockErosion.h"

std::map<int, GridF> RockErosion::precomputedAttackMasks;
std::map<int, GridF> RockErosion::precomputedAttackMasks2D;

RockErosion::RockErosion()
{

}
RockErosion::RockErosion(int size, float maxStrength)
    : size(size), maxStrength(maxStrength)
{
    /*
    attackMask = GridF(size, size, size);
    float radius = size / 2.f;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                float t_i = (i - radius) / radius;
                float t_j = (j - radius) / radius;
                float t_k = (k - radius) / radius;

                float dist = std::min(std::sqrt(t_i*t_i + t_j*t_j + t_k*t_k), 1.f);
                attackMask.at(i, j, k) = -(1 - dist) * this->maxStrength;
            }
        }
    }*/
//    this->attackMask = this->createPrecomputedAttackMask(this->size) * this->maxStrength;
}

void RockErosion::Apply(std::shared_ptr<VoxelGrid> grid, const Vector3& pos, bool addingMatterMode) {
    GridF erosionMatrix(grid->getDimensions());
    grid->applyModification(this->computeErosionMatrix(erosionMatrix, pos, addingMatterMode));
}

GridF& RockErosion::computeErosionMatrix(GridF& blankMatrix, const Vector3& pos, bool addingMatterMode, bool useMax)
{
//    GridF realMask = this->attackMask * this->maxStrength;
//    return this->computeErosionMatrix(blankMatrix, attackMask, pos, addingMatterMode, useMax);
    return this->computeErosionMatrix(blankMatrix, this->createPrecomputedAttackMask(this->size) * this->maxStrength, pos, addingMatterMode, useMax);
}

GridF& RockErosion::computeErosionMatrix(GridF& blankMatrix, GridF modifs, const Vector3& pos, bool addingMatterMode, bool useMax)
{
    return this->computeErosionMatrix(blankMatrix, modifs, pos, addingMatterMode, Vector3(modifs.sizeX/2.f, modifs.sizeY/2.f, modifs.sizeZ/2.f), useMax);
}
GridF& RockErosion::computeErosionMatrix(GridF& blankMatrix, GridF modifs, const Vector3& pos, bool addingMatterMode, const Vector3& anchor, bool useMax)
{
    if (useMax) {
        if (addingMatterMode)
            blankMatrix.max(modifs, pos - anchor);
        else
            blankMatrix.min(modifs, pos - anchor);
    } else {
        auto mat = modifs * (addingMatterMode ? -1.f : 1.f);
        blankMatrix.add(mat, pos - anchor);
//        blankMatrix.add(modifs * (addingMatterMode ? -1.f : 1.f), pos - anchor, true);
    }
    return blankMatrix;
}

GridF &RockErosion::computeErosionMatrix2D(GridF &blankMatrix, const Vector3& pos, bool addingMatterMode, bool useMax)
{
    auto mask = this->createPrecomputedAttackMask2D(this->size);
    mask *= this->maxStrength;
    return this->computeErosionMatrix(blankMatrix, mask, pos.xy(), addingMatterMode, useMax);
}

GridF& RockErosion::createPrecomputedAttackMask(int size)
{
    if (RockErosion::precomputedAttackMasks.count(size) == 0) {
        GridF attackMask(size, size, size);
        float radius = size / 2.f;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    float t_i = (i - radius) / radius;
                    float t_j = (j - radius) / radius;
                    float t_k = (k - radius) / radius;

                    float dist = std::min(std::sqrt(t_i*t_i + t_j*t_j + t_k*t_k), 1.f);
                    attackMask.at(i, j, k) = -(1 - dist);
                }
            }
        }
        attackMask.raiseErrorOnBadCoord = false,
        RockErosion::precomputedAttackMasks[size] = -attackMask / attackMask.sum();
    }
    return RockErosion::precomputedAttackMasks[size];
}

GridF& RockErosion::createPrecomputedAttackMask2D(int size)
{
    if (RockErosion::precomputedAttackMasks2D.count(size) == 0) {
        GridF attackMask(size, size);
        float radius = size / 2.f;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    float t_i = (i - radius) / radius;
                    float t_j = (j - radius) / radius;
                    float t_k = (k - radius) / radius;

                    float dist = std::min(std::sqrt(t_i*t_i + t_j*t_j + t_k*t_k), 1.f);
                    attackMask.at(i, j) += -(1 - dist);
                }
            }
        }
        attackMask.raiseErrorOnBadCoord = false,
        RockErosion::precomputedAttackMasks2D[size] = -(attackMask * .5f) / attackMask.sum();
    }
    return RockErosion::precomputedAttackMasks2D[size];
}
