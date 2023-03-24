#include "TerrainModification/RockErosion.h"

std::map<int, Matrix3<float>> RockErosion::precomputedAttackMasks;
std::map<int, Matrix3<float>> RockErosion::precomputedAttackMasks2D;

RockErosion::RockErosion()
{

}
RockErosion::RockErosion(int size, float maxStrength)
    : size(size), maxStrength(maxStrength)
{
    /*
    attackMask = Matrix3<float>(size, size, size);
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

void RockErosion::Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode) {
    Matrix3<float> erosionMatrix(grid->getDimensions());
    grid->applyModification(this->computeErosionMatrix(erosionMatrix, pos, addingMatterMode));
}

Matrix3<float>& RockErosion::computeErosionMatrix(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode, bool useMax)
{
//    Matrix3<float> realMask = this->attackMask * this->maxStrength;
//    return this->computeErosionMatrix(blankMatrix, attackMask, pos, addingMatterMode, useMax);
    return this->computeErosionMatrix(blankMatrix, this->createPrecomputedAttackMask(this->size) * this->maxStrength, pos, addingMatterMode, useMax);
}

Matrix3<float>& RockErosion::computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float> modifs, Vector3 pos, bool addingMatterMode, bool useMax)
{
    return this->computeErosionMatrix(blankMatrix, modifs, pos, addingMatterMode, Vector3(modifs.sizeX/2.f, modifs.sizeY/2.f, modifs.sizeZ/2.f), useMax);
}
Matrix3<float>& RockErosion::computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float> modifs, Vector3 pos, bool addingMatterMode, Vector3 anchor, bool useMax)
{
    if (useMax) {
        if (addingMatterMode)
            blankMatrix.max(modifs, pos - anchor);
        else
            blankMatrix.min(modifs, pos - anchor);
    } else {
        blankMatrix.add(modifs * (addingMatterMode ? -1.f : 1.f), pos - anchor);
//        blankMatrix.add(modifs * (addingMatterMode ? -1.f : 1.f), pos - anchor, true);
    }
    return blankMatrix;
}

Matrix3<float> &RockErosion::computeErosionMatrix2D(Matrix3<float> &blankMatrix, Vector3 pos, bool addingMatterMode, bool useMax)
{
    return this->computeErosionMatrix(blankMatrix, this->createPrecomputedAttackMask2D(this->size) * this->maxStrength, pos, addingMatterMode, useMax);
}

Matrix3<float>& RockErosion::createPrecomputedAttackMask(int size)
{
    if (RockErosion::precomputedAttackMasks.count(size) == 0) {
        Matrix3<float> attackMask(size, size, size);
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
        RockErosion::precomputedAttackMasks[size] = attackMask;
    }
    return RockErosion::precomputedAttackMasks[size];
}

Matrix3<float>& RockErosion::createPrecomputedAttackMask2D(int size)
{
    if (RockErosion::precomputedAttackMasks2D.count(size) == 0) {
        Matrix3<float> attackMask(size, size);
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
        RockErosion::precomputedAttackMasks2D[size] = attackMask * .5f;
    }
    return RockErosion::precomputedAttackMasks2D[size];
}
