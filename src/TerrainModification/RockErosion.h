#ifndef ROCKEROSION_H
#define ROCKEROSION_H

#include "DataStructure/Voxel.h"
#include "DataStructure/Matrix3.h"

class RockErosion
{
public:
    RockErosion();
    RockErosion(int size, float maxStrength);

    void Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode = false);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float> modifs, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float> modifs, Vector3 pos, bool addingMatterMode,
                                         Vector3 anchor, bool useMax = false);

    Matrix3<float>& computeErosionMatrix2D(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    static Matrix3<float>& createPrecomputedAttackMask(int size);
    static Matrix3<float>& createPrecomputedAttackMask2D(int size);

    int size;
    float maxStrength;
    Matrix3<float> attackMask;

protected:
    static std::map<int, Matrix3<float>> precomputedAttackMasks;
    static std::map<int, Matrix3<float>> precomputedAttackMasks2D;
};

#endif // ROCKEROSION_H
