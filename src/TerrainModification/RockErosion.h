#ifndef ROCKEROSION_H
#define ROCKEROSION_H

#include "DataStructure/Voxel.h"
#include "DataStructure/Matrix3.h"

class RockErosion
{
public:
    RockErosion();
    RockErosion(int size, float maxStrength);

    void Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode = false, bool applyRemeshing = true);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Vector3 pos, bool addingMatterMode = false);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float>& modifs, Vector3 pos, bool addingMatterMode = false);
    Matrix3<float>& computeErosionMatrix(Matrix3<float>& blankMatrix, Matrix3<float>& modifs, Vector3 pos, bool addingMatterMode,
                                         Vector3 anchor);

    int size;
    float maxStrength;
    Matrix3<float> attackMask;
};

#endif // ROCKEROSION_H
