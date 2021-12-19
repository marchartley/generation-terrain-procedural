#ifndef ROCKEROSION_H
#define ROCKEROSION_H

#include "Voxel.h"
#include "Matrix3.h"

class RockErosion
{
public:
    RockErosion();
    RockErosion(int size, float maxStrength);

    void Apply(std::shared_ptr<VoxelGrid> grid, Vector3 pos, bool addingMatterMode = false, bool applyRemeshing = true);

    int size;
    float maxStrength;
    Matrix3<float> attackMask;
};

#endif // ROCKEROSION_H
