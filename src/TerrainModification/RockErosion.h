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
    GridF& computeErosionMatrix(GridF& blankMatrix, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    GridF& computeErosionMatrix(GridF& blankMatrix, GridF modifs, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    GridF& computeErosionMatrix(GridF& blankMatrix, GridF modifs, Vector3 pos, bool addingMatterMode,
                                         Vector3 anchor, bool useMax = false);

    GridF& computeErosionMatrix2D(GridF& blankMatrix, Vector3 pos, bool addingMatterMode = false, bool useMax = false);
    static GridF& createPrecomputedAttackMask(int size);
    static GridF& createPrecomputedAttackMask2D(int size);

    int size;
    float maxStrength;
    GridF attackMask;

protected:
    static std::map<int, GridF> precomputedAttackMasks;
    static std::map<int, GridF> precomputedAttackMasks2D;
};

#endif // ROCKEROSION_H
