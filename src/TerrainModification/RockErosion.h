#ifndef ROCKEROSION_H
#define ROCKEROSION_H

#include "DataStructure/Matrix3.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/ImplicitPatch.h"

class RockErosion
{
public:
    RockErosion();
    RockErosion(int size, float maxStrength);

    void Apply(std::shared_ptr<VoxelGrid> grid, const Vector3& pos, bool addingMatterMode = false);
    GridF& computeErosionMatrix(GridF& blankMatrix, const Vector3& pos, bool addingMatterMode = false, bool useMax = false);
    GridF& computeErosionMatrix(GridF& blankMatrix, GridF modifs, const Vector3& pos, bool addingMatterMode = false, bool useMax = false);
    GridF& computeErosionMatrix(GridF& blankMatrix, GridF modifs, const Vector3& pos, bool addingMatterMode,
                                         const Vector3& anchor, bool useMax = false);

    GridF& computeErosionMatrix2D(GridF& blankMatrix, const Vector3& pos, bool addingMatterMode = false, bool useMax = false);
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
