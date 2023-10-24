#ifndef SPHEROIDALWEATHERING_H
#define SPHEROIDALWEATHERING_H

#include "TerrainGen/VoxelGrid.h"

class SpheroidalWeathering
{
public:
    SpheroidalWeathering(std::shared_ptr<VoxelGrid> voxelGrid = nullptr);

    void applyErosion();

    void _erode(GridF& voxels, const Vector3& pos);
    void _deposition(GridF& voxels);

    std::shared_ptr<VoxelGrid> voxelGrid;

    GridF decimation;
    GridF voxels;
};

#endif // SPHEROIDALWEATHERING_H
