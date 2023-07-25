#ifndef SPHEROIDALWEATHERING_H
#define SPHEROIDALWEATHERING_H

#include "TerrainGen/VoxelGrid.h"

class SpheroidalWeathering
{
public:
    SpheroidalWeathering(std::shared_ptr<VoxelGrid> voxelGrid = nullptr);

    void applyErosion();

    void _erode(Matrix3<float>& voxels, const Vector3& pos);
    void _deposition(Matrix3<float>& voxels);

    std::shared_ptr<VoxelGrid> voxelGrid;

    Matrix3<float> decimation;
    Matrix3<float> voxels;
};

#endif // SPHEROIDALWEATHERING_H
