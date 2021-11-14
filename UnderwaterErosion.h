#ifndef UNDERWATEREROSION_H
#define UNDERWATEREROSION_H

#include "VoxelGrid.h"

class UnderwaterErosion
{
public:
    UnderwaterErosion();
    UnderwaterErosion(VoxelGrid* grid, int maxRockSize, float maxRockStrength, int rockAmount);

    std::vector<std::vector<Vector3>> Apply();
    std::vector<std::vector<Vector3>> Apply(Vector3 startingPoint);

    VoxelGrid* grid;
    int maxRockSize, rockAmount;
    float maxRockStrength;
};

#endif // UNDERWATEREROSION_H
