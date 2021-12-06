#ifndef UNDERWATEREROSION_H
#define UNDERWATEREROSION_H

#include "VoxelGrid.h"

class UnderwaterErosion
{
public:
    UnderwaterErosion();
    UnderwaterErosion(VoxelGrid* grid, int maxRockSize, float maxRockStrength, int rockAmount);

//    std::vector<std::vector<Vector3>> Apply(int avoidMatter = -1);
    std::vector<std::vector<Vector3>> Apply(Vector3* startingPoint = nullptr, Vector3* originalDirection = nullptr, int avoidMatter = -1);
    std::vector<Vector3> CreateTunnel(Vector3* startingPoint = nullptr, Vector3* endingPoint = nullptr, int numberPoints = 2, bool addingMatter = false);


    VoxelGrid* grid;
    int maxRockSize, rockAmount;
    float maxRockStrength;
};

#endif // UNDERWATEREROSION_H
