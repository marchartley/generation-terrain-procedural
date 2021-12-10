#ifndef ROCKEROSION_H
#define ROCKEROSION_H

#include "Voxel.h"

class RockErosion
{
public:
    RockErosion();
    RockErosion(int size, float maxStrength);

    void Apply(std::shared_ptr<Voxel> v, bool addingMatterMode = false, bool applyRemeshing = true);

    int size;
    float maxStrength;
    std::vector<std::vector<std::vector<float>>> attackMask;
};

#endif // ROCKEROSION_H
