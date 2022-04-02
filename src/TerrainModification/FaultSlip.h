#ifndef FAULTSLIP_H
#define FAULTSLIP_H

#include "DataStructure/Vector3.h"
#include "TerrainGen/VoxelGrid.h"

class FaultSlip
{
public:
    FaultSlip();

    void Apply(std::shared_ptr<VoxelGrid> grid, bool applyRemeshing = true);

    Vector3 firstPointFault;
    Vector3 secondPointFault;
    Vector3 slippingDirection;
    float slippingDistance;

    bool positiveSideFalling;
};

#endif // FAULTSLIP_H
