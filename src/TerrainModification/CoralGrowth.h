#ifndef CORALGROWTH_H
#define CORALGROWTH_H

#include "DataStructure/Matrix3.h"

class CoralGrowth
{
public:
    CoralGrowth();

    void step();

    GridF volume;
    GridF coralArea;
    GridF highErosionArea;
    GridF highDepositArea;

    Vector3 waterflow = Vector3(1, 0, 0);
    Vector3 terrainSize = Vector3(30, 30, 30);
    float waterHeight = 25.f;
};

#endif // CORALGROWTH_H
