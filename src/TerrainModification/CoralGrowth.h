#ifndef CORALGROWTH_H
#define CORALGROWTH_H

#include "DataStructure/Matrix3.h"

class CoralGrowth
{
public:
    CoralGrowth();

    void step();

    Matrix3<float> volume;
    Matrix3<float> coralArea;
    Matrix3<float> highErosionArea;
    Matrix3<float> highDepositArea;

    Vector3 waterflow = Vector3(1, 0, 0);
    Vector3 terrainSize = Vector3(30, 30, 30);
    float waterHeight = 25.f;
};

#endif // CORALGROWTH_H
