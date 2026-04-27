#ifndef CURVE1D_H
#define CURVE1D_H

#include "Curve.h"

class Curve1D
{
public:
    Curve1D();
    Curve1D(const Curve& curve);
    ~Curve1D();

    float get(float x) const;

protected:
    Curve* curve;
};

#endif // CURVE1D_H
