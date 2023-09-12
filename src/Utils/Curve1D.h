#ifndef CURVE1D_H
#define CURVE1D_H

#include "BSpline.h"

class Curve1D : public BSpline
{
public:
    Curve1D();
//    Curve1D(float mini, float maxi, int numberOfPoints);
    Curve1D(const std::vector<Vector3>& points);

    float get(float x) const;
};

#endif // CURVE1D_H
