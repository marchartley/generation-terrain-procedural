#ifndef SHAPECURVE_H
#define SHAPECURVE_H

#include "BSpline.h"

class ShapeCurve : public BSpline
{
public:
    ShapeCurve();
    ShapeCurve(std::vector<Vector3> points);
    ShapeCurve(BSpline path);

    bool inside(Vector3 pos);
    float estimateDistanceFrom(Vector3 pos);
};

#endif // SHAPECURVE_H