#ifndef SHAPECURVE_H
#define SHAPECURVE_H

#include "BSpline.h"

class ShapeCurve : public BSpline
{
public:
    ShapeCurve();
    ShapeCurve(std::vector<Vector3> points);
    ShapeCurve(BSpline path);

    bool inside(Vector3 pos, bool useNativeShape = true);
    float estimateDistanceFrom(Vector3 pos);
    float computeArea();

    ShapeCurve intersect(ShapeCurve other);

    Vector3 planeNormal();
    ShapeCurve grow(float increase);
    ShapeCurve shrink(float decrease);

    ShapeCurve& removeDuplicates();

    std::vector<Vector3> randomPointsInside(int numberOfPoints = 1);
};

#endif // SHAPECURVE_H
