#ifndef RADIALSHAPE_H
#define RADIALSHAPE_H

#include "ShapeCurve.h"

class RadialShape : public ShapeCurve
{
public:
    RadialShape();
    RadialShape(float radius, int number_points);

    float estimateDistanceFrom(const Vector3 &pos) const;
    bool containsXY(const Vector3 &pos, bool useNativeShape = false, int increaseAccuracy = false) const;
};

#endif // RADIALSHAPE_H
