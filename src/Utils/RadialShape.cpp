#include "RadialShape.h"

RadialShape::RadialShape() {}

RadialShape::RadialShape(float radius, int number_points)
    : ShapeCurve({Vector3(0, radius), Vector3(1.f, radius)})
{
    this->resamplePoints(number_points);
}

float RadialShape::estimateDistanceFrom(const Vector3 &pos) const
{
    float t = pos.x;
    return abs(getPoint(t).y - pos.y);
}

bool RadialShape::containsXY(const Vector3 &pos, bool useNativeShape, int increaseAccuracy) const
{
    float t = pos.x;
    return getPoint(t).y > pos.y;
}
