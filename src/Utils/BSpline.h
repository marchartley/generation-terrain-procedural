#ifndef BSPLINE_H
#define BSPLINE_H

#include <vector>
#include "DataStructure/Vector3.h"

class BSpline
{
public:
    BSpline();
    BSpline(int numberOfPoints);
    BSpline(std::vector<Vector3> points);

    std::vector<Vector3> getPath(float resolution);
    Vector3 getPoint(float x);
    Vector3 getPoint(float x, Vector3 a, Vector3 b);
    float length();

    Vector3 getCatmullPoint(float x);

    std::vector<Vector3> points;

};

#endif // BSPLINE_H
