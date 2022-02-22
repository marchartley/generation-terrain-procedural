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
    BSpline(std::vector<BSpline> subsplines);

    std::vector<Vector3> getPath(float resolution);
    Vector3 getPoint(float x);
    Vector3 getPoint(float x, Vector3 a, Vector3 b);
    Vector3 getDerivative(float x);
    float estimateClosestTime(Vector3 pos);
    Vector3 estimateClosestPos(Vector3 pos);
    float length();

    BSpline& close();

    Vector3 getCatmullPoint(float x);

    BSpline simplifyByRamerDouglasPeucker(float epsilon, BSpline subspline = BSpline());

    std::vector<Vector3> points;

};

#endif // BSPLINE_H
