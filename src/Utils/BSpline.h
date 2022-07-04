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

    std::vector<Vector3> getPath(int numberOfPoints);
    Vector3 getPoint(float x);
    Vector3 getPoint(float x, Vector3 a, Vector3 b);
    Vector3 getDerivative(float x);
    Vector3 getSecondDerivative(float x);
    float estimateClosestTime(Vector3 pos);
    Vector3 estimateClosestPos(Vector3 pos);
    float estimateDistanceFrom(Vector3 pos);
    float length();

    Vector3 getFrenetDirection(float x);
    Vector3 getFrenetNormal(float x);
    Vector3 getFrenetBinormal(float x);

    Vector3 getCenterCircle(float x);
    Vector3 getDirection(float x);
    Vector3 getNormal(float x);
    Vector3 getBinormal(float x);
    float getCurvature(float x);

    Vector3 center();

    BSpline& close();

    Vector3 getCatmullPoint(float x);

    BSpline simplifyByRamerDouglasPeucker(float epsilon, BSpline subspline = BSpline());

    std::tuple<Vector3, Vector3> AABBox();
    Vector3 containingBoxSize();

    BSpline& grow(float increase);

    BSpline computeConvexHull();

    BSpline& removeDuplicates();

    std::string toString() const;

    friend std::ostream& operator<<(std::ostream& io, const BSpline& s);
    friend std::ostream& operator<<(std::ostream& io, std::shared_ptr<BSpline> s);

    std::vector<Vector3> points;
    bool closed = false;
};

#endif // BSPLINE_H
