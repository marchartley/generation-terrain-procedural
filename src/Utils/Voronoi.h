#ifndef VORONOI_H
#define VORONOI_H


#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"
#include "Utils/ShapeCurve.h"
#include <vector>

class Voronoi
{
public:
    Voronoi();
    Voronoi(int numberRandomPoints, Vector3 maxBoundarie);
    Voronoi(int numberRandomPoints, Vector3 minBoundarie, Vector3 maxBoundarie);
    Voronoi(int numberRandomPoints, ShapeCurve boundingShape);
    Voronoi(std::vector<Vector3> pointset);
    Voronoi(std::vector<Vector3> pointset, Vector3 maxBoundarie);
    Voronoi(std::vector<Vector3> pointset, Vector3 minBoundarie, Vector3 maxBoundarie);
    std::vector<BSpline> solve(int numberOfRelaxations = 0);
    std::vector<BSpline> relax(int numberOfRelaxations = 1);


public:
    std::vector<Vector3> pointset;
    std::vector<std::vector<Vector3>> intersectionPoints;
    std::vector<std::vector<int>> neighbors;
    std::vector<BSpline> areas;
    BSpline boundingShape;
    Vector3 minBoundarie, maxBoundarie;
};

#endif // VORONOI_H
