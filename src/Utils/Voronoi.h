#ifndef VORONOI_H
#define VORONOI_H


#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"
#include <vector>

class Voronoi
{
public:
    Voronoi();
    Voronoi(int numberRandomPoints, Vector3 minBoundarie, Vector3 maxBoundarie);
    Voronoi(std::vector<Vector3> pointset);
    std::vector<BSpline> solve();


public:
    std::vector<Vector3> pointset;
    std::vector<std::vector<Vector3>> intersectionPoints;
    BSpline boundingShape;
};

#endif // VORONOI_H
