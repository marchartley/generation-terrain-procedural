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
    Voronoi(int numberRandomPoints, const Vector3& maxBoundarie);
    Voronoi(int numberRandomPoints, const Vector3& minBoundarie, const Vector3& maxBoundarie);
    Voronoi(int numberRandomPoints, ShapeCurve boundingShape);
    Voronoi(std::vector<Vector3> pointset);
    Voronoi(std::vector<Vector3> pointset, const Vector3& maxBoundarie);
    Voronoi(std::vector<Vector3> pointset, const Vector3& minBoundarie, const Vector3& maxBoundarie);
    Voronoi(std::vector<Vector3> pointset, ShapeCurve boundingShape);
    std::vector<ShapeCurve> solve(bool randomizeUntilAllPointsAreSet = true, int numberOfRelaxations = 10);
    std::vector<ShapeCurve> relax(int numberOfRelaxations = 1);


public:
    std::vector<Vector3> pointset;
    std::vector<std::vector<Vector3>> intersectionPoints;
    std::vector<std::vector<int>> neighbors;
    std::vector<ShapeCurve> areas;
    ShapeCurve boundingShape;
    Vector3 minBoundarie, maxBoundarie;
};

#endif // VORONOI_H
