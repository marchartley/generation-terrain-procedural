#ifndef VORONOI_H
#define VORONOI_H


#include "DataStructure/Vector3.h"
#include "Curves/CatmullRomSpline.h"
#include "Curves/Contour.h"
#include <vector>

#include "Graph/Graph.h"

class Voronoi
{
public:
    Voronoi();
    Voronoi(int numberRandomPoints, const Vector3& maxBoundarie);
    Voronoi(int numberRandomPoints, const Vector3& minBoundarie, const Vector3& maxBoundarie);
    Voronoi(int numberRandomPoints, Contour boundingShape);
    Voronoi(std::vector<Vector3> pointset);
    Voronoi(std::vector<Vector3> pointset, const Vector3& maxBoundarie);
    Voronoi(std::vector<Vector3> pointset, const Vector3& minBoundarie, const Vector3& maxBoundarie);
    Voronoi(std::vector<Vector3> pointset, Contour boundingShape);
    std::vector<Contour> solve(bool randomizeUntilAllPointsAreSet = true, int numberOfRelaxations = 10);
    std::vector<Contour> relax(int numberOfRelaxations = 1);

    std::vector<std::vector<Vector3> > computeIntersectionPoints();

    Graph toGraph();


public:
    std::vector<Vector3> pointset;
    std::vector<std::vector<Vector3>> intersectionPoints;
    std::vector<std::vector<int>> neighbors;
    std::vector<Contour> areas;
    Contour boundingShape;
    Vector3 minBoundarie, maxBoundarie;
};

#endif // VORONOI_H
