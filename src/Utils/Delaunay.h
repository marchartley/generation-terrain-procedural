#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "Utils/Voronoi.h"
#include "Graph/Graph.h"

class Delaunay
{
public:
    Delaunay();
    Delaunay(const Voronoi& voronoi);

    Delaunay& fromVoronoi(const Voronoi& voronoi);


    std::vector<Vector3> points;
    Graph graph;
};

#endif // DELAUNAY_H
