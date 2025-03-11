#ifndef MATRIX3GRAPH_H
#define MATRIX3GRAPH_H

#include "Graph/Graph.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

class Matrix3Graph : public Graph
{
public:
    Matrix3Graph();
    Matrix3Graph(GridI matrix);

    Matrix3Graph& initFromBinary(GridI matrix);

    Matrix3Graph& computeSurface();
    Matrix3Graph& computeSurface(GridI matrix);

    Matrix3Graph& randomizeEdges(float randomFactor = 1.f, bool allowNegatives = false);

    std::vector<Vector3> shortestPath(const Vector3 &start, const Vector3& end);

    GridI originalMatrix;
};

#endif // MATRIX3GRAPH_H
