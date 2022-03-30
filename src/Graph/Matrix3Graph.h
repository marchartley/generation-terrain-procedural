#ifndef MATRIX3GRAPH_H
#define MATRIX3GRAPH_H

#include "Graph/Graph.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

class Matrix3Graph : public Graph<int>
{
public:
    Matrix3Graph();
    Matrix3Graph(Matrix3<int> matrix);

    Matrix3Graph& initFromBinary(Matrix3<int> matrix);

    Matrix3Graph& computeSurface();
    Matrix3Graph& computeSurface(Matrix3<int> matrix);

    std::vector<Vector3> shortestPath(Vector3 start, Vector3 end);

    Matrix3<int> originalMatrix;
};

#endif // MATRIX3GRAPH_H
