#ifndef GRAPH_H
#define GRAPH_H

#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"

template <class T>
class Graph
{
public:
    Graph();

    std::vector<std::shared_ptr<GraphNode<T>>> nodes;

    Matrix3<int> nodesIndices;
    Matrix3<int> connectionMatrix;
    Matrix3<int> precedenceMatrix;
    Matrix3<float> adjencyMatrix;
};

template<class T>
Graph<T>::Graph() {

}
#endif // GRAPH_H
