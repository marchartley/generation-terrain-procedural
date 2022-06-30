#ifndef GRAPH_H
#define GRAPH_H

#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"
#include "Utils/ConstraintsSolver.h"

template <class T>
class Graph
{
public:
    Graph();

    void forceDrivenPositioning();

    std::vector<std::shared_ptr<GraphNode<T>>> nodes;

    Matrix3<int> nodesIndices;
    Matrix3<int> connectionMatrix;
    Matrix3<int> precedenceMatrix;
    Matrix3<float> adjencyMatrix;
};

template<class T>
Graph<T>::Graph() {

}

template<class T>
void Graph<T>::forceDrivenPositioning()
{
    int maxTriesForAdjencyConstraint = 1000;
    ConstraintsSolver solver;

    for (auto& n : nodes) {
        solver.addItem(new Vector3(n->pos));
    }

    for (int iNode = 0; iNode < nodes.size(); iNode++) {
        for (int jNode = iNode + 1; jNode < nodes.size(); jNode++) {
            if (connectionMatrix.at(iNode, jNode) == 1)
                solver.addDistanceConstraint(iNode, jNode, adjencyMatrix.at(iNode, jNode));
        }
    }
    std::map<int, Vector3> positions;
    for (int iteration = 0; iteration < maxTriesForAdjencyConstraint; iteration++) {
        positions = solver.solve();

        // Now check for crossing adjencies
        // Wait... how... Maybe using the Voronoi algorithm
        break;
    }
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i]->pos = positions[i];
    }
}
#endif // GRAPH_H
