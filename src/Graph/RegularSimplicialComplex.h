#ifndef REGULARSIMPLICIALCOMPLEX_H
#define REGULARSIMPLICIALCOMPLEX_H

#include "Graph/Graph.h"

class RegularSimplicialComplex : public Graph<int>
{
public:
    RegularSimplicialComplex();
    RegularSimplicialComplex(int sizeX, int sizeY);

    std::shared_ptr<GraphNode<int>> getNode(int x, int y);
    std::shared_ptr<GraphNode<int>> getNode(Vector3 pos);

    void display();

    int maxNodesValues();

    bool checkConnection(Vector3 posA, Vector3 posB);
    bool checkConnection(std::shared_ptr<GraphNode<int>> nodeA, std::shared_ptr<GraphNode<int>> nodeB);
    void removeUnavailableLinks();

    int sizeX, sizeY;
};

#endif // REGULARSIMPLICIALCOMPLEX_H
