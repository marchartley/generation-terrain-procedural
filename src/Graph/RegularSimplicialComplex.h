#ifndef REGULARSIMPLICIALCOMPLEX_H
#define REGULARSIMPLICIALCOMPLEX_H

#include "Graph/Graph.h"

class RegularSimplicialComplex : public Graph
{
public:
    RegularSimplicialComplex();
    RegularSimplicialComplex(int sizeX, int sizeY);

    std::shared_ptr<GraphNodeTemplate<int>> getNode(int x, int y);
    std::shared_ptr<GraphNodeTemplate<int>> getNode(const Vector3& pos);

    void display();

    int maxNodesValues();

    bool checkConnection(const Vector3& posA, const Vector3& posB);
    bool checkConnection(std::shared_ptr<GraphNodeTemplate<int>> nodeA, std::shared_ptr<GraphNodeTemplate<int>> nodeB);
    void removeUnavailableLinks();

    int sizeX, sizeY;
};

#endif // REGULARSIMPLICIALCOMPLEX_H
