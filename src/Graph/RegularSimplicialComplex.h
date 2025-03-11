#ifndef REGULARSIMPLICIALCOMPLEX_H
#define REGULARSIMPLICIALCOMPLEX_H

#include "Graph/Graph.h"

class RegularSimplicialComplex : public Graph
{
public:
    RegularSimplicialComplex();
    RegularSimplicialComplex(int sizeX, int sizeY);

    GraphNode* getNode(int x, int y);
    GraphNode* getNode(const Vector3& pos);

    void display();

    int maxNodesValues();

    bool checkConnection(const Vector3& posA, const Vector3& posB);
    bool checkConnection(GraphNode* nodeA, GraphNode* nodeB);
    void removeUnavailableLinks();

    int sizeX, sizeY;
};

#endif // REGULARSIMPLICIALCOMPLEX_H
