#ifndef TOPOMAP_H
#define TOPOMAP_H

#include <boost/graph/graph_concepts.hpp>
#include <vector>
#include "Graph/Graph.h"

class CombinMap;

class Brin {
public:
    Brin();
    Brin(Brin* beta2);
    Brin(Brin* beta2, Brin* beta1);

    Brin* previous();
    std::vector<Brin*> getFace();
    std::vector<Brin*> getOppositeFace();
    std::pair<std::vector<Brin*>, std::vector<Brin*>> getTwoFaces();
    std::vector<std::vector<Brin*>> getAllFaces();
    std::vector<Brin*> getAllBrins();

    std::pair<Brin*, Brin*> subdivide();
    std::vector<Brin*> orbit();
    std::vector<Brin*> inversedOrbit();

    void affectSource(int source);
    int getSource();
    void affectDest(int dest);
    int getDest();
    void affectSourceAndDest(int source, int dest);

    bool isValid();
    int degree();

    std::string toString();

public:
    Brin* beta1 = nullptr;
    Brin* beta2 = nullptr;
    int index = -1;
    int affectedSource = -1;
    int affectedDest = -1;

    static int currentBrinIndex;
};

class CombinMap
{
    friend class Brin;
public:
    CombinMap(Brin* root = nullptr);

    std::vector<Brin*> allBrins();
    // Maybe add the option to directly add a node ID to the new edges
    std::vector<Brin*> addFace(int numberOfEdges = 3, std::vector<Brin*> adjacentBrins = std::vector<Brin*>(), std::vector<int> newIndices = std::vector<int>());
    std::vector<std::vector<Brin*>> allFaces(bool sortBySize = false);
    std::vector<std::vector<Brin*>> allOrbits();
    std::vector<std::vector<Brin*>> allIngoingOrbits();
    std::vector<std::vector<Brin*>> allOutgoingOrbits();
    std::vector<Brin*> exteriorFace(Brin* oneExteriorBrin = nullptr);

    std::vector<std::vector<Brin*>> getHoles(int minHoleSize = 4);
    void triangulate(Brin* oneExteriorBrin = nullptr);

    // At least one of the two to remove...
    std::vector<std::pair<int, int>> getUnorientedEdges();
//    std::vector<std::pair<int, int>> getCorrectionUnorientedEdges();

    std::vector<int> affectNodes();

    void collapseDuplicateNodes();

    CombinMap dual(Brin* oneExteriorBrin = nullptr, Brin* oneEdgeTowardInfinity = nullptr);

    // ToNetworkX()
    // ToGraph()
    // ToGeometricGraph()

    void debug();

    void clean();
    void collapse(Brin* brinToMerge);
    int addNeutralComponents();

    std::string toString();

    Graph<int> toGraph();


public:
    Brin* root = nullptr;
};



#endif // TOPOMAP_H
