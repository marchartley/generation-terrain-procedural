#ifndef ADJENCYSOLVER_H
#define ADJENCYSOLVER_H

#include "Graph/Graph.h"
#include "Biomes/BiomeUtils.h"
#include "Utils/Voronoi.h"

class AdjencySolver
{
public:
    AdjencySolver();

    std::vector<int> solveGeomToTopo(Voronoi& diagram, std::vector<std::string> restrictedBiomesNames, std::vector<std::pair<std::string, std::string>> impossibleAdjency);

    void solveTopoToGeom(std::vector<BiomeInstance> instancesToPlace, std::vector<std::pair<std::string, std::string>> impossibleAdjency, ShapeCurve availableSpace);

    GraphTemplate<std::string> graph;
    GraphTemplate<std::string> adjency;
};

#endif // ADJENCYSOLVER_H
