#ifndef ADJENCYSOLVER_H
#define ADJENCYSOLVER_H

#include "Graph/Graph.h"
#include "Biomes/BiomeUtils.h"
#include "Utils/Voronoi.h"

class AdjencySolver
{
public:
    AdjencySolver();

    std::vector<int> solve(Voronoi& diagram, std::vector<std::string> restrictedBiomesNames, std::vector<std::pair<std::string, std::string>> impossibleAdgency);

    Graph<std::string> graph;
    Graph<std::string> adjency;
};

#endif // ADJENCYSOLVER_H
