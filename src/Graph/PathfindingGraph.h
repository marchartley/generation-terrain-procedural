#ifndef PATHFINDINGGRAPH_H
#define PATHFINDINGGRAPH_H


#include <vector>
#include "DataStructure/Matrix3.h"
#include "Graph/Graph.h"

class Pathfinding
{
public:
    template<class T>
    static std::pair<float, std::vector<int>> ShortestPathFrom(int source, int dest, std::vector<std::shared_ptr<GraphNodeTemplate<T>>>& nodes, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<float, std::vector<int>> ShortestPathFrom(int source, int dest, GridF& adjencyMap, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFrom(int source, GridF& adjencyMap);
    static std::pair<GridF, GridI> ShortestPathFrom(GridF& adjencyMap);

    static std::vector<int> getPath(int dest, std::vector<int> prec);

    template<class T>
    static std::pair<float, std::vector<int>> AStar(std::map<T, std::shared_ptr<GraphNodeTemplate<T>>>& nodes, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
//    static std::pair<float, std::vector<int>> AStar(GridF& adjencyMap, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<float, std::vector<int>> AStar(Graph& graph, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> Djikstra(Graph& graph, int source);
    static std::pair<std::vector<float>, std::vector<int>> BellmanFord(GridF& adjencyMap, int source);
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFasterAlgorithm(GridF& adjencyMap, int source, bool useSmallLabelFirst = false, bool useLargeLabelLast = false);
    static bool containsNegativeCycle(const std::vector<float>& distanceVector, GridF& adjencyMap);
    static std::pair<GridF, GridI> FloydWarshall(GridF& adjencyMap);
    static std::pair<GridF, GridI> FloydWarshallImproved(GridF& adjencyMap);
    static std::pair<GridF, GridI> Johnson(GridF& adjencyMap);

};

template<class T>
std::pair<float, std::vector<int> > Pathfinding::AStar(std::map<T, std::shared_ptr<GraphNodeTemplate<T>>>& nodes, int source, int dest, std::function<float(int)> heuristicFunction)
{
    int n = nodes.size();
    std::vector<int> prec(n, -1);
    std::vector<float> gScore(n, std::numeric_limits<float>::max());
    std::vector<float> fScore(n, std::numeric_limits<float>::max());
    std::vector<int> openSet = {source};
    std::vector<int> closedSet;

    gScore[source] = 0;
    fScore[source] = heuristicFunction(source);

    while (!openSet.empty())
    {
        int current = -1;
        float minFScore = std::numeric_limits<float>::max();
        for (const int& q : openSet) {
            if (fScore[q] < minFScore) {
                current = q;
                minFScore = fScore[q];
            }
        }
        if (current == dest) {
            float totalDist = 0.0;
            int next = current;
            while (current != source) {
                current = prec[current];
                totalDist += nodes[current]->getNeighborDistanceByIndex(next);
                next = current;
            }
            return std::make_pair(totalDist, prec);
        }
        // Remove this node from the open set
        openSet.erase(std::find(openSet.begin(), openSet.end(), current));

        for (const auto& n : nodes[current]->neighbors) {
            int u = std::get<0>(n)->index;
            float possibleGScore = gScore[current] + std::get<1>(n) + heuristicFunction(u); //nodes[current]->getNeighborDistanceByIndex(u);
            if (possibleGScore < gScore[u]) {
                prec[u] = current;
                gScore[u] = possibleGScore;
                fScore[u] = possibleGScore + heuristicFunction(u); // + nodes[current]->getNeighborDistanceByIndex(u);
                if (std::find(openSet.begin(), openSet.end(), u) == openSet.end())
                    openSet.push_back(u);
            }
        }
    }
    return std::make_pair(std::numeric_limits<float>::max(), prec);
}

template<class T>
std::pair<float, std::vector<int>> Pathfinding::ShortestPathFrom(int source, int dest, std::vector<std::shared_ptr<GraphNodeTemplate<T>>>& nodes, std::function<float(int)> heuristicFunction)
{
    return Pathfinding::AStar(nodes, source, dest, heuristicFunction);
}

#endif // PATHFINDINGGRAPH_H
