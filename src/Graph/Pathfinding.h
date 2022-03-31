#ifndef PATHFINDING_H
#define PATHFINDING_H

#include <vector>
#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"

class Pathfinding
{
public:
    template<class T>
    static std::pair<float, std::vector<int>> ShortestPathFrom(int source, int dest, std::vector<std::shared_ptr<GraphNode<T>>>& nodes, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<float, std::vector<int>> ShortestPathFrom(int source, int dest, Matrix3<float>& adjencyMap, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFrom(int source, Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> ShortestPathFrom(Matrix3<float>& adjencyMap);

    static std::vector<int> getPath(int dest, std::vector<int> prec);

    template<class T>
    static std::pair<float, std::vector<int>> AStar(std::vector<std::shared_ptr<GraphNode<T>>>& nodes, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<float, std::vector<int>> AStar(Matrix3<float>& adjencyMap, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> Djikstra(Matrix3<float>& adjencyMap, int source);
    static std::pair<std::vector<float>, std::vector<int>> BellmanFord(Matrix3<float>& adjencyMap, int source);
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFasterAlgorithm(Matrix3<float>& adjencyMap, int source, bool useSmallLabelFirst = false, bool useLargeLabelLast = false);
    static bool containsNegativeCycle(const std::vector<float>& distanceVector, Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> FloydWarshall(Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> FloydWarshallImproved(Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> Johnson(Matrix3<float>& adjencyMap);

};

template<class T>
std::pair<float, std::vector<int> > Pathfinding::AStar(std::vector<std::shared_ptr<GraphNode<T>>>& nodes, int source, int dest, std::function<float(int)> heuristicFunction)
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
std::pair<float, std::vector<int>> Pathfinding::ShortestPathFrom(int source, int dest, std::vector<std::shared_ptr<GraphNode<T>>>& nodes, std::function<float(int)> heuristicFunction)
{
    return Pathfinding::AStar(nodes, source, dest, heuristicFunction);
}
#endif // PATHFINDING_H
