#ifndef PATHFINDING_H
#define PATHFINDING_H

#include <vector>
#include "DataStructure/Matrix3.h"

class Pathfinding
{
public:
    static std::pair<float, std::vector<int>> ShortestPathFrom(int source, int dest, Matrix3<float>& adjencyMap, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFrom(int source, Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> ShortestPathFrom(Matrix3<float>& adjencyMap);

    static std::vector<int> getPath(int dest, std::vector<int> prec);

    static std::pair<float, std::vector<int>> AStar(Matrix3<float>& adjencyMap, int source, int dest, std::function<float(int)> heuristicFunction = [](int)-> float {return 0.f;});
    static std::pair<std::vector<float>, std::vector<int>> Djikstra(Matrix3<float>& adjencyMap, int source);
    static std::pair<std::vector<float>, std::vector<int>> BellmanFord(Matrix3<float>& adjencyMap, int source);
    static std::pair<std::vector<float>, std::vector<int>> ShortestPathFasterAlgorithm(Matrix3<float>& adjencyMap, int source, bool useSmallLabelFirst = false, bool useLargeLabelLast = false);
    static bool containsNegativeCycle(const std::vector<float>& distanceVector, Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> FloydWarshall(Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> FloydWarshallImproved(Matrix3<float>& adjencyMap);
    static std::pair<Matrix3<float>, Matrix3<int>> Johnson(Matrix3<float>& adjencyMap);

};

#endif // PATHFINDING_H
