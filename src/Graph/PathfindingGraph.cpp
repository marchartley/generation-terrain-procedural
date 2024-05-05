#include "PathfindingGraph.h"

std::pair<float, std::vector<int> > Pathfinding::ShortestPathFrom(int source, int dest, Graph &graph, std::function<float (int)> heuristicFunction)
{
    return AStar(graph, source, dest, heuristicFunction);
}
std::pair<std::vector<float>, std::vector<int> > Pathfinding::ShortestPathFrom(int source, Graph &graph)
{
    return ShortestPathFasterAlgorithm(graph, source);
}
std::pair<GridF, GridI > Pathfinding::ShortestPathFrom(Graph& graph)
{
    return FloydWarshallImproved(graph);
}

std::vector<int> Pathfinding::getPath(int dest, std::vector<int> prec)
{
    std::vector<int> path;
    int current = dest;
    while (current != -1) {
        path.insert(path.begin(), current);
        current = prec[current];
    }
    return path;
}

std::vector<int> Pathfinding::getPath(int source, int dest, GridI prec)
{
    std::vector<int> path;
    int current = source;
    while (current != -1) {
        path.push_back(current);
        current = prec(dest, current);
    }
    return path;
}

std::pair<float, std::vector<int> > Pathfinding::AStar(Graph &graph, int source, int dest, std::function<float(int)> heuristicFunction)
{
//    return Pathfinding::AStar(graph, source, dest, heuristicFunction);
//    /*
    int n = graph.size();
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
                totalDist += graph.getDistance(graph[current], graph[next]); // graph[current]->getNeighborDistanceByIndex(next);
                next = current;
            }
            return std::make_pair(totalDist, prec);
        }
        // Remove this node from the open set
        openSet.erase(std::find(openSet.begin(), openSet.end(), current));

        for (int u = 0; u < n; u++) {
            float distance = graph.getDistance(graph[current], graph[u]) /*graph[current]->getNeighborDistanceByIndex(u)*/;
            float possibleGScore = gScore[current] + distance + heuristicFunction(u);
            if (possibleGScore < gScore[u]) {
                prec[u] = current;
                gScore[u] = possibleGScore;
                fScore[u] = possibleGScore + graph.getDistance(graph[current], graph[u]) /*graph[current]->getNeighborDistanceByIndex(u)*/ + heuristicFunction(u);
                if (std::find(openSet.begin(), openSet.end(), u) == openSet.end())
                    openSet.push_back(u);
            }
        }
    }
    return std::make_pair(std::numeric_limits<float>::max(), prec);
//    */
}

std::pair<std::vector<float>, std::vector<int> > Pathfinding::Djikstra(Graph& graph, int source)
{
    int n = graph.size();
    std::vector<float> distances(n, std::numeric_limits<float>::max());
    std::vector<int> prec(n, -1);
    distances[source] = 0.0;

    std::vector<int> Q; // Stores all the nodes
    for (int i = 0; i < n; i++)
        Q.push_back(i);
    while (!Q.empty()) {
        // Find the closest node in the active list
        float minDist = std::numeric_limits<float>::max();
        int minNodeNumber = -1;
        int indexToRemoveFromQ = -1;
        for (size_t i = 0; i < Q.size(); i++) {
//            if (Q[i] == source) continue;
            if (distances[Q[i]] < minDist) {
                minDist = distances[Q[i]];
                minNodeNumber = Q[i];
                indexToRemoveFromQ = i;
            }
        }
        if (indexToRemoveFromQ == -1) { // We have seen all the nodes accessible from the source, meaning the graph is not connex, stop here
            return std::make_pair(distances, prec);
        }
        Q.erase(Q.begin() + indexToRemoveFromQ);
        for (int i = 0; i < n; i++) {
            if (distances[i] > distances[minNodeNumber] + graph.getDistance(graph[minNodeNumber], graph[i]) /*graph[minNodeNumber]->getNeighborDistanceByIndex(i)*/) {
                distances[i] = distances[minNodeNumber] + graph.getDistance(graph[minNodeNumber], graph[i]); //graph[minNodeNumber]->getNeighborDistanceByIndex(i);
                prec[i] = minNodeNumber;
            }
        }
    }
    return std::make_pair(distances, prec);
}

std::pair<std::vector<float>, std::vector<int> > Pathfinding::BellmanFord(Graph& graph, int source)
{
    int n = graph.size();
    std::vector<float> distances(n, std::numeric_limits<float>::max());
    std::vector<int> prec(n, -1);
    distances[source] = 0.0;
    bool changeMade = false;
    for (int k = 0; k < n; k++) {
        changeMade = false;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (distances[j] > distances[i] + graph.getDistance(graph[i], graph[j]) /*graph[i]->getNeighborDistanceByIndex(j)*/) {
                    distances[j] = distances[i] + graph.getDistance(graph[i], graph[j]); //graph[i]->getNeighborDistanceByIndex(j);
                    prec[j] = i;
                    changeMade = true;
                }
            }
        }
        if (!changeMade)
            break;
    }
    return std::make_pair(distances, prec);
}

std::pair<std::vector<float>, std::vector<int> > Pathfinding::ShortestPathFasterAlgorithm(Graph& graph, int source, bool useSmallLabelFirst, bool useLargeLabelLast)
{
    std::function<void(std::vector<int>& Q, const std::vector<float>& dist)> SLF = [&](std::vector<int>& Q, const std::vector<float>& dist) -> void {
        if (Q.empty()) return;
        if (dist[Q[Q.size()-1]] < dist[Q[0]]) {
            int u = Q[Q.size()-1]; Q.pop_back();
            Q.insert(Q.begin(), u);
        }
    };
    std::function<void(std::vector<int>& Q, const std::vector<float>& dist)> LLL = [&](std::vector<int>& Q, const std::vector<float>& dist) -> void {
        if (Q.empty()) return;
        float x = 0.f;
        for (const int& q : Q)
            x += dist[q];
        x /= (float)Q.size();
        while (dist[Q[0]] > x) {
            int u = Q[0]; Q.erase(Q.begin());
            Q.push_back(u);
        }
    };
    int n = graph.size();
    std::vector<float> distances(n, std::numeric_limits<float>::max());
    std::vector<int> prec(n, -1);
    distances[source] = 0.0;
    std::vector<int> Q = {source};
    while (!Q.empty()) {
        int u = Q[0];
        Q.erase(Q.begin());
        for (int v = 0; v < n; v++) {
            if (distances[u] + graph.getDistance(graph[u], graph[v]) /*graph[u]->getNeighborDistanceByIndex(v)*/ < distances[v]) {
                distances[v] = distances[u] + graph.getDistance(graph[u], graph[v]); //graph[u]->getNeighborDistanceByIndex(v);
                prec[v] = u;
                if (std::find(Q.begin(), Q.end(), v) == Q.end())
                    Q.push_back(v);
            }
        }
        if (useSmallLabelFirst)
            SLF(Q, distances);
        if (useLargeLabelLast)
            LLL(Q, distances);
    }
    return std::make_pair(distances, prec);
}

bool Pathfinding::containsNegativeCycle(const std::vector<float> &distanceVector, GridF& adjencyMap)
{
    int n = int(distanceVector.size());
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (distanceVector[i] > distanceVector[j] + adjencyMap.at(j, i))
                return true;
    return false;
}

std::pair<GridF, GridI> Pathfinding::FloydWarshall(Graph& graph)
{
    GridF W = graph.getAdjacencyMatrix();
    int n = W.sizeX;
    GridI prec(n, n, 1, -1);
    prec.iterate([&](int u, int v, int _) {
        if (W(u, v) < 100000) prec(u, v) = u;
        if (u == v) prec(u, v) = -1;
    });
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (W.at(i, j) > W.at(i, k) + W.at(k, j)) {
                    W.at(i, j) = W.at(i, k) + W.at(k, j);
                    prec(i, j) = prec(k, j);
                }
                /*if (W.at(i, j) < W.at(i, k) + W.at(k, j)) {
                    W.at(i, j) = W.at(i, j);
                    prec.at(i, j) = j;
                } else {
                    W.at(i, j) = W.at(i, k) + W.at(k, j);
                    prec.at(i, j) = k;
                }*/
            }
        }
    }
    return std::make_pair(W, prec);
}

std::pair<GridF, GridI > Pathfinding::FloydWarshallImproved(Graph& graph)
{
    GridF W = graph.getAdjacencyMatrix();
    int n = W.sizeX;
    GridI prec(n, n, 1, -1);
    prec.iterate([&](int u, int v, int _) {
        if (W(u, v) < 100000) prec(u, v) = u;
        if (u == v) prec(u, v) = -1;
    });
    std::vector<std::vector<int>> in(n), out(n);
    std::vector<int> allK(n);
    for (int i = 0; i < n; i++) allK[i] = i;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && W.at(i, j) < std::numeric_limits<float>::max()) {
                out[i].push_back(j);
                in[j].push_back(i);
            }
        }
    }
    while (!allK.empty())
    {
        int k = std::numeric_limits<int>::max(), minK = std::numeric_limits<int>::max();
        // Look for the k which minimize in[k] * out[k], for all k in allK
        for (const int& _k : allK) {
            if (minK > int(in[_k].size() * out[_k].size())) {
                minK = int(in[_k].size() * out[_k].size());
                k = _k;
            }
        }
        allK.erase(std::find(allK.begin(), allK.end(), k));

        for (const int& i : in[k]) {
            for (const int& j : out[k]) {
                if (W.at(i, j) > W.at(i, k) + W.at(k, j)) {
//                    float x = W.at(i, j);
                    if (W.at(i, j) > std::numeric_limits<float>::max()/2.f) {
                        out[i].push_back(j);
                        in[j].push_back(i);
                    }
                    W.at(i, j) = W.at(i, k) + W.at(k, j);
                    prec.at(i, j) = prec(k, j);//k;
                }
            }
        }
    }
    return std::make_pair(W, prec);
}

std::pair<GridF, GridI > Pathfinding::Johnson(Graph& graph)
{
    GridF W = graph.getAdjacencyMatrix();
    int n = W.sizeX;
    auto newNodePtr = graph.addNode(graph.size());
    for (auto& node : graph.nodes) {
        graph.addConnection(node, newNodePtr, 0.f);
        graph.addConnection(newNodePtr, node, 0.f);
    }
    GridI prec(n, n, 1, -1);
    int newNode = n + 1;

    // G' = G with an extra node connected to everybody
    GridF newAdjencyMap(W);
    // Add a column and a row to create the links between the nodes
    newAdjencyMap.insertRow(n, 0, 0.0);
    newAdjencyMap.insertRow(n, 1, 0.0);

    std::vector<float> BellmanDistances;
    std::vector<int> BellmanPrec;
    // Run Bellman-Ford from this new node
    std::tie(BellmanDistances, BellmanPrec) = BellmanFord(graph, newNode - 1);

    if (containsNegativeCycle(BellmanDistances, newAdjencyMap)) {
        throw std::invalid_argument("There is a negative cycle in the graph, sorry mate...");
    }
    // Remove the negative values, thanks to that first Bellman pass
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            W.at(i, j) = W.at(i, j) + (BellmanDistances[i] - BellmanDistances[j]);
        }
    }
    for (auto& node : graph.nodes)
        graph.removeConnection(node, newNodePtr);
    graph.removeNode(newNodePtr);
    // Run Djikstra from each node as source
    for (int i = 0; i < n; i++) {
        std::vector<float> DjikstraDistances;
        std::vector<int> DjikstraPrec;
        std::tie(DjikstraDistances, DjikstraPrec) = Djikstra(graph, i);
        for(int x = 0; x < n; x++) {
            W.at(i, x) = DjikstraDistances[x];
            prec.at(i, x) = DjikstraPrec[x];
        }
    }
    // Put back the negative values, now that Djikstra is done
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            W.at(i, j) = W.at(i, j) - (BellmanDistances[i] - BellmanDistances[j]);
        }
    }
    return std::make_pair(W, prec);
}
