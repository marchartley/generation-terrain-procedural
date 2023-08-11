#include "Matrix3Graph.h"
#include "Graph/Pathfinding.h"

Matrix3Graph::Matrix3Graph()
{
    /*
    std::vector<std::shared_ptr<GraphNode<T>>> nodes;

    GridI nodesIndices;
    GridI connectionMatrix;
    GridI precedenceMatrix;
    GridF adjencyMatrix;
    */
}

Matrix3Graph::Matrix3Graph(GridI matrix)
{
    this->initFromBinary(matrix);
}

Matrix3Graph& Matrix3Graph::initFromBinary(GridI matrix)
{
    this->originalMatrix = matrix; // Save it for later use
    this->nodes.clear();
    for (size_t i = 0; i < matrix.size(); i++) {
        GraphNode<int> node(matrix[i], matrix.getCoordAsVector3(i),i);
        this->nodes.push_back(std::make_shared<GraphNode<int>>(node));
    }
    // this->connectionMatrix = GridI(this->nodes.size(), this->nodes.size());
    // this->adjencyMatrix = GridF(this->nodes.size(), this->nodes.size(), 1, std::numeric_limits<float>::max());

    for (int x = 0; x < matrix.sizeX; x++) {
        for (int y = 0; y < matrix.sizeY; y++) {
            for (int z = 0; z < matrix.sizeZ; z++) {
                Vector3 nodePos(x, y, z);
                int node_id = matrix.getIndex(x, y, z);
                int node_label = this->nodes[node_id]->value;

                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            Vector3 dt(dx, dy, dz);
                            Vector3 neighborPos = nodePos + dt;
                            if (matrix.checkCoord(neighborPos)) {
                                int neighbor_id = matrix.getIndex(neighborPos);
                                int neighbor_label = this->nodes[neighbor_id]->value;
                                if (neighbor_label == node_label) {
                                    this->nodes[node_id]->addNeighbor(this->nodes[neighbor_id], dt.norm());
                                }
                                /*if (neighbor_label == node_label) {
                                    this->connectionMatrix.at(node_id, neighbor_id) = 1;
                                    this->connectionMatrix.at(neighbor_id, node_id) = 1;
                                    this->adjencyMatrix.at(node_id, neighbor_id) = dt.norm();
                                    this->adjencyMatrix.at(neighbor_id, node_id) = dt.norm();
                                } else {
                                    this->connectionMatrix.at(node_id, neighbor_id) = 0;
                                    this->connectionMatrix.at(neighbor_id, node_id) = 0;
                                }*/
                            }
                        }
                    }
                }

            }
        }
    }
    return *this;
}

Matrix3Graph& Matrix3Graph::computeSurface()
{
    return this->computeSurface(this->originalMatrix);
    /*
    for (auto& node : this->nodes) {
        for (int i = node->neighbors.size() - 1; i >= 0; i--) {
            if (std::get<1>(node->neighbors[i]) > 1.01f) {
                node->neighbors.erase(node->neighbors.begin() + i);
            }
        }
    }*/
    /*
    for (size_t i = 0; i < this->adjencyMatrix.size(); i++) {
        if (this->connectionMatrix[i] == 1 && this->adjencyMatrix[i] > 1.01f) {
            this->connectionMatrix[i] = 0;
            this->adjencyMatrix[i] = std::numeric_limits<float>::max();
        }
    }*/
    return *this;
}

Matrix3Graph& Matrix3Graph::computeSurface(GridI matrix)
{
    GridF distMap = matrix.toDistanceMap();
    GridI binaryMap(distMap.sizeX, distMap.sizeY, distMap.sizeZ);
    for (size_t i = 0; i < binaryMap.size(); i++) {
        if (std::abs(distMap[i] - 1.f) < 0.01f)
            binaryMap[i] = 1;
    }
    return this->initFromBinary(binaryMap);

}

Matrix3Graph& Matrix3Graph::randomizeEdges(float randomFactor, bool allowNegatives)
{
    for (auto& node : this->nodes) {
        for (size_t i = 0; i < node->neighbors.size(); i++) {
            float previous_distance = std::get<1>(node->neighbors[i]);
            float rand = 1 + (random_gen::generate(-1.f, 1.f) * randomFactor);
            float new_distance = previous_distance * rand;
            if (!allowNegatives && new_distance < 0)
                new_distance = 0;
            node->neighbors[i] = std::make_pair(std::get<0>(node->neighbors[i]), new_distance);
        }
    }
    return *this;
}

std::vector<Vector3> Matrix3Graph::shortestPath(Vector3 start, Vector3 end)
{
    float closestDistToStart = std::numeric_limits<float>::max();
    float closestDistToEnd = std::numeric_limits<float>::max();
    Vector3 newStart;
    Vector3 newEnd;
    for (auto& node : this->nodes) {
        Vector3 pos = node->pos;
        if ((pos - start).norm2() < closestDistToStart) {
            closestDistToStart = (pos - start).norm2();
            newStart = pos;
        }
        if ((pos - end).norm2() < closestDistToEnd) {
            closestDistToEnd = (pos - end).norm2();
            newEnd = pos;
        }
    }
    start = newStart;
    end = newEnd;
    float totalDistance;
    std::vector<int> path;
    std::tie(totalDistance, path) = Pathfinding::ShortestPathFrom(this->originalMatrix.getIndex(start), this->originalMatrix.getIndex(end),
                                  this->nodes/*, [&](int index) -> float {
        return (this->nodes[index]->pos - end).norm2();
    }*/);

    path = Pathfinding::getPath(this->originalMatrix.getIndex(end), path);
    std::vector<Vector3> returnedPath(path.size());
    for (size_t i = 0; i < path.size(); i++)
        returnedPath[i] = this->nodes[path[i]]->pos;
    return returnedPath;
}
