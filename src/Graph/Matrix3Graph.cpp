#include "Matrix3Graph.h"
#include "Graph/Pathfinding.h"

Matrix3Graph::Matrix3Graph()
{
    /*
    std::vector<std::shared_ptr<GraphNode<T>>> nodes;

    Matrix3<int> nodesIndices;
    Matrix3<int> connectionMatrix;
    Matrix3<int> precedenceMatrix;
    Matrix3<float> adjencyMatrix;
    */
}

Matrix3Graph::Matrix3Graph(Matrix3<int> matrix)
{
    this->initFromBinary(matrix);
}

Matrix3Graph& Matrix3Graph::initFromBinary(Matrix3<int> matrix)
{
    this->originalMatrix = matrix; // Save it for later use
    this->nodes.clear();
    for (size_t i = 0; i < matrix.size(); i++) {
        GraphNode<int> node(matrix[i], matrix.getCoordAsVector3(i),i);
        this->nodes.push_back(std::make_shared<GraphNode<int>>(node));
    }
    // this->connectionMatrix = Matrix3<int>(this->nodes.size(), this->nodes.size());
    // this->adjencyMatrix = Matrix3<float>(this->nodes.size(), this->nodes.size(), 1, std::numeric_limits<float>::max());

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

Matrix3Graph& Matrix3Graph::computeSurface(Matrix3<int> matrix)
{
    Matrix3<float> distMap = matrix.toDistanceMap();
    Matrix3<int> binaryMap(distMap.sizeX, distMap.sizeY, distMap.sizeZ);
    for (size_t i = 0; i < binaryMap.size(); i++) {
        if (std::abs(distMap[i] - 1.f) < 0.01f)
            binaryMap[i] = 1;
    }
    return this->initFromBinary(binaryMap);

}

std::vector<Vector3> Matrix3Graph::shortestPath(Vector3 start, Vector3 end)
{
    float totalDistance;
    std::vector<int> path;
    std::tie(totalDistance, path) = Pathfinding::ShortestPathFrom(this->originalMatrix.getIndex(start), this->originalMatrix.getIndex(end),
                                  this->nodes, [&](int index) -> float {
        return (this->nodes[index]->pos - end).norm2();
    });

    path = Pathfinding::getPath(this->originalMatrix.getIndex(end), path);
    std::vector<Vector3> returnedPath(path.size());
    for (size_t i = 0; i < path.size(); i++)
        returnedPath[i] = this->nodes[path[i]]->pos;
    return returnedPath;
}
