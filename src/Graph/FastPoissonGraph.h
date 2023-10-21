#ifndef FASTPOISSONGRAPH_H
#define FASTPOISSONGRAPH_H

template<class T>
class FastPoissonGraph;

#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"
#include "DataStructure/Vector3.h"
#include "Utils/Globals.h"
#include "Graph/Pathfinding.h"
#include "Graph/Graph.h"

template<class T>
class FastPoissonGraph : public Graph<T>
{
public:
    FastPoissonGraph();
    FastPoissonGraph(int sizeX, int sizeY, int sizeZ, float radius = 1.f, int max_tries = 30);
    template<class U>
    FastPoissonGraph(Matrix3<U>& available_space_matrix, float radius = 1.f, int max_tries = 30);

    void createEdges(int maxNumberNeighbors, float maxDistanceToNeighbor, bool justUpdateConnectionMatrix = false);

//protected:
    template<class U>
    void initNodes(Matrix3<U>& available_space_matrix, float radius, int max_tries);

    int sizeX, sizeY, sizeZ;
    float cellWidth;

    float poissonRadius;
    float maxNeighborRadius;
    int poissonNeighbors;
};


template<class T>
FastPoissonGraph<T>::FastPoissonGraph()
{
}

template<class T>
FastPoissonGraph<T>::FastPoissonGraph(int sizeX, int sizeY, int sizeZ, float radius, int max_tries)
{
    Matrix3<short int> available_space_matrix(sizeX, sizeY, sizeZ, 1);
    initNodes(available_space_matrix, radius, max_tries);
}


template<class T> template<class U>
FastPoissonGraph<T>::FastPoissonGraph(Matrix3<U> &available_space_matrix, float radius, int max_tries)
{
    initNodes(available_space_matrix, radius, max_tries);
}

template<class T> template<class U>
void FastPoissonGraph<T>::initNodes(Matrix3<U> &_available_space_matrix, float radius, int max_tries)
{
    if (_available_space_matrix.getDimensions().minComp() <= radius)
        radius = _available_space_matrix.getDimensions().minComp();
    this->poissonRadius = radius;
    // Fllowing Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions. SIGGRAPH sketches, 10(1).
    // Step 0 : Create a n-dim grid with cells of size r/sqrt(n) => here n = 3
    this->cellWidth = radius / std::sqrt(3);
    int prevSizeX = _available_space_matrix.sizeX;
    int prevSizeY = _available_space_matrix.sizeY;
    int prevSizeZ = _available_space_matrix.sizeZ;
    Matrix3<U> available_space_matrix = _available_space_matrix.resize(std::floor(prevSizeX / cellWidth),
                                                                       std::floor(prevSizeY / cellWidth),
                                                                       std::floor(prevSizeZ / cellWidth)/*,
                                                                       RESIZE_MODE::MIN_VAL*/);
//    available_space_matrix = available_space_matrix.rounded();
    available_space_matrix.raiseErrorOnBadCoord = false;

    radius = 1.0;
    float sqr_diameter = radius * radius;
    this->sizeX = available_space_matrix.sizeX;
    this->sizeY = available_space_matrix.sizeY;
    this->sizeZ = available_space_matrix.sizeZ;
    this->nodesIndices = GridI(this->sizeX, this->sizeY, this->sizeZ, -1); // Fill all indices to -1 as default values
    Vector3 sizeVec(this->sizeX, this->sizeY, this->sizeZ);

    std::vector<std::shared_ptr<GraphNode<T>>> allNodes;
    std::vector<std::shared_ptr<GraphNode<T>>> activeList;
    // Step 1 : Create initial sample, register it and put it in the active list
    Vector3 randomStart;
    int tries = 0;
    do {
        randomStart = Vector3(random_gen::generate() * sizeX, random_gen::generate() * sizeY, random_gen::generate() * sizeZ);
        tries ++;
        if (tries > 100000)
            return; // Error, not enough space to place nodes
    } while (!available_space_matrix.at(randomStart));
    std::shared_ptr<GraphNode<T>> firstNode = std::make_shared<GraphNode<T>>(T(), randomStart, 0);
    activeList.push_back(firstNode);
    allNodes.push_back(firstNode);
    int current_node_index = 1;

    // Step 2 : process all samples recursively
    int steps = 0;
    while (!activeList.empty()) {
        steps++;
        int randomIndex = random_gen::generate(0, activeList.size());
        std::shared_ptr<GraphNode<T>> currentNode = activeList[randomIndex];
        bool atLeastOneNodeAdded = false;
        for (int i = 0; i < max_tries && !atLeastOneNodeAdded; i++) {
            Vector3 pos = (currentNode->pos + Vector3::random() * radius * random_gen::generate(1.0, 2.0)).abs();
            if (!this->nodesIndices.checkCoord(pos)) // If it's out of frame, don't keep it
                continue;
            if (this->nodesIndices.at(pos) > -1) // If the cell is occupied, we know we are too close to another sample
                continue;
            if (!available_space_matrix.at(pos))
                continue;

            bool keep = true;
            Vector3 indicePos = pos.rounded();
            for (int x = indicePos.x - 1; x <= indicePos.x + 1 && keep; x++) {
                for (int y = indicePos.y - 1; y <= indicePos.y + 1 && keep; y++) {
                    for (int z = indicePos.z - 1; z <= indicePos.z + 1 && keep; z++) {
                        if (!this->nodesIndices.checkCoord(x, y, z))
                            continue;
                        if (x == indicePos.x && y == indicePos.y && z == indicePos.z)
                            continue;
                        if (this->nodesIndices.at(x, y, z) != -1) {
                            if ((pos - allNodes[this->nodesIndices.at(x, y, z)]->pos).norm2() < sqr_diameter) {
                                keep = false;
                            } else {

                            }
                        }
                    }
                }
            }
            if (keep) {
                std::shared_ptr<GraphNode<T>> newSample = std::make_shared<GraphNode<T>>(T(), pos, current_node_index);
                newSample->privateVector = pos;
                activeList.push_back(newSample);
                allNodes.push_back(newSample);
                this->nodesIndices.at(pos) = current_node_index;
                current_node_index++;
                atLeastOneNodeAdded = true;
            } else {

            }
        }
        if (!atLeastOneNodeAdded) {
            currentNode->privateIndex = this->nodes.size();
            currentNode->privateVector = currentNode->pos;
            this->nodes.push_back(currentNode);
            activeList.erase(activeList.begin() + randomIndex);
        }
    }

    // At this point, put back all positions in their initial reference and update the indices
    for (auto& val : this->nodes) {
        val->pos /= sizeVec;
        val->pos *= Vector3(_available_space_matrix.sizeX, _available_space_matrix.sizeY, _available_space_matrix.sizeZ);

        this->nodesIndices.at(val->privateVector) = val->privateIndex;
        val->index = val->privateIndex;
    }
}

template<class T>
void FastPoissonGraph<T>::createEdges(int maxNumberNeighbors, float maxDistanceToNeighbor, bool justUpdateConnectionMatrix)
{
    this->maxNeighborRadius = maxDistanceToNeighbor;
    this->poissonNeighbors = maxNumberNeighbors;
    float sqr_dist = maxNeighborRadius * maxNeighborRadius;
    this->connectionMatrix = GridI(this->nodes.size(), this->nodes.size(), 1, 0);
    this->adjencyMatrix = GridF(this->nodes.size(), this->nodes.size(), 1, std::numeric_limits<float>::max());
    this->precedenceMatrix = GridI(this->nodes.size(), this->nodes.size(), 1, -1);

    int cellsToCheckOnX = std::ceil(this->maxNeighborRadius / this->poissonRadius);
    int cellsToCheckOnY = std::ceil(this->maxNeighborRadius / this->poissonRadius);
    int cellsToCheckOnZ = std::ceil(this->maxNeighborRadius / this->poissonRadius);

    for (size_t i = 0; i < this->nodes.size(); i++) {
        auto sample = this->nodes[i];
        std::vector<std::shared_ptr<GraphNode<T>>> candidates;
        for (size_t j = 0; j < this->nodes.size(); j++){
            if (i == j) continue;
            if ((sample->pos - this->nodes[j]->pos).norm2() < sqr_dist)
                candidates.push_back(this->nodes[j]);
        }
        /*
        for (int x = sample->privateVector.x-cellsToCheckOnX; x < sample->privateVector.x+cellsToCheckOnX+1; x++) {
            for (int y = sample->privateVector.y-cellsToCheckOnY; y < sample->privateVector.y+cellsToCheckOnY+1; y++) {
                for (int z = sample->privateVector.z-cellsToCheckOnZ; z < sample->privateVector.z+cellsToCheckOnZ+1; z++) {
                    int index = this->nodesIndices.checkCoord(x, y, z) ? this->nodesIndices.at(x, y, z) : -1;
//                    if (index == sample->index)
//                        continue;

                    // First, check if there is a node in here
                    if(index != -1) {
                        // Now check with real precision
                        if ((this->nodes[index]->pos - sample->pos).norm2() < sqr_dist)
                            candidates.push_back(this->nodes[index]);
                    }
                }
            }
        }*/
        // Sort the neighbors by the distance to the current node
        std::sort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) -> bool {
            return (a->pos - sample->pos) < (b->pos - sample->pos);
        });

        for (int iNeighbor = 0; iNeighbor < std::min(this->poissonNeighbors+1, int(candidates.size())); iNeighbor++) {
            // Just a one-way connection for now
            this->connectionMatrix.at(sample->privateIndex, candidates[iNeighbor]->privateIndex) = 1;
            if (!justUpdateConnectionMatrix)
                this->adjencyMatrix.at(sample->privateIndex, candidates[iNeighbor]->privateIndex) = (candidates[iNeighbor]->pos - sample->pos).norm();
        }
    }
}
#endif // FASTPOISSONGRAPH_H
