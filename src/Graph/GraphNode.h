#ifndef GRAPHNODE_H
#define GRAPHNODE_H

#include "DataStructure/Vector3.h"
template <class T>
class GraphNodeTemplate
{
public:
    GraphNodeTemplate();
    GraphNodeTemplate(T value);
    GraphNodeTemplate(T value, const Vector3& pos, int index = 0);

    void addNeighbor(GraphNodeTemplate<T>* neighbor);
    void addNeighbor(GraphNodeTemplate<T>* neighbor, float distance);

    void removeNeighbor(GraphNodeTemplate<T>* neighbor);

    bool hasNeighbor(GraphNodeTemplate<T>* neighbor);

    float getNeighborDistanceByIndex(int index);

    void cleanNeighborhood();

    Vector3 pos;
    T value;
    int index;

    Vector3 privateVector;
    int privateIndex;

    std::vector<std::pair<GraphNodeTemplate<T>*, float>> neighbors;
};

template<class T>
GraphNodeTemplate<T>::GraphNodeTemplate() : GraphNodeTemplate(T(), Vector3(), 0)
{

}
template<class T>
GraphNodeTemplate<T>::GraphNodeTemplate(T value) : GraphNodeTemplate(value, Vector3(), 0)
{

}
template<class T>
GraphNodeTemplate<T>::GraphNodeTemplate(T value, const Vector3& pos, int index)
    : value(value), pos(pos), index(index)
{

}

template<class T>
void GraphNodeTemplate<T>::addNeighbor(GraphNodeTemplate<T>* neighbor)
{
    this->neighbors.push_back(std::make_pair(neighbor, (this->pos - neighbor->pos).norm()));
}

template<class T>
void GraphNodeTemplate<T>::addNeighbor(GraphNodeTemplate<T>* neighbor, float distance)
{
    this->neighbors.push_back(std::make_pair(neighbor, distance));
}

template<class T>
void GraphNodeTemplate<T>::removeNeighbor(GraphNodeTemplate<T>* neighbor)
{
    for (int i = this->neighbors.size() - 1; i >= 0; i--)
        if (this->neighbors[i].first == neighbor)
            this->neighbors.erase(this->neighbors.begin() + i);
//    this->neighbors.erase(std::find(this->neighbors.begin(), this->neighbors.end(), neighbor));
}

template<class T>
bool GraphNodeTemplate<T>::hasNeighbor(GraphNodeTemplate<T>* neighbor)
{
    for (auto& [n, w] : this->neighbors)
        if (n == neighbor)
            return true;
    return false;
}

template<class T>
float GraphNodeTemplate<T>::getNeighborDistanceByIndex(int index)
{
    float distance = 0.f;
    bool neighborFound = false;
    for (auto& neighbor : this->neighbors) {
        if (std::get<0>(neighbor)->index == index) {
            neighborFound = true;
            distance += std::get<1>(neighbor);
        }
    }
    if (neighborFound)
        return distance;
    return std::numeric_limits<float>::max();
}

template<class T>
void GraphNodeTemplate<T>::cleanNeighborhood()
{
    std::vector<GraphNodeTemplate<T>*> nodes;
    std::vector<float> distances;
    for (int i = this->neighbors.size() - 1; i >= 0; i--) {
        float dist;
        GraphNodeTemplate<T>* neighbor;
        std::tie(neighbor, dist) = this->neighbors[i];

        auto foundIndex = std::find(nodes.begin(), nodes.end(), neighbor);
        if (foundIndex != std::end(nodes)) {
            distances[std::distance(nodes.begin(), foundIndex)] += dist;
        } else {
            nodes.push_back(neighbor);
            distances.push_back(dist);
        }
    }

    this->neighbors.resize(nodes.size());
    for (int i = 0; i < nodes.size(); i++) {
        this->neighbors[i] = {nodes[i], distances[i]};
    }
}

typedef GraphNodeTemplate<int> GraphNode;

#endif // GRAPHNODE_H
