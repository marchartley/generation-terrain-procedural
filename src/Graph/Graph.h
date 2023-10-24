#ifndef GRAPH_H
#define GRAPH_H

#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"
#include "Utils/ConstraintsSolver.h"

#include "Graphics/DisplayGraphics.h"

template <class T>
class Graph
{
public:
    Graph();

    Graph<T>& circularLayout(bool randomOrder = false);
    Graph<T> &forceDrivenPositioning(bool startWithCircularLayout = true);

    std::shared_ptr<GraphNode<T>> addNode(int node);
    std::shared_ptr<GraphNode<T>> addNode(std::shared_ptr<GraphNode<T>> newNode);
    std::vector<std::shared_ptr<GraphNode<T>> > addNodes(std::vector<int> nodes);
    std::vector<std::shared_ptr<GraphNode<T> > > addNodes(std::vector<std::shared_ptr<GraphNode<T>>> newNodes);

    void addConnection(int nodeA, int nodeB, float weight = 1.f);
    void addConnection(std::shared_ptr<GraphNode<T>> nodeA, std::shared_ptr<GraphNode<T>> nodeB, float weight = 1.f);
    void removeConnection(int nodeA, int nodeB);
    void removeConnection(std::shared_ptr<GraphNode<T>> nodeA, std::shared_ptr<GraphNode<T>> nodeB);

    std::shared_ptr<GraphNode<T>> findNodeByID(int ID);

    std::vector<std::shared_ptr<GraphNode<T>>> nodes;

    GridI nodesIndices;
    GridI connectionMatrix;
    GridI precedenceMatrix;
    GridF adjencyMatrix;

    void draw();
};

template<class T>
Graph<T>::Graph() {

}

template<class T>
Graph<T> &Graph<T>::circularLayout(bool randomOrder)
{
    std::vector<size_t> indices(this->nodes.size());
    for (size_t i = 0; i < indices.size(); i++) indices[i] = i;
    if (randomOrder)
        std::shuffle(indices.begin(), indices.end(), random_gen::random_generator);
    for (size_t i = 0; i < this->nodes.size(); i++) {
        this->nodes[indices[i]]->pos = Vector3(std::cos(2.f * PI * float(i)/float(this->nodes.size())),
                                      std::sin(2.f * PI * float(i)/float(this->nodes.size())));
    }
    return *this;
}

template<class T>
Graph<T>& Graph<T>::forceDrivenPositioning(bool startWithCircularLayout)
{
    if (startWithCircularLayout)
        this->circularLayout();
    int maxTriesForAdjencyConstraint = 1000;
    ConstraintsSolver solver;
    solver.numberIterations = maxTriesForAdjencyConstraint;
    solver.stoppingEpsilon = 1.f;

    for (auto& n : nodes) {
        solver.addItem(new Vector3(n->pos));
    }

    std::map<int, Vector3> initalPositions = solver.solve(false, 0.f, 0.f);

    for (int iNode = 0; iNode < nodes.size(); iNode++) {
        for (int jNode = iNode + 1; jNode < nodes.size(); jNode++) {
            if (connectionMatrix.at(iNode, jNode) == 1)
                solver.addDistanceConstraint(iNode, jNode, adjencyMatrix.at(iNode, jNode));
        }
    }
    std::map<int, Vector3> positions;
    for (int iteration = 0; iteration < maxTriesForAdjencyConstraint; iteration++) {
        positions = solver.solve();

        // Now check for crossing adjencies
        // Wait... how... Maybe using the Voronoi algorithm
        break;
    }
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i]->pos = positions[i];
    }
    return *this;
}

template<class T>
std::shared_ptr<GraphNode<T>> Graph<T>::addNode(int node)
{
    return this->addNode(std::make_shared<GraphNode<T>>(T(), Vector3(), node));
}

template<class T>
std::shared_ptr<GraphNode<T>> Graph<T>::addNode(std::shared_ptr<GraphNode<T> > newNode)
{
    this->addNodes({newNode});
    return newNode;
}

template<class T>
std::vector<std::shared_ptr<GraphNode<T> > > Graph<T>::addNodes(std::vector<int> nodes)
{
    std::vector<std::shared_ptr<GraphNode<T>>> newNodes;
    for (auto& node : nodes)
        newNodes.push_back(std::make_shared<GraphNode<T>>(T(), Vector3(), node));
    return this->addNodes(newNodes);
}

template<class T>
std::vector<std::shared_ptr<GraphNode<T> > > Graph<T>::addNodes(std::vector<std::shared_ptr<GraphNode<T> > > newNodes)
{
    int numberOfNodesToAdd = newNodes.size();
    GridI nodesIndices = GridI(this->nodesIndices.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
    nodesIndices.paste(this->nodesIndices);
    this->nodesIndices = nodesIndices;

    GridI connectionMatrix = GridI(this->connectionMatrix.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
    connectionMatrix.paste(this->connectionMatrix);
    this->connectionMatrix = connectionMatrix;

    GridI precedenceMatrix = GridI(this->precedenceMatrix.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
    precedenceMatrix.paste(this->precedenceMatrix);
    this->precedenceMatrix = precedenceMatrix;

    GridF adjencyMatrix = GridI(this->adjencyMatrix.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
    adjencyMatrix.paste(this->adjencyMatrix);
    this->adjencyMatrix = adjencyMatrix;

    this->nodes.insert(this->nodes.end(), newNodes.begin(), newNodes.end());
    return newNodes;
}

template<class T>
void Graph<T>::addConnection(int nodeA, int nodeB, float weight)
{
    return this->addConnection(this->findNodeByID(nodeA), this->findNodeByID(nodeB), weight);
}

template<class T>
void Graph<T>::addConnection(std::shared_ptr<GraphNode<T> > nodeA, std::shared_ptr<GraphNode<T> > nodeB, float weight)
{
    int indexA = -1;
    int indexB = -1;
    for (size_t i = 0; i < this->nodes.size(); i++) {
        if (this->nodes[i] == nodeA)
            indexA = i;
        if (this->nodes[i] == nodeB)
            indexB = i;
        if (indexA > -1 && indexB > -1)
            break;
    }

    nodeA->addNeighbor(nodeB, weight);
    nodeB->addNeighbor(nodeA, weight);

    this->connectionMatrix.at(indexA, indexB) = 1;
    this->connectionMatrix.at(indexB, indexA) = 1;
    this->adjencyMatrix.at(indexA, indexB) = weight;
    this->adjencyMatrix.at(indexB, indexA) = weight;
}

template<class T>
void Graph<T>::removeConnection(int nodeA, int nodeB)
{
    return this->removeConnection(this->findNodeByID(nodeA), this->findNodeByID(nodeB));
}

template<class T>
void Graph<T>::removeConnection(std::shared_ptr<GraphNode<T> > nodeA, std::shared_ptr<GraphNode<T> > nodeB)
{
    int indexA = -1;
    int indexB = -1;
    for (size_t i = 0; i < this->nodes.size(); i++) {
        if (this->nodes[i] == nodeA)
            indexA = i;
        if (this->nodes[i] == nodeB)
            indexB = i;
        if (indexA > -1 && indexB > -1)
            break;
    }

    nodeA->removeNeighbor(nodeB);
    nodeB->removeNeighbor(nodeA);

    this->connectionMatrix.at(indexA, indexB) = 0;
    this->connectionMatrix.at(indexB, indexA) = 0;
    this->adjencyMatrix.at(indexA, indexB) = 0.f;
    this->adjencyMatrix.at(indexB, indexA) = 0.f;
}

template<class T>
std::shared_ptr<GraphNode<T> > Graph<T>::findNodeByID(int ID)
{
    for (auto& node : this->nodes) {
        if (node->index == ID)
            return node;
    }
    return nullptr;
}

template<class T>
void Graph<T>::draw()
{
    /*
     * TODO : Install QT Charts to have access to the plotter
    Plotter plt;
    std::vector<Vector3> nodesPos;
    std::vector<std::string> labels;
    std::vector<std::vector<Vector3>> linksPos;
    for (auto& node : this->nodes) {
        nodesPos.push_back(node->pos);
        labels.push_back(std::to_string(node->index));

        for (auto& [neighbor, weight] : node->neighbors) {
//            linksPos.push_back({{node->pos.x, node->pos.y}, {neighbor->pos.x, neighbor->pos.y}});
            plt.addPlot({node->pos, neighbor->pos});
        }
    }
    plt.addScatter(nodesPos,  "", labels);
    plt.exec();
    */
}
#endif // GRAPH_H
