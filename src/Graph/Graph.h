#ifndef GRAPH_H
#define GRAPH_H

#include "DataStructure/Matrix3.h"
#include "Graph/GraphNode.h"
#include "Utils/ConstraintsSolver.h"

#include "Graphics/DisplayGraphics.h"

template <class T>
class GraphTemplate
{
public:
    GraphTemplate(bool useMatrices = true);
    ~GraphTemplate();

    GraphTemplate<T>& circularLayout(bool randomOrder = false);
    GraphTemplate<T> &forceDrivenPositioning(bool startWithCircularLayout = true);

    GraphNodeTemplate<T>* addNode(int node);
    GraphNodeTemplate<T>* addNode(GraphNodeTemplate<T>* newNode);
    void removeNode(int node);
    void removeNode(GraphNodeTemplate<T>* newNode);
//    std::vector<std::shared_ptr<GraphNode<T>> > addNodes(std::vector<int> nodes);
//    std::vector<std::shared_ptr<GraphNode<T> > > addNodes(std::vector<std::shared_ptr<GraphNode<T>>> newNodes);

    std::map<int, GraphNodeTemplate<T>* > addNodes(std::vector<int> nodes);
    std::map<int, GraphNodeTemplate<T>* > addNodes(std::vector<GraphNodeTemplate<T>*> newNodes);
    void removeNodes(std::vector<int> nodes);
    void removeNodes(std::vector<GraphNodeTemplate<T>*> newNodes);


    void addConnection(int nodeA, int nodeB, float weight = 1.f);
    void addConnection(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB, float weight = 1.f);
    void removeConnection(int nodeA, int nodeB);
    void removeConnection(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB);
    void setConnection(int nodeA, int nodeB, float weight);
    void setConnection(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB, float weight);

    bool connected(int nodeA, int nodeB) const;
    bool connected(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB) const;

    GraphNodeTemplate<T>* findNodeByID(int ID) const;

    GraphNodeTemplate<T>* operator[](T index) const;

    void initAllNodesValues(T value) const;
    bool empty() const;
    size_t size() const;

    GridF getAdjacencyMatrix() const;

//    std::vector<std::shared_ptr<GraphNode<T>>> nodes;
    std::map<int, GraphNodeTemplate<T>*> nodes;

    GridI connectionMatrix;
    GridF adjencyMatrix;

    void draw();

    bool useMatrices;
};

template<class T>
GraphTemplate<T>::GraphTemplate(bool useMatrices)
    : useMatrices(useMatrices)
{
    this->useMatrices = false; // FORCE
}

template<class T>
GraphTemplate<T>::~GraphTemplate()
{
    for (auto& [ID, node] : nodes) {
//        delete node;
    }
    nodes.clear();
}

template<class T>
GraphTemplate<T> &GraphTemplate<T>::circularLayout(bool randomOrder)
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
GraphTemplate<T>& GraphTemplate<T>::forceDrivenPositioning(bool startWithCircularLayout)
{
    if (startWithCircularLayout)
        this->circularLayout();
    int maxTriesForAdjencyConstraint = 1000;
    ConstraintsSolver solver;
    solver.numberIterations = maxTriesForAdjencyConstraint;
    solver.stoppingEpsilon = 1.f;

    for (auto& [ID, n] : nodes) {
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
GraphNodeTemplate<T>* GraphTemplate<T>::addNode(int node)
{
    return this->addNode(new GraphNodeTemplate<T>(T(), Vector3(), node));
}

template<class T>
GraphNodeTemplate<T>* GraphTemplate<T>::addNode(GraphNodeTemplate<T>* newNode)
{
    this->addNodes({newNode});
    return newNode;
}

template<class T>
void GraphTemplate<T>::removeNode(int node)
{
    if (!this->nodes.count(node)) return;
    this->nodes.erase(node);
}

template<class T>
void GraphTemplate<T>::removeNode(GraphNodeTemplate<T>* newNode)
{
    this->removeNode(newNode->index);
}

template<class T>
std::map<int, GraphNodeTemplate<T>* > GraphTemplate<T>::addNodes(std::vector<int> nodes)
{
    std::vector<GraphNodeTemplate<T>*> newNodes;
    for (auto& node : nodes)
        newNodes.push_back(new GraphNodeTemplate<T>(T(), Vector3(), node));
    return this->addNodes(newNodes);
}

template<class T>
std::map<int, GraphNodeTemplate<T>* > GraphTemplate<T>::addNodes(std::vector<GraphNodeTemplate<T>* > newNodes)
{
//    this->nodes.insert(this->nodes.end(), newNodes.begin(), newNodes.end());
    std::map<int, GraphNodeTemplate<T>*> newNodesMap;
    for (const auto& node : newNodes) {
        newNodesMap[node->index] = node;
    }
    this->nodes.insert(newNodesMap.begin(), newNodesMap.end());

    if (useMatrices) {
        int numberOfNodesToAdd = newNodes.size();
//        GridI nodesIndices = GridI(this->nodesIndices.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
//        nodesIndices.paste(this->nodesIndices);
//        this->nodesIndices = nodesIndices;

        GridI connectionMatrix = GridI(this->connectionMatrix.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
        connectionMatrix.paste(this->connectionMatrix);
        this->connectionMatrix = connectionMatrix;

        GridF adjencyMatrix = GridI(this->adjencyMatrix.getDimensions() + Vector3(numberOfNodesToAdd, numberOfNodesToAdd, 0));
        adjencyMatrix.paste(this->adjencyMatrix);
        this->adjencyMatrix = adjencyMatrix;
    }

    return newNodesMap;
}

template<class T>
void GraphTemplate<T>::removeNodes(std::vector<int> nodes)
{
    for (auto node : nodes)
        this->removeNode(node);
}

template<class T>
void GraphTemplate<T>::removeNodes(std::vector<GraphNodeTemplate<T>* > newNodes)
{
    for (auto node : nodes)
        this->removeNode(node);
}

template<class T>
void GraphTemplate<T>::addConnection(int nodeA, int nodeB, float weight)
{
    return this->addConnection(this->findNodeByID(nodeA), this->findNodeByID(nodeB), weight);
}

template<class T>
void GraphTemplate<T>::addConnection(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB, float weight)
{
    nodeA->addNeighbor(nodeB, weight);
//    nodeB->addNeighbor(nodeA, weight);

    if (useMatrices) {
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
        this->connectionMatrix.at(indexA, indexB) = 1;
        this->connectionMatrix.at(indexB, indexA) = 1;
        this->adjencyMatrix.at(indexA, indexB) = weight;
        this->adjencyMatrix.at(indexB, indexA) = weight;
    }
}

template<class T>
void GraphTemplate<T>::removeConnection(int nodeA, int nodeB)
{
    return this->removeConnection(this->findNodeByID(nodeA), this->findNodeByID(nodeB));
}

template<class T>
void GraphTemplate<T>::removeConnection(GraphNodeTemplate<T>* nodeA, GraphNodeTemplate<T>* nodeB)
{
    nodeA->removeNeighbor(nodeB);
    nodeB->removeNeighbor(nodeA);

    if (useMatrices) {
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
        this->connectionMatrix.at(indexA, indexB) = 0;
        this->connectionMatrix.at(indexB, indexA) = 0;
        this->adjencyMatrix.at(indexA, indexB) = 0.f;
        this->adjencyMatrix.at(indexB, indexA) = 0.f;
    }
}

template<class T>
void GraphTemplate<T>::setConnection(int nodeA, int nodeB, float weight)
{
    return this->setConnection(findNodeByID(nodeA), findNodeByID(nodeB), weight);
}

template<class T>
void GraphTemplate<T>::setConnection(GraphNodeTemplate<T> *nodeA, GraphNodeTemplate<T> *nodeB, float weight)
{
    nodeA->removeNeighbor(nodeB);
    nodeA->addNeighbor(nodeB, weight); // Not best, but hey!
}

template<class T>
bool GraphTemplate<T>::connected(GraphNodeTemplate<T> *nodeA, GraphNodeTemplate<T> *nodeB) const
{
    return nodeA->hasNeighbor(nodeB);
}

template<class T>
bool GraphTemplate<T>::connected(int nodeA, int nodeB) const
{
    return this->connected(findNodeByID(nodeA), findNodeByID(nodeB));
}

template<class T>
GraphNodeTemplate<T>* GraphTemplate<T>::findNodeByID(int ID) const
{
    if (this->nodes.count(ID)) return this->nodes.at(ID);
    /*for (auto& node : this->nodes) {
        if (node->index == ID)
            return node;
    }*/
    return nullptr;
}

template<class T>
GraphNodeTemplate<T>* GraphTemplate<T>::operator[](T index) const
{
    if (!nodes.count(index)) return nullptr;
    return this->nodes.at(index);
}

template<class T>
void GraphTemplate<T>::initAllNodesValues(T value) const
{
    for (auto& [ID, node] : nodes) {
        node->value = value;
    }
}

template<class T>
bool GraphTemplate<T>::empty() const
{
    return this->nodes.empty();
}

template<class T>
size_t GraphTemplate<T>::size() const
{
    return this->nodes.size();
}

template<class T>
GridF GraphTemplate<T>::getAdjacencyMatrix() const
{
    if (this->useMatrices) return this->adjencyMatrix;

    GridF matrix(this->size(), this->size(), 1);
    for (auto& [IDA, nodeA] : nodes) {
        for (auto& [IDB, nodeB] : nodes) {
            if (IDA == IDB) matrix(IDA, IDB) = 0;
            else matrix(IDA, IDB) = nodeA->getNeighborDistanceByIndex(IDB);
        }
    }
    return matrix;
}

template<class T>
void GraphTemplate<T>::draw()
{
    std::vector<Vector3> nodesPos;
    std::vector<std::string> labels;
    std::vector<std::vector<Vector3>> linksPos;
    for (auto& [ID, node] : this->nodes) {
        nodesPos.push_back(node->pos);
        labels.push_back(std::to_string(node->index));

        for (auto& [neighbor, weight] : node->neighbors) {
//            linksPos.push_back({{node->pos.x, node->pos.y}, {neighbor->pos.x, neighbor->pos.y}});
            Plotter::get()->addPlot({node->pos, neighbor->pos});
        }
    }
    Plotter::get()->addScatter(nodesPos,  "", labels);
    Plotter::get()->show();
}

typedef GraphTemplate<int> Graph;

#endif // GRAPH_H
