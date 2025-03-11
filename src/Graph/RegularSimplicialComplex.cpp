#include "RegularSimplicialComplex.h"

RegularSimplicialComplex::RegularSimplicialComplex()
    : RegularSimplicialComplex(10, 10)
{

}

RegularSimplicialComplex::RegularSimplicialComplex(int sizeX, int sizeY)
    : GraphTemplate(), sizeX(sizeX), sizeY(sizeY)
{
    std::vector<int> newNodes(sizeX * sizeY);
    for (int i = 0; i < sizeX * sizeY; i++) newNodes[i] = i;
    this->addNodes(newNodes);

    /*for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
//            auto node = this->addNode(x * sizeX + y);
            auto node = this->getNode(x, y);
            node->value = -1;
            node->pos = Vector3(x, y);
        }
    }*/
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            auto node = this->getNode(x, y);
            node->value = -1;
            node->pos = Vector3(x, y);
            if (x > 0) {
                this->addConnection(this->getNode(x - 1, y), node);
                this->addConnection(this->getNode(y, x - 1), node);
            }
            if (y < sizeY - 1 && x > 0) {
                this->addConnection(this->getNode(x - 1, y + 1), node);
                this->addConnection(this->getNode(y + 1, x - 1), node);
            }
            if (y < sizeY - 1) {
                this->addConnection(this->getNode(x, y + 1), node);
                this->addConnection(this->getNode(y + 1, x), node);
            }
        }
    }
}

GraphNode* RegularSimplicialComplex::getNode(int x, int y)
{
    if (x < 0 || x >= this->sizeX || y < 0 || y >= this->sizeY)
        return nullptr;
    return this->findNodeByID(x * this->sizeX + y);
}

GraphNode* RegularSimplicialComplex::getNode(const Vector3& pos)
{
    return this->getNode(pos.x, pos.y);
}

void RegularSimplicialComplex::display()
{
    float maxNodeValue = this->maxNodesValues() + 1.f;
    std::vector<std::vector<Vector3>> plots;
    std::vector<Vector3> scatter;
    std::vector<std::string> scatter_labels;
    std::vector<QColor> plots_colors;
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            auto myNode = this->getNode(x, y);
            scatter.push_back(myNode->pos);
            scatter_labels.push_back("(" + std::to_string(x) + ", " + std::to_string(y) + ")");

            if (x > 0) {
                auto secondNode = this->getNode(x - 1, y);
                if (myNode->hasNeighbor(secondNode)) {
                    plots.push_back({myNode->pos, secondNode->pos});
                    if (myNode->value < 0 || secondNode->value < 0 || myNode->value != secondNode->value) {
                        plots_colors.push_back(Qt::gray); // plots_colors.push_back(PlotColor::GRAY);
                    } else {
                        plots_colors.push_back(QColor::fromHsvF((float)myNode->value / maxNodeValue, 1.f, 1.f)); // plots_colors.push_back((PlotColor)(PlotColor::RED + myNode->value));
                    }
                }
            }
            if (y < sizeY - 1 && x > 0) {
                auto secondNode = this->getNode(x - 1, y + 1);
                if (myNode->hasNeighbor(secondNode)) {
                    plots.push_back({myNode->pos, secondNode->pos});
                    if (myNode->value < 0 || secondNode->value < 0 || myNode->value != secondNode->value) {
                        plots_colors.push_back(Qt::gray); // plots_colors.push_back(PlotColor::GRAY);
                    } else {
                        plots_colors.push_back(QColor::fromHsvF((float)myNode->value / maxNodeValue, 1.f, 1.f)); // plots_colors.push_back((PlotColor)(PlotColor::RED + myNode->value));
                    }
                }
            }
            if (y < sizeY - 1) {
                auto secondNode = this->getNode(x, y + 1);
                if (myNode->hasNeighbor(secondNode)) {
                    plots.push_back({myNode->pos, secondNode->pos});
                    if (myNode->value < 0 || secondNode->value < 0 || myNode->value != secondNode->value) {
                        plots_colors.push_back(Qt::gray); // plots_colors.push_back(PlotColor::GRAY);
                    } else {
                        plots_colors.push_back(QColor::fromHsvF((float)myNode->value / maxNodeValue, 1.f, 1.f)); // plots_colors.push_back((PlotColor)(PlotColor::RED + myNode->value));
                    }
                }
            }
        }
    }

    Plotter::get()->reset();
    for (size_t iPlot = 0; iPlot < plots.size(); iPlot++)
        Plotter::get()->addPlot(plots[iPlot], "", plots_colors[iPlot]);
//    plt.addScatter(scatter, "", scatter_labels, {PlotColor::GRAY});

    Plotter::get()->show();
}

int RegularSimplicialComplex::maxNodesValues()
{
    int max = -1;
    for (const auto& node : this->nodes)
        max = std::max(max, node->value);
    return max;
}

bool RegularSimplicialComplex::checkConnection(const Vector3& posA, const Vector3& posB)
{
    auto nodeA = this->getNode(posA);
    auto nodeB = this->getNode(posB);
    if (nodeA == nullptr || nodeB == nullptr)
        return false;
    else
        return this->checkConnection(nodeA, nodeB);
}

bool RegularSimplicialComplex::checkConnection(GraphNode* nodeA, GraphNode* nodeB)
{
    return nodeA->hasNeighbor(nodeB);
}

void RegularSimplicialComplex::removeUnavailableLinks()
{
    // Remove connections if two neighbors have different values
    for (size_t iNode = 0; iNode < this->nodes.size(); iNode++) {
        auto& node = this->nodes[iNode];
        if (node->value != -1) {
            for (int iNeighbor = node->neighbors.size() - 1; iNeighbor >= 0; iNeighbor--) {
                auto& [neighbor, weight] = node->neighbors[iNeighbor];
                if (neighbor->value != node->value) {
                    this->removeConnection(node, neighbor);
                }
            }
        }
    }

    // Remove connections if a node doesn't have 2 consecutive links available.
    // Because in this case, that means that no path can be drawn without englobing another.
    /// This condition is only available because we have a regular grid here.
    std::vector<Vector3> neighborOffsets = {
        Vector3(-1,  0, 0),
        Vector3( 0, -1, 0),
        Vector3( 1, -1, 0),
        Vector3( 1,  0, 0),
        Vector3( 0,  1, 0),
        Vector3(-1,  1, 0)
    };
    for (size_t iNode = 0; iNode < this->nodes.size(); iNode++) {
        auto& node = this->nodes[iNode];
        if (node->value == -1) {
            Vector3 pos = node->pos;
            for (size_t iOffset = 0; iOffset < neighborOffsets.size(); iOffset++) {
                // Does it have a neighbor here ?
                if (this->checkConnection(pos, pos + neighborOffsets[iOffset])) {
                    if (!this->checkConnection(pos, pos + neighborOffsets[(iOffset - 1 + neighborOffsets.size()) % neighborOffsets.size()]) &&
                            !this->checkConnection(pos, pos + neighborOffsets[(iOffset + 1 + neighborOffsets.size()) % neighborOffsets.size()]))
                        this->removeConnection(node, this->getNode(pos + neighborOffsets[iOffset]));
                }
            }
        }
    }
}
