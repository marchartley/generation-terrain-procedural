#include "Delaunay.h"

Delaunay::Delaunay()
{

}

Delaunay::Delaunay(const Voronoi& voronoi)
{
    this->fromVoronoi(voronoi);
}

Delaunay& Delaunay::fromVoronoi(const Voronoi& voronoi)
{
    auto points = voronoi.pointset;
    std::vector<int> nodeIds(points.size());
    std::vector<std::shared_ptr<GraphNodeTemplate<int>>> nodes(points.size());
    for (int i = 0; i < nodeIds.size(); i++) {
        nodeIds[i] = i;
        nodes[i] = std::make_shared<GraphNodeTemplate<int>>(i, points[i], i);
    }

    this->graph = Graph();
    graph.addNodes(nodes);

    for (int iNode = 0; iNode < nodes.size(); iNode++) {
        for (auto& neighbor : voronoi.neighbors[iNode]) {
            graph.addConnection(nodes[iNode], nodes[neighbor]);
        }
    }

    return *this;
}
