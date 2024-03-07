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
    std::vector<GraphNode*> nodes(points.size());
    for (int i = 0; i < nodeIds.size(); i++) {
        nodeIds[i] = i;
        nodes[i] = new GraphNode(i, points[i], i);
    }

    this->graph = Graph(false);
    graph.addNodes(nodes);

    for (int iNode = 0; iNode < nodes.size(); iNode++) {
        for (auto& neighbor : voronoi.neighbors[iNode]) {
            float dist = (nodes[iNode]->pos - nodes[neighbor]->pos).norm();
            graph.addConnection(nodes[iNode], nodes[neighbor], dist);
            graph.addConnection(nodes[neighbor], nodes[iNode], dist);
        }
    }

    return *this;
}
