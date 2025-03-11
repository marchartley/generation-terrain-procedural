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

std::vector<std::vector<Vector3> > Delaunay::getTriangles() const
{
    std::set<std::set<GraphNode*>> triangles;
    std::map<GraphNode*, std::set<GraphNode*>> neighborhoods;
    for (int iNode = 0; iNode < graph.size(); iNode++) {
        auto& n = graph.nodes[iNode];
        for (auto& [nn, dist] : n->neighbors) {
            //if (nn->index < n->index)
                neighborhoods[n].insert(nn);
        }
    }

    for (int iNode = 0; iNode < graph.size(); iNode++) {
        auto& n = graph.nodes[iNode];
        for (auto& nn : neighborhoods[n]) {
            for (auto& nnn : neighborhoods[nn]) {
                if (neighborhoods[n].count(nnn)) {
                    triangles.insert({n, nn, nnn});
                }
            }
        }
    }

    std::vector<std::vector<Vector3>> spatialTriangles;
    for (auto& triangle : triangles) {
        std::vector<Vector3> points;
        for (auto& n : triangle) {
            points.push_back(n->pos);
        }
        spatialTriangles.push_back(points);
    }
    return spatialTriangles;
}
