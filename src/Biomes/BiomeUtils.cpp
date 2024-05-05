#include "BiomeUtils.h"

int biomeID(std::string biomeName) {
    auto find = std::find(biomesNames.begin(), biomesNames.end(), toLower(biomeName));
    if (find != biomesNames.end())
        return std::distance(biomesNames.begin(), find);
    return -1;
}

std::string biomeName(int ID) {
    return biomesNames[ID];
}

FastPoissonGraph<int> generateHugeBiomesGraphe(std::vector<int> desiredBiomes, Graph adjencyGraph) {
    GridI available(100, 100, 1, 1); // All available
    int maxNeighbors = 6;
    int maxNeighboringTries = 20;
    FastPoissonGraph<int> graph(available);

    if (graph.empty()) return graph;
    graph.initAllNodesValues(-1);

    std::set<GraphNode*> pendingNodes;
//    pendingNodes.insert(graph.nodes.begin());
    pendingNodes.insert(graph.nodes.front());
    while (!pendingNodes.empty()) {
        auto& current = *pendingNodes.begin();
        pendingNodes.erase(pendingNodes.begin());

        std::set<GraphNode*> nodesToAddToPending;

        // If it's an unknown node type, give it a type
        if (current->value == -1) {
            current->value = desiredBiomes[int(random_gen::generate(0, desiredBiomes.size()))];
        }

        // If it isn't saturated in neighbors yet, add some randomly
        if (int(current->neighbors.size()) < maxNeighbors)
        {
            for (int i = 0; i < maxNeighboringTries && int(current->neighbors.size()) < maxNeighbors; i++) {
                int neighborIndex = int(random_gen::generate(0, graph.nodes.size()));
                auto& neighbor = graph.nodes[neighborIndex];

                if (neighbor->value != -1) {
                    // Don't use this node if it's neighboring count is already saturated
                    if (int(neighbor->neighbors.size()) >= maxNeighbors)
                        continue;
                    // Don't use it if the biomes are not compatible
                    if (adjencyGraph.connectionMatrix.at(current->value, neighbor->value) == 0)
                        continue;
                } else {
                    // Find all possible biomes ID that can take the neighbor
                    std::vector<int> possibleBiomes;
                    for (int iBiome = 0; iBiome < adjencyGraph.connectionMatrix.sizeX; iBiome++) {
                        if (adjencyGraph.connectionMatrix.at(current->value, iBiome) != 0)
                            possibleBiomes.push_back(iBiome);
                    }
                    // If the biome have no possible neighbor, it is an isolated node in the adjency graph
                    if (possibleBiomes.empty())
                        break;

                    neighbor->value = possibleBiomes[int(random_gen::generate(0, possibleBiomes.size()))];
                }
                /// TODO : here it is possible to specify a distance value, might be useful to add
                /// a value depending on the desired biome size and apply a solver to position everything
                current->addNeighbor(neighbor, 1.0);
                neighbor->addNeighbor(current, 1.0);
                nodesToAddToPending.insert(neighbor);
            }
        }
        pendingNodes.insert(nodesToAddToPending.begin(), nodesToAddToPending.end());
    }
    graph.forceDrivenPositioning();
    Vector3 minPosition = Vector3::max();
    for (auto& n : graph.nodes) {
        minPosition.x = std::min(minPosition.x, n->pos.x);
        minPosition.y = std::min(minPosition.y, n->pos.y);
        minPosition.z = std::min(minPosition.z, n->pos.z);
    }
    for (auto& n : graph.nodes)
        n->pos -= minPosition;

    return graph;
}

Graph subsetToFitMostBiomes(Graph graph, std::vector<std::string> biomesNames)
{
    Graph returnedGraph;

    ShapeCurve biomesHull;
    for (size_t iNode = 0; iNode < graph.nodes.size(); iNode++) {
        biomesHull.points.push_back(graph.nodes[iNode]->pos);
    }
    biomesHull = biomesHull.computeConvexHull();
    Vector3 AABBoxMin, AABBoxMax;
    std::tie(AABBoxMin, AABBoxMax) = biomesHull.AABBox();
    Vector3 desiredAABBoxMin = AABBoxMin + (AABBoxMax - AABBoxMin) * .3f;
    Vector3 desiredAABBoxMax = AABBoxMin + (AABBoxMax - AABBoxMin) * .6f;

    for (size_t iNode = 0; iNode < graph.nodes.size(); iNode++) {
        Vector3 pos = graph.nodes[iNode]->pos;
        if (desiredAABBoxMin.x < pos.x && pos.x < desiredAABBoxMax.x && desiredAABBoxMin.y < pos.y && pos.y < desiredAABBoxMax.y)
//            returnedGraph.nodes.push_back(graph.nodes[iNode]);
            returnedGraph.nodes[iNode] = graph.nodes[iNode];
    }
    return returnedGraph;
}

std::shared_ptr<BiomeInstance> recursivelyCreateBiomeInstance(nlohmann::json json_content, const Vector3& biomePosition, ShapeCurve area) {
    std::string biomeClass = json_content.at("class").get<std::string>();
    // Should be able to retrieve the parameters of the biome...
    std::shared_ptr<BiomeInstance> instance = std::make_shared<BiomeInstance>(BiomeInstance::fromClass(biomeClass));
    instance->position = biomePosition;
    instance->area = area;
    auto children = json_content.at("children");
    Voronoi diagram(children.size(), area);
    std::vector<ShapeCurve> subarea_borders = diagram.solve();
    for (size_t i = 0; i < children.size(); i++) {
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiomeInstance(children[i], diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}

std::shared_ptr<BiomeInstance> generateBiome(Graph biomeGraph)
{

}
