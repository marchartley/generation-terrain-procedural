#include "KarstPathsGeneration.h"

#include "Graph/FastPoissonGraph.h"
#include "Graph/Pathfinding.h"

WaterHeight::WaterHeight(float height, float max_dist, float effect)
    : height(height), effect(effect), max_distance_of_effect(max_dist)
{
}

float WaterHeight::cost(const Vector3& node) const
{
    float dist = std::abs(node.z - this->height);
    if (dist < this->max_distance_of_effect) {
        return (dist / max_distance_of_effect) * this->effect;
    } else {
        return 1.f;
    }
}

PorositySphere::PorositySphere(const Vector3& pos, float radius, float effect)
    : pos(pos), radius(radius), effect(effect)
{
}

float PorositySphere::cost(const Vector3& node) const
{
    if ((node - this->pos).norm2() < this->radius*this->radius) {
        return (node - this->pos).norm() / this->radius;
    } else {
        return 0.f;
    }
}

FractureDirection::FractureDirection(const Vector3& direction, float effect)
    : direction(direction), effect(effect)
{
}

float FractureDirection::cost(const Vector3& nodeA, const Vector3& nodeB) const
{
    Vector3 dirBA = (nodeA - nodeB).normalized();
    float cost = this->effect * (1 - (std::pow(this->direction.normalized().dot(dirBA), 2)));
    if (cost < 0) {
        std::cout << this->effect << "\n" << this->direction.normalized().dot(dirBA) << "\n";
        std::cout << this->direction.normalized() << "\n" << dirBA << "\n";
        std::cout << nodeA << " " << nodeB << std::endl;
    }
    return cost;
}

KarstPathsGeneration::KarstPathsGeneration()
//    : KarstPathsGeneration(GridI(0, 0, 0), Vector3(0, 0, 0))
{

}

KarstPathsGeneration::KarstPathsGeneration(GridI availablityMap, const Vector3& karst_dimensions, float poissonDistance, float gamma)
    : gamma(gamma)
{
    this->karst_available_matrix = availablityMap.resize(karst_dimensions, RESIZE_MODE::MIN_VAL);
    this->porosityMap = GridF(karst_dimensions);
    this->graph = FastPoissonGraph<int>(this->karst_available_matrix, poissonDistance, 30);
    this->tortuosityOffsets = std::vector<Vector3>(this->graph.nodes.size());
}

void KarstPathsGeneration::resetSpecialNodes()
{
    this->specialNodes.clear();
}

void KarstPathsGeneration::setSpecialNode(int node, KARST_NODE_TYPE mode)
{
    this->specialNodes.push_back(std::make_tuple(node, mode));
}

float KarstPathsGeneration::computeCost(int nodeA, int nodeB)
{
    float distanceCost = 0.f;
    float fractureCost = 0.f;
    float waterCost = 0.f;
    float porosityCost = 0.f;

    if (distanceWeight > 0) {
        distanceCost = 1.f;
        distanceCost = (this->getNodePos(nodeA) - this->getNodePos(nodeB)).norm();
        distanceCost *= distanceWeight;
    }
    if (fractureWeight > 0) {
        fractureCost = 1.f;
        for (const FractureDirection& frac : this->fracturesDirections)
            fractureCost += frac.cost(this->getNodePos(nodeA), this->getNodePos(nodeB));
        fractureCost *= fractureWeight;
    }
    if (porosityWeight > 0) {
        porosityCost = 1.f;
        for (const PorositySphere& poro : this->porositySpheres)
            porosityCost += poro.cost(this->getNodePos(nodeA));
        this->porosityMap.raiseErrorOnBadCoord = false; // Errors can rise if we manually bring a node outside the map
        porosityCost += this->porosityMap.at(this->getNodePos(nodeA));
        this->porosityMap.raiseErrorOnBadCoord = true;
        porosityCost *= porosityWeight;
    }
    if (waterHeightWeight > 0) {
        waterCost = 1.f;
        for (const WaterHeight& waterHeight : this->waterHeights)
            waterCost += KarstPathsGeneration::smoothStep(waterHeight.cost(this->getNodePos(nodeA)));
        waterCost *= waterHeightWeight;
    }

//    float totalCost = (distanceCost * fractureCost * porosityCost * waterCost);
    float totalCost = (distanceCost + fractureCost + porosityCost + waterCost);
    return totalCost;
}

void KarstPathsGeneration::addWaterHeight(WaterHeight waterHeight)
{
    this->waterHeights.push_back(waterHeight);
}

void KarstPathsGeneration::addPorositySphere(PorositySphere porositySphere)
{
    this->porositySpheres.push_back(porositySphere);
}

void KarstPathsGeneration::addFractureDirection(FractureDirection fractureDirection)
{
    this->fracturesDirections.push_back(fractureDirection);
}

void KarstPathsGeneration::addPorosityMap(GridF porosity)
{
    this->porosityMap += porosity.resize(this->porosityMap.sizeX, this->porosityMap.sizeY, this->porosityMap.sizeZ, LINEAR);
}

void KarstPathsGeneration::createEdges(int maxNumberNeighbors, float maxNeighboringDistance)
{
    this->graph.createEdges(maxNumberNeighbors, maxNeighboringDistance, true); // The "true" means the cost are not precomputed, just the connection matrix
    for (int i = 0; i < this->graph.connectionMatrix.sizeX; i++) {
        for (int j = 0; j < this->graph.connectionMatrix.sizeY; j++) {
            if (this->graph.connectionMatrix.at(i, j)) {
                this->graph.adjencyMatrix.at(i, j) = std::pow(this->computeCost(i, j), this->gamma);
            }
        }
    }
}

void KarstPathsGeneration::computeAllPathsBetweenSpecialNodes(int uniqueNodeToRecompute)
{
    /*
    NORMAL = 0,
    FORCED_POINT = 1,
    DEAD_END = 2,
    SOURCE = 3,
    EXIT = 4*/

    std::vector<std::tuple<int, KARST_NODE_TYPE>> specialNodesToCompute = this->specialNodes;
    if (uniqueNodeToRecompute != -1)
        for (auto& n : this->specialNodes)
            if (std::get<0>(n) == uniqueNodeToRecompute)
                specialNodesToCompute = std::vector<std::tuple<int, KARST_NODE_TYPE>>({n});

    this->finalPaths.clear();
    this->allPaths = Matrix3<std::vector<int>>(this->specialNodes.size(), this->specialNodes.size());
    GridF pathDistances(this->specialNodes.size(), this->specialNodes.size(), 1, std::numeric_limits<float>::max());
    // Compute all paths
    std::vector<float> distances;
    std::vector<int> prec;
    for(size_t i = 0; i < specialNodesToCompute.size(); i++) {
        std::cout << 100 * i / (float)(specialNodesToCompute.size() - 1) << "%" << std::flush;
        if (std::get<1>(specialNodesToCompute[i]) == KARST_NODE_TYPE::DEAD_END)
            continue;
        int nodeA = std::get<0>(specialNodesToCompute[i]);
        auto start = std::chrono::system_clock::now();
        std::tie(distances, prec) = Pathfinding::ShortestPathFrom(nodeA, this->graph.adjencyMatrix);

        for(size_t j = 0; j < this->specialNodes.size(); j++) {
//            if (i == j) continue;
            int nodeB = std::get<0>(this->specialNodes[j]);
            std::vector<int> path = Pathfinding::getPath(nodeB, prec);
            if (this->allPaths.at(j, i).size() == 0 || distances[nodeB] < pathDistances.at(j, i)) // distances[j]/(float)path.size() < pathDistances.at(j, i)/(float)this->allPaths.at(j, i).size())
            {
                this->allPaths.at(i, j) = path;
                this->allPaths.at(j, i) = path; // Delete the opposite path if this one is shorter
                pathDistances.at(i, j) = distances[nodeB]; //(float)path.size();
                pathDistances.at(j, i) = distances[nodeB]; //(float)path.size();
            }
        }
        std::cout << " (" << std::chrono::duration<float>(std::chrono::system_clock::now() - start).count() << "s)" << std::endl;
    }

    float min = pathDistances.min(), max = pathDistances.max();
//    std::cout << min << " " << max << std::endl;
//    exit(0);
    pathDistances = (pathDistances - min) / (max - min);
    // Pruning the paths like in Paris, A., Peytavie, A., Guérin, E., Collon, P., & Galin, E. (2021). Synthesizing Geologically Coherent Cave Networks.
    // Meaning using a Relative Neighboring Graph by using path(i, k) + path(k, j) instead of path(i, j) if the distances at the power gamma is better
    this->finalPaths.clear();
    for (int i = 0; i < this->allPaths.sizeX; i++) { // Starting node
        for (int j = 0; j < this->allPaths.sizeY; j++) { // Ending node
            if (this->allPaths.at(i, j).empty()) continue;
            bool keepPath = true;
            if (true || this->gamma != 1.f)
            {
                float gammaDistanceIJ = std::pow(pathDistances.at(i, j), this->gamma);
                for (int k = 0; k < this->allPaths.sizeX; k++) { // Itermediate node
                    if (i == k || j == k) continue;
//                    if (this->allPaths.at(i, k).empty() || this->allPaths.at(k, j).empty()) continue;
                    float gammaDistanceIK = std::pow(pathDistances.at(i, k), this->gamma);
                    float gammaDistanceKJ = std::pow(pathDistances.at(k, j), this->gamma);
                    if (gammaDistanceIJ > gammaDistanceIK + gammaDistanceKJ)
                    {
                        keepPath = false;
                        break;
                    }
                }
            }
            if (keepPath)
            {
                finalPaths.push_back(this->allPaths.at(i, j));
            }
        }
    }
}

void KarstPathsGeneration::updateTortuosity(float distanceToDisplace, std::vector<int> ignoreForSpecificNodes)
{
    // If the offset vectors are already computed, juste use these
    if (this->tortuosityOffsets.size() != this->graph.nodes.size())
        this->tortuosityOffsets = std::vector<Vector3>(this->graph.nodes.size());

    for (size_t i = 0; i < this->tortuosityOffsets.size(); i++) {
        if (std::find(ignoreForSpecificNodes.begin(), ignoreForSpecificNodes.end(), int(i)) != ignoreForSpecificNodes.end())
            continue; // If we wnat to ignore this node, just continue

        // If the tortuosity offset has already been defined (not first call to this function), just use it
        if (this->tortuosityOffsets[i].norm() > 0) {
            this->tortuosityOffsets[i] = this->tortuosityOffsets[i].normalized() * distanceToDisplace;
        } else {
            this->tortuosityOffsets[i] = Vector3::random() * distanceToDisplace;
        }
        // If the point went outside because of this distance, find a new tortuosity vector
        while(!this->porosityMap.checkCoord(this->getNodePos(i))) {
            this->tortuosityOffsets[i] = Vector3::random() * distanceToDisplace;
        }
    }
}

Vector3 KarstPathsGeneration::getNodePos(int node)
{
    // Don't make it vary if it's the source node
    for (auto& n : this->specialNodes)
        if (std::get<0>(n) == node && std::get<1>(n) == KARST_NODE_TYPE::SOURCE)
            return this->graph.nodes[node]->pos;

    return this->graph.nodes[node]->pos + this->tortuosityOffsets[node];
}
