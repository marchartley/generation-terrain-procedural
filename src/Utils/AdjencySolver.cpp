#include "AdjencySolver.h"
#include "Utils/Utils.h"

typedef int ID;

AdjencySolver::AdjencySolver()
{

}


void completeAdjenciesWithRandomNeighbors(std::vector<ID>& bestCombinationYet, std::vector<std::vector<ID>>& neighbors, std::vector<std::vector<ID>>& possibleIdPerNeighborId) {
    size_t numberOfPossibleBiomes = bestCombinationYet.size();
    size_t numberOfNodes = bestCombinationYet.size();

    // Working on having random values for unset biomes
    std::vector<ID> unsetBiomes;
    std::vector<ID> unusedIDs;
    std::vector<ID> usedIDs(numberOfPossibleBiomes, -1);
    for (size_t i = 0; i < numberOfNodes; i++) {
        if (bestCombinationYet[i] == -1) {
            unsetBiomes.push_back(i);
        } else {
            usedIDs[bestCombinationYet[i]] = 1;
        }
    }
    for (int i = 0; i < numberOfPossibleBiomes; i++)
        if (usedIDs[i] == -1)
            unusedIDs.push_back(i);
    /*
     for (size_t iNeighbor : neighbors[index]) {
        availableBiomePerNode[iNeighbor] = vectorIntersection(availableBiomePerNode[iNeighbor], convertSetToVector(neighborsPossibilities));
     }
            */
    for (ID unsetID : unsetBiomes) {
        std::vector<ID> unsetNeighbors = neighbors[unsetID];
        std::shuffle(unsetNeighbors.begin(), unsetNeighbors.end(), random_gen::random_generator);
        int newBiomeID = -1;
        while (newBiomeID == -1) {
            std::vector<ID> availableBiomes = unusedIDs;
                    /*(numberOfNodes);
            for (size_t i = 0; i < numberOfNodes; i++) availableBiomes[i] = i;*/

            for (size_t iNeighbor : unsetNeighbors) {
               std::vector<ID> nextAvailableBiomesList = vectorIntersection(availableBiomes, possibleIdPerNeighborId[iNeighbor]);
               if (nextAvailableBiomesList.empty()) {
                   newBiomeID = availableBiomes[(int)random_gen::generate(0, availableBiomes.size())];
                   break;
               }
               availableBiomes = nextAvailableBiomesList;
            }
            if (!availableBiomes.empty()) {
                newBiomeID = availableBiomes[(int)random_gen::generate(0, availableBiomes.size())];
                break;
            }
        }
        bestCombinationYet[unsetID] = newBiomeID;
        unusedIDs.erase(std::find(unusedIDs.begin(), unusedIDs.end(), newBiomeID));
    }
}

float computeEntropy(std::vector<ID>& availableClasses) {
    return availableClasses.size();
}

std::vector<ID> simplifiedWFC(std::vector<std::string> restrictedBiomesNames, std::vector<std::pair<ID, ID>>& possibleAdjency, std::vector<std::vector<int>> neighbors, std::ostream& out) {
    size_t numberOfNodes = restrictedBiomesNames.size();
    // We want to set classes, taking into account duplicates
    std::vector<ID> remainingClassesToSet = std::vector<ID>(restrictedBiomesNames.size());
    for (size_t i = 0; i < remainingClassesToSet.size(); i++) remainingClassesToSet[i] = i;
    // Start with all nodes completely free
    std::vector<std::vector<ID>> availableBiomePerNode(numberOfNodes, remainingClassesToSet);
    std::vector<ID> biomesAlreadySet(numberOfNodes, -1); // Empty string => not set

    while (true) { // This must be changed later!!
        float minEntropy = std::numeric_limits<float>::max();
        for (size_t i = 0; i < numberOfNodes; i++) {
            if (!availableBiomePerNode[i].empty())
                minEntropy = std::min(minEntropy, computeEntropy(availableBiomePerNode[i]));
        }
        if (minEntropy == std::numeric_limits<float>::max()) { // All biomes are already set I guess
            out << "No nodes to check anymore" << std::endl;
            break;
        }
        // Get all nodes with the lower entropy and select one randomly
        std::vector<size_t> candidatesIndices;
        for (size_t i = 0; i < numberOfNodes; i++) {
            if (availableBiomePerNode[i].size() == minEntropy)
                candidatesIndices.push_back(i);
        }
        std::shuffle(candidatesIndices.begin(), candidatesIndices.end(), random_gen::random_generator);
        ID index = candidatesIndices.front();

        // Increase the possibility that we fall on the same class than a neighbor, or one of the compatible ones
        for (size_t iNeighbor : neighbors[index]) {
            if (biomesAlreadySet[iNeighbor] != -1) {
                for (int _ = 0; _ < 100; _++) {
//                    availableBiomePerNode[index].push_back(biomesAlreadySet[iNeighbor]);
                }
            } else {
//                availableBiomePerNode[index] = vectorUnion(availableBiomePerNode[index], availableBiomePerNode[iNeighbor]);
            }
        }
        // Set this biome to a random class and remove all other possibilities
        ID selectedClass = availableBiomePerNode[index][(int)random_gen::generate(0, availableBiomePerNode[index].size())];
        biomesAlreadySet[index] = selectedClass;
        availableBiomePerNode[index].clear();

        // Remove the affected classID to all the nodes
        for (auto& possibilities : availableBiomePerNode) {
            auto found = std::find(possibilities.begin(), possibilities.end(), selectedClass);
            if (found != possibilities.end())
                possibilities.erase(found);
        }

        // Now list all the classes the neighbors can have
        std::set<ID> neighborsPossibilities;
        for (const auto& biomePair : possibleAdjency) {
            if (biomePair.first == selectedClass)
                neighborsPossibilities.insert(biomePair.second);
        }

        // Remove possibilities of each neighbor that don't match the new class
        for (size_t iNeighbor : neighbors[index]) {
            availableBiomePerNode[iNeighbor] = vectorIntersection(availableBiomePerNode[iNeighbor], convertSetToVector(neighborsPossibilities));
        }
    }

    // If at least one biome is unset, repeat the function
    bool allBiomesSet = true;
    int nbValidBiomes = 0;
    for (auto biomeClass : biomesAlreadySet) {
        nbValidBiomes++;
        if (biomeClass == -1) {
            allBiomesSet = false;
            nbValidBiomes --;
//                break;
        }
    }
    if (!allBiomesSet)
        out << "So close, we had " << nbValidBiomes << "/" << numberOfNodes << " valid nodes..." << std::endl;
    else
        out << "All nodes are valid" << std::endl;

    return biomesAlreadySet;
}

void computePossibleAndImpossibleNeighborings(std::vector<std::string> restrictedBiomesNames,
                                              std::vector<std::pair<std::string, std::string> > impossibleAdjency,
                                              std::vector<std::vector<ID>>& possibleIdPerNeighborId,
                                              std::vector<std::vector<ID>>& impossibleIdPerNeighborId,
                                              std::vector<std::pair<ID, ID>>& possibleAdjency)
{
    possibleAdjency.clear();
    possibleIdPerNeighborId.resize(restrictedBiomesNames.size());
    impossibleIdPerNeighborId.resize(restrictedBiomesNames.size());
    for (ID biome1 = 0; biome1 < restrictedBiomesNames.size(); biome1++) {
        for (ID biome2 = 0; biome2 < restrictedBiomesNames.size(); biome2++) {
//            if (biome1 == biome2) continue; // Don't allow to have 2 identical biomes instance
            std::pair<std::string, std::string> biomePair = std::make_pair(restrictedBiomesNames[biome1], restrictedBiomesNames[biome2]);
            if (!(std::find(impossibleAdjency.begin(), impossibleAdjency.end(), biomePair) != impossibleAdjency.end())) {
                possibleAdjency.push_back(std::make_pair(biome1, biome2));
                possibleIdPerNeighborId[biome1].push_back(biome2);
            } else {
                impossibleIdPerNeighborId[biome1].push_back(biome2);
            }
        }
    }
}
std::vector<int> AdjencySolver::solveGeomToTopo(Voronoi& diagram, std::vector<std::string> restrictedBiomesNames, std::vector<std::pair<std::string, std::string> > impossibleAdjency)
{
    std::ostream& out = std::cout;
    // List all the biome names
//    std::vector<std::string> allPossibleBiomesNames = removeDuplicatesFromVector(restrictedBiomesNames);
    std::vector<ID> biomesAlreadySet;

    size_t numberOfNodes = diagram.pointset.size();
    std::vector<std::pair<ID, ID>> possibleAdjency;
    std::vector<std::vector<ID>> possibleIdPerNeighborId;
    std::vector<std::vector<ID>> impossibleIdPerNeighborId;
    computePossibleAndImpossibleNeighborings(restrictedBiomesNames, impossibleAdjency, possibleIdPerNeighborId, impossibleIdPerNeighborId, possibleAdjency);
    // Get the adjency map resulting from the Voronoi diagram
    std::vector<std::vector<int>> neighbors = diagram.neighbors;
    if (numberOfNodes <= 1)
    {
        std::vector<int> returnValue = std::vector<int>(numberOfNodes);
        for (size_t i = 0; i < returnValue.size(); i++) returnValue[i] = i;
        return returnValue;
    }

    bool allBiomesSet = false;
    std::vector<ID> bestCombinationYet = std::vector<ID>(numberOfNodes, -1);
    size_t maxIterations = numberOfNodes * 2; // * numberOfNodes * 1;

    int previousNbValidBiomes;

    out << "Starting ordering" << std::endl;
    // Repeat the process a few times to find a viable solution
    for (size_t iteration = 0; iteration < maxIterations && !allBiomesSet; iteration++) {
        out << "Iteration #" << iteration+1 << "/" << maxIterations << "..." << std::endl;
        biomesAlreadySet = simplifiedWFC(restrictedBiomesNames, possibleAdjency, neighbors, out);

        allBiomesSet = true;
        int nbValidBiomes = 0;
        for (auto biomeClass : biomesAlreadySet) {
            nbValidBiomes++;
            if (biomeClass == -1) {
                allBiomesSet = false;
                nbValidBiomes --;
            }
        }
        previousNbValidBiomes = 0;
        for (auto biomeClass : bestCombinationYet) {
            previousNbValidBiomes += 1 * (biomeClass == -1);
        }
        if (previousNbValidBiomes >= nbValidBiomes)
            bestCombinationYet = biomesAlreadySet;
    }
    completeAdjenciesWithRandomNeighbors(bestCombinationYet, neighbors, possibleIdPerNeighborId);

    out << "Done. ";
    if (allBiomesSet)
        out << "All OK" << std::endl;
    else
        out << previousNbValidBiomes << " nodes have been set relatively randomly" << std::endl;
    // Transform the IDs to the real classname
    std::vector<std::string> finalClassNames(numberOfNodes);
    for (size_t i = 0; i < finalClassNames.size(); i++) {
        finalClassNames[i] = restrictedBiomesNames[bestCombinationYet[i]];
    }
//    return finalClassNames;
    return bestCombinationYet;
}

void AdjencySolver::solveTopoToGeom(std::vector<BiomeInstance> instancesToPlace, std::vector<std::pair<std::string, std::string> > impossibleAdjency, ShapeCurve availableSpace)
{
    size_t numberOfNodes = instancesToPlace.size();

    // For n < 2, no need for more computation
    if (numberOfNodes == 0) {
        return;
    } else if (numberOfNodes == 1) {
        instancesToPlace[0].position = availableSpace.randomPointsInside()[0];
        instancesToPlace[0].area = availableSpace;
        return;
    }

    std::vector<std::string> restrictedBiomesNames;
    for (auto& biome : instancesToPlace)
        restrictedBiomesNames.push_back(biome.classname);

    // Create a fake Voronoi diagram
    Voronoi fakeDiagram = Voronoi(numberOfNodes, Vector3(1000, 1000)); // Use huge space, just to be sure that there are no floating decimal errors
    fakeDiagram.solve(false, 1);
    // Use the "geom. -> topo" solver to get a partial adjency graph
    std::ostringstream out;
    std::vector<ID> biomesAlreadySet;

    std::vector<std::pair<ID, ID>> possibleAdjency;
    std::vector<std::vector<ID>> possibleIdPerNeighborId;
    std::vector<std::vector<ID>> impossibleIdPerNeighborId;
    computePossibleAndImpossibleNeighborings(restrictedBiomesNames, impossibleAdjency, possibleIdPerNeighborId, impossibleIdPerNeighborId, possibleAdjency);

    // Create the (empty) topology graph
    Graph<ID> adjencyGraph;
    for (size_t i = 0; i < numberOfNodes; i++)
        adjencyGraph.nodes.push_back(std::make_shared<GraphNode<ID>>(i));

    // Get the adjency map resulting from the Voronoi diagram
    std::vector<std::vector<int>> neighbors = fakeDiagram.neighbors;

    bool allBiomesSet = false;
    std::vector<ID> bestCombinationYet = std::vector<ID>(numberOfNodes, -1);
    size_t maxIterations = numberOfNodes;

    int previousNbValidBiomes;

    // Repeat the process a few times to find a viable solution
    for (size_t iteration = 0; iteration < maxIterations && !allBiomesSet; iteration++) {
        biomesAlreadySet = simplifiedWFC(restrictedBiomesNames, possibleAdjency, neighbors, out);

        allBiomesSet = true;
        int nbValidBiomes = 0;
        for (auto biomeClass : biomesAlreadySet) {
            nbValidBiomes++;
            if (biomeClass == -1) {
                allBiomesSet = false;
                nbValidBiomes --;
            }
        }
        previousNbValidBiomes = 0;
        for (auto biomeClass : bestCombinationYet) {
            previousNbValidBiomes += 1 * (biomeClass == -1);
        }
        if (previousNbValidBiomes > nbValidBiomes)
            bestCombinationYet = biomesAlreadySet;
    }

    for (size_t iArea = 0; iArea < numberOfNodes; iArea++) {
        ID biomeID = bestCombinationYet[iArea];
        std::vector<int> biomeNeighbors = neighbors[iArea];
        for (int neighbor : biomeNeighbors) {
            ID neighborBiomeID = bestCombinationYet[neighbor];
            if (neighborBiomeID == -1) continue;
            adjencyGraph.nodes[biomeID]->addNeighbor(adjencyGraph.nodes[neighborBiomeID], std::max(0.f, instancesToPlace[biomeID].idealSize) + std::max(0.f, instancesToPlace[neighborBiomeID].idealSize));
        }
    }

    // Complete the missing biomes by adding edges randomly
    std::vector<ID> unsetBiomes;
    std::vector<ID> unusedIDs;
    std::vector<ID> usedIDs(numberOfNodes, -1);
    for (size_t i = 0; i < numberOfNodes; i++) {
        if (bestCombinationYet[i] == -1) {
            unsetBiomes.push_back(i);
        } else {
            usedIDs[bestCombinationYet[i]] = 1;
        }
    }
    for (size_t i = 0; i < numberOfNodes; i++)
        if (usedIDs[i] == -1)
            unusedIDs.push_back(i);

    for (ID unsetID : unsetBiomes) {
        std::vector<ID> newPossibleNeighbors = possibleIdPerNeighborId[unsetID];
        std::shuffle(newPossibleNeighbors.begin(), newPossibleNeighbors.end(), random_gen::random_generator);
        ID neighborID = newPossibleNeighbors[0];
        adjencyGraph.nodes[unsetID]->addNeighbor(adjencyGraph.nodes[neighborID], std::max(0.f, instancesToPlace[unsetID].idealSize) + std::max(0.f, instancesToPlace[neighborID].idealSize));
        /*
        std::vector<ID> unsetNeighbors = neighbors[unsetID];
        std::shuffle(unsetNeighbors.begin(), unsetNeighbors.end(), random_gen::random_generator);
        int newBiomeID = -1;
        while (newBiomeID == -1) {
            std::vector<ID> availableBiomes = unusedIDs;

            for (size_t iNeighbor : unsetNeighbors) {
               std::vector<ID> nextAvailableBiomesList = vectorIntersection(availableBiomes, possibleIdPerNeighborId[iNeighbor]);
               if (nextAvailableBiomesList.empty()) {
                   newBiomeID = availableBiomes[(int)random_gen::generate(0, availableBiomes.size())];
                   break;
               }
               availableBiomes = nextAvailableBiomesList;
            }
            if (!availableBiomes.empty()) {
                newBiomeID = availableBiomes[(int)random_gen::generate(0, availableBiomes.size())];
                break;
            }
        }
        bestCombinationYet[unsetID] = newBiomeID;
        unusedIDs.erase(std::find(unusedIDs.begin(), unusedIDs.end(), newBiomeID));*/
    }
    return;

    // Generate a new geometry from it
    //
}
