#include "AdjencySolver.h"
#include "Utils/Utils.h"

AdjencySolver::AdjencySolver()
{

}


std::vector<int> AdjencySolver::solve(Voronoi& diagram, std::vector<std::string> restrictedBiomesNames, std::vector<std::pair<std::string, std::string> > impossibleAdgency)
{
    std::ostringstream out;
    typedef int ID;
    // List all the biome names
//    std::vector<std::string> allPossibleBiomesNames = removeDuplicatesFromVector(restrictedBiomesNames);
    std::vector<ID> biomesAlreadySet;

    size_t numberOfNodes = diagram.pointset.size();
    size_t numberOfPossibleBiomes = restrictedBiomesNames.size();
    std::vector<std::pair<ID, ID>> possibleAdgency;
    std::vector<std::vector<ID>> possibleIdPerNeighborId(numberOfPossibleBiomes);
    std::vector<std::vector<ID>> impossibleIdPerNeighborId(numberOfPossibleBiomes);
    for (ID biome1 = 0; biome1 < restrictedBiomesNames.size(); biome1++) {
        for (ID biome2 = 0; biome2 < restrictedBiomesNames.size(); biome2++) {
//            if (biome1 == biome2) continue; // Don't allow to have 2 identical biomes instance
            std::pair<std::string, std::string> biomePair = std::make_pair(restrictedBiomesNames[biome1], restrictedBiomesNames[biome2]);
            if (!(std::find(impossibleAdgency.begin(), impossibleAdgency.end(), biomePair) != impossibleAdgency.end())) {
                possibleAdgency.push_back(std::make_pair(biome1, biome2));
                possibleIdPerNeighborId[biome1].push_back(biome2);
            } else {
                impossibleIdPerNeighborId[biome1].push_back(biome2);
            }
        }
    }
    // Get the adgency map resulting from the Voronoi diagram
    std::vector<std::vector<int>> neighbors = diagram.neighbors;
    if (numberOfNodes <= 1)
    {
        std::vector<int> returnValue = std::vector<int>(numberOfNodes);
        for (size_t i = 0; i < returnValue.size(); i++) returnValue[i] = i;
        return returnValue;
    }

    bool allBiomesSet = false;
    std::vector<ID> bestCombinationYet = std::vector<ID>(numberOfNodes, -1);
    size_t maxIterations = numberOfNodes * numberOfNodes;

    out << "Starting ordering" << std::endl;
    // Repeat the process a few times to find a viable solution
    for (size_t iteration = 0; iteration < maxIterations && !allBiomesSet; iteration++) {
        out << "Iteration #" << iteration+1 << "/" << maxIterations << "..." << std::endl;
        // We want to set classes, taking into account duplicates
        std::vector<ID> remainingClassesToSet = std::vector<ID>(restrictedBiomesNames.size());
        for (size_t i = 0; i < remainingClassesToSet.size(); i++) remainingClassesToSet[i] = i;
        // Start with all nodes completely free
        std::vector<std::vector<ID>> availableBiomePerNode(numberOfNodes, remainingClassesToSet);
        biomesAlreadySet = std::vector<ID>(numberOfNodes, -1); // Empty string => not set

        while (true) { // This must be changed later!!
            size_t minEntropy = std::numeric_limits<size_t>::max();
            for (size_t i = 0; i < numberOfNodes; i++) {
                if (!availableBiomePerNode[i].empty())
                    minEntropy = std::min(minEntropy, availableBiomePerNode[i].size());
            }
            if (minEntropy == std::numeric_limits<size_t>::max()) { // All biomes are already set I guess
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
            for (const auto& biomePair : possibleAdgency) {
                if (biomePair.first == selectedClass)
                    neighborsPossibilities.insert(biomePair.second);
            }

            // Remove possibilities of each neighbor that don't match the new class
            for (size_t iNeighbor : neighbors[index]) {
                availableBiomePerNode[iNeighbor] = vectorIntersection(availableBiomePerNode[iNeighbor], convertSetToVector(neighborsPossibilities));
            }
        }

        // If at least one biome is unset, repeat the function
        allBiomesSet = true;
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

        int previousNbValidBiomes = 0;
        for (auto biomeClass : bestCombinationYet) {
            previousNbValidBiomes += 1 * (biomeClass == -1);
        }
        if (previousNbValidBiomes > nbValidBiomes)
            bestCombinationYet = biomesAlreadySet;
    }
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

    out << "Done. ";
    if (allBiomesSet)
        out << "All OK" << std::endl;
    else
        out << unusedIDs.size() << " nodes have been set relatively randomly" << std::endl;
    // Transform the IDs to the real classname
    std::vector<std::string> finalClassNames(numberOfNodes);
    for (size_t i = 0; i < finalClassNames.size(); i++) {
        finalClassNames[i] = restrictedBiomesNames[bestCombinationYet[i]];
    }
//    return finalClassNames;
    return bestCombinationYet;
}
