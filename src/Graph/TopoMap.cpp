#include "TopoMap.h"
#include "Utils/Utils.h"

int Brin::currentBrinIndex = -1;

Brin::Brin() : Brin(nullptr)
{}

Brin::Brin(Brin *beta2) : Brin(beta2, nullptr)
{}

Brin::Brin(Brin *beta2, Brin *beta1)
    : beta1(beta1), beta2(beta2)
{
    Brin::currentBrinIndex ++;
    this->index = Brin::currentBrinIndex;
}

Brin *Brin::previous()
{
    return this->getFace().back();
}

std::vector<Brin *> Brin::getFace()
{
    std::vector<Brin*> listOfBrins({this});
    Brin* current = this->beta1;
    while (current != this) { // && std::find(listOfBrins.begin(), listOfBrins.end(), current) == listOfBrins.end()) {
        listOfBrins.push_back(current);
//        if (current == current->beta1)
//            return listOfBrins;
        current = current->beta1;
    }
    return listOfBrins;
}

std::vector<Brin *> Brin::getOppositeFace()
{
    return this->beta2->getFace();
}

std::pair<std::vector<Brin *>, std::vector<Brin *> > Brin::getTwoFaces()
{
    return {this->getFace(), this->getOppositeFace()};
}

std::vector<std::vector<Brin *> > Brin::getAllFaces()
{
    std::vector<std::vector<Brin*>> faces;
    auto unseenBrins = this->getAllBrins();
    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();

        auto face = current->getFace();
        for (auto& brin : face) {
            if (isIn(brin, unseenBrins)) {
                unseenBrins.erase(std::find(unseenBrins.begin(), unseenBrins.end(), brin));
            }
        }
        if (face.size() > 1)
            faces.push_back(face);
    }
    return faces;
}

std::vector<Brin *> Brin::getAllBrins()
{
    std::vector<Brin*> seenBrins;
    std::vector<Brin*> unseenBrins({this});

    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();
        if (isIn(current, seenBrins))
            continue;
        seenBrins.push_back(current);
        if (!isIn(current->beta2, seenBrins))
            unseenBrins.push_back(current->beta2);
        if (!isIn(current->beta1, seenBrins))
            unseenBrins.push_back(current->beta1);
    }
    return seenBrins;
}

std::pair<Brin *, Brin *> Brin::subdivide()
{
    Brin* nextBrin = this->beta1;
    Brin* twinBrin = this->beta2;
    Brin* prevTwin = twinBrin->previous();

    this->beta1 = new Brin(nullptr, nextBrin);
    prevTwin->beta1 = new Brin(this->beta1, twinBrin);

    this->beta1->beta2 = prevTwin->beta1;

    return {this, this->beta1};
}

std::vector<Brin *> Brin::orbit()
{
//    if (this->beta1 == this->beta1->beta1)
//        return {this};
    std::vector<Brin*> brins({this});
    Brin* current = this->beta2->beta1;
    while (current != this) {
        brins.push_back(current);
//        if (!current->isValid() || !current->beta2->isValid()) {
//            return {};
//        }
        current = current->beta2->beta1;
    }
    return brins;
}

std::vector<Brin *> Brin::inversedOrbit()
{
    return this->beta2->orbit();
}

void Brin::affectSource(int source)
{
    for (auto& brin : this->orbit()) {
        brin->beta2->affectedDest = source;
        brin->affectedSource = source;
    }
}

int Brin::getSource()
{
    if (this->affectedSource <= -1) {
        // Source unknown, check the neighbors
        for (auto& brin : this->orbit()) {
            if (brin->affectedSource > -1)
                this->affectedSource = brin->affectedSource;
            else if (brin->beta2->affectedDest > -1)
                this->affectedSource = brin->beta2->affectedDest;
        }
    }
    return this->affectedSource;
}

void Brin::affectDest(int dest)
{
    for (auto& brin : this->inversedOrbit()) {
        brin->beta2->affectedDest = dest;
        brin->affectedSource = dest;
    }
}

int Brin::getDest()
{
    if (this->affectedDest <= -1) {
        // Source unknown, check the neighbors
        for (auto& brin : this->inversedOrbit()) {
            if (brin->beta2->affectedDest > -1)
                this->affectedDest = brin->beta2->affectedDest;
            else if (brin->affectedSource > -1)
                this->affectedDest = brin->affectedSource;
        }
    }
    return this->affectedDest;
}

void Brin::affectSourceAndDest(int source, int dest)
{
    this->affectSource(source);
    this->affectDest(dest);
}

bool Brin::isValid()
{
    return this->beta1 != nullptr && this->beta2 != nullptr;
}

int Brin::degree()
{
    return this->inversedOrbit().size();
}

std::string Brin::toString()
{
    bool valid = this->isValid();
    return "Brin #" + std::to_string(this->index) + ", "
            "twin = " + (this->beta2 ? std::to_string(this->beta2->index) : "none") + " "
            "next = " + (this->beta1 ? std::to_string(this->beta1->index) : "none") + " "
            "(" + (valid ? std::to_string(this->getSource()) : "--") + " -> " + (valid ? std::to_string(this->getDest()) : "--") + ")";
}

CombinMap::CombinMap(Brin* root) : root(root)
{

}

std::vector<Brin *> CombinMap::allBrins()
{
    if (this->root != nullptr)
        return this->root->getAllBrins();
    return {};
}

std::vector<Brin *> CombinMap::addFace(int numberOfEdges, std::vector<Brin *> adjacentBrins, std::vector<int> newIndices)
{
    if (numberOfEdges < 3)
        throw std::invalid_argument("We need at least 3 edges to create a new face. (Only " + std::to_string(numberOfEdges) + " edges asked)");

    int nbEdgesToCreate = numberOfEdges - adjacentBrins.size();
    if (newIndices.empty())
        newIndices = std::vector<int>(nbEdgesToCreate - (adjacentBrins.empty() ? 0 : 1), -1); // Fill nodes indices with "-1"
    else if (int(newIndices.size()) != nbEdgesToCreate - (adjacentBrins.empty() ? 0 : 1))
        throw std::invalid_argument("We need to assign 0 or " + std::to_string(nbEdgesToCreate - (adjacentBrins.empty() ? 0 : 1)) + " nodes ID. (" + std::to_string(newIndices.size()) + " has been given)");
    std::vector<Brin*> newBrins(2 * nbEdgesToCreate);
    for (size_t i = 0; i < newBrins.size(); i++) {
        newBrins[i] = new Brin();
    }

    Brin* entryBrin;
    Brin* exitBrin;
    if (this->root == nullptr) {
        entryBrin = newBrins.front();
        exitBrin = newBrins.back();
        this->root = newBrins[1]; // Brin#1 is on the outside, so it's not so bad... (arbitrary choice)
    } else if (!adjacentBrins.empty()) {
        entryBrin = adjacentBrins.front()->previous();
        exitBrin = adjacentBrins.back()->beta1;
    } else {
        throw std::invalid_argument("Either you are creating a new graph, either you define at least one edge to paste the "
                                    "new face to. ");
    }

    // Twins and successions
    for (size_t i = 0; i < newBrins.size() / 2; i++) {
        newBrins[2 * i]->beta2 = newBrins[2 * i + 1];
        newBrins[2 * i + 1]->beta2 = newBrins[2 * i];

        if (i == 0) {
            // Special case for the first pair of new brins
            newBrins[1]->beta1 = exitBrin;
            newBrins[0]->beta1 = (newBrins.size() > 2 ? newBrins[2] : adjacentBrins.front());
        } else if (i == newBrins.size()/2 - 1) {
            // Special case for the last pair of new brins
            newBrins[2 * i]->beta1 = (adjacentBrins.empty() ? entryBrin : adjacentBrins.front());
            newBrins[2 * i + 1]->beta1 = newBrins[2 * i - 1];
        } else {
            // Normal case
            newBrins[2 * i]->beta1 = newBrins[2 * (i+1)];
            newBrins[2 * i + 1]->beta1 = newBrins[2 * i - 1];
        }
    }
    for (size_t i = 0; i < newIndices.size(); i++) {
        // Assign nodes indices
        newBrins[2 * i]->affectDest(newIndices[i]);
    }

    if (!adjacentBrins.empty()) {
        adjacentBrins.back()->beta1 = newBrins.front();
        entryBrin->beta1 = newBrins.back();
    }
    return newBrins;
}

std::vector<std::vector<Brin *> > CombinMap::allFaces(bool sortBySize)
{
    if (this->root != nullptr) {
        auto faces = this->root->getAllFaces();
        if (sortBySize) {
            std::sort(faces.begin(), faces.end(),
                      [](const auto& a, const auto& b) -> bool {
                return a.size() < b.size();
            });
        }
        return faces;
    }
    return {};
}

std::vector<std::vector<Brin *> > CombinMap::allOrbits()
{
    std::vector<std::vector<Brin*>> orbits;
    std::vector<Brin*> unseenBrins({this->root});
    std::vector<Brin*> seenBrins;

    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();
        if (isIn(current, seenBrins))
            continue;
        auto orbit = current->inversedOrbit();
        orbits.push_back(orbit);

        for (auto& brin : orbit) {
            seenBrins.push_back(brin->beta2);
            if (!isIn(brin, seenBrins))
                unseenBrins.push_back(brin);
        }
    }
    return orbits;
}

std::vector<std::vector<Brin *> > CombinMap::allIngoingOrbits()
{
    auto allOutgoingOrbits = this->allOutgoingOrbits();
    for (auto& orbit : allOutgoingOrbits)
        for(auto& brin : orbit)
            brin = brin->beta2; // Just inverse all the results
    return allOutgoingOrbits;
}

std::vector<std::vector<Brin *> > CombinMap::allOutgoingOrbits()
{
    return this->allOrbits();
}

std::vector<Brin *> CombinMap::exteriorFace(Brin *oneExteriorBrin)
{
    /// Needs to be improved... not a viable solution
    if (oneExteriorBrin != nullptr)
        return oneExteriorBrin->getFace();
    auto faces = this->allFaces(true);
    if (faces.empty())
        return {};
    return faces.back();
}

std::vector<std::vector<Brin *> > CombinMap::getHoles(int minHoleSize)
{
    auto faces = this->allFaces(true);
    while (!faces.empty() && faces.front().size() < minHoleSize)
        faces.erase(faces.begin());
    return faces;
}

void CombinMap::triangulate(Brin *oneExteriorBrin)
{
    auto faces = this->getHoles();
    auto exterior = this->exteriorFace(oneExteriorBrin);

    for (auto& face : faces) {
        if (isIn(face[0], exterior))
            continue; // Don't triangulate the exterior (?)
        Brin* addedEdge = face[0];
        for (size_t i = 2; i < face.size() - 1; i++) {
            auto newEdges = this->addFace(3, {addedEdge, face[i - 1]});
            addedEdge = (isIn(newEdges[0], addedEdge->inversedOrbit()) ? newEdges[0] : newEdges[1]);
        }
    }
}

std::vector<std::pair<int, int> > CombinMap::getUnorientedEdges()
{
    this->affectNodes();

    std::vector<std::pair<int, int>> unorientedEdges;
    auto unseenBrins = this->allBrins();
    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();
        // No need to check existance, the beta2 MUST be in the list
        unseenBrins.erase(std::find(unseenBrins.begin(), unseenBrins.end(), current->beta2));
        auto source = current->getSource();
        auto dest = current->getDest();
        auto pairA = std::make_pair(source, dest);
        unorientedEdges.push_back(pairA);
    }
    return unorientedEdges;
}
/*
std::vector<std::pair<int, int> > CombinMap::getCorrectionUnorientedEdges()
{
    // Add a check that there is a unique edge between node pairs
    this->affectNodes();

    std::vector<std::pair<int, int>> unorientedEdges;
    auto unseenBrins = this->allBrins();
    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();
        // No need to check existance, the beta2 MUST be in the list
        unseenBrins.erase(std::find(unseenBrins.begin(), unseenBrins.end(), current->beta2));
        auto source = current->getSource();
        auto dest = current->getDest();
        auto pairA = std::make_pair(source, dest);
        auto pairB = std::make_pair(dest, source);
        if (!isIn(pairA, unorientedEdges) && !isIn(pairB, unorientedEdges))
            unorientedEdges.push_back(pairA);
    }
    return unorientedEdges;
}
*/
std::vector<int> CombinMap::affectNodes()
{
    auto allEdges = this->allBrins();
    int nodeIndex = 0;
    for (auto& brin : allEdges)
        nodeIndex = std::max({nodeIndex, brin->getSource() + 1, brin->getDest() + 1});

    for (auto& brin : allEdges) {
        if (brin->getSource() == -1) {
            brin->affectSource(nodeIndex);
            nodeIndex++;
        }
    }
    std::vector<int> allNodesIndices(nodeIndex);
    for (size_t i = 0; i < allNodesIndices.size(); i++)
        allNodesIndices[i] = i;
    return allNodesIndices;
}

void CombinMap::collapseDuplicateNodes()
{
    bool needCollapsing = true;
    while (needCollapsing) {
        needCollapsing = false;
        auto faces = this->allFaces();
        for (auto& face : faces) {
            for (auto& brin : face) {
                if (brin->getSource() == brin->getDest()) {
                    needCollapsing = true;
                    this->collapse(brin);
                    break;
                } else if (brin->getSource() == brin->beta1->getDest()) {
                    needCollapsing = true;
                    auto newEdges = this->addFace(3, {brin, brin->beta1});
                    this->collapse(newEdges[0]);
                    break;
                }
            }
            if (needCollapsing)
                break;
        }
    }
}

CombinMap CombinMap::dual(Brin *oneExteriorBrin, Brin *oneEdgeTowardInfinity)
{
    std::vector<Brin*> myExterior;
    int infinityNodeIndex = -1;
    if (oneEdgeTowardInfinity == nullptr) {
        myExterior = this->exteriorFace(oneExteriorBrin);
        int infinityNodeIndex = 0; // Will contains the Index to the "Infinite node"
        for (auto& brin : this->allBrins())
            infinityNodeIndex = std::max({infinityNodeIndex, brin->getSource() + 1, brin->getDest() + 1});

        auto newEdges = this->addFace(3, {myExterior[0]}, {infinityNodeIndex}); // Create a triangle from the exterior
        for (size_t i = myExterior.size() - 1; i >= 2; i--) {
            // Iterate around the exterior, in counter-clockwise to use the newly created edges
            newEdges = this->addFace(3, {myExterior[i], myExterior[i]->beta1});
        }

    } else {
        // Care about the dual of dual later (oneEdgeTowardInfinity != nullptr)
        infinityNodeIndex = oneEdgeTowardInfinity->getDest();
    }

    this->debug();

    // Create new topo
    CombinMap dual;

    auto myBrins = this->allBrins();
    auto myFaces = this->allFaces();
    auto myOrbits = this->allIngoingOrbits();

    std::map<Brin*, int> edgeToFace;
    // Can be accelerated by removing the found brins...
    for (size_t iFace = 0; iFace < myFaces.size(); iFace++) {
        for (size_t iBrin = 0; iBrin < myBrins.size(); iBrin++) {
            if (isIn(myBrins[iBrin], myFaces[iFace])) {
                edgeToFace[myBrins[iBrin]] = iFace;
                std::cout << "Brin " << myBrins[iBrin]->index << " -> Face " << iFace << std::endl;
            }
        }
    }
//    std::vector<CombinMap> maps;
    std::vector<std::vector<int>> newFaces;
    for (size_t iOrbit = 0; iOrbit < myOrbits.size(); iOrbit++) {
        std::vector<int> faceIds;
        for (auto& brin : myOrbits[iOrbit]) {
            faceIds.push_back(edgeToFace[brin]);
        }
        newFaces.push_back(faceIds);
//        CombinMap m;
//        m.addFace(faceIds.size(), {}, faceIds);
//        maps.push_back(m);
    }

    // Merge all partial graphs toghether
    // To try later

    while (!newFaces.empty()) {
        auto face = newFaces.back();
        newFaces.pop_back();
        if (face.size() < 3) {
            face.push_back(1000 * face[0] + face[1]);
        }
        bool faceIntegrated = false;

        if (dual.root == nullptr) {
            dual.addFace(face.size(), {}, face);
            faceIntegrated = true;
        } else {
            auto currentListOfBrins = dual.allBrins();
            Brin* attachedBrin = nullptr;
            for (auto& brin : currentListOfBrins) {
                auto source = brin->getSource();
                auto dest = brin->getDest();
                /// TODO : Check for adjacent source and dest, then attach it.
                /// After that, collapse as much as possible the graph. Repeat.
                if (isIn(source, face) && isIn(dest, face)) {
                    auto sourcePos = std::find(face.begin(), face.end(), source);
                    auto destPos = std::find(face.begin(), face.end(), dest);
                    auto distance = std::distance(sourcePos, destPos);
                    // Check if they are next to each other on the face
                    if (distance == -1 || distance == -int(face.size() - 1)) {
                        attachedBrin = brin;
                        break;
                    }
                }
            }
            if (attachedBrin != nullptr) {
                faceIntegrated = true;
                auto source = attachedBrin->getSource();
                auto dest = attachedBrin->getDest();
                std::vector<int> newLabels;
                for (auto it = std::find(face.begin(), face.end(), source) + 1; it < face.end(); it++)
                    newLabels.push_back(*it);
                for (auto it = face.begin(); it < std::find(face.begin(), face.end(), dest); it++)
                    newLabels.push_back(*it);
                dual.addFace(face.size(), {attachedBrin->beta2}, newLabels);

            }
        }
        if (!faceIntegrated) {
            newFaces.insert(newFaces.begin(), face); // Put it back in the array, to be treated later
        } else {
        }
    }
    dual.collapseDuplicateNodes();
    return dual;
}

void CombinMap::debug()
{
    /// Todo, quickly!
    std::cout << this->toString() << std::endl;
}

void CombinMap::clean()
{
    std::vector<Brin*> seenBrins;
    std::vector<Brin*> unseenBrins({this->root});

    while (!unseenBrins.empty()) {
        Brin* current = unseenBrins.back();
        unseenBrins.pop_back();
        if (isIn(current, seenBrins))
            continue;
        seenBrins.push_back(current);

        // Find double-edges
        if (current == current->beta1->beta1) {
            Brin* twinA = current->beta2;
            Brin* twinB = current->beta1->beta2;
            // "Current" is going to get deleted
            if (isIn(this->root, {current, current->beta1}))
                this->root = twinA;
            twinA->beta2 = twinB;
            twinB->beta2 = twinA;
        } else {
            // All clear
            if (!isIn(current->beta2, seenBrins))
                unseenBrins.push_back(current->beta2);
            if (!isIn(current->beta1, seenBrins))
                unseenBrins.push_back(current->beta1);
        }
    }
}

void CombinMap::collapse(Brin *brinToMerge)
{
    if (brinToMerge->degree() == 1 && brinToMerge->beta2->degree() == 1) {
        // Just destroy the graph
        this->root = nullptr;
        return;
    } else if (brinToMerge->beta2->degree() == 1) {
        // Special case

        // Source node becomes Dest node
        brinToMerge->affectDest(brinToMerge->getSource());

        brinToMerge->beta2->previous()->beta1 = brinToMerge->beta1;
        if (this->root == brinToMerge || this->root == brinToMerge->beta2)
            this->root = brinToMerge->beta1;
    } else if (brinToMerge->degree() == 1) {
        // Special case, just a repetition of the above case, without the relabelling of the node
        brinToMerge->previous()->beta1 = brinToMerge->beta2->beta1;
        if (this->root == brinToMerge || this->root == brinToMerge->beta2)
            this->root = brinToMerge->beta2->beta1;
    } else {
        // Normal case

        // Source node becomes Dest node
        brinToMerge->affectDest(brinToMerge->getSource());

        Brin* previousA = brinToMerge->previous();
        Brin* previousB = brinToMerge->beta2->previous();
        previousA->beta1 = brinToMerge->beta1;
        previousB->beta1 = brinToMerge->beta2->beta1;

        if (this->root == brinToMerge || this->root == brinToMerge->beta2)
            this->root = previousA;
    }
    this->clean();
}

int CombinMap::addNeutralComponents()
{
    return -1;
}

std::string CombinMap::toString()
{
    if (this->root == nullptr)
        return "Empty map";
    std::string out = "Graph: \n";
    auto brins = this->allBrins();
    std::sort(brins.begin(), brins.end(), [](const auto& a, const auto& b) -> bool {
        return a->index < b->index;
    });
    for (auto& brin : brins)
        out += brin->toString() + (this->root == brin ? " (root) " : "") + "\n";
    return out;
}

Graph<int> CombinMap::toGraph()
{
    auto unorientedEdges = this->getUnorientedEdges();
//    std::vector<int> nodesIDs;
    Graph<int> g;
    for (auto [startIndex, endIndex] : unorientedEdges) {
        if (g.findNodeByID((int)startIndex)) {
        } else {
            g.addNode(startIndex);
        }
        if (g.findNodeByID((int)endIndex)) {
        } else {
            g.addNode(endIndex);
        }
        g.addConnection(startIndex, endIndex);
    }
    return g;
}
