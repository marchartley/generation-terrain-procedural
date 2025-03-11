#include "WaveFunctionCollapse.h"

#include "Utils/Utils.h"
#include <set>

WaveFunctionCollapse::WaveFunctionCollapse()
{

}

WaveFunctionCollapse::WaveFunctionCollapse(std::vector<GridI > constraints)
    : constraints(constraints)
{
    this->tilesDim = constraints[0].getDimensions();
    this->tilesAnchor = tilesDim / 2.0;
}

void WaveFunctionCollapse::run(int sizeX, int sizeY, int sizeZ)
{
    std::cout << "Started" << std::endl;
    do {
        std::cout << "Reseting the algorithm" << std::endl;
        this->finalMap = GridI(sizeX, sizeY, sizeZ, -1);
        std::vector<int> allPossible(this->constraints.size());
        for (size_t i = 0; i < constraints.size(); i++)
            allPossible[i] = i;
        this->availablePatternsOnCell = Matrix3<std::vector<int>>(sizeX, sizeY, sizeZ, allPossible);
        this->modifiedCell = GridI(sizeX, sizeY, sizeZ, 0);

        for (size_t i = 0; i < finalMap.size(); i++) {
            Vector3 pos = finalMap.getCoordAsVector3(i);
            if (pos.x == 0 || pos.x == finalMap.sizeX - 1 || pos.y == 0 || pos.y == finalMap.sizeY - 1) {
                finalMap.at(pos) = int(pos.x) % 2; //0;
            }
        }

        while(this->step()) {
            modifiedCell.reset();
            std::cout << this->finalMap.displayValues() << "\n";
            std::cout << "Next step" << std::endl;
        }
    } while (finalMap.min() == -1);
    std::cout << "Ended" << std::endl;
    std::cout << this->finalMap.displayValues() << "\n";
    return;
}

bool WaveFunctionCollapse::step()
{
    Vector3 lowEntroIndex = getLowestEntropyCellIndex();
    if (!lowEntroIndex.isValid() || availablePatternsOnCell.at(lowEntroIndex).size() <= 1)
        return false; // Either all done or it has failed

    this->finalMap.at(lowEntroIndex) = constraints[availablePatternsOnCell.at(lowEntroIndex)[int(random_gen::generate(0, availablePatternsOnCell.at(lowEntroIndex).size()))]].at(tilesAnchor);

    return this->propagate(lowEntroIndex, false);
}

Vector3 WaveFunctionCollapse::getLowestEntropyCellIndex()
{
    int minEntropyValue = -1;
    std::vector<int> minEntropyIndices;
    for (size_t i = 0; i < availablePatternsOnCell.size(); i++) {
        if (finalMap.at(i) == -1 && (minEntropyValue == -1 || int(availablePatternsOnCell[i].size()) < minEntropyValue)) {
            minEntropyValue = availablePatternsOnCell[i].size();
        }
    }
    if (minEntropyValue < 1)
        return Vector3(false);

    for (size_t i = 0; i < availablePatternsOnCell.size(); i++) {
        auto& patterns = availablePatternsOnCell.at(i);
        if (finalMap.at(i) == -1 && patterns.size() == minEntropyValue)
            minEntropyIndices.push_back(i);
    }
    int minEntropyIndex = minEntropyIndices[int(random_gen::generate(0, minEntropyIndices.size()))];
    return availablePatternsOnCell.getCoordAsVector3(minEntropyIndex);
}

/*std::set<int> convertToSet(std::vector<int> v)
{
    return std::set<int>(v.begin(), v.end());
}
std::vector<int> convertToVector(std::set<int> s) {
    return std::vector<int>(s.begin(), s.end());
}*/
bool WaveFunctionCollapse::propagate(const Vector3& from, bool forceNeighborsPropagation)
{
    if (!finalMap.checkCoord(from)) return true;


    int myLabel = this->finalMap.at(from);

    bool cellChanged = forceNeighborsPropagation;
    if (myLabel == -1) {
        std::set<int> possibleStates = convertVectorToSet(availablePatternsOnCell.at(from));

        for (int x = from.x - 1; x <= from.x + 1; x++) {
            for (int y = from.y - 1; y <= from.y + 1; y++) {
                for (int z = from.z - 1; z <= from.z + 1; z++) {
                    Vector3 neighborPos(x, y, z);
                    if (neighborPos == from || !finalMap.checkCoord(neighborPos)) continue;

                    Vector3 neighborRelativePos = neighborPos - from;
                    std::set<int> possibleNeighborLabel;
                    if (finalMap.at(neighborPos) != -1)
                        possibleNeighborLabel.insert(finalMap.at(neighborPos));
                    else {
                        for (int neighborStateID : availablePatternsOnCell.at(neighborPos)) {
                            possibleNeighborLabel.insert(constraints[neighborStateID].at(tilesAnchor));
                        }
                    }

                    std::set<int> tmpStates = possibleStates;
                    for (int stateID : tmpStates) {
                        bool keepPossibleStates = false;
                        for (int label : possibleNeighborLabel) {
                            if (constraints[stateID].at(neighborRelativePos + tilesAnchor) == label && std::find(possibleStates.begin(), possibleStates.end(), stateID) != possibleStates.end())
                                keepPossibleStates = true;
                        }
                        if (!keepPossibleStates)
                            possibleStates.erase(std::find(possibleStates.begin(), possibleStates.end(), stateID));
                    }
                    if (tmpStates.size() != possibleStates.size())
                    {
                        this->modifiedCell.at(from) = 1;
                        this->availablePatternsOnCell.at(from) = convertSetToVector(possibleStates);
                        cellChanged = true;
                    }
                }
            }
        }
    } else {
//        return;
    }
    bool returnValue = true;
    if (cellChanged) {
        if (availablePatternsOnCell.at(from).size() == 1)
            finalMap.at(from) = constraints[availablePatternsOnCell.at(from)[0]].at(tilesAnchor);
        else if (availablePatternsOnCell.at(from).size() == 0)
            return false;
        for (int x = from.x - 1; x <= from.x + 1; x++) {
            for (int y = from.y - 1; y <= from.y + 1; y++) {
                for (int z = from.z - 1; z <= from.z + 1; z++) {
                    if (from != Vector3(x, y, z))
                        returnValue &= propagate(Vector3(x, y, z));
                }
            }
        }
    }
    return returnValue;

    /*
    // Label is now set yet
    if (myLabel == -1) {
        // Check all neighbors
        for (int x = from.x - 2; x <= from.x + 2; x++) {
            for (int y = from.y - 2; y <= from.y + 2; y++) {
                for (int z = from.z - 2; z <= from.z + 2; z++) {
                    Vector3 neighborPos(x, y, z);
                    // If neighbor didn't change since step, just ignore it
                    if (modifiedCell.at(neighborPos) == 0) continue;

                    std::vector<int> possibleStates;

                    // Check all his possible state against all ours
                    for (size_t iState = 0; iState < availablePatternsOnCell.at(from).size(); iState ++) {
                        GridI& checkedState = constraints[availablePatternsOnCell.at(from)[iState]];
                        bool stopTestingThisState = false;
                        // Neighbor is set, just check his state
                        if (finalMap.at(neighborPos) != -1) {
                            // If he's not fitting, ignore this state
                            if (this->constraints[availablePatternsOnCell.at(from)[iState]].at(neighborPos - (from - Vector3(1, 1, 1))) != finalMap.at(neighborPos))
                                stopTestingThisState = true;
                        } else {
                            Vector3 neighborOffset = neighborPos - from;
                            for (int _x = -1; _x <= 1 && !stopTestingThisState; _x++) {
                                for (int _y = -1; _y <= 1 && !stopTestingThisState; _y++) {
                                    for (int _z = -1; _z <= 1 && !stopTestingThisState; _z++) {
                                        Vector3 conflictedPos = from + Vector3(_x, _y, _z);
                                        for (size_t nState = 0; nState < availablePatternsOnCell.at(neighborPos).size() && !stopTestingThisState; nState ++) {
                                            if (!checkedState.checkCoord(conflictedPos - (neighborPos - Vector3(1, 1, 1)))) break;
                                            int stateValue = constraints[availablePatternsOnCell.at(neighborPos)[nState]].at(conflictedPos - (neighborPos - Vector3(1, 1, 1)));
                                            if (stateValue != checkedState.at(Vector3(_x, _y, _z) + Vector3(1, 1, 1)))
                                                stopTestingThisState = true;
                                        }
                                    }
                                }
                            }
                        }

                        if (!stopTestingThisState)
                            possibleStates.push_back(availablePatternsOnCell.at(from)[iState]);
                    }
                }
            }
        }
    } else {
        return;
    }*/
}

std::vector<GridI > WaveFunctionCollapse::createLabelsFromImage(GridI image, const Vector3& tilesDim)
{
    std::vector<GridI > patterns;
    Vector3 halfDimLow = Vector3(1, 1, 0);
    Vector3 halfDimHigh = Vector3(1, 1, 0);
    for (int x = halfDimLow.x; x < image.sizeX - halfDimHigh.x; x++) {
        for (int y = halfDimLow.y; y < image.sizeY - halfDimHigh.y; y++) {
            for (int z = halfDimLow.z; z < image.sizeZ - halfDimHigh.z; z++) {
                patterns.push_back(image.subset(x-1, x + tilesDim.x-1, y-1, y + tilesDim.y-1, z, z + tilesDim.z));
            }
        }
    }
    return patterns;
}
