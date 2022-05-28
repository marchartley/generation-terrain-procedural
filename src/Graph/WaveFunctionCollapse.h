#ifndef WAVEFUNCTIONCOLLAPSE_H
#define WAVEFUNCTIONCOLLAPSE_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

class WaveFunctionCollapse;


class WaveFunctionCollapse
{
public:
    WaveFunctionCollapse();
    WaveFunctionCollapse(std::vector<Matrix3<int>> constraints);
    void run(int sizeX, int sizeY, int sizeZ = 1);
    bool step();

    Vector3 getLowestEntropyCellIndex();
    bool propagate(Vector3 from, bool forceNeighborsPropagation = false);

    static std::vector<Matrix3<int>> createLabelsFromImage(Matrix3<int> image, Vector3 tilesDim = Vector3(3, 3, 1));

    Vector3 tilesDim;
    Vector3 tilesAnchor;
    std::vector<Matrix3<int>> constraints;
    Matrix3<int> finalMap;
    Matrix3<int> modifiedCell;
    Matrix3<std::vector<int>> availablePatternsOnCell;
};

#endif // WAVEFUNCTIONCOLLAPSE_H
