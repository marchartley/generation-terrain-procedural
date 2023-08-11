#ifndef WAVEFUNCTIONCOLLAPSE_H
#define WAVEFUNCTIONCOLLAPSE_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

class WaveFunctionCollapse;


class WaveFunctionCollapse
{
public:
    WaveFunctionCollapse();
    WaveFunctionCollapse(std::vector<GridI> constraints);
    void run(int sizeX, int sizeY, int sizeZ = 1);
    bool step();

    Vector3 getLowestEntropyCellIndex();
    bool propagate(Vector3 from, bool forceNeighborsPropagation = false);

    static std::vector<GridI> createLabelsFromImage(GridI image, Vector3 tilesDim = Vector3(3, 3, 1));

    Vector3 tilesDim;
    Vector3 tilesAnchor;
    std::vector<GridI> constraints;
    GridI finalMap;
    GridI modifiedCell;
    Matrix3<std::vector<int>> availablePatternsOnCell;
};

#endif // WAVEFUNCTIONCOLLAPSE_H
