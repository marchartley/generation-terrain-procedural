#ifndef CONSTRAINTSSOLVER_H
#define CONSTRAINTSSOLVER_H

#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"
#include "TerrainGen/VoxelGrid.h"
#include "Utils/BSpline.h"

class ConstraintsSolver
{
public:
    ConstraintsSolver();

    int addItem(Vector3* point);
    int addItem(BSpline* curve);

    std::map<int, Vector3> solveWithVoxelGrid(std::shared_ptr<VoxelGrid> mainGrid);
    std::map<int, Vector3> solve(bool checkPossible = true, float deltaMoveForHigherDistances = 1.f, float deltaMoveForLowerDistances = 1.f);

    void addDistanceConstraint(int object1, int object2, float distance);
    void addNormalConstraint(int object, float minAngle, float maxAngle);

    bool checkFeasibility(bool verbose = false);

protected:
    void addConstraintSlot();

    std::map<int, Vector3*> pointsConstrainted;
    std::map<int, BSpline*> curvesConstrainted;

    Matrix3<int> distanceConstraintsAvailable;
    Matrix3<float> distanceConstraints;
    std::map<int, std::tuple<float, float>> objectsNormalConstraints;
};

#endif // CONSTRAINTSSOLVER_H
