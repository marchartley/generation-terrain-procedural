#include "ConstraintsSolver.h"

ConstraintsSolver::ConstraintsSolver()
    : distanceConstraintsAvailable(GridI(0, 0, 0)),
      distanceConstraints(GridF(0, 0, 0))
{

}

int ConstraintsSolver::addItem(Vector3 *point)
{
    int ID = distanceConstraints.sizeX;
    this->pointsConstrainted.insert(std::pair<int, Vector3*>(ID, point));
    this->addConstraintSlot();
    this->addNormalConstraint(ID, -10, 10);
    return ID;
}

int ConstraintsSolver::addItem(BSpline *curve)
{
    int ID = distanceConstraints.sizeX;
    this->curvesConstrainted.insert(std::pair<int, BSpline*>(ID, curve));
    this->addConstraintSlot();
    this->addNormalConstraint(ID, -10, 10);
    return ID;
}

std::map<int, Vector3> ConstraintsSolver::solveWithVoxelGrid(std::shared_ptr<VoxelGrid> mainGrid)
{
    auto startTime = std::chrono::system_clock::now();
    int number_elements = distanceConstraints.sizeX;
    float epsilon = 1e-5;
    float deltaMove = .1f;
    float deltaMoveDamping = 1.f;
    int maxTries = 10000;
    std::map<int, Vector3> returnedPositions;
    float minError = std::numeric_limits<float>::max();
    float error = 0;

    std::vector<float> errorCounter;
    std::vector<int> stepCounter;

    if (!checkFeasibility()) {
        std::cout << "No config possible! : \n";
        checkFeasibility(true);
        return returnedPositions;
    }

    // If an object has a "normal constraint", extract the available voxels first
    GridF voxels = mainGrid->getVoxelValues();
    GridV3 gradient = voxels.gradient() * -1.f;
    GridF borders_f = voxels.binarize().toDistanceMap();
    GridI   borders(borders_f.getDimensions(), 0);
    borders.raiseErrorOnBadCoord = false;
    for (size_t i = 0; i < borders.size(); i++)
        borders[i] = (borders_f[i] == 1 ? 1 : 0);
    GridF angles(voxels.getDimensions(), -100.f); // Init to -100rad, update only the borders
    angles.raiseErrorOnBadCoord = false;
    angles.defaultValueOnBadCoord = -100.f;
    Vector3 verticalVector(0, 0, 1);
    std::vector<int> borders_IDs;
    for (size_t i = 0; i < voxels.size(); i++) {
        // Don't allow any element to get on the borders of the terrain (it could be the bottom of the map, for example)
        Vector3 vox = voxels.getCoordAsVector3(i);
        if (vox.x == 0 || vox.x == voxels.sizeX - 1 || vox.y == 0 || vox.y == voxels.sizeY - 1 || vox.z == 0 || vox.z == voxels.sizeZ - 1)
            borders[i] = 0;

        if (borders[i] == 1) {
            // On the border, compute the vertical angle of the gradient
            angles[i] = std::acos(verticalVector.dot(gradient[i].normalized()));
            // Save the borders IDs to gain speed... maybe...
            borders_IDs.push_back(i);
        }
    }
    // Matrix used if a position is placed in the air, to guide it to the closest border
    GridI reconnectionMatrix_i = borders;
    for (auto& val : reconnectionMatrix_i)
        val = 1 - val; // Switch 0 to 1
    GridV3 reconnectionMatrix = reconnectionMatrix_i.toDistanceMap().gradient() * -1.f;
    int iteration = 0;
    for (iteration = 0; iteration < 350 && minError > 1.0 * number_elements; iteration++) {
        // Shuffle the border voxels, allowing to assign positions to elements easily
        std::shuffle(borders_IDs.begin(), borders_IDs.end(), random_gen::random_generator);

        std::vector<int> availableBorderID = borders_IDs;
        std::vector<Vector3> elements_position(number_elements);
        for (int i = 0; i < number_elements; i++) {
            float minSlope = -1.f, maxSlope = 2*PI; // more than extem values, just to loose all errors of imprecision (especially for PI...)
            if (objectsNormalConstraints.find(i) != objectsNormalConstraints.end()) {
                // "Normal-cosntrainted" elements should stay in a slope
                std::tie(minSlope, maxSlope) = objectsNormalConstraints.find(i)->second;
            }
            // Check if a position is OK for each element
            for (size_t iBorder = 0; iBorder < availableBorderID.size(); iBorder++) {
                if (minSlope <= angles[availableBorderID[iBorder]] && angles[availableBorderID[iBorder]] <= maxSlope) {
                    // If it's OK, assign the position and remove the position for the other elements
                    elements_position[i] = angles.getCoordAsVector3(availableBorderID[iBorder]);
                    availableBorderID.erase(availableBorderID.begin() + iBorder);
                    break;
                }
            }
        }


        int tries = 0;
        error = 0;
        while (true) { // Change this condition, one day
            bool atLeastOneDisplacementToDo = false;

            GridF currentDistance(number_elements, number_elements, 1, 0.f);
            GridV3 nodesPairs(number_elements, number_elements);
            std::vector<Vector3> moves(number_elements);
            for (int i = 0; i < number_elements; i++) {
                for (int j = 0; j < number_elements; j++) {
                    nodesPairs.at(i, j) = elements_position[j] - elements_position[i];
                    currentDistance.at(i, j) = nodesPairs.at(i, j).norm();
                }
            }

            for (int i = 0; i < number_elements; i++) {
                Vector3 displacement;
                int divisor = 0;
                for (int j = 0; j < number_elements; j++) {
                    if (distanceConstraintsAvailable.at(i, j) == 1) {
                        float distanceToSatisfaction = (currentDistance.at(i, j) - distanceConstraints.at(i, j));
                        displacement += nodesPairs.at(i, j).normalized() * distanceToSatisfaction;
                        error += std::abs(distanceToSatisfaction);
                        divisor ++;
                    }
                }
                if (divisor > 0)
                    displacement /= (float)divisor;
                displacement *= deltaMove;

                float minSlope, maxSlope;
                std::tie(minSlope, maxSlope) = objectsNormalConstraints.find(i)->second;
                bool displacementApproved = false;

                // If the dispacement is moving the element in the air or in an incorrect normal,
                // bring it back to a border and reduce the distance until the ground normal is good
                int displacementApprovalTries = 0;
                while (!displacementApproved) {
                    float displacementDist = displacement.norm();
                    // Reduce the length by 1 if the normal is incompatible
                    displacementDist = displacement.norm();
                    if (minSlope > angles.at(elements_position[i] + displacement) ||
                        angles.at(elements_position[i] + displacement) > maxSlope)
                    {
                        displacement.setMag(--displacementDist);
                    }
                    // Bring back to border
                    int reconnectionTries = 0;
                    while (borders.at(elements_position[i] + displacement) == 0) {
                        Vector3 reconnection;
                        if (reconnectionMatrix.checkCoord(elements_position[i] + displacement)) {
                            reconnection = reconnectionMatrix.at(elements_position[i] + displacement);
                        } else {
                            // If it's outside, bring it back inside
                            reconnection = (elements_position[i] + displacement) - reconnectionMatrix.getMirrorPosition(elements_position[i] + displacement);
                        }
                        displacement += reconnection;
                        if (reconnectionTries > displacementDist)
                            break;
                        if (displacementDist > 100000 || displacementDist != displacementDist) { // displacement goes to inf or NaN
                            displacement = Vector3();
                            displacementDist = 0;
                        }
                        reconnectionTries++;
                    }

                    displacementApprovalTries++;
                    if (displacementApprovalTries > 1000)
                    {
                        displacement = Vector3();
                        break;
                    }

                    // If displacementDist < 0, stop trying
                    displacementApproved = (borders.at(elements_position[i] + displacement) == 1 &&
                                            minSlope <= angles.at(elements_position[i] + displacement) &&
                                            angles.at(elements_position[i] + displacement) <= maxSlope) || displacementDist <= 0;
                }
                moves[i] = displacement;

                elements_position[i] += displacement;

                if (displacement.norm2() > epsilon) atLeastOneDisplacementToDo = true;
            }
            tries++;
            if (!atLeastOneDisplacementToDo || tries >= maxTries)
                break;
            deltaMove *= deltaMoveDamping;
        }
        if (error < minError) {
            minError = error;
            for(int i = 0; i < int(elements_position.size()); i++)
                returnedPositions[i] = elements_position[i];
        }
        errorCounter.push_back(error);
        stepCounter.push_back(tries);
        std::cout << "Iter #" << iteration << ": error = " << error << " after " << tries << " steps" << std::endl;
    }
    auto endTime = std::chrono::system_clock::now();
    std::cout << "Error: " << minError << " on iteration #" << iteration << "\nTime: " << std::chrono::duration<float>(endTime - startTime).count() << "s\nDetails:";
//    for (int i = 0; i < iteration; i++) {
//        std::cout << "\n- iter #" << (i+1) << ": Error = " << errorCounter[i] << " after " << stepCounter[i] << " steps";
//    }
    std::cout << std::endl;

    return returnedPositions;
}

std::map<int, Vector3> ConstraintsSolver::solve(bool checkPossible, float deltaMoveForHigherDistances, float deltaMoveForLowerDistances)
{
    auto startTime = std::chrono::system_clock::now();
    int number_elements = distanceConstraints.sizeX;
    float epsilon = 1e-5;
    float deltaMove = .1f;
    float deltaMoveDamping = 1.f;

    float deltaMoveIfDistanceIsGreater = deltaMoveForHigherDistances;
    float deltaMoveIfDistanceIsLower = deltaMoveForLowerDistances;

    int maxTries = 10000;
    std::map<int, Vector3> returnedPositions;
    float minError = std::numeric_limits<float>::max();
    float error = 0;

    std::vector<float> errorCounter;
    std::vector<int> stepCounter;

    if (checkPossible && !checkFeasibility()) {
        std::cout << "No config possible! : \n";
        checkFeasibility(true);
        return returnedPositions;
    }

    int iteration = 0;
    for (iteration = 0; iteration < this->numberIterations && minError > this->stoppingEpsilon; iteration++) {

        std::vector<Vector3> elements_position(number_elements);
        for (int i = 0; i < number_elements; i++) {
            elements_position[i] = *this->pointsConstrainted[i];
        }
        if (iteration > 0) {
            std::shuffle(elements_position.begin(), elements_position.end(), random_gen::random_generator);
        }

        int tries = 0;
        error = 0;
        while (true) { // Change this condition, one day
            bool atLeastOneDisplacementToDo = false;

            GridF currentDistance(number_elements, number_elements, 1, 0.f);
            GridV3 nodesPairs(number_elements, number_elements);
            std::vector<Vector3> moves(number_elements);
            for (int i = 0; i < number_elements; i++) {
                for (int j = 0; j < number_elements; j++) {
                    nodesPairs.at(i, j) = elements_position[j] - elements_position[i];
                    currentDistance.at(i, j) = nodesPairs.at(i, j).norm();
                }
            }

            for (int i = 0; i < number_elements; i++) {
                Vector3 displacement;
                int divisor = 0;
                for (int j = 0; j < number_elements; j++) {
                    if (distanceConstraintsAvailable.at(i, j) == 1) {
                        float distanceToSatisfaction = (currentDistance.at(i, j) - distanceConstraints.at(i, j));
                        displacement += nodesPairs.at(i, j).normalized() * distanceToSatisfaction * (distanceToSatisfaction > 0 ? deltaMoveIfDistanceIsGreater : deltaMoveIfDistanceIsLower);
                        error += std::abs(distanceToSatisfaction);
                        divisor ++;
                    }
                }
                if (divisor > 0)
                    displacement /= (float)divisor;
                displacement *= deltaMove;

                moves[i] = displacement;

                elements_position[i] += displacement;

                if (displacement.norm2() > epsilon) atLeastOneDisplacementToDo = true;
            }
            tries++;
            if (!atLeastOneDisplacementToDo || tries >= maxTries)
                break;
            deltaMove *= deltaMoveDamping;
        }
        if (error < minError) {
            minError = error;
            for(int i = 0; i < int(elements_position.size()); i++)
                returnedPositions[i] = elements_position[i];
        }
        errorCounter.push_back(error);
        stepCounter.push_back(tries);
//        std::cout << "Iter #" << iteration << ": error = " << error << " after " << tries << " steps" << std::endl;
    }
    auto endTime = std::chrono::system_clock::now();
//    std::cout << "Error: " << minError << " on iteration #" << iteration << "\nTime: " << std::chrono::duration<float>(endTime - startTime).count() << "s\nDetails:";
//    for (int i = 0; i < iteration; i++) {
//        std::cout << "\n- iter #" << (i+1) << ": Error = " << errorCounter[i] << " after " << stepCounter[i] << " steps";
//    }
    std::cout << std::endl;

    return returnedPositions;
}

void ConstraintsSolver::addDistanceConstraint(int object1, int object2, float distance)
{
    distanceConstraintsAvailable.at(object1, object2) = 1;
    distanceConstraints.at(object1, object2) = distance;
    distanceConstraintsAvailable.at(object2, object1) = 1;
    distanceConstraints.at(object2, object1) = distance;
}

void ConstraintsSolver::addNormalConstraint(int object, float minAngle, float maxAngle)
{
    objectsNormalConstraints[object] = std::tuple<float, float>(minAngle, maxAngle);
}

std::string id_to_char(int id) { return std::to_string('A' + id); }
bool ConstraintsSolver::checkFeasibility(bool verbose)
{
    // For every element, check if AB + BC >= AC
    int number_elements = distanceConstraints.sizeX;
    for (int i = 0; i < number_elements - 1; i++) {
        for (int j = i + 2; j < number_elements; j++) {
            // If there is no constraints between AB, BC or AC, no need to check AC
            if (distanceConstraintsAvailable.at(i  , j-1) == 0 ||
                    distanceConstraintsAvailable.at(j-1, j ) == 0 ||
                    distanceConstraintsAvailable.at(i  , j  ) == 0)
                continue;
            float distanceAB = distanceConstraints.at(i  , j-1);
            float distanceBC = distanceConstraints.at(j-1, j  );
            float distanceAC = distanceConstraints.at(i  , j  );

            if (verbose)
                std::cout << id_to_char(i) << id_to_char(j-1) << "(" << distanceAB << ") + " << id_to_char(j-1) << id_to_char(j) << "(" << distanceBC << ") >= " << id_to_char(i) << id_to_char(j) << "(" << distanceAC << ")? " << (distanceAB + distanceBC >= distanceAC ? "Yes" : "No") << "\n";
            if (distanceAB + distanceBC < distanceAC)
                return false;
        }
    }
    return true;
}

void ConstraintsSolver::addConstraintSlot()
{
    if (distanceConstraintsAvailable.empty()) {
        distanceConstraintsAvailable = GridI(1, 1, 1);
        distanceConstraints = GridF(1, 1, 1);
    } else {
        distanceConstraintsAvailable.insertRow(distanceConstraintsAvailable.sizeX, 0);
        distanceConstraintsAvailable.insertRow(distanceConstraintsAvailable.sizeY, 1);
        distanceConstraints.insertRow(distanceConstraints.sizeX, 0);
        distanceConstraints.insertRow(distanceConstraints.sizeY, 1);
    }
}
