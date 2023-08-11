#ifndef KARSTHOLE_H
#define KARSTHOLE_H

#include "Karst/KarstHoleProfile.h"
#include "DataStructure/Matrix3.h"

class KarstHole
{
public:
    KarstHole(float width = 1.f, float height = 1.f, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    KarstHole(Vector3 start, Vector3 end, float width = 1.f, float height = 1.f, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    KarstHole(BSpline fullPath, float width = 1.f, float height = 1.f, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);

    KarstHoleProfile interpolate(float t, float previousAcceptedTime = -1.f, float nextAcceptedTime = -1.f);

    std::vector<std::vector<Vector3>> generateMesh();
    std::vector<std::vector<Vector3>> computeClosingMesh(std::vector<Vector3>& vertices);

    std::tuple<GridF, Vector3> generateMask(std::vector<std::vector<Vector3>> precomputedTriangles = std::vector<std::vector<Vector3>>());

    KarstHoleProfile startingProfile;
    KarstHoleProfile endingProfile;
    BSpline path;
    float width, height;
    int number_of_points = 10;
    int number_of_intermediates = 2;

    std::vector<std::tuple<Vector3, float>> vertexCylinders;
    std::vector<std::tuple<Vector3, Vector3>> cylinders;
    std::map<int, int> vertexGroups;
};

#endif // KARSTHOLE_H
