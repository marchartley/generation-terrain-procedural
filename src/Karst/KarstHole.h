#ifndef KARSTHOLE_H
#define KARSTHOLE_H

#include "Karst/KarstHoleProfile.h"
#include "DataStructure/Matrix3.h"

class KarstHole
{
public:
    KarstHole(float size = 1.f);
    KarstHole(Vector3 start, Vector3 end, float size = 1.f);
    KarstHole(BSpline fullPath, float size = 1.f);

    KarstHoleProfile interpolate(float t, float previousAcceptedTime = -1.f, float nextAcceptedTime = -1.f);

    std::vector<std::vector<Vector3>> generateMesh();
    std::vector<std::vector<Vector3>> computeClosingMesh(std::vector<Vector3>& vertices);

    std::tuple<Matrix3<float>, Vector3> generateMask(std::vector<std::vector<Vector3>> precomputedTriangles = std::vector<std::vector<Vector3>>());

    static int segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3);

    KarstHoleProfile startingProfile;
    KarstHoleProfile endingProfile;
    BSpline path;
    float size;
    int number_of_points = 10;
    int number_of_intermediates = 2;

    std::vector<std::tuple<Vector3, float>> vertexCylinders;
    std::vector<std::tuple<Vector3, Vector3>> cylinders;
    std::map<int, int> vertexGroups;
};

#endif // KARSTHOLE_H
