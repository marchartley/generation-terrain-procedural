#ifndef KARSTHOLE_H
#define KARSTHOLE_H

#include "Karst/KarstHoleProfile.h"
#include "DataStructure/Matrix3.h"

class KarstHole
{
public:
    KarstHole();
    KarstHole(Vector3 start, Vector3 end);
    KarstHole(BSpline fullPath);

    KarstHoleProfile interpolate(float t);

    std::tuple<Matrix3<float>, Vector3> generateMask();

    bool segmentToTriangleCollision(Vector3 s1, Vector3 s2, Vector3 t1, Vector3 t2, Vector3 t3);

    KarstHoleProfile startingProfile;
    KarstHoleProfile endingProfile;
    BSpline path;
};

#endif // KARSTHOLE_H
