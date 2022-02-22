#ifndef KARSTHOLEPROFILE_H
#define KARSTHOLEPROFILE_H

#include <vector>
#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"

enum KarstHolePredefinedShapes
{
    TUBE = 0,
    SOLUBLE_BED = 1,
    PASSAGE = 2,
    KEYHOLE = 3,
    CANYON = 4,
};

class KarstHoleProfile
{
public:
    KarstHoleProfile();
    KarstHoleProfile(KarstHolePredefinedShapes shape);
    KarstHoleProfile(BSpline shape);
    KarstHoleProfile(std::vector<Vector3> shape);

    KarstHoleProfile& setNumberOfVertices(int vertice_count);
    KarstHoleProfile& setSize(float sizeX, float sizeY);

    static BSpline createTubeProfile();
    static BSpline createSolubleBedProfile();
    static BSpline createPassageProfile();
    static BSpline createKeyholeProfile();
    static BSpline createCanyonProfile();

    BSpline vertices;
};

#endif // KARSTHOLEPROFILE_H
