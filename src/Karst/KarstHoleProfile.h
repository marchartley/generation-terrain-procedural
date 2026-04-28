#ifndef KARSTHOLEPROFILE_H
#define KARSTHOLEPROFILE_H

#include <vector>
#include "DataStructure/Vector3.h"
#include "Curves/CatmullRomSpline.h"

enum KarstHolePredefinedShapes
{
    TUBE = 0,
    SOLUBLE_BED = 1,
    PASSAGE = 2,
    KEYHOLE = 3,
    CANYON = 4,
    CRACK = 5,
    STAR = 6
};

class KarstHoleProfile
{
public:
    KarstHoleProfile();
    KarstHoleProfile(KarstHolePredefinedShapes shape);
    KarstHoleProfile(CatmullRomSpline shape);
    KarstHoleProfile(std::vector<Vector3> shape);

    KarstHoleProfile& rotateTowardVector(CatmullRomSpline path, float t);
    KarstHoleProfile& translate(const Vector3& new_pos, bool verbose = false);
    KarstHoleProfile& scale(float scale_factor, bool verbose = false);
    KarstHoleProfile interpolate(KarstHoleProfile other, CatmullRomSpline path, float t, float previousAcceptedTime = -1.f, float nextAcceptedTime = -1.f);
    std::pair<KarstHoleProfile, std::vector<std::vector<Vector3>>> interpolateAndGetMesh(KarstHoleProfile other, CatmullRomSpline path, float t);

    KarstHoleProfile& rotateIndicesUntilBestFitWith(KarstHoleProfile& otherProfile, int numberOfPointsUsed);

    KarstHoleProfile& setNumberOfVertices(int vertice_count);
    KarstHoleProfile& setSize(float sizeX, float sizeY);

    std::vector<std::vector<int>> computeTrianglesIndices(const std::vector<Vector3>& points);

    static CatmullRomSpline createTubeProfile();
    static CatmullRomSpline createSolubleBedProfile();
    static CatmullRomSpline createPassageProfile();
    static CatmullRomSpline createKeyholeProfile();
    static CatmullRomSpline createCanyonProfile();
    static CatmullRomSpline createCrackProfile();
    static CatmullRomSpline createStarProfile();

    CatmullRomSpline vertices;
    Vector3 scaling = Vector3(1.f, 1.f, 1.f);
};

#endif // KARSTHOLEPROFILE_H
