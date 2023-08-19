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
    CRACK = 5,
    STAR = 6
};

class KarstHoleProfile
{
public:
    KarstHoleProfile();
    KarstHoleProfile(KarstHolePredefinedShapes shape);
    KarstHoleProfile(BSpline shape);
    KarstHoleProfile(std::vector<Vector3> shape);

    KarstHoleProfile& rotateTowardVector(BSpline path, float t);
    KarstHoleProfile& translate(const Vector3& new_pos, bool verbose = false);
    KarstHoleProfile& scale(float scale_factor, bool verbose = false);
    KarstHoleProfile interpolate(KarstHoleProfile other, BSpline path, float t, float previousAcceptedTime = -1.f, float nextAcceptedTime = -1.f);
    std::pair<KarstHoleProfile, std::vector<std::vector<Vector3>>> interpolateAndGetMesh(KarstHoleProfile other, BSpline path, float t);

    KarstHoleProfile& rotateIndicesUntilBestFitWith(KarstHoleProfile& otherProfile, int numberOfPointsUsed);

    KarstHoleProfile& setNumberOfVertices(int vertice_count);
    KarstHoleProfile& setSize(float sizeX, float sizeY);

    std::vector<std::vector<int>> computeTrianglesIndices(const std::vector<Vector3>& points);

    static BSpline createTubeProfile();
    static BSpline createSolubleBedProfile();
    static BSpline createPassageProfile();
    static BSpline createKeyholeProfile();
    static BSpline createCanyonProfile();
    static BSpline createCrackProfile();
    static BSpline createStarProfile();

    BSpline vertices;
    Vector3 scaling = Vector3(1.f, 1.f, 1.f);
};

#include "Utils/json.h"
nlohmann::json karst_profile_to_json(KarstHoleProfile profile);
KarstHoleProfile json_to_karst_profile(nlohmann::json json);
#endif // KARSTHOLEPROFILE_H
