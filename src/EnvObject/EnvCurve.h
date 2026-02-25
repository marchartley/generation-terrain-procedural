#ifndef ENVCURVE_H
#define ENVCURVE_H

class EnvCurve;

#include "EnvObject/EnvObject.h"

class EnvCurve : public EnvObject {
public:
    EnvCurve();

    // static EnvCurve* fromJSON(nlohmann::json content);

    BSpline curve;
    float width;
    float length;

    std::vector<Kelvinlet*> startingPointKelvinlets;
    std::vector<Kelvinlet*> endingPointKelvinlets;
    std::vector<Kelvinlet*> curveKelvinlets;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvCurve* clone();
    // static EnvCurve* instantiate(const std::string& objectName);

    virtual bool placeInTerrain(const Vector3& seedPosition);
    virtual bool placeInTerrain(const BSpline& seedCurve);

    virtual void improvePositionning(float steps);

    virtual void recomputeEvaluationPoints();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual GridV3& computeFlowModification(GridV3& waterFlow);
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    // virtual GridF createHeightfield();

    virtual EnvCurve& translate(const Vector3& translation);
    void updateCurve(const BSpline &newCurve);

    // virtual nlohmann::json toJSON() const;

    enum CURVE_FOLLOW { GRADIENTS, ISOVALUE, SKELETON };

    CURVE_FOLLOW curveFollow = CURVE_FOLLOW::GRADIENTS;


    virtual bool isCurve() const { return true; }
};

#endif // ENVCURVE_H
