#ifndef ENVCURVE_H
#define ENVCURVE_H

#include "EnvObject/EnvObject.h"

class EnvCurve : public EnvObject {
public:
    EnvCurve();

    static EnvCurve* fromJSON(nlohmann::json content);

    BSpline curve;
    float width;
    float length;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvCurve* clone();
    static EnvCurve* instantiate(std::string objectName);

    virtual void recomputeEvaluationPoints();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    // virtual GridF createHeightfield();

    virtual EnvCurve& translate(const Vector3& translation);
    void updateCurve(const BSpline &newCurve);

    virtual nlohmann::json toJSON() const;

    enum CURVE_FOLLOW { GRADIENTS, ISOVALUE, SKELETON };

    CURVE_FOLLOW curveFollow = CURVE_FOLLOW::GRADIENTS;
};

#endif // ENVCURVE_H
