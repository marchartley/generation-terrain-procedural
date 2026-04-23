#ifndef ENVCURVE_H
#define ENVCURVE_H

class EnvCurve;
class EnvCurveInstance;

#include "EnvObject/EnvObject.h"

class EnvCurve : public EnvObject {
public:
    EnvCurve();

    virtual EnvObjectInstance* instantiate();

    virtual EnvCurve* clone() const;

    virtual void clearKelvinlets();

    // static EnvCurve* fromJSON(nlohmann::json content);
    float width;
    float length;

    std::vector<Kelvinlet*> startingPointKelvinlets;
    std::vector<Kelvinlet*> endingPointKelvinlets;
    std::vector<KelvinletCurve*> curveKelvinlets;

    // virtual nlohmann::json toJSON() const;

    enum CURVE_FOLLOW { GRADIENTS, ISOVALUE, SKELETON };

    CURVE_FOLLOW curveFollow = CURVE_FOLLOW::GRADIENTS;


    virtual bool isCurve() const { return true; }
};

class EnvCurveInstance : public EnvObjectInstance {
public:
    EnvCurveInstance();
    EnvCurveInstance(EnvCurve* definition);
    virtual ~EnvCurveInstance() {}

    EnvCurve* getDefinition() const { return dynamic_cast<EnvCurve*>(definition); }

    BSpline curve;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvCurveInstance* clone();
    // static EnvCurveInstance* instantiate(const std::string& objectName);

    virtual bool placeInTerrain(const Vector3& seedPosition);
    virtual bool placeInTerrain(const BSpline& seedCurve);

    virtual void improvePositionning(float steps);

    virtual void recomputeEvaluationPoints();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual GridV3& computeFlowModification(GridV3& waterFlow, float scale = 1.f);
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    // virtual GridF createHeightfield();

    virtual EnvCurveInstance& translate(const Vector3& translation);
    void updateCurve(const BSpline &newCurve);
};

#endif // ENVCURVE_H
