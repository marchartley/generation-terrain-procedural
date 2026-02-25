#ifndef ENVPOINT_H
#define ENVPOINT_H

class EnvPoint;

#include "EnvObject/EnvObject.h"

class EnvPoint : public EnvObject {
public:
    EnvPoint();

    // static EnvPoint* fromJSON(nlohmann::json content);

    Vector3 position;
    float radius;
    float height = 10.f;

    std::vector<Kelvinlet*> mainKelvinlets;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvPoint* clone();
    // static EnvPoint* instantiate(const std::string& objectName);

    virtual bool placeInTerrain(const Vector3& seedPosition);
    virtual bool placeInTerrain(const BSpline& seedCurve);

    virtual void improvePositionning(float maxDistance);

    virtual void recomputeEvaluationPoints();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual GridV3& computeFlowModification(GridV3& waterFlow);
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    // virtual GridF createHeightfield();

    virtual EnvPoint& translate(const Vector3& translation);

    virtual bool isPoint() const { return true; }
};

#endif // ENVPOINT_H
