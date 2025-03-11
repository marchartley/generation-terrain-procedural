#ifndef ENVPOINT_H
#define ENVPOINT_H

#include "EnvObject/EnvObject.h"

class EnvPoint : public EnvObject {
public:
    EnvPoint();

    static EnvPoint* fromJSON(nlohmann::json content);

    Vector3 position;
    float radius;
    float height = 10.f;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvPoint* clone();
    static EnvPoint* instantiate(std::string objectName);

    virtual void recomputeEvaluationPoints();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    // virtual GridF createHeightfield();

    virtual EnvPoint& translate(const Vector3& translation);

    virtual nlohmann::json toJSON() const;

//    float estimateShadowing(const GridV3 &flow, const Vector3& pos);
};

#endif // ENVPOINT_H
