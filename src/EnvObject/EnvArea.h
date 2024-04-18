#ifndef ENVAREA_H
#define ENVAREA_H

#include "EnvObject/EnvObject.h"

class EnvArea : public EnvObject {
public:
    EnvArea();

    static EnvObject* fromJSON(nlohmann::json content);

    ShapeCurve curve;
    float width;
    float length;
    float height;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvArea* clone();

    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);
    virtual void applyDepositionOnDeath();

    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    virtual GridF createHeightfield() const;

    virtual EnvObject& translate(const Vector3& translation);
    void updateCurve(const BSpline& newCurve);

    virtual nlohmann::json toJSON() const;
};

#endif // ENVAREA_H
