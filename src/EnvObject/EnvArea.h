#ifndef ENVAREA_H
#define ENVAREA_H

class EnvArea;

#include "EnvObject/EnvObject.h"

class EnvArea : public EnvObject {
public:
    EnvArea();

    // static EnvArea* fromJSON(nlohmann::json content);

    ShapeCurve curve;
    float width;
    float length;
    float flowAttenuation;

    std::vector<Kelvinlet*> curveKelvinlets;

    bool evaluateInside = false; // Ture = evaluation points inside, false = evaluation points on borders

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvArea* clone();
    // static EnvArea* instantiate(const std::string& objectName);

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

    virtual EnvArea& translate(const Vector3& translation);
    void updateCurve(const BSpline& newCurve);

    virtual bool isArea() const { return true; }
};

#endif // ENVAREA_H
