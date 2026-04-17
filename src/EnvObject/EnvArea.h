#ifndef ENVAREA_H
#define ENVAREA_H

class EnvArea;
class EnvAreaInstance;

#include "EnvObject/EnvObject.h"

class EnvArea : public EnvObject {
public:
    EnvArea();

    virtual EnvObjectInstance* instantiate();
    virtual EnvArea* clone() const;

    float width;
    float length;
    float flowAttenuation;

    std::vector<Kelvinlet*> curveKelvinlets;

    bool evaluateInside = false; // Ture = evaluation points inside, false = evaluation points on borders


    virtual bool isArea() const { return true; }
};


class EnvAreaInstance : public EnvObjectInstance {
public:
    EnvAreaInstance();
    EnvAreaInstance(EnvArea* definition);
    virtual ~EnvAreaInstance() {}

    EnvArea* getDefinition() const { return dynamic_cast<EnvArea*>(definition); }

    ShapeCurve curve;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvAreaInstance* clone();
    // static EnvAreaInstance* instantiate(const std::string& objectName);

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

    virtual EnvAreaInstance& translate(const Vector3& translation);
    void updateCurve(const BSpline& newCurve);
};

#endif // ENVAREA_H
