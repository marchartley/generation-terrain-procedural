#ifndef ENVPOINT_H
#define ENVPOINT_H

class EnvPoint;
class EnvPointInstance;

#include "EnvObject/EnvObject.h"

class EnvPoint : public EnvObject {
public:
    EnvPoint();

    float radius;
    float height;

    std::vector<Kelvinlet*> mainKelvinlets;

    virtual EnvPoint* clone() const;
    virtual EnvObjectInstance* instantiate();

    virtual bool isPoint() const { return true; }
};


class EnvPointInstance : public EnvObjectInstance {
public:
    EnvPointInstance();
    EnvPointInstance(EnvPoint* definition);
    virtual ~EnvPointInstance() {}

    EnvPoint* getDefinition() const { return dynamic_cast<EnvPoint*>(definition); }

    Vector3 position;

    virtual float getSqrDistance(const Vector3& position);
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvPointInstance* clone();
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

    virtual EnvPointInstance& translate(const Vector3& translation);
};

#endif // ENVPOINT_H
