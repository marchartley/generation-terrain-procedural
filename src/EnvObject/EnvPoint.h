#ifndef ENVPOINT_H
#define ENVPOINT_H

#include "EnvObject/EnvObject.h"

class EnvPoint : public EnvObject {
public:
    EnvPoint();

    static EnvObject* fromJSON(nlohmann::json content);

    Vector3 position;
    float radius;
    float height = 10.f;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "");
    virtual Vector3 getVector(const Vector3& position, std::string complement = "");
    virtual Vector3 getNormal(const Vector3& position);
    virtual Vector3 getDirection(const Vector3& position);
    virtual Vector3 getProperty(const Vector3& position, std::string prop) const;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvPoint* clone();
    virtual void applySandDeposit();
    virtual void applySandAbsorption();
    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch(ImplicitPrimitive *previousPrimitive = nullptr);
    virtual GridF createHeightfield() const;

    virtual EnvObject& translate(const Vector3& translation);
};

#endif // ENVPOINT_H