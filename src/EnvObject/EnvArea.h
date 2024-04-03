#ifndef ENVAREA_H
#define ENVAREA_H

#include "EnvObject/EnvObject.h"

class EnvArea : public EnvObject {
public:
    EnvArea();

    static EnvObject* fromJSON(nlohmann::json content);

    ShapeCurve area;
    float width;
    float height;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "");
    virtual Vector3 getVector(const Vector3& position, std::string complement = "");
    virtual Vector3 getNormal(const Vector3& position);
    virtual Vector3 getDirection(const Vector3& position);
    virtual Vector3 getProperty(const Vector3& position, std::string prop) const;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const;
    virtual EnvArea* clone();

    /*virtual void applySandDeposit();
    virtual void applySandAbsorption();
    virtual void applyPolypDeposit();
    virtual void applyPolypAbsorption();*/
    virtual void applyDeposition(EnvMaterial& material);
    virtual void applyAbsorption(EnvMaterial& material);

    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch(const GridF& heights, ImplicitPrimitive *previousPrimitive = nullptr);
    virtual GridF createHeightfield() const;

    virtual EnvObject& translate(const Vector3& translation);
};

#endif // ENVAREA_H
