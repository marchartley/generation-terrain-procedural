#ifndef ENVOBJECT_H
#define ENVOBJECT_H

#include "Utils/ShapeCurve.h"
#include "Utils/json.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

#include "TerrainGen/ImplicitPatch.h"

class EnvPoint;
class EnvCurve;
class EnvArea;


class EnvObject
{
public:
    EnvObject();
    virtual ~EnvObject();

    static void readFile(std::string filename);

    static EnvObject* fromJSON(nlohmann::json content);


    std::string name;
    std::function<float(Vector3)> fittingFunction;
    Vector3 flowEffect;
    float sandEffect;
    Vector3 inputDimensions;

    TerrainTypes material;
    ImplicitPatch::PredefinedShapes implicitShape;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "") = 0;
    virtual Vector3 getVector(const Vector3& position, std::string complement = "") = 0;
    virtual EnvObject* clone() = 0;


    virtual void applySandDeposit() = 0;
    virtual std::pair<GridV3, GridF> computeFlowModification() = 0;


    static std::function<float(Vector3)> parseFittingFunction(std::string formula);


    static GridV3 flowfield;
    static GridV3 terrainNormals;
    static GridF sandDeposit;

    static float flowImpactFactor;

    static std::map<std::string, EnvObject*> availableObjects;
    static std::vector<EnvObject*> instantiatedObjects;

    static std::pair<std::string, std::string> extractNameAndComplement(std::string variable);
    static std::pair<float, EnvObject*> getSqrDistanceTo(std::string objectName, const Vector3& position);
    static std::pair<Vector3, EnvObject*> getVectorOf(std::string objectName, const Vector3& position);

    static EnvObject* instantiate(std::string objectName);
    static void removeAllObjects();
    virtual ImplicitPatch* createImplicitPatch() = 0;

    static void applyEffects();
};

class EnvPoint : public EnvObject {
public:
    EnvPoint();

    static EnvObject* fromJSON(nlohmann::json content);

    Vector3 position;
    float radius;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "");
    virtual Vector3 getVector(const Vector3& position, std::string complement = "");
    virtual EnvPoint* clone();
    virtual void applySandDeposit();
    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch();
};

class EnvCurve : public EnvObject {
public:
    EnvCurve();

    static EnvObject* fromJSON(nlohmann::json content);

    BSpline curve;
    float width;
    float length;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "");
    virtual Vector3 getVector(const Vector3& position, std::string complement = "");
    virtual EnvCurve* clone();
    virtual void applySandDeposit();
    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch();
};

class EnvArea : public EnvObject {
public:
    EnvArea();

    static EnvObject* fromJSON(nlohmann::json content);

    ShapeCurve area;
    float width;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "");
    virtual Vector3 getVector(const Vector3& position, std::string complement = "");
    virtual EnvArea* clone();
    virtual void applySandDeposit();
    virtual std::pair<GridV3, GridF> computeFlowModification();
    virtual ImplicitPatch* createImplicitPatch();
};

#endif // ENVOBJECT_H
