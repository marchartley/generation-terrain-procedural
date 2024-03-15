#ifndef ENVOBJECT_H
#define ENVOBJECT_H

#include "Utils/ShapeCurve.h"
#include "Utils/json.h"
#include "DataStructure/Matrix3.h"
#include "DataStructure/Vector3.h"

#include "TerrainGen/Heightmap.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/ImplicitPatch.h"
#include "Graph/FastPoissonGraph.h"

#include "EnvMaterial.h"

class EnvPoint;
class EnvCurve;
class EnvArea;

class EnvObject;

typedef GraphTemplate<EnvObject*> GraphObj;
typedef GraphNodeTemplate<EnvObject*> GraphNodeObj;

class EnvObject
{
public:
    EnvObject();
    virtual ~EnvObject();

    static void readFile(std::string filename);

    static EnvObject* fromJSON(nlohmann::json content);


    std::string name;
    std::string s_FittingFunction;
    std::function<float(Vector3)> fittingFunction;
    Vector3 flowEffect;
//    float sandEffect;
//    float polypEffect;
//    std::map<std::string, float> materialEffect;
    std::map<std::string, float> materialDepositionRate;
    std::map<std::string, float> materialAbsorptionRate;
    Vector3 inputDimensions;
    float age = 0.f;
    std::map<std::string, float> needsForGrowth;
    std::map<std::string, float> currentSatisfaction;

    TerrainTypes material;
    ImplicitPatch::PredefinedShapes implicitShape;
    int ID = -1;
    float growingState = 0.f;

    virtual float getSqrDistance(const Vector3& position, std::string complement = "") = 0;
    virtual Vector3 getVector(const Vector3& position, std::string complement = "") = 0;
    virtual Vector3 getNormal(const Vector3& position) = 0;
    virtual Vector3 getDirection(const Vector3& position) = 0;
    virtual Vector3 getProperty(const Vector3& position, std::string prop) const = 0;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const = 0;
    virtual EnvObject* clone() = 0;


    virtual float computeGrowingState();

    /*virtual void applySandDeposit() = 0;
    virtual void applySandAbsorption() = 0;
    virtual void applyPolypDeposit() = 0;
    virtual void applyPolypAbsorption() = 0;*/
    virtual void applyDeposition(EnvMaterial& material) = 0;
    virtual void applyAbsorption(EnvMaterial& material) = 0;
    virtual std::pair<GridV3, GridF> computeFlowModification() = 0;


    static std::function<float(Vector3)> parseFittingFunction(std::string formula, std::string currentObject);


    static GridV3 flowfield;
    static GridV3 initialFlowfield;
    static GridV3 terrainNormals;
//    static GridF sandDeposit;
//    static GridF polypDeposit;
//    static std::map<std::string, GridF> materialDeposit;

    static std::map<std::string, EnvMaterial> materials;

    static float flowImpactFactor;

    static std::map<std::string, EnvObject*> availableObjects;
    static std::vector<EnvObject*> instantiatedObjects;

    static std::pair<std::string, std::string> extractNameAndComplement(std::string variable);
    static std::pair<float, EnvObject*> getSqrDistanceTo(std::string objectName, const Vector3& position);
    static std::pair<Vector3, EnvObject*> getVectorOf(std::string objectName, const Vector3& position);

    static EnvObject* instantiate(std::string objectName);
    static void removeAllObjects();
    virtual ImplicitPatch* createImplicitPatch(ImplicitPrimitive *previousPrimitive = nullptr) = 0;
    virtual GridF createHeightfield() const = 0;

    static void applyEffects();
    static void updateSedimentation();
    static void updateFlowfield();
    static void beImpactedByEvents();

    float evaluate(const Vector3& position);

    virtual EnvObject& translate(const Vector3& translation) {}

    static int currentMaxID;

    static std::map<std::string, GridV3> allVectorProperties;
    static std::map<std::string, GridF> allScalarProperties;
    static void precomputeTerrainProperties(const Heightmap& heightmap);
    static void recomputeTerrainPropertiesForObject(const Heightmap& heightmap, std::string objectName);
    static void recomputeFlowAndSandProperties(const Heightmap &heightmap);

    static void reset();

    static GraphObj sceneToGraph();
};

#endif // ENVOBJECT_H
