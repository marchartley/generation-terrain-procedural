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

#include "DataStructure/Kelvinlet.h"

#include "EnvObject/PositionOptimizer.h"
#include "EnvObject/EnvScenario.h"

class EnvPoint;
class EnvCurve;
class EnvArea;

class EnvObject;

typedef GraphTemplate<EnvObject*> GraphObj;
typedef GraphNodeTemplate<EnvObject*> GraphNodeObj;

typedef std::pair<std::map<std::string, float>, std::map<std::string, float>> MaterialsTransformation;


class EnvObject
{
public:
    EnvObject();
    virtual ~EnvObject();

    static void readEnvObjectsFile(std::string filename);
    static void readEnvObjectsFileContent(std::string content);

    static void readEnvMaterialsFile(std::string filename);
    static void readEnvMaterialsFileContent(std::string content);

    static void readEnvMaterialsTransformationsFile(std::string filename);
    static void readEnvMaterialsTransformationsFileContent(std::string content);

    static void readScenarioFile(std::string filename);
    static void readScenarioFileContent(std::string content);

    static EnvObject* fromJSON(nlohmann::json content);

    virtual nlohmann::json toJSON() const;


    std::string name;
    std::string s_FittingFunction;
    std::function<float(Vector3)> fittingFunction;
    Vector3 flowEffect;
    std::map<std::string, float> materialDepositionRate;
    std::map<std::string, float> materialAbsorptionRate;
    std::map<std::string, float> materialDepositionOnDeath;
    Vector3 inputDimensions;
    float age = 0.f;
    std::map<std::string, float> needsForGrowth;
    std::map<std::string, float> currentSatisfaction;
    float fittingScoreAtCreation = -1.f;
    Vector3 evaluationPosition;

    TerrainTypes material;
    ImplicitPatch::PredefinedShapes implicitShape;
    int ID = -1;
    int spawnTime = 0;

    virtual float getSqrDistance(const Vector3& position) = 0;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const = 0;

    virtual EnvObject* clone() = 0;
    virtual float computeGrowingState();
    virtual float computeGrowingState2();
    virtual void applyDeposition(EnvMaterial& material) = 0;
    virtual void applyAbsorption(EnvMaterial& material) = 0;
    virtual void applyDepositionOnDeath() = 0;
    virtual std::pair<GridV3, GridF> computeFlowModification() = 0;
    virtual ImplicitPatch* createImplicitPatch(const GridF& height, ImplicitPrimitive *previousPrimitive = nullptr) = 0;
    virtual GridF createHeightfield() const = 0;
    virtual EnvObject& translate(const Vector3& translation) = 0;
    float evaluate(const Vector3& position);

    void die();

    bool createdManually = false;


    static std::function<float(Vector3)> parseFittingFunction(std::string formula, std::string currentObject, bool removeSelfInstances = false, EnvObject* myObject = nullptr);


    static GridV3 flowfield;
    static GridV3 initialFlowfield;
    static GridV3 terrainNormals;

    static std::map<std::string, EnvMaterial> materials;

    static float flowImpactFactor;

    static std::map<std::string, EnvObject*> availableObjects;
    static std::vector<EnvObject*> instantiatedObjects;

    static EnvObject* findClosest(std::string objectName, const Vector3& pos);

    static EnvObject* instantiate(std::string objectName);
    static void removeObject(EnvObject* obj);
    static void removeAllObjects();
    static bool applyEffects(const GridF& heights);
    static bool updateSedimentation(const GridF& heights);
    static void applyMaterialsTransformations();
    static void updateFlowfield();
    static void beImpactedByEvents();

    static int currentMaxID;

    static std::map<std::string, GridV3> allVectorProperties;
    static std::map<std::string, GridF> allScalarProperties;
    static void precomputeTerrainProperties(const Heightmap& heightmap);
    static void recomputeTerrainPropertiesForObject(std::string objectName);
    static void recomputeFlowAndSandProperties(const Heightmap &heightmap);

    static void reset();

    static GraphObj sceneToGraph();

    static std::vector<MaterialsTransformation> transformationRules;

    static int currentTime;

    GridV3 _cachedFlowModif;

    static Scenario scenario;
};

#endif // ENVOBJECT_H
