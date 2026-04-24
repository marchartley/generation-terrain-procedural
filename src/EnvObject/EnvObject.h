#ifndef ENVOBJECT_H
#define ENVOBJECT_H

#include "Curves/ShapeCurve.h"
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
// #include "EnvObject/EnvScenario.h"

class EnvPoint;
class EnvCurve;
class EnvArea;

class EnvObject;


class EnvPointInstance;
class EnvCurveInstance;
class EnvAreaInstance;

class EnvObjectInstance;

class EnvironmentalScene;
// #include "EnvObject/EnvironmentalScene.h"


struct DepositionRate {
    float radius = 0;
    float rate = 0;
};

struct AbsorptionRate {
    float radius = 0;
    float rate = 0;
};


struct EnvPointKelvinletAttachment {
    float relativeDistanceToCenter;
    float angleFromCurrents;
    KelvinletPoint* kelvinlet;
};

struct EnvObjectKelvinletAttachment {
    enum AnchorPoint { START, END, CURVE };

    AnchorPoint anchor;
    EnvPointKelvinletAttachment* kelvinletFromPoint;
    KelvinletCurve* kelvinletFromCurve;
};


class EnvObject
{
public:
    EnvObject();
    virtual ~EnvObject();
    static std::function<float(const Vector3&)> parseFittingFunction(std::string formula, std::string currentObject, EnvironmentalScene *scene, bool removeSelfInstances = false, EnvObjectInstance* myObject = nullptr);

    float height;

    void updateFittingFunction();


    std::string name;
    std::string s_FittingFunction;
    std::function<float(const Vector3&)> fittingFunction;
    std::string s_FitnessFunction;
    std::function<float(const Vector3&)> fitnessFunction;
    // Vector3 flowEffect;
    std::map<std::string, DepositionRate> materialDepositionRate;
    std::map<std::string, AbsorptionRate> materialAbsorptionRate;
    std::map<std::string, DepositionRate> materialDepositionOnDeath;
    float age = 0.f;
    std::map<std::string, float> needsForGrowth;
    std::map<std::string, float> currentSatisfaction;
    float fitnessScoreAtCreation = -1.f;
    std::vector<Vector3> evaluationPositions;
    float minScore = 0.f;

    TerrainTypes material;
    ImplicitPatch::PredefinedShapes implicitShape;

    SnakeSegmentationParameters* snakeParameters;
    SnakeImageFieldImplicit* snakeField;

    virtual EnvObjectInstance* instantiate() = 0;


    virtual bool isPoint() const { return false; }
    virtual bool isCurve() const { return false; }
    virtual bool isArea() const { return false; }

    virtual EnvObject* clone() const = 0;

    virtual void clearKelvinlets() = 0;


    enum HeightmapFrom {
        SURFACE, GROUND, WATER
    };
    HeightmapFrom heightFrom = SURFACE;



    EnvironmentalScene* scene;

    /*

    int ID = -1;
    int spawnTime = 0;

    SnakeSegmentation snake;
    // bool snakeDefined = false;





    virtual float getSqrDistance(const Vector3& position) = 0;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const = 0;

    virtual void recomputeEvaluationPoints() = 0;

    virtual EnvObject* clone() = 0;
    virtual float computeGrowingState();
    virtual float computeGrowingState2();
    virtual void applyDeposition(EnvMaterial& material) = 0;
    virtual void applyAbsorption(EnvMaterial& material) = 0;
    virtual void applyDepositionOnDeath() = 0;
    virtual GridV3& computeFlowModification(GridV3& waterFlow) = 0;
    virtual ImplicitPatch* createImplicitPatch(const GridF& height, ImplicitPrimitive *previousPrimitive = nullptr) = 0;
    virtual GridF createHeightfield();
    virtual EnvObject& translate(const Vector3& translation) = 0;
    float evaluate(const Vector3& position);
    float evaluate();

    virtual bool placeInTerrain(const Vector3& seedPosition) = 0;
    virtual bool placeInTerrain(const BSpline& seedCurve) = 0;

    virtual void improvePositionning(float stepsOrDistance) = 0;

    void die();

    bool premature = false;

    bool createdManually = false;
    bool geometryNeedsUpdate = true;

    enum HeightmapFrom {
        SURFACE, GROUND, WATER
    };
    HeightmapFrom heightFrom = SURFACE;

    GridV3 _cachedFlowModif;
    ImplicitPatch* _patch = nullptr;
    GridF _cachedHeightfield;

    GridF _cachedAbsorptionDepositionField;

    Vector3 storedOrientation = Vector3::invalid;



    EnvironmentalScene* scene;
    */
};




class EnvObjectInstance
{
public:
    EnvObjectInstance();
    EnvObjectInstance(EnvObject *definition);
    virtual ~EnvObjectInstance() {}

    float age = 0.f;
    float fitnessScoreAtCreation = -1.f;
    std::vector<Vector3> evaluationPositions;
    int ID = -1;
    int spawnTime = 0;

    EnvObject* definition;

    virtual EnvObject* getDefinition() const { return dynamic_cast<EnvObject*>(definition); }

    SnakeSegmentation snake; //

    virtual float getSqrDistance(const Vector3& position) = 0;
    virtual std::map<std::string, Vector3> getAllProperties(const Vector3& position) const = 0;

    virtual void recomputeEvaluationPoints() = 0;

    virtual EnvObjectInstance* clone() = 0;
    virtual float computeGrowingState();
    virtual float computeGrowingState2();
    virtual void applyDeposition(EnvMaterial& material) = 0;
    virtual void applyAbsorption(EnvMaterial& material) = 0;
    virtual void applyDepositionOnDeath() = 0;
    virtual GridV3& computeFlowModification(GridV3& waterFlow, float scale = 1.f) = 0;
    virtual ImplicitPatch* createImplicitPatch(const GridF& height, ImplicitPrimitive *previousPrimitive = nullptr) = 0;
    virtual GridF createHeightfield();
    virtual EnvObjectInstance& translate(const Vector3& translation) = 0;
    float evaluate(const Vector3& position);
    float evaluate();

    virtual bool placeInTerrain(const Vector3& seedPosition) = 0;
    virtual bool placeInTerrain(const BSpline& seedCurve) = 0;

    virtual void improvePositionning(float stepsOrDistance) = 0;

    void die();

    bool premature = false;

    bool createdManually = false;
    bool geometryNeedsUpdate = true;

    GridV3 _cachedFlowModif;
    ImplicitPatch* _patch = nullptr;
    GridF _cachedHeightfield;

    GridF _cachedAbsorptionDepositionField;

    Vector3 storedOrientation = Vector3::invalid;

    EnvironmentalScene* scene;
};

#endif // ENVOBJECT_H
