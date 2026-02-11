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
// #include "EnvObject/EnvScenario.h"

class EnvPoint;
class EnvCurve;
class EnvArea;

class EnvObject;

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

struct FlowEffect {

};


class EnvObject
{
public:
    EnvObject();
    virtual ~EnvObject();

    // static EnvObject* fromJSON(nlohmann::json content);

    // virtual nlohmann::json toJSON() const;

    static std::function<float(const Vector3&)> parseFittingFunction(std::string formula, std::string currentObject, EnvironmentalScene *scene, bool removeSelfInstances = false, EnvObject* myObject = nullptr);

    float height;


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
    // Vector3 evaluationPosition;

    TerrainTypes material;
    ImplicitPatch::PredefinedShapes implicitShape;
    int ID = -1;
    int spawnTime = 0;

    SnakeSegmentationImplicit snake;
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
    virtual GridV3 computeFlowModification() = 0;
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

    Vector3 storedOrientation = Vector3::invalid();

    EnvironmentalScene* scene;
};

#endif // ENVOBJECT_H
