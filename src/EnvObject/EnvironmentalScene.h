#ifndef ENVIRONMENTALSCENE_H
#define ENVIRONMENTALSCENE_H

class EnvironmentalScene;

#include "EnvObject/EnvObject.h"
#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"
#include "EnvObject/EnvScenario.h"
#include "EnvObject/EnvMaterial.h"

// typedef GraphTemplate<EnvObject*> GraphObj;
// typedef GraphNodeTemplate<EnvObject*> GraphNodeObj;

typedef std::pair<std::map<std::string, float>, std::map<std::string, float>> MaterialsTransformation;

class EnvironmentalScene
{
public:
    EnvironmentalScene();

    GridV3 initFlow(bool force = false);

    void readEnvObjectsFile(std::string filename);
    void readEnvObjectsFileContent(std::string content);

    void readEnvMaterialsFile(std::string filename);
    void readEnvMaterialsFileContent(std::string content);

    void readEnvMaterialsTransformationsFile(std::string filename);
    void readEnvMaterialsTransformationsFileContent(std::string content);

    void readScenarioFile(std::string filename);
    void readScenarioFileContent(std::string content);

    std::vector<std::string> getMaterialsToUpdate() const;


    GridV3 flowfield;
    GridV3 initialFlowfield;
    GridV3 terrainNormals;

    std::map<std::string, EnvMaterial> materials;

    float flowImpactFactor;

    std::map<std::string, EnvObject*> availableObjects;
    std::vector<EnvObject*> instantiatedObjects;

    EnvObject* findClosest(std::string objectName, const Vector3& pos);

    EnvObject* instantiate(std::string objectName);
    void removeObject(EnvObject* obj);
    void removeAllObjects();
    // bool applyEffects(const GridF& heights, const GridV3 &userFlow = GridV3());
    // bool updateSedimentation(const GridF& heights);
    std::vector<std::string> updateSedimentationKnowingFluidsAndGradients(const GridF& heights, const GridV3& heightsGradients, const GridV3& smoothFluids, std::vector<std::string> unstableMaterials);
    void stabilizeMaterials(const GridF& heights, int maxIterations = 40); // 40 is enough iterations to find a good stability usually, without taking too much time
    void applyMaterialsTransformations();
    const GridV3& updateFlowfield(const GridV3& userFlow = GridV3(), const GridV3& simulationFlow = GridV3());
    void beImpactedByEvents();

    int currentMaxID;

    std::map<std::string, GridV3> allVectorProperties;
    std::map<std::string, GridF> allScalarProperties;
    void precomputeTerrainProperties(const GridF &heightmap, float waterLevel, float maxHeight);
    void recomputeTerrainPropertiesForObject(std::string objectName);
    void recomputeFlowAndSandProperties(const GridF &heightmap, float waterLevel, float maxHeight);
    void recomputeFlow();

    GridF getHeightmap(const GridF &initialHeightmap, float absoluteWaterLevel, float flowErosionFactor = 0.f, bool displayGrooves = false);

    void reset();

    // GraphObj sceneToGraph();

    std::vector<MaterialsTransformation> transformationRules;

    int currentTime;

    Scenario scenario;
};

#endif // ENVIRONMENTALSCENE_H
