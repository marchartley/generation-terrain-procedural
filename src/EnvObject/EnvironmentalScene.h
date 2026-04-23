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

    GridV3 &initFlow(bool force = false);

    nlohmann::ordered_json readEnvObjectsFile(const std::string& filename);
    nlohmann::ordered_json readEnvObjectsFileContent(const std::string& content);
    void updateEnvObjectsFileContent(const std::string& filename);

    nlohmann::ordered_json readEnvMaterialsFile(const std::string& filename);
    nlohmann::ordered_json readEnvMaterialsFileContent(const std::string& content);
    void updateEnvMaterialsFileContent(const std::string& filename);

    void readEnvMaterialsTransformationsFile(const std::string& filename);
    void readEnvMaterialsTransformationsFileContent(const std::string& content);
    void updateEnvMaterialsTransformationFileContent(const std::string& filename);

    nlohmann::ordered_json readScenarioFile(const std::string& filename);
    nlohmann::ordered_json readScenarioFileContent(const std::string& content);
    void updateScenarioFileContent(const std::string& filename);

    std::vector<std::string> getMaterialsToUpdate() const;


    GridV3 flowfield;
    GridV3 initialFlowfield;

    std::map<std::string, EnvMaterial> materials;

    std::map<std::string, EnvObject*> availableObjects;
    std::vector<EnvObjectInstance*> instantiatedObjects;

    EnvObjectInstance* findClosest(const std::string& objectName, const Vector3& pos);

    EnvObjectInstance* instantiate(const std::string& objectName);
    void removeObject(EnvObjectInstance *obj);
    void removeAllObjects();
    // bool applyEffects(const GridF& heights, const GridV3 &userFlow = GridV3());
    // bool updateSedimentation(const GridF& heights);
    EnvMaterial& stabilizeOneMaterial(const GridV3& heightsGradients, const GridV3& smoothFluids, EnvMaterial& material, int maxIterations = 40);
    // std::vector<std::string> updateAllMaterialsOnce(const GridV3& heightsGradients, const GridV3& smoothFluids, const std::vector<std::string>& unstableMaterials);
    void stabilizeMaterials(const GridF& heights, int maxIterations = 40); // 40 is enough iterations to find a good stability usually, without taking too much time
    void applyMaterialsTransformations();
    const GridV3& updateFlowfield(const GridV3& userFlow = GridV3(), const GridV3& simulationFlow = GridV3());
    void beImpactedByEvents();

    int currentMaxID;

    std::map<std::string, GridV3> allVectorProperties;
    std::map<std::string, GridF> allScalarProperties;
    void precomputeTerrainProperties(const GridF &heightmap, float waterLevel, float maxHeight);
    void recomputeTerrainPropertiesForObject(const std::string& objectName);
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
