#ifndef TERRAINGENERATIONINTERFACE_H
#define TERRAINGENERATIONINTERFACE_H

class TerrainGenerationInterface;
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "Interface/Viewer.h"

class TerrainGenerationInterface : public ActionInterface
{
    Q_OBJECT
public:
    TerrainGenerationInterface(QWidget *parent = nullptr);

    void display(const Vector3& camPos = Vector3(false));
    void displayWaterLevel();

    void createTerrainFromNoise(int nx, int ny, int nz, bool noise2D = false, float noiseStrength = 1.f, float frequency = 1.f, float lacunarity = 2.f, float noise_shifting = 0.0);
    void createTerrainFromFile(std::string filename, std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>());
    void createTerrainFromBiomes(nlohmann::json json_content);
    void createTerrainFromImplicitPatches(nlohmann::json json_content);
    void saveTerrain(std::string filename, Vector3 dimensions = Vector3(false));

    void reloadShaders();

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    void prepareShader(bool reload = false);

    void setWaterLevel(float newLevel);
    void setAmbiantOcclusion(float newValue);

    void updateDisplayedView(const Vector3& newVoxelGridOffset, float newVoxelGridScaling);
    void reloadTerrain(std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>());

    void afterTerrainUpdated();

    void heightmapToVoxels();
    void heightmapToLayers();
    void heightmapToImplicit();
    void heightmapToAll();

    void voxelsToHeightmap();
    void voxelsToLayers();
    void voxelsToImplicit();
    void voxelsToAll();

    void layersToVoxels();
    void layersToHeightmap();
    void layersToImplicit();
    void layersToAll();

    void implicitToVoxels();
    void implicitToLayers();
    void implicitToHeightmap();
    void implicitToAll();

    void openMapUI();
    void saveMapUI();

    void reinforceVoxels();

    void saveErosionDepositionTextureMasks(std::string savingFolder, std::string savingName);
    void saveErosionDepositionTextureMasksOnMultiple();

    void changeDisplayToComparativeMode(bool toComparative);
    void setHeightFactor(float newHeightFactor);

    void changeDisplayDepthMode(bool display);
    void changeDisplayShadowsMode(bool display);

    void updateScalarFieldToDisplay(const GridF& scalarField, float min = 0.f, float max = 1.f);

public:
    float minIsoLevel = -1000.0;
    float maxIsoLevel =  1000.0;

    bool biomeGenerationNeeded = false;
    nlohmann::json biomeGenerationModelData;

    std::map<std::string, int> colorTexturesIndex;
    std::map<std::string, int> normalTexturesIndex;
    std::map<std::string, int> displacementTexturesIndex;

    void setVisu(MapMode _mapMode, SmoothingAlgorithm _smoothingAlgorithm, bool _displayParticles);
    MapMode mapMode = MapMode::GRID_MODE; // MapMode::VOXEL_MODE;
    SmoothingAlgorithm smoothingAlgorithm = SmoothingAlgorithm::MARCHING_CUBES;
//    bool displayParticles = false;

    std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces;
    std::string mapSavingFolder = "saved_maps/";

protected:
    Mesh marchingCubeMesh;
    GLuint dataFieldTex;
    Mesh heightmapMesh;
    GLuint heightmapFieldTex, biomeFieldTex;
    GLuint edgeTableTex, triTableTex;
    GLuint biomesAndDensitiesTex;

    Mesh layersMesh;
    Mesh implicitMesh;

    GLuint allBiomesColorTextures, allBiomesNormalTextures, allBiomesDisplacementTextures;

    Vector3 voxelGridOffset = Vector3(0, 0, 0);
    float voxelGridScaling = 1.f;

//    std::vector<Vector3> randomParticlesPositions;

    std::string lastLoadedMap = "";


    std::chrono::system_clock::time_point startingTime;
    size_t voxelsPreviousHistoryIndex = -1;
    size_t layersPreviousHistoryIndex = -1;

    GridF initialTerrainValues;

public:
    Mesh waterLevelMesh;
    float waterLevel = 0.0f;
    float ambiantOcclusionFactor = 0.f;
    bool displayAsComparativeMode = false;
    float heightFactor = 1.f;
    bool displayDepth = false;
    bool displayShadows = false;

    GridF scalarFieldToDisplay = GridF(1, 1, 1, 0.5f); // Default to "nothing interesting"
};

#endif // TERRAINGENERATIONINTERFACE_H
