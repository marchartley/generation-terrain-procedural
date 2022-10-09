#ifndef TERRAINGENERATIONINTERFACE_H
#define TERRAINGENERATIONINTERFACE_H

class TerrainGenerationInterface;
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Grid.h"
#include "Interface/Viewer.h"

class TerrainGenerationInterface : public ActionInterface
{
    Q_OBJECT
public:
    TerrainGenerationInterface(QWidget *parent = nullptr);

    void display(MapMode mapMode, SmoothingAlgorithm smoothingAlgorithm, bool displayParticles);

    void createTerrainFromNoise(int nx, int ny, int nz, float blockSize, float noise_shifting = 0.0);
    void createTerrainFromFile(std::string filename, std::map<std::string, std::shared_ptr<ActionInterface>> actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>());
    void createTerrainFromBiomes(nlohmann::json json_content);
    void saveTerrain(std::string filename);

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    void prepareShader();

    void updateDisplayedView(Vector3 newVoxelGridOffset, float newVoxelGridScaling);


public:
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Grid> heightmapGrid;
    float minIsoLevel = -1000.0;
    float maxIsoLevel =  1000.0;

    bool biomeGenerationNeeded = false;
    nlohmann::json biomeGenerationModelData;

    std::map<std::string, int> colorTexturesIndex;
    std::map<std::string, int> normalTexturesIndex;
    std::map<std::string, int> displacementTexturesIndex;

protected:
    void regenerateRocksAndParticles();

    Mesh marchingCubeMesh;
    GLuint dataFieldTex;
    Mesh heightmapMesh;
    GLuint heightmapFieldTex, biomeFieldTex;
    GLuint edgeTableTex, triTableTex;

    GLuint allBiomesColorTextures, allBiomesNormalTextures, allBiomesDisplacementTextures;

    // TODO : Transform this into a "particle" system
    std::vector<Mesh> possibleRocks;
    std::vector<std::tuple<int, Vector3, float>> rocksIndicesAndPositionAndSize;
    int numberOfRocksDisplayed = 0;
    std::vector<Mesh> possibleCorals;
    std::vector<std::tuple<int, Vector3, float>> coralsIndicesAndPositionAndSize;

    Vector3 voxelGridOffset = Vector3(0, 0, 0);
    float voxelGridScaling = 1.f;

    Mesh particlesMesh;
    std::vector<Vector3> randomParticlesPositions;


    std::chrono::system_clock::time_point startingTime;
    size_t previousHistoryIndex = -1;
};

#endif // TERRAINGENERATIONINTERFACE_H
