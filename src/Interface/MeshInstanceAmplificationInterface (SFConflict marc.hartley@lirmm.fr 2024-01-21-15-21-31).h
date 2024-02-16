#ifndef MESHINSTANCEAMPLIFICATIONINTERFACE_H
#define MESHINSTANCEAMPLIFICATIONINTERFACE_H

#include "ActionInterface.h"

struct InstantiationMeshOption {
    InstantiationMeshOption(std::string name, std::string folderName, std::pair<float, float> minMaxSizes, std::vector<float> color, const Vector3& translation = Vector3())
        : name(name), minMaxSizes(minMaxSizes), color(color), folderName(folderName), requiredTranslation(translation)
    {}
    InstantiationMeshOption(std::string name, std::pair<float, float> minMaxSizes, std::vector<float> color, const Vector3& translation = Vector3())
        : InstantiationMeshOption(name, name, minMaxSizes, color, translation)
    {}
    std::string name;
    std::string folderName;
    bool displayed = false;
    int numberDisplayed = 10000;
    int numberOfLoadedMesh = -1;
    std::vector<Mesh> possibleMeshes;
    std::vector<std::tuple<int, Vector3, float>> indicesAndPositionsAndSizes;
    std::pair<float, float> minMaxSizes = {10.f, 15.f};
    std::vector<float> color = {.0f, 1.f, .5f, 1.f};
    Vector3 requiredTranslation;
};

class MeshInstanceAmplificationInterface : public ActionInterface
{
    Q_OBJECT
public:
    MeshInstanceAmplificationInterface(QWidget* parent = nullptr);

    void display(const Vector3& camPos = Vector3(false));
    void replay(nlohmann::json action); // Nothing to replay
    void reloadShaders();


    QLayout* createGUI();

    std::vector<AABBox> getAvailablePositionsForMaterial(TerrainTypes target);
    std::vector<AABBox> getCoralAvailablePositions();
    std::vector<AABBox> getRocksAvailablePositions();
    std::vector<std::pair<Vector3, float> > getPositionsFor(std::string type);

public Q_SLOTS:
    void setCoralsDisplayed(bool display);
    void setRocksDisplayed(bool display);
    void setDisplayingType(InstantiationMeshOption &options, bool display);

    void afterTerrainUpdated();
    void regenerateRocksPositions();
    void regenerateAllTypePositions();

public:
    std::vector<InstantiationMeshOption> meshesOptions;
    /*std::map<std::string, bool> displayType;
    std::map<std::string, int> numberOfTypeDisplayed;
    std::map<std::string, int> numberOfLoadedModel;
    std::map<std::string, std::vector<Mesh>> possibleMeshType;
    std::map<std::string, std::vector<std::tuple<int, Vector3, float>>> indicesAndPositionAndSize;*/

    bool displayCorals = false;
    bool displayRocks = false;

    int numberOfRocksDisplayed = 1000;
    int numberOfCoralsDisplayed = 1000;

    int numberOfLoadedRocks = -1; // Set to -1 to get all
    int numberOfLoadedCorals = -1;
    std::vector<Mesh> possibleRocks;
    std::vector<std::tuple<int, Vector3, float>> rocksIndicesAndPositionAndSize;
    std::vector<Mesh> possibleCorals;
    std::vector<std::tuple<int, Vector3, float>> coralsIndicesAndPositionAndSize;

    size_t previousHistoryIndex = 0;
};

#endif // MESHINSTANCEAMPLIFICATIONINTERFACE_H
