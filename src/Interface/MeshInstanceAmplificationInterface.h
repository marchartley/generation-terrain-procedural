#ifndef MESHINSTANCEAMPLIFICATIONINTERFACE_H
#define MESHINSTANCEAMPLIFICATIONINTERFACE_H

#include "ActionInterface.h"

class MeshInstanceAmplificationInterface : public ActionInterface
{
    Q_OBJECT
public:
    MeshInstanceAmplificationInterface(QWidget* parent = nullptr);

    void display();
    void replay(nlohmann::json action); // Nothing to replay
    void reloadShaders();


    QLayout* createGUI();

    std::vector<std::pair<Vector3, Vector3>> getAvailablePositionsForMaterial(TerrainTypes target);
    std::vector<std::pair<Vector3, Vector3>> getCoralAvailablePositions();
    std::vector<std::pair<Vector3, Vector3>> getRocksAvailablePositions();

public Q_SLOTS:
    void setCoralsDisplayed(bool display);
    void setRocksDisplayed(bool display);

    void afterTerrainUpdated();
    void regenerateRocksPositions();

public:
    bool displayCorals = false;
    bool displayRocks = false;

    int numberOfRocksDisplayed = 1000;
    int numberOfCoralsDisplayed = 1000;

    int numberOfLoadedRocks = 0; // Set to -1 to get all
    int numberOfLoadedCorals = 0;
    std::vector<Mesh> possibleRocks;
    std::vector<std::tuple<int, Vector3, float>> rocksIndicesAndPositionAndSize;
    std::vector<Mesh> possibleCorals;
    std::vector<std::tuple<int, Vector3, float>> coralsIndicesAndPositionAndSize;

    size_t previousHistoryIndex = 0;
};

#endif // MESHINSTANCEAMPLIFICATIONINTERFACE_H
