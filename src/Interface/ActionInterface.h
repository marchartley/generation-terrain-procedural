#ifndef ACTIONINTERFACE_H
#define ACTIONINTERFACE_H

#include "Interface/CustomInteractiveObject.h"
#include "TerrainModification/TerrainAction.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"

class ActionInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    ActionInterface(std::string actionTypeName, QWidget* parent = nullptr);
//    ~ActionInterface();

    virtual void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid) {
        this->voxelGrid = voxelGrid;
    }
    virtual void affectHeightmap(std::shared_ptr<Heightmap> heightmap) {
        this->heightmap = heightmap;
    }
    virtual void affectLayerGrid(std::shared_ptr<LayerBasedGrid> layerGrid) {
        this->layerGrid = layerGrid;
    }
    virtual void affectImplicitTerrain(std::shared_ptr<ImplicitNaryOperator> implicitPatch) {
        this->implicitTerrain = implicitPatch;
    }

    virtual void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr) {
        this->affectHeightmap(heightmap);
        this->affectVoxelGrid(voxelGrid);
        this->affectLayerGrid(layerGrid);
        this->affectImplicitTerrain(implicitPatch);
    }

    virtual QLayout* createGUI() {
        return nullptr;
    }
    virtual void display(const Vector3& camPos = Vector3(false)) {

    }
    virtual void reloadShaders() {

    }

    void addTerrainAction(nlohmann::json actionParameters) {
        jsonActionsHistory->push_back(nlohmann::json({
                                                         {"type", this->actionType},
                                                         {"parameters", actionParameters}
                                                     }));
//        saveAllActions();
    }

    void affectSavingFile(std::shared_ptr<std::vector<nlohmann::json>> history, std::shared_ptr<std::fstream> file, std::string filename) {
        this->jsonActionsHistory = history;
        this->savingFile = file;
        this->savingFilename = filename;
    }

    void saveAllActions(std::string filename = "") {
        savingFile->close();
        savingFile->open((filename.empty() ? savingFilename : filename), std::fstream::in | std::fstream::out | std::fstream::trunc);
        *savingFile << nlohmann::json({{"actions", *jsonActionsHistory}}).dump(1, '\t');
        savingFile->flush();
    }

    bool isConcerned(nlohmann::json& action) { return action.contains("type") && action.at("type").get<std::string>() == this->actionType; }

    virtual void replay(nlohmann::json action) = 0;
    std::string actionType;
    std::string savingFilename;
    std::shared_ptr<std::fstream> savingFile;
    std::shared_ptr<std::vector<nlohmann::json>> jsonActionsHistory;

    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Heightmap> heightmap;
    std::shared_ptr<LayerBasedGrid> layerGrid;
    std::shared_ptr<ImplicitNaryOperator> implicitTerrain = nullptr;

Q_SIGNALS:
    void updated();
    void terrainUpdated();

public Q_SLOTS:
    virtual void afterTerrainUpdated() {

    }
    virtual void mouseClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model) {

    }
    virtual void mouseDoubleClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model) {

    }
    virtual void mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model) {

    }
};

#endif // ACTIONINTERFACE_H
