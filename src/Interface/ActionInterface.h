#ifndef ACTIONINTERFACE_H
#define ACTIONINTERFACE_H

class ActionInterface;

#include "Interface/CustomInteractiveObject.h"
#include "TerrainModification/TerrainAction.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "Interface/Viewer.h"

class ActionInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    ActionInterface(std::string actionTypeName,
                    std::string interfaceName,
                    std::string interfaceType,
                    std::string mainActionDescription,
                    std::string mainActionButtonLogo,
                    QWidget* parent = nullptr);
//    ~ActionInterface();

    virtual void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    virtual void affectHeightmap(std::shared_ptr<Heightmap> heightmap);
    virtual void affectLayerGrid(std::shared_ptr<LayerBasedGrid> layerGrid);
    virtual void affectImplicitTerrain(std::shared_ptr<ImplicitNaryOperator> implicitPatch);

    virtual void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    virtual QLayout* createGUI();
    virtual void display(const Vector3& camPos = Vector3(false));
    virtual void reloadShaders();

    void addTerrainAction(nlohmann::json actionParameters);

    void affectSavingFile(std::shared_ptr<std::vector<nlohmann::json>> history, std::shared_ptr<std::fstream> file, std::string filename);

    void saveAllActions(std::string filename = "");

    bool isConcerned(nlohmann::json& action);

    void log(const std::string& message, bool verbose = true);
    void error(const std::string& message, bool verbose = true);

    virtual void replay(nlohmann::json action) = 0;

    std::string actionType;
    std::string interfaceName;
    std::string interfaceType;
    std::string mainActionDescription;
    std::string mainActionButtonLogo;


    std::string savingFilename;
    std::shared_ptr<std::fstream> savingFile;
    std::shared_ptr<std::vector<nlohmann::json>> jsonActionsHistory;

    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Heightmap> heightmap;
    std::shared_ptr<LayerBasedGrid> layerGrid;
    std::shared_ptr<ImplicitNaryOperator> implicitTerrain = nullptr;

    Viewer* viewer = nullptr;

    std::shared_ptr<ActionInterface> findOtherInterface(std::string name) const;

Q_SIGNALS:
    void updated();
    void terrainUpdated();
    void waterLevelChanged(float newLevel);

public Q_SLOTS:
    virtual void afterTerrainUpdated();
    virtual void afterWaterLevelChanged();
    virtual void mouseClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    virtual void mouseDoubleClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    virtual void mouseReleasedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    virtual void mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model);
};

#endif // ACTIONINTERFACE_H
