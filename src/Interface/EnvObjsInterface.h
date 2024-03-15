#ifndef ENVOBJSINTERFACE_H
#define ENVOBJSINTERFACE_H

#include <QWidget>
#include "Interface/ActionInterface.h"
#include "Interface/HierarchicalListWidget.h"

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

class EnvObjsInterface : public ActionInterface
{
    Q_OBJECT
public:
    EnvObjsInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

    std::tuple<GridF, GridV3> extractErosionDataOnTerrain();

    void createEnvObjectsFromImplicitTerrain();
    void setDefinitionFile(std::string filename);

    EnvObject* instantiateObjectAtBestPosition(std::string objectName, const Vector3& position, const GridF& score);

public Q_SLOTS:
    void show();
    void hide();
    virtual void afterTerrainUpdated();
    virtual void afterWaterLevelChanged();

    virtual void mouseClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);

    void instantiateObject();
    void instantiateSpecific(std::string objectName);
    void fromGanUI();

    void recomputeErosionValues();

    void updateEnvironmentFromEnvObjects(bool updateImplicitTerrain = false, bool emitUpdateSignal = true);
    void onlyUpdateFlowAndSandFromEnvObjects();

    void displayProbas(std::string objectName);
    void displaySedimentsDistrib();
    void displayPolypDistrib();
    void displayFlowfieldAsImage();

    void updateObjectsList();

    void updateObjectsListSelection(QListWidgetItem* newSelectionItem);

    void hotReloadFile();

    void evaluateAndDisplayCustomCostFormula(std::string formula) const;

    BSpline computeNewObjectsCurveAtPosition(const Vector3& seedPosition, const GridV3 &gradients, float directionLength, float widthMaxLength);
    ShapeCurve computeNewObjectsShapeAtPosition(const Vector3& seedPosition, const GridV3 &gradients, float directionLength, float widthMaxLength);

    void runPerformanceTest();

public:
    Mesh velocitiesMesh;
    Mesh highErosionsMesh;
    Mesh highDepositionMesh;
    Mesh objectsMesh;

    HierarchicalListWidget* objectsListWidget = nullptr;

    bool displayVelocities = true;
    bool displayHighErosions = true;
    bool displaySediments = true;
    bool displayHighCurrents = true;
    bool waitAtEachFrame = false;


    GridF erosionGrid;
    GridV3 velocitiesGrid;

    std::string primitiveDefinitionFile;
    QDateTime lastTimeFileHasBeenModified;

    std::map<EnvObject*, ImplicitPatch*> implicitPatchesFromObjects;
    Implicit2DNary* rootPatch;

    EnvObject* currentSelection = nullptr;
};


#endif // ENVOBJSINTERFACE_H
