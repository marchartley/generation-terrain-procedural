#ifndef ENVOBJSINTERFACE_H
#define ENVOBJSINTERFACE_H

#include <QWidget>
#include "Interface/ActionInterface.h"
#include "Interface/HierarchicalListWidget.h"
#include "Utils/HotreloadFile.h"

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
    void setMaterialsDefinitionFile(std::string filename);
    void setDefinitionFile(std::string filename);
    void setTransformationsFile(std::string filename);

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
    void destroyEnvObject(EnvObject* object);

    void displayProbas(std::string objectName);
//    void displaySedimentsDistrib();
//    void displayPolypDistrib();
    void displayMaterialDistrib(std::string materialName);
    void displayFlowfieldAsImage();

    void updateObjectsList();

    void updateObjectsListSelection(QListWidgetItem* newSelectionItem);
    void updateSelectionMesh();

    void updateObjectsDefinitions(const std::string& newDefinition);
    void updateMaterialsDefinitions(const std::string& newDefinition);
    void updateMaterialsTransformationsDefinitions(const std::string& newDefinition);

    // void hotReloadFile();

    void evaluateAndDisplayCustomCostFormula(std::string formula) const;

    BSpline computeNewObjectsCurveAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float widthMaxLength);
    ShapeCurve computeNewObjectsShapeAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float widthMaxLength);

    void runPerformanceTest();

    void resetScene();

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

    std::string currentlyPreviewedObject;


    GridF erosionGrid;
    GridV3 velocitiesGrid;

//    std::string primitiveDefinitionFile;
//    std::string materialsDefinitionFile;

    HotreloadFile primitiveDefinitionFile;
    HotreloadFile materialsDefinitionFile;
    HotreloadFile transformationsFile;
//    QDateTime lastTimeFileHasBeenModified;

    std::map<EnvObject*, ImplicitPatch*> implicitPatchesFromObjects;
    ImplicitNaryOperator* rootPatch;

    EnvObject* currentSelection = nullptr;

    std::string previousFileContent = "";
    std::string previousMaterialsFileContent = "";
    std::string previousMaterialsTransformationsFileContent = "";
};

BSpline followIsovalue(const GridF &values, const GridV3& gradients, const Vector3& startPoint, float maxDist);
BSpline followGradient(const GridV3 gradients, const Vector3& startPoint, float maxDist, bool followInverse = false) ;
std::vector<Vector3> findCandidatesPositions(const Vector3& startPosition, const Vector3& direction, float angle, float radius, int nbCandidates);
std::vector<BSpline> getCandidatesPaths(const GridV3& gradients, const std::vector<Vector3>& positions, float directionLength);
BSpline getBestCandidatesPath(const GridF &score, const BSpline& initialPath, const std::vector<BSpline>& paths);

#endif // ENVOBJSINTERFACE_H
