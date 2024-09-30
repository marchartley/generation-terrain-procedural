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
    void setScenarioFile(std::string filename);

    EnvObject* instantiateObjectAtBestPosition(std::string objectName, Vector3 position, const GridF& score);
    EnvObject* instantiateObjectUsingSpline(std::string objectName, const BSpline& spline);

public Q_SLOTS:
    void show();
    void hide();
    virtual void afterTerrainUpdated();
    virtual void afterWaterLevelChanged();

    virtual void mouseClickedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    virtual void mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model);
    virtual void mouseReleasedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    virtual void keyPressEvent(QKeyEvent* event);

public:
    EnvObject* instantiateSpecific(std::string objectName, GridF score = GridF(), bool waitForFullyGrown = true, bool updateScreen = false);
    EnvObject* fakeInstantiate(std::string objectName, GridF score = GridF());

    bool checkIfObjectShouldDie(EnvObject* obj, float limitFactorForDying = .2f);

    void fromGanUI();

    void recomputeErosionValues();

    void runNextStep();
    void runScenario();

    void updateEnvironmentFromEnvObjects(bool updateImplicitTerrain = false, bool emitUpdateSignal = true, bool killObjectsIfPossible = true);
    void onlyUpdateFlowAndSandFromEnvObjects();
    void destroyEnvObject(EnvObject* object, bool applyDying = true, bool recomputeTerrainPropertiesForObject = true);

    void displayProbas(std::string objectName);
    void displayMaterialDistrib(std::string materialName);
    void displayFlowfieldAsImage();

    void manualModificationOfFocusArea();
    void manualModificationOfFlowfield();

    void updateObjectsList();

    void updateObjectsListSelection(QListWidgetItem* __newSelectionItem = nullptr);
    void updateSelectionMesh();
    void updateNewObjectMesh();

    void updateObjectsDefinitions(const std::string& newDefinition);
    void updateMaterialsDefinitions(const std::string& newDefinition);
    void updateMaterialsTransformationsDefinitions(const std::string& newDefinition);
    void updateScenarioDefinition(const std::string& newDefinition);

    void evaluateAndDisplayCustomFitnessFormula(std::string formula);
    void evaluateAndDisplayCustomFittingFormula(std::string formula);

    BSpline computeNewObjectsCurveAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float widthMaxLength, bool followIsolevel = false);
    ShapeCurve computeNewObjectsShapeAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength);
    ShapeCurve computeNewObjectsShapeAtPositionForceCircle(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength);
    ShapeCurve computeNewObjectsShapeAtPositionForceCircleOptimizedArea(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float targetArea);

    void runPerformanceTest();

    void resetScene();

    void loadScene(std::string filename);
    void saveScene(std::string filename);

    GridV3 renderFocusArea() const;
    GridV3 renderFlowfield() const;

    void previewCurrentEnvObjectPlacement(Vector3 position);
    void previewFocusAreaEdition(Vector3 mousePos, bool addingFocus);
    void previewFlowEdition(Vector3 mousePos, Vector3 brushDir);

    void showAllElementsOnPlotter() const;

    void addObjectsHeightmaps();
    void flowErosionSimulation();

    void startNewObjectCreation();
    void addPointOnNewObjectCreation(const Vector3& position, bool addPoint = true, float removeRadius = 2.f);
    void endNewObjectCreation();

    void startDraggingObject(const Vector3& position, bool singleVertexMoved);
    void moveDraggedObject(const Vector3& position);
    void endDraggingObject(bool destroyObjects);

    std::string getCurrentObjectName() const;

    void updateVectorFieldVisu();

    StatsValues displayStatsForObjectCreation(std::string objectName, int nbSamples = 10);

public:
    Mesh velocitiesMesh;
    Mesh highErosionsMesh;
    Mesh highDepositionMesh;
    Mesh selectedObjectsMesh;
    Mesh newObjectMesh;

    GridV3 userFlowField;
    GridV3 simulationFlowField;
    bool displayFlow = false;

    HierarchicalListWidget* objectsListWidget = nullptr;

    Vector3 draggingPoint = Vector3(false);
    Vector3 draggingFullObject = Vector3(false);
    Vector3 draggingHasBeenApplied = Vector3(false);

    bool displayVelocities = true;
    bool displayHighErosions = true;
    bool displaySediments = true;
    bool displayHighCurrents = true;
    bool waitAtEachFrame = false;
    bool displayGrooves = false;

    float flowErosionFactor = 0.f;

    // std::string currentlyPreviewedObject;

    bool materialSimulationStable = false;


    GridF erosionGrid;
    GridV3 velocitiesGrid;

    HotreloadFile primitiveDefinitionFile;
    HotreloadFile materialsDefinitionFile;
    HotreloadFile transformationsFile;
    HotreloadFile scenarioFile;

    std::map<EnvObject*, ImplicitPatch*> implicitPatchesFromObjects;
    ImplicitNaryOperator* rootPatch;
    // Implicit2DNary* rootPatch;

    std::vector<EnvObject*> currentSelections;

    std::string previousFileContent = "";
    std::string previousMaterialsFileContent = "";
    std::string previousMaterialsTransformationsFileContent = "";
    std::string previousScenarioFileContent = "";

    GridF focusedArea;

    bool focusAreaEditing = false;
    ComboboxElement* objectCombobox;

    bool flowfieldEditing = false;

    BSpline objectSkeletonCreation;
    bool manuallyCreatingObject = false;

    bool forceScenarioInterruption = true;

    bool fluidSimulationIsStable = false;

    GridF initialHeightmap;
    GridF subsidedHeightmap;

    bool displayDepositionOnHeightmap = true;
};

BSpline followIsovalue(const GridF &values, const GridV3& gradients, const Vector3& startPoint, float maxDist);
BSpline followGradient(const GridV3 gradients, const Vector3& startPoint, float maxDist, bool followInverse = false) ;
std::vector<Vector3> findCandidatesPositions(const Vector3& startPosition, const Vector3& direction, float angle, float radius, int nbCandidates);
std::vector<BSpline> getCandidatesPaths(const GridV3& gradients, const std::vector<Vector3>& positions, float directionLength);
BSpline getBestCandidatesPath(const GridF &score, const BSpline& initialPath, const std::vector<BSpline>& paths);

#endif // ENVOBJSINTERFACE_H
