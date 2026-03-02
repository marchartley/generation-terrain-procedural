#ifndef ENVOBJSINTERFACE_H
#define ENVOBJSINTERFACE_H

#include <QWidget>
#include "Interface/ActionInterface.h"
#include "GUIElements/HierarchicalListWidget.h"
#include "Utils/HotreloadFile.h"

// #include "EnvObject/EnvPoint.h"
// #include "EnvObject/EnvCurve.h"
// #include "EnvObject/EnvArea.h"
#include "EnvObject/EnvironmentalScene.h"

class EnvObjsInterface : public ActionInterface
{
    Q_OBJECT
public:
    EnvObjsInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3::invalid);

    void replay(nlohmann::json action);

    QLayout* createGUI();

    std::tuple<GridF, GridV3> extractErosionDataOnTerrain();

    void createEnvObjectsFromImplicitTerrain();
    void setMaterialsDefinitionFile(const std::string& filename);
    void setDefinitionFile(const std::string& filename);
    void setTransformationsFile(const std::string& filename);
    void setScenarioFile(const std::string& filename);

    EnvObjectInstance* instantiateObjectAtBestPosition(const std::string& objectName, Vector3 position, const GridF& score);
    EnvObjectInstance* instantiateObjectAtBestPositionWithoutScoreMap(const std::string& objectName, Vector3 position, const Vector3 &maxPos);
    EnvObjectInstance* instantiateObjectUsingSpline(const std::string& objectName, const BSpline& spline);

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
    EnvObjectInstance* instantiateSpecific(const std::string& objectName, const Vector3& targetPosition = Vector3::invalid, const GridF& score = GridF(), bool waitForFullyGrown = true, bool updateScreen = false, bool updateEnvironmentDirectly = true);
    EnvObjectInstance* fakeInstantiate(const std::string& objectName, const GridF& score = GridF());

    bool checkIfObjectShouldDie(EnvObjectInstance* obj, float limitFactorForDying = .2f);

    void fromGanUI();

    void recomputeErosionValues();

    void runNextStep();
    void runScenario();

    void updateEnvironmentFromEnvObjects(bool updateImplicitTerrain = false, bool emitUpdateSignal = true, bool killObjectsIfPossible = true);
    void updateUntilStabilization();
    void destroyEnvObject(EnvObjectInstance* object, bool applyDying = true, bool recomputeTerrainPropertiesForObject = true);

    void displayProbas(const std::string& objectName);
    void displayMaterialDistrib(const std::string& materialName);

    void manualModificationOfFocusArea();
    void manualModificationOfFlowfield();
    void resetFlowfield();

    void updateObjectsList();

    void updateObjectsListSelection(QListWidgetItem* __newSelectionItem = nullptr);
    void updateSelectionMesh();
    void updateNewObjectMesh();

    void updateObjectsDefinitions(const std::string& newDefinition);
    void updateMaterialsDefinitions(const std::string& newDefinition);
    void updateMaterialsTransformationsDefinitions(const std::string& newDefinition);
    void updateScenarioDefinition(const std::string& newDefinition);

    void evaluateAndDisplayCustomFitnessFormula(const std::string& formula);
    void evaluateAndDisplayCustomFittingFormula(const std::string& formula);
    void evaluateAndDisplayCustomFitnessAndFittingFormula(const std::string& fitnessFuncFormula, std::string fittingFuncFormula);

    BSpline computeNewObjectsCurveAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float widthMaxLength, bool followIsolevel = false);
    ShapeCurve computeNewObjectsShapeAtPosition(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength);
    ShapeCurve computeNewObjectsShapeAtPositionForceCircle(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength);
    ShapeCurve computeNewObjectsShapeAtPositionForceCircleOptimizedArea(const Vector3& seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float targetArea);

    void runPerformanceTest();

    void resetScene();

    void loadScene(const std::string& filename);
    void saveScene(const std::string& filename);

    void previewCurrentEnvObjectPlacement(const Vector3& position);
    void previewFlowEdition(const Vector3& mousePos, const Vector3& brushDir);
    void previewMaterialEdition(const Vector3& position, bool addingMaterial);

    void showAllElementsOnPlotter();

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

    void saveForRenders();

    StatsValues displayStatsForObjectCreation(const std::string& objectName, int nbSamples = 10);

    GridV3 computeUserKelvinletField() const;

public:
    Mesh velocitiesMesh;
    Mesh highErosionsMesh;
    Mesh highDepositionMesh;
    Mesh selectedObjectsMesh;
    Mesh newObjectMesh;

    GridV3 userFlowField;
    GridV3 simulationFlowField;
    std::vector<Kelvinlet*> userKelvinlets;
    std::string KelvinletChoice = "grab";
    bool displayFlow = false;

    HierarchicalListUI* objectsListWidget = nullptr;

    Vector3 draggingPoint = Vector3::invalid;
    Vector3 draggingFullObject = Vector3::invalid;
    Vector3 draggingHasBeenApplied = Vector3::invalid;


    Vector3 kelvinletDraggingPoint = Vector3::invalid;
    Vector3 kelvinletDraggingFullObject = Vector3::invalid;
    Vector3 kelvinletDraggingHasBeenApplied = Vector3::invalid;

    bool displayVelocities = true;
    bool displayHighErosions = true;
    bool displaySediments = true;
    bool displayHighCurrents = true;
    bool waitAtEachFrame = false;
    bool displayGrooves = false;

    float flowErosionFactor = 0.f;

    // std::string currentlyPreviewedObject;

    bool materialSimulationStable = false;

    std::string currentMaterialEdited = "";


    GridF erosionGrid;
    GridV3 velocitiesGrid;

    HotreloadFile primitiveDefinitionFile;
    HotreloadFile materialsDefinitionFile;
    HotreloadFile transformationsFile;
    HotreloadFile scenarioFile;

    std::map<EnvObjectInstance*, ImplicitPatch*> implicitPatchesFromObjects;
    ImplicitNaryOperator* rootPatch;
    // Implicit2DNary* rootPatch;

    std::vector<EnvObjectInstance*> currentSelections;

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
    bool previewingObjectInPlotter = false;

    bool forceScenarioInterruption = true;

    bool fluidSimulationIsStable = false;

    GridF initialHeightmap;
    GridF subsidedHeightmap;




    GridF groundConstraintedHeights;
    GridF waterConstraintedHeights;
    GridF surfaceHeights;

    bool displayDepositionOnHeightmap = true;

    std::string testedFitnessFunction = "";
    std::string testedFittingFunction = "";



    std::shared_ptr<EnvironmentalScene> scene;
};

BSpline followIsovalue(const GridF &values, const GridV3& gradients, const Vector3& startPoint, float maxDist);
BSpline followGradient(const GridV3 gradients, const Vector3& startPoint, float maxDist, bool followInverse = false) ;
std::vector<Vector3> findCandidatesPositions(const Vector3& startPosition, const Vector3& direction, float angle, float radius, int nbCandidates);
std::vector<BSpline> getCandidatesPaths(const GridV3& gradients, const std::vector<Vector3>& positions, float directionLength);
BSpline getBestCandidatesPath(const GridF &score, const BSpline& initialPath, const std::vector<BSpline>& paths);

#endif // ENVOBJSINTERFACE_H
