#ifndef PRIMITIVEPATCHESINTERFACE_H
#define PRIMITIVEPATCHESINTERFACE_H

class ImplicitPatch;
class PatchReplacementDialog;

#include "ActionInterface.h"
#include "Interface/HierarchicalListWidget.h"
#include "Interface/ControlPoint.h"


class PrimitivePatchesInterface : public ActionInterface
{
    Q_OBJECT
public:
    PrimitivePatchesInterface(QWidget *parent = nullptr);

    void display();
    void reloadShaders();

    void replay(nlohmann::json action);

//    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid);

    QLayout* createGUI();

    void show();
    void hide();

    void keyPressEvent(QKeyEvent *e);

    ImplicitPatch* selectedPatch();

public Q_SLOTS:
    void mouseMovedOnMapEvent(Vector3 newPosition, TerrainModel *model);
    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event, TerrainModel *model);
    void createPatchWithOperation(Vector3 pos);

    void setSelectedShape(ImplicitPatch::PredefinedShapes newShape, Vector3 newPosition = Vector3());
    void setCurrentOperation(ImplicitPatch::CompositionFunction newOperation);
    void setCurrentPositioning(ImplicitPatch::PositionalLabel positioning);

    void resetPatch();

    void updateMapWithCurrentPatch();

    void setSelectedWidth(float newVal);
    void setSelectedHeight(float newVal);
    void setSelectedDepth(float newVal);
    void setSelectedSigma(float newVal);
    void setSelectedBlendingFactor(float newVal);

    void addNoiseOnSelectedPatch();
    void addDistortionOnSelectedPatch();
    void addSpreadOnSelectedPatch();

    void savePatchesAsFile(std::string filename);
    void loadPatchesFromFile(std::string filename);
    void hotReloadFile();

    void updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous);
    void openPrimitiveModificationDialog(QListWidgetItem* item);
    void modifyPrimitiveHierarchy(int ID_to_move, int relatedID, HIERARCHY_TYPE relation, QDropEvent* event);

    void openFileForNewPatch();

    void setMainFilename(std::string newMainFilename);
    void findAllSubfiles();

    void loadTransformationRules();

    void addParametricPoint(Vector3 point);
    void displayParametricCurve();

    void displayPatchesTree();

protected:
    ImplicitPatch* createPatchFromParameters(Vector3 position, ImplicitPatch* replacedPatch = nullptr);
    ImplicitPatch* createOperationPatchFromParameters(ImplicitPatch* composableA = nullptr, ImplicitPatch* composableB = nullptr, ImplicitOperator* replacedPatch = nullptr);
    ImplicitPatch* findPrimitiveById(int ID);
    ImplicitPatch* naiveApproachToGetParent(ImplicitPatch* child);

    QDateTime lastTimeFileHasBeenModified;

    ImplicitPatch::PredefinedShapes currentShapeSelected = ImplicitPatch::PredefinedShapes::Sphere;
    ImplicitPatch::CompositionFunction currentOperation = ImplicitPatch::CompositionFunction::BLEND;
    ImplicitPatch::PositionalLabel currentPositioning = ImplicitPatch::PositionalLabel::ABOVE;
//    float selectedDensity = 1.5f;
    TerrainTypes selectedTerrainType = TerrainTypes::ROCK;

    float selectedWidth = 20.f;
    float selectedDepth = 20.f;
    float selectedHeight = 20.f;
    float selectedSigma = 5.f;
    float selectedBlendingFactor = 2.f;
    bool applyIntersection = false;

    Mesh previewMesh;
    Mesh patchAABBoxMesh;

    ImplicitPatch* mainPatch;
    ImplicitPatch* desiredPatchFromFile = nullptr;
    std::string desiredPatchFilename = "";

//    Patch3D* currentPatch;
    std::vector<ImplicitPatch*> storedPatches;

//    std::function<float(Vector3)> currentSelectionFunction;

    Vector3 currentPos;

    Vector3 functionSize = Vector3(1, 1, 1);
    void updateFunctionSize();


    HierarchicalListWidget* primitiveSelectionGui;
    void updatePrimitiveList();
    void cleanPatch(ImplicitPatch* patch);

    friend class PatchReplacementDialog;

    std::unique_ptr<ControlPoint> primitiveControlPoint;

    ImplicitPatch* currentlyManipulatedPatch = nullptr;
    ImplicitPatch* currentlySelectedPatch = nullptr;

    std::string mainFilename;
    std::vector<std::string> allSubfiles;
    std::string rulesFilename = "saved_maps/transformation_rules.txt";
    bool enableHotReloading = false;

    void displayDebuggingVoxels();
    Matrix3<float> debuggingVoxels;
    Mesh debuggingVoxelsMesh;
    bool debugMeshDisplayed = false;

    int nbPrimitives = 0;
    int nbUnaryOperators = 0;
    int nbBinaryOperators = 0;

    BSpline parametricCurve;
    Mesh parametricCurveMesh;
};








class PatchReplacementDialog : public QDialog {
    Q_OBJECT
public:
    PatchReplacementDialog(PrimitivePatchesInterface *caller = nullptr, ImplicitPatch* patchToModify = nullptr);


public Q_SLOTS:
    void open();
    void cancel();
    void confirm();

public:
    QPushButton* cancelButton;
    QPushButton* validButton;
    PrimitivePatchesInterface* caller = nullptr;
    ImplicitPatch* patch = nullptr;
};

#endif // PRIMITIVEPATCHESINTERFACE_H
