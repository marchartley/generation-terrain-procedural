#ifndef PRIMITIVEPATCHESINTERFACE_H
#define PRIMITIVEPATCHESINTERFACE_H

class ImplicitPatch;
class PatchReplacementDialog;

#include "ActionInterface.h"
#include "Interface/HierarchicalListWidget.h"
#include "Interface/ControlPoint.h"

/*
enum PRIMITIVE_SHAPE {
    SPHERE = 0,
    BLOCK = 1,
    GAUSSIAN = 2
};
enum PATCH_OPERATION {
    STACK = 0,
    BLEND = 1,
    REPLACE = 2
};*/

class PrimitivePatchesInterface : public ActionInterface
{
    Q_OBJECT
public:
    PrimitivePatchesInterface(QWidget *parent = nullptr);

    void display();
    void reloadShaders();

    void replay(nlohmann::json action);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    QLayout* createGUI();

    void show();
    void hide();

public Q_SLOTS:
    void mouseMovedOnMap(Vector3 newPosition);
    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event);

    void setSelectedShape(ImplicitPatch::PredefinedShapes newShape, Vector3 newPosition = Vector3());
    void setCurrentOperation(ImplicitPatch::CompositionFunction newOperation);

    void resetPatch();

    void updateMapWithCurrentPatch();

    void setSelectedWidth(float newVal);
    void setSelectedHeight(float newVal);
    void setSelectedDepth(float newVal);
    void setSelectedSigma(float newVal);
    void setSelectedBlendingFactor(float newVal);

    void savePatchesAsFile(std::string filename);
    void loadPatchesFromFile(std::string filename);
    void hotReloadFile();

    void updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous);
    void openPrimitiveModificationDialog(QListWidgetItem* item);
    void modifyPrimitiveHierarchy(int ID_to_move, int relatedID, HIERARCHY_TYPE relation, QDropEvent* event);

protected:
    ImplicitPatch* createPatchFromParameters(Vector3 position, ImplicitPatch* replacedPatch = nullptr);
    ImplicitPatch* createOperationPatchFromParameters(ImplicitPatch* composableA = nullptr, ImplicitPatch* composableB = nullptr, ImplicitPatch* replacedPatch = nullptr);
    ImplicitPatch* findPrimitiveById(int ID);
    ImplicitPatch* naiveApproachToGetParent(ImplicitPatch* child);

    QDateTime lastTimeFileHasBeenModified;

    ImplicitPatch::PredefinedShapes currentShapeSelected = ImplicitPatch::PredefinedShapes::Sphere;
    ImplicitPatch::CompositionFunction currentOperation = ImplicitPatch::CompositionFunction::BLEND;
    float selectedDensity = 1.5f;

    float selectedWidth = 20.f;
    float selectedDepth = 20.f;
    float selectedHeight = 20.f;
    float selectedSigma = 5.f;
    float selectedBlendingFactor = 1.f;

    Mesh previewMesh;
    Mesh patchAABBoxMesh;

    ImplicitPatch* mainPatch;

    Patch3D* currentPatch;
//    std::vector<Patch3D*> storedPatches;
    std::vector<ImplicitPatch*> storedPatches;

    std::function<float(Vector3)> currentSelectionFunction;

    std::function<float(Vector3)> blockFunction(float height);
    std::function<float(Vector3)> sphereFunction(float radius);
    std::function<float(Vector3)> gaussianFunction(float sigma, float height);

    Vector3 currentPos;

    Vector3 functionSize = Vector3(1, 1, 1);
    void updateFunctionSize();


    HierarchicalListWidget* primitiveSelectionGui;
    void updatePrimitiveList();
    void cleanPatch(ImplicitPatch* patch);

    friend class PatchReplacementDialog;

    std::unique_ptr<ControlPoint> primitiveControlPoint;

    ImplicitPatch* currentlyManipulatedPatch = nullptr;

    std::string mainFilename;
    bool enableHotReloading = false;

    void displayDebuggingVoxels();
    Matrix3<float> debuggingVoxels;
    Mesh debuggingVoxelsMesh;
    bool debugMeshDisplayed = false;
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
    HierarchicalListWidget* allPrimitivesList;
    QPushButton* cancelButton;
    QPushButton* validButton;
    PrimitivePatchesInterface* caller = nullptr;
//    int selectedPrimitiveIndex = -1;
    ImplicitPatch* patch = nullptr;
};

#endif // PRIMITIVEPATCHESINTERFACE_H
