#ifndef PRIMITIVEPATCHESINTERFACE_H
#define PRIMITIVEPATCHESINTERFACE_H

#include "ActionInterface.h"

enum PRIMITIVE_SHAPE {
    SPHERE = 0,
    BLOCK = 1,
    GAUSSIAN = 2
};
enum PATCH_OPERATION {
    STACK = 0,
    BLEND = 1,
    REPLACE = 2
};

class PrimitivePatchesInterface : public ActionInterface
{
    Q_OBJECT
public:
    PrimitivePatchesInterface(QWidget *parent = nullptr);

    void display();
    void reloadShaders();

    void replay(nlohmann::json action);

    QLayout* createGUI();

    void show();
    void hide();

public Q_SLOTS:
    void mouseMovedOnMap(Vector3 newPosition);
    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event);

    void setSelectedShape(PRIMITIVE_SHAPE newShape, Vector3 newPosition = Vector3());
    void setCurrentOperation(PATCH_OPERATION newOperation);

    void resetPatch();

protected:
    PRIMITIVE_SHAPE currentShapeSelected = PRIMITIVE_SHAPE::SPHERE;
    PATCH_OPERATION currentOperation = PATCH_OPERATION::STACK;
    float selectedDensity = 1.5f;

    float selectedWidth = 10.f;
    float selectedDepth = 10.f;
    float selectedHeight = 10.f;
    float selectedSigma = 10.f;

    Mesh previewMesh;

    Patch3D* currentPatch;
    std::vector<Patch3D*> storedPatches;

    std::function<float(Vector3)> currentSelectionFunction;

    std::function<float(Vector3)> blockFunction(float height);
    std::function<float(Vector3)> sphereFunction(float radius);
    std::function<float(Vector3)> gaussianFunction(float sigma, float height);

    Vector3 currentPos;

    Vector3 functionSize = Vector3(1, 1, 1);
    void updateFunctionSize();
};

#endif // PRIMITIVEPATCHESINTERFACE_H
