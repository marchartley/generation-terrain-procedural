#ifndef SPACECOLONIZATIONINTERFACE_H
#define SPACECOLONIZATIONINTERFACE_H

class SpaceColonizationInterface;
#include "Interface/ControlPoint.h"
#include "Karst/KarstPathsGeneration.h"
#include "Graphics/Mesh.h"
#include "Interface/FancySlider.h"
#include "Utils/BSpline.h"
#include <QWidget>
#include "TreeColonisation/TreeColonisation.h"
#include "TerrainGen/VoxelGrid.h"
#include "Interface/VisitingCamera.h"
#include "Interface/PathCameraConstraint.h"

using namespace TreeColonisationAlgo;

class SpaceColonizationInterface : public CustomInteractiveObject
{
    Q_OBJECT
public:
    SpaceColonizationInterface(QWidget *parent = nullptr);
    ~SpaceColonizationInterface();

    void display();

//    bool isHidden;
    void hide();
    void show();

    void keyPressEvent(QKeyEvent* event);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    std::shared_ptr<VoxelGrid> voxelGrid;

    std::unique_ptr<ControlPoint> startingPoint;
    std::vector<std::unique_ptr<ControlPoint>> controlPoints;
    std::vector<BSpline> karstPaths;
    Mesh pathsMeshes;

    QHBoxLayout* createGUI();
    QHBoxLayout* spaceColonizationLayout;

Q_SIGNALS:
    void useAsMainCamera(qglviewer::Camera* cam, bool useMyCamera);
    void karstPathUpdated();

public Q_SLOTS:
    void initSpaceColonizer();
    void computeKarst();
    void updateKarstPath();
    void createKarst(bool usingSpheres = false);

public:
    float karstWidth = 10.f;
    TreeColonisation* colonizer = nullptr;
    VisitingCamera* visitingCamera = nullptr;
    PathCameraConstraint *cameraConstraint;
};

#endif // SPACECOLONIZATIONINTERFACE_H
