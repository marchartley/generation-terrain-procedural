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
#include <QGLViewer/camera.h>
#include "Interface/PathCameraConstraint.h"

using namespace TreeColonisationAlgo;

class SpaceColonizationInterface : public QWidget
{
    Q_OBJECT
public:
    SpaceColonizationInterface();
    ~SpaceColonizationInterface();

    void display();

    bool isHidden;

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    std::shared_ptr<VoxelGrid> voxelGrid;

    ControlPoint* startingPoint;
    std::vector<ControlPoint*> controlPoints;
    std::vector<BSpline> karstPaths;
    Mesh pathsMeshes;

    QHBoxLayout* createGUI();
    QHBoxLayout* spaceColonizationLayout;

Q_SIGNALS:
    void useAsMainCamera(qglviewer::Camera* cam, bool useMyCamera);

public Q_SLOTS:
    void initSpaceColonizer();
    void computeKarst();
    void updateKarstPath();
    void createKarst();

public:
    TreeColonisation* colonizer = nullptr;
    qglviewer::Camera* visitingCamera = nullptr;
    PathCameraConstraint *cameraConstraint;
};

#endif // SPACECOLONIZATIONINTERFACE_H
