#ifndef KARSTPATHGENERATIONINTERFACE_H
#define KARSTPATHGENERATIONINTERFACE_H

class KarstPathGenerationInterface;
#include "Interface/ControlPoint.h"
#include "Interface/InteractiveVector.h"
#include "Interface/Slider3D.h"
#include "Karst/KarstPathsGeneration.h"
#include "Graphics/Mesh.h"
#include "Interface/FancySlider.h"
#include "Utils/BSpline.h"
#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "Interface/PathCameraConstraint.h"
#include <QGLViewer/manipulatedCameraFrame.h>

class KarstPathGenerationInterface : public QWidget
{
    Q_OBJECT
public:
    KarstPathGenerationInterface();
//    KarstPathGenerationInterface(KarstPathsGeneration karstCreator, Vector3 AABBoxMinPos, Vector3 AABBoxMaxPos);

    void display();

    InteractiveVector *fractureVector;
    Slider3D *waterHeightSlider;
    Mesh waterLevelMesh;

    bool isHidden;

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);
    std::shared_ptr<VoxelGrid> voxelGrid;

    std::vector<BSpline> karstPaths;

    Mesh pathsMeshes;

    QHBoxLayout* createGUI();
    QHBoxLayout* karstCreationLayout;

Q_SIGNALS:
    void useAsMainCamera(qglviewer::Camera* cam, bool useMyCamera);

public Q_SLOTS:
    void updateFracture(Vector3 newFractureDir = Vector3());
    void updateWaterHeight(float newHeight = 0.f);

    void computeKarst();
    void updateKarstPath();
    void createKarst();

public:
    Vector3 AABBoxMinPos;
    Vector3 AABBoxMaxPos;

    KarstPathsGeneration* karstCreator = nullptr;
    ControlPoint *sourceControlPoint;
    qglviewer::Camera* visitingCamera = nullptr;
    PathCameraConstraint *cameraConstraint;
};

#endif // KARSTPATHGENERATIONINTERFACE_H
