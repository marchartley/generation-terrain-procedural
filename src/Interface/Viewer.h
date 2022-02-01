#ifndef VIEWER_H
#define VIEWER_H

class Viewer;

#include "TerrainGen/Grid.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/LayerBasedGrid.h"
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <qmessagebox.h>
#include "TerrainModification/RockErosion.h"
#include "Graphics/Shader.h"
#include "Graphics/DebugShader.h"
#include "Graphics/Sphere.h"
#include <QObject>
#include "Karst/KarstPathsGeneration.h"
#include "Utils/BSpline.h"
#include "TreeColonisation/TreeColonisation.h"
#include "Interface/ControlPoint.h"
#include "Interface/KarstPathGenerationInterface.h"
#include "Interface/SpaceColonizationInterface.h"

enum MapMode {
    GRID_MODE  = 0b001,
    VOXEL_MODE = 0b010,
    LAYER_MODE = 0b100,
};
enum ViewerMode {
    FILL_MODE  = 0b001,
    WIRE_MODE  = 0b010,
    NO_DISPLAY = 0b100,
};
enum SmoothingAlgorithm {
    NONE            = 0b001,
    MARCHING_CUBES  = 0b010,
    DUAL_CONTOURING = 0b100
};
enum DebugMeshesNames {
    ROCK_TRAILS  = 0b000001,
    FAILED_ROCKS = 0b000010,
    FLOW_TRAILS  = 0b000100,
    TUNNEL_PATHS = 0b001000,
    KARST_PATHS  = 0b010000,
    SPACE_COLONI = 0b100000
};

class Viewer : public QGLViewer {
    Q_OBJECT
public:
    Viewer(QWidget *parent = nullptr);
    Viewer(std::shared_ptr<Grid> grid, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, MapMode map = VOXEL_MODE, ViewerMode mode = FILL_MODE, QWidget *parent = nullptr);
    Viewer(std::shared_ptr<Grid> g, QWidget *parent = nullptr);
    Viewer(std::shared_ptr<VoxelGrid> g, QWidget *parent = nullptr);
    ~Viewer();


public Q_SLOTS:
    void setErosionRocksSize(int newSize) { this->erosionSize = newSize;}
    void setErosionRocksStrength(float newStrength) { this->erosionStrength = newStrength;}
    void setErosionRocksQuantity(int newQtt) { this->erosionQtt = newQtt;}
    void erodeMap(bool sendFromCam = false);

    void recomputeFlowfield();

    void setManualErosionRocksSize(int newSize);
    void setManualErosionRocksStrength(float newStrength) { this->manualErosionStrength = newStrength;}
    void setAddingMatterMode(bool addingMode) { this->addingMatterMode = addingMode; }
    void throwRock();

    void setCurvesErosionSize(int newSize) { this->curvesErosionSize = newSize; }
    void setCurvesErosionStrength(float newStrength) { this->curvesErosionStrength = newStrength;}
    void setCurvesErosionAddingMatterMode(bool addingMode) { this->addingCurvesErosionMatterMode = addingMode; }
    void addCurvesControlPoint(Vector3 pos, bool justUpdatePath = false);
    void setCurvesErosionConstructionMode(bool isConstructing) {this->curvesErosionConstructionMode = isConstructing; }
    void createTunnel(bool removingMatter = true);

    bool createGlobalGravity();
    bool createSandGravity();

    bool startStopRecording();

    void setViewerMode(ViewerMode newMode) { this->viewerMode = newMode; update(); }
    void setMapMode(MapMode newMode) { this->mapMode = newMode; update(); }
    void setSmoothingAlgorithm(SmoothingAlgorithm newAlgo) { this->algorithm = newAlgo;
                                                             voxelGrid->displayWithMarchingCubes = this->algorithm == MARCHING_CUBES;
                                                             voxelGrid->createMesh();
                                                            update();}
    void setErosionFlowfieldFactor(float newVal) { this->erosionFlowfieldFactor = newVal; }
    void setErosionFlowfieldRandomness(float newVal) { this->erosionFlowfieldRandomness = newVal; }

    void setLoD(int newLoD) { this->LoD = newLoD; }
    void computeLoD();

    void swapCamera(qglviewer::Camera* altCamera, bool useAltCamera);
    void frameInterpolated();

    void loadMapUI();
    void saveMapUI();

//    void saveScreenshot();

//protected:
public:
    virtual void init();
    virtual void draw();

    void mousePressEvent(QMouseEvent* e);
    void keyPressEvent(QKeyEvent *e);

    void mouseMoveEvent(QMouseEvent* e);

    void closeEvent(QCloseEvent* e);

    void animate();

    ViewerMode viewerMode;
    MapMode mapMode;
    SmoothingAlgorithm algorithm = MARCHING_CUBES;

//private:
    std::shared_ptr<Grid> grid;
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<LayerBasedGrid> layerGrid;
    bool display_vertices = true;
    qglviewer::Vec selectedPoint, orig, dir;

    int erosionSize = 1;
    float erosionStrength = .0;
    int erosionQtt = 1000;
    float erosionFlowfieldFactor = 1.0;
    float erosionFlowfieldRandomness = 0.05;

    int manualErosionSize = 15;
    float manualErosionStrength = .5;

    int curvesErosionSize = 15;
    float curvesErosionStrength = .5;
    bool addingCurvesErosionMatterMode = true;

//    KarstPathsGeneration karstPathCreator;
//    std::vector<BSpline> karstPaths;
    std::shared_ptr<KarstPathGenerationInterface> karstPathInterface = nullptr;
    std::shared_ptr<SpaceColonizationInterface> spaceColonizationInterface = nullptr;

//    TreeColonisationAlgo::TreeColonisation spaceColonizer;
//    std::vector<BSpline> spaceColonizerPaths;

    int LoD = 0;


    int getMaxLoDAvailable() { return this->voxelGrid->getMaxLoD(); }

    bool addingMatterMode = true;
    RockErosion matter_adder;
    std::vector<std::vector<Vector3>> lastRocksLaunched;
    std::vector<std::vector<Vector3>> lastFailedRocksLaunched;
    bool curvesErosionConstructionMode = false;
    std::vector<Vector3> currentTunnelPoints;
    std::vector<Vector3> tunnelPath;
    bool mouseInWorld = false;
    Vector3 mousePosWorld;
    QPoint mousePos;

    Vector3 minVoxelsShown();
    Vector3 maxVoxelsShown();

    float minSliceMapX = 0.f;
    float maxSliceMapX = 1.f;
    float minSliceMapY = 0.f;
    float maxSliceMapY = 1.f;
    float minSliceMapZ = 0.f;
    float maxSliceMapZ = 1.f;

    std::shared_ptr<Shader> shader;
    PositionalLight light;

    std::string screenshotFolder;
    bool isTakingScreenshots = false;
    int screenshotIndex = 0;

    bool applyLetItFall = false;
    bool applyLetSandFall = false;

    bool checkMouseOnVoxel();

    unsigned int frame_num = 0;

//    std::vector<ControlPoint> grabbers;
    ControlPoint *mainGrabber;

    void updateFlowfieldDebugMesh();

    std::map<DebugMeshesNames, Mesh> debugMeshes;
    std::map<DebugMeshesNames, std::vector<ControlPoint*>> debugControlPoints;

    std::string mapSavingFilename = "map1.data";
    std::string mapSavingFolder;

    qglviewer::Camera *mainCamera;
    qglviewer::Camera *alternativeCamera;
    bool usingMainCamera = true;
};


#endif // VIEWER_H
