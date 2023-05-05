#ifndef VIEWER_H
#define VIEWER_H

class Viewer;


enum MapMode {
    GRID_MODE       = 0b0001,
    VOXEL_MODE      = 0b0010,
    LAYER_MODE      = 0b0100,
    IMPLICIT_MODE   = 0b1000
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

#include "TerrainGen/Heightmap.h"
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
#include "Interface/ActionInterface.h"


class Viewer : public QGLViewer {
    Q_OBJECT
public:
    Viewer(QWidget *parent = nullptr);
    Viewer(std::shared_ptr<Heightmap> grid, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPath, MapMode map = VOXEL_MODE, ViewerMode mode = FILL_MODE, QWidget *parent = nullptr);
    Viewer(std::shared_ptr<Heightmap> g, QWidget *parent = nullptr);
    Viewer(std::shared_ptr<VoxelGrid> g, QWidget *parent = nullptr);
    ~Viewer();

Q_SIGNALS:
    void mouseClickOnMap(Vector3 mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    void mouseDoubleClickedOnMap(Vector3 mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model);
    void mouseMovedOnMap(Vector3 mouseWorldPosition, TerrainModel* model);

public Q_SLOTS:
    bool startRecording(std::string folderUsed = "");
    bool stopRecording();
    bool startStopRecording();

    void setViewerMode(ViewerMode newMode) { this->viewerMode = newMode; update(); }
    void setMapMode(MapMode newMode) { this->mapMode = newMode; update(); }
    void setSmoothingAlgorithm(SmoothingAlgorithm newAlgo) { this->algorithm = newAlgo;
//                                                             voxelGrid->displayWithMarchingCubes = this->algorithm == MARCHING_CUBES;
//                                                             voxelGrid->createMesh();
                                                            update();}

    void swapCamera(qglviewer::Camera* altCamera, bool useAltCamera);

    void clipViewTemporarily(Vector3 direction, Vector3 center, bool active);

    void drawingProcess();
    void reloadAllShaders();

    void setupViewFromFile(std::string filename);
    void saveViewToFile(std::string filename);

    void screenshot();

//protected:
public:
    virtual void init();
    virtual void draw();

    TerrainModel* getCurrentTerrainModel();

    bool eventFilter(QObject* obj, QEvent* event);

    void mousePressEvent(QMouseEvent* e);
    void keyPressEvent(QKeyEvent *e);
    void keyReleaseEvent(QKeyEvent* e);

    void mouseMoveEvent(QMouseEvent* e);
    void mouseDoubleClickEvent(QMouseEvent* e);

    void closeEvent(QCloseEvent* e);

    void animate();

    ViewerMode viewerMode;
    MapMode mapMode;
    SmoothingAlgorithm algorithm = MARCHING_CUBES;

//private:
    std::shared_ptr<Heightmap> heightmap;
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<LayerBasedGrid> layerGrid;
    std::shared_ptr<ImplicitNaryOperator> implicitTerrain;
    bool display_vertices = true;
    qglviewer::Vec selectedPoint, orig, dir;

    std::map<std::string, std::shared_ptr<ActionInterface>> interfaces;

    int current_frame = 0;

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
    int voxelsSmoothedOnBorders = 1;

    PositionalLight light;

    std::string main_screenshotFolder;
    std::string screenshotFolder;
    bool isTakingScreenshots = false;
    int screenshotIndex = 0;

    bool checkMouseOnVoxel();

    unsigned int frame_num = 0;

    ControlPoint *mainGrabber;

    bool displayParticles = false;
    float fogNear = 5.f; //1000.f;
    float fogFar = 30.f; //5000.f;
    bool usingSpotlight = false;

    std::string mapSavingFilename = "map1.data";
    std::string mapSavingFolder;

    qglviewer::Camera *mainCamera;
    qglviewer::Camera *alternativeCamera;
    bool usingMainCamera = true;

    bool temporaryClipPlaneActivated = false;
    Vector3 clipPlanePosition;
    Vector3 clipPlaneDirection;

    std::shared_ptr<Shader> raymarchingShader;
    ShaderUBO sceneUBO;
    Mesh raymarchingQuad;
};


#endif // VIEWER_H
