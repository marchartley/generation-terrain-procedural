#ifndef VIEWER_H
#define VIEWER_H

#include "Grid.h"
#include "VoxelGrid.h"
#include "LayerBasedGrid.h"
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <qmessagebox.h>
#include "RockErosion.h"
#include "Shader.h"
#include <QObject>

enum MapMode {
    GRID_MODE  = 0b001,
    VOXEL_MODE = 0b010,
    LAYER_MODE = 0b100
};
enum ViewerMode {
    FILL_MODE  = 0b01,
    WIRE_MODE  = 0b10
};
enum SmoothingAlgorithm {
    NONE            = 0b001,
    MARCHING_CUBES  = 0b010,
    DUAL_CONTOURING = 0b100
};

class Viewer : public QGLViewer {
    Q_OBJECT
public:
    Viewer(QWidget *parent = nullptr);
    Viewer(Grid* grid, VoxelGrid* voxelGrid, LayerBasedGrid* layerGrid, MapMode map = VOXEL_MODE, ViewerMode mode = FILL_MODE, QWidget *parent = nullptr);
    Viewer(Grid* g, QWidget *parent = nullptr);
    Viewer(VoxelGrid* g, QWidget *parent = nullptr);
    ~Viewer();


public Q_SLOTS:
    void setErosionRocksSize(int newSize) { this->erosionSize = newSize;}
    void setErosionRocksStrength(float newStrength) { this->erosionStrength = newStrength;}
    void setErosionRocksQuantity(int newQtt) { this->erosionQtt = newQtt;}
    void erodeMap(bool sendFromCam = false);

    void recomputeFlowfield();

    void setManualErosionRocksSize(int newSize) { this->manualErosionSize = newSize;}
    void setManualErosionRocksStrength(float newStrength) { this->manualErosionStrength = newStrength;}
    void setAddingMatterMode(bool addingMode) { this->addingMatterMode = addingMode; }
    void throwRock();

    void setCurvesErosionSize(int newSize) { this->curvesErosionSize = newSize; }
    void setCurvesErosionStrength(float newStrength) { this->curvesErosionStrength = newStrength;}
    void setCurvesErosionAddingMatterMode(bool addingMode) { this->addingCurvesErosionMatterMode = addingMode; }
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
    Grid* grid;
    VoxelGrid* voxelGrid;
    LayerBasedGrid* layerGrid;
    bool display_vertices = true;
    qglviewer::Vec selectedPoint, orig, dir;

    bool displayRockTrajectories = false;
    bool displayFailedRockTrajectories = false;
    bool displayFlowfield = false;
    bool displayTunnelsPath = false;

    int erosionSize = 15;
    float erosionStrength = .5;
    int erosionQtt = 10;
    float erosionFlowfieldFactor = 0.1;
    float erosionFlowfieldRandomness = 0.05;

    int manualErosionSize = 15;
    float manualErosionStrength = .5;

    int curvesErosionSize = 15;
    float curvesErosionStrength = .5;
    bool addingCurvesErosionMatterMode = true;

    int LoD = 0;


    int getMaxLoDAvailable() { return this->voxelGrid->getMaxLoD(); }

    float*** selectionShape;
    bool addingMatterMode = true;
    RockErosion matter_adder;
    std::vector<std::vector<Vector3>> lastRocksLaunched;
    std::vector<std::vector<Vector3>> lastFailedRocksLaunched;
    std::vector<Vector3> tunnelPath;
    bool mouseInWorld = false;
    Vector3 mousePosWorld;
    QPoint mousePos;

    Shader shader;
    PositionalLight light;
    GLuint rocksVBO;

    std::string screenshotFolder;
    bool isTakingScreenshots = false;
    int screenshotIndex = 0;

    bool applyLetItFall = false;
    bool applyLetSandFall = false;

    bool checkMouseOnVoxel();

    unsigned int frame_num = 0;

    Mesh rocksMeshes;
    Mesh failedRocksMeshes;
    Mesh flowDebugMeshes;
    Mesh tunnelsMesh;
};


#endif // VIEWER_H
