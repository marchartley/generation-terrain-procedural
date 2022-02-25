#include "Utils/Globals.h"
#include "Interface/Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>
#include <chrono>
#include "TerrainModification/UnderwaterErosion.h"
#include "DataStructure/Matrix.h"
#include "Utils/Utils.h"


//using namespace qglviewer;
//using namespace std;

Viewer::Viewer(QWidget *parent):
    Viewer(
        std::shared_ptr<Grid>(nullptr), // new Grid(100, 100, 40, 1.0),
        std::make_shared<VoxelGrid>(3, 3, 3, 1.0, .30),
        std::shared_ptr<LayerBasedGrid>(nullptr), // new LayerBasedGrid(10, 10, 50),
        VOXEL_MODE,
        FILL_MODE,
        parent
        )
{
}
Viewer::Viewer(std::shared_ptr<Grid> grid, std::shared_ptr<VoxelGrid> voxelGrid,
               std::shared_ptr<LayerBasedGrid> layerGrid, MapMode map,
               ViewerMode mode, QWidget *parent)
    : QGLViewer(parent), viewerMode(mode), mapMode(map), grid(grid), voxelGrid(voxelGrid), layerGrid(layerGrid)
{
    this->mainCamera = this->camera();
}
Viewer::Viewer(std::shared_ptr<Grid> g, QWidget *parent)
    : Viewer(g, nullptr, nullptr, GRID_MODE, FILL_MODE, parent) {

}
Viewer::Viewer(std::shared_ptr<VoxelGrid> g, QWidget *parent)
    : Viewer(nullptr, g, nullptr, VOXEL_MODE, FILL_MODE, parent) {

}
Viewer::~Viewer()
{
}

void Viewer::init() {
    restoreStateFromFile();
    setSceneRadius(500.0);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    GlobalsGL::generateBuffers();

    this->camera()->setType(qglviewer::Camera::PERSPECTIVE);

    setTextIsEnabled(true);
    setMouseTracking(true);

    const char* vShader_grid = ":/src/Shaders/grid_vertex_shader_blinn_phong.glsl";
    const char* fShader_grid = ":/src/Shaders/grid_fragment_shader_blinn_phong.glsl";
    const char* vShader_voxels = ":/src/Shaders/voxels_vertex_shader_blinn_phong.glsl";
    const char* fShader_voxels = ":/src/Shaders/voxels_fragment_shader_blinn_phong.glsl";
//    const char* vShader_layer = ":/src/Shaders/layer_based_vertex_shader.glsl";
//    const char* fShader_layer = ":/src/Shaders/layer_based_fragment_shader.glsl";
    const char* vNoShader = ":/src/Shaders/no_vertex_shader.glsl";
    const char* fNoShader = ":/src/Shaders/no_fragment_shader.glsl";

    glEnable              ( GL_DEBUG_OUTPUT );
    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 );

    Shader::default_shader = std::make_shared<Shader>(vNoShader, fNoShader);
    ControlPoint::base_shader = std::make_shared<Shader>(vNoShader, fNoShader);

    this->debugMeshes[ROCK_TRAILS] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->debugMeshes[FAILED_ROCKS] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->debugMeshes[FLOW_TRAILS] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->debugMeshes[TUNNEL_PATHS] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->debugMeshes[KARST_PATHS] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->debugMeshes[SPACE_COLONI] = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);

    this->debugControlPoints[ROCK_TRAILS] = std::vector<ControlPoint*>();
    this->debugControlPoints[FAILED_ROCKS] = std::vector<ControlPoint*>();
    this->debugControlPoints[FLOW_TRAILS] = std::vector<ControlPoint*>();
    this->debugControlPoints[TUNNEL_PATHS] = std::vector<ControlPoint*>();
    this->debugControlPoints[KARST_PATHS] = std::vector<ControlPoint*>();
    this->debugControlPoints[SPACE_COLONI] = std::vector<ControlPoint*>();

    this->debugMeshes[ROCK_TRAILS].shader->setVector("color", std::vector<float>({86/255.f, 176/255.f, 12/255.f, .5f}));
    this->debugMeshes[FAILED_ROCKS].shader->setVector("color", std::vector<float>({176/255.f, 72/255.f, 12/255.f, .4f}));
    this->debugMeshes[FLOW_TRAILS].shader->setVector("color", std::vector<float>({143/255.f, 212/255.f, 255/255.f, .5f}));
    this->debugMeshes[TUNNEL_PATHS].shader->setVector("color", std::vector<float>({152/255.f, 94/255.f, 209/255.f, 1.0}));
    this->debugMeshes[KARST_PATHS].shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.0}));
    this->debugMeshes[SPACE_COLONI].shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.0}));
    ControlPoint::base_shader->setVector("color", std::vector<float>({160/255.f, 5/255.f, 0/255.f, 1.f}));
    this->mainGrabber = new ControlPoint(Vector3(), 1.f, ACTIVE, false);

    // Don't compute the indices for this meshes, there's no chance any two vertex are the same
    for (auto& debugMesh : this->debugMeshes)
        std::get<1>(debugMesh).useIndices = false;

    this->matter_adder = RockErosion(this->erosionSize, 1.0);

    this->light = PositionalLight(
                new float[4]{.5, .5, .5, 1.},
                new float[4]{.2, .2, .2, 1.},
                new float[4]{.5, .5, .5, 1.},
                Vector3(0.0, 0.0, 100.0)
                );

    this->setAnimationPeriod(0);

    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);

#ifdef _WIN32
    this->screenshotFolder = "C:/codes/Qt/generation-terrain-procedural/screenshots/";
    this->mapSavingFolder = "C:/codes/Qt/generation-terrain-procedural/saved_maps/";
#elif linux
    this->screenshotFolder = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/screenshots/";
    this->mapSavingFolder = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/saved_maps/";
#endif
    if(!makedir(this->screenshotFolder)) {
        std::cerr << "Not possible to create folder " << this->screenshotFolder << std::endl;
        exit(-1);
    }
    if(!makedir(this->mapSavingFolder)) {
        std::cerr << "Not possible to create folder " << this->mapSavingFolder << std::endl;
        exit(-1);
    }
    if (this->voxelGrid != nullptr) {
        this->screenshotFolder += std::string(s_time) + "__" + voxelGrid->toShortString() + "/";
//        // this->displayMessage(QString::fromStdString(std::string("Screenshots will be saved in folder ") + std::string(this->screenshotFolder)));
    }

    if (grid != nullptr) {
        this->grid->createMesh();
        this->grid->mesh.shader = std::make_shared<Shader>(vShader_grid, fShader_grid);
    }
    if (layerGrid != nullptr) {
        this->layerGrid->createMesh();
        this->layerGrid->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
    }
    if (voxelGrid != nullptr) {
        voxelGrid->retrieveMap(this->mapSavingFolder + "cube.data");
        voxelGrid->fromIsoData();
        voxelGrid->displayWithMarchingCubes = (this->algorithm == MARCHING_CUBES);
        this->voxelGrid->createMesh();
        for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks)
            vc->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
        this->setSceneCenter(qglviewer::Vec(voxelGrid->blockSize * voxelGrid->sizeX/2, voxelGrid->blockSize * voxelGrid->sizeY/2, voxelGrid->blockSize * voxelGrid->sizeZ/2));
//        this->initSpaceColonizer();
        updateFlowfieldDebugMesh();
    }

    QObject::connect(this->spaceColonizationInterface.get(), &SpaceColonizationInterface::useAsMainCamera, this, &Viewer::swapCamera);
    QObject::connect(this->karstPathInterface.get(), &KarstPathGenerationInterface::useAsMainCamera, this, &Viewer::swapCamera);

    Mesh::setShaderToAllMeshesWithoutShader(*Shader::default_shader);
//    startAnimation();
    QGLViewer::init();
}

void Viewer::draw() {
    this->frame_num ++;
    glClear(GL_DEPTH_BUFFER_BIT);
    if (this->viewerMode == ViewerMode::WIRE_MODE)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);


    Material ground_material(
                    new float[4]{.48, .16, .04, 1.},
                    new float[4]{.60, .20, .08, 1.},
                    new float[4]{.62, .56, .37, 1.},
                    51.2f
                    );
    Material grass_material(
                    new float[4]{.28, .90, .00, 1.},
                    new float[4]{.32, .80, .00, 1.},
                    new float[4]{.62, .56, .37, 1.},
                    51.2f
                    );
//    this->light.position = Vector3(100.0 * std::cos(this->frame_num / (float)10), 100.0 * std::sin(this->frame_num / (float)10), 0.0);
    float globalAmbiant[4] = {.10, .10, .10, 1.0};

    Shader::applyToAllShaders([&](std::shared_ptr<Shader> shader) -> void {
        shader->setMatrix("proj_matrix", pMatrix);
        shader->setMatrix("mv_matrix", mvMatrix);
        shader->setPositionalLight("light", this->light);
        shader->setMaterial("ground_material", ground_material);
        shader->setMaterial("grass_material", grass_material);
        shader->setVector("globalAmbiant", globalAmbiant, 4);
        shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
        shader->setBool("display_light_source", true);
    });

    if (this->viewerMode != NO_DISPLAY)
    {
        if (this->mapMode == GRID_MODE) {
            if (this->grid == nullptr) {
                std::cerr << "No grid to display" << std::endl;
            } else {
                this->grid->mesh.shader->setMatrix("proj_matrix", pMatrix);
                this->grid->mesh.shader->setMatrix("mv_matrix", mvMatrix);
                this->grid->mesh.shader->setPositionalLight("light", this->light);
                this->grid->mesh.shader->setMaterial("ground_material", ground_material);
                this->grid->mesh.shader->setMaterial("grass_material", grass_material);
                this->grid->mesh.shader->setVector("globalAmbiant", globalAmbiant, 4);
                this->grid->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
                this->grid->mesh.shader->setBool("display_light_source", true);
                this->grid->display(true);
            }
        }
        else if (this->mapMode == VOXEL_MODE) {
            if (this->voxelGrid == nullptr) {
                std::cerr << "No voxel grid to display" << std::endl;
            } else {
                for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks) {
                    vc->mesh.shader->setMatrix("proj_matrix", pMatrix);
                    vc->mesh.shader->setMatrix("mv_matrix", mvMatrix);
                    vc->mesh.shader->setPositionalLight("light", this->light);
                    vc->mesh.shader->setMaterial("ground_material", ground_material);
                    vc->mesh.shader->setMaterial("grass_material", grass_material);
                    vc->mesh.shader->setVector("globalAmbiant", globalAmbiant, 4);
                    vc->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
                    vc->mesh.shader->setBool("display_light_source", true);
                    vc->mesh.shader->setVector("min_vertice_positions", minVoxelsShown());
                    vc->mesh.shader->setVector("max_vertice_positions", maxVoxelsShown());
                }
                this->voxelGrid->display();
            }
        }
        else if (this->mapMode == LAYER_MODE) {
            if (this->layerGrid == nullptr) {
                std::cerr << "No layer based grid to display" << std::endl;
            } else {
                this->layerGrid->mesh.shader->setPositionalLight("light", this->light);
                this->layerGrid->mesh.shader->setMaterial("ground_material", ground_material);
                this->layerGrid->mesh.shader->setMaterial("grass_material", grass_material);
                this->layerGrid->mesh.shader->setVector("globalAmbiant", globalAmbiant, 4);
                this->layerGrid->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
                this->layerGrid->mesh.shader->setMatrix("mv_matrix", mvMatrix);
                this->layerGrid->mesh.shader->setMatrix("proj_matrix", pMatrix);
                this->layerGrid->mesh.shader->setBool("display_light_source", true);
//                this->layerGrid->mesh.shader->setFloat("offsetX", 0.0);
//                this->layerGrid->mesh.shader->setFloat("offsetY", 0.0);
                this->layerGrid->display();
            }
        }
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    GLfloat previousLineWidth[1];
    glGetFloatv(GL_LINE_WIDTH, previousLineWidth);
    glLineWidth(5.f);
    for (auto& debugMesh : this->debugMeshes) {
        std::get<1>(debugMesh).shader->setMatrix("proj_matrix", pMatrix);
        std::get<1>(debugMesh).shader->setMatrix("mv_matrix", mvMatrix);
        std::get<1>(debugMesh).shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
        std::get<1>(debugMesh).display(GL_LINES, 3.f);
    }
    for (auto& controlPointsArray : this->debugControlPoints) {
        std::vector<ControlPoint*> grabbers = std::get<1>(controlPointsArray);
        for (auto& grabber : std::get<1>(controlPointsArray)) {
            grabber->mesh.shader->setMatrix("proj_matrix", pMatrix);
            grabber->mesh.shader->setMatrix("mv_matrix", mvMatrix);
            grabber->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
            grabber->mesh.isDisplayed = this->debugMeshes[std::get<0>(controlPointsArray)].isDisplayed;
            grabber->display();
        }
    }
    //    glBlendFunc(GL_SRC_ALPHA, GL_ZERO);
    glLineWidth(previousLineWidth[0]);

    this->mainGrabber->mesh.shader->setMatrix("proj_matrix", pMatrix);
    this->mainGrabber->mesh.shader->setMatrix("mv_matrix", mvMatrix);
    this->mainGrabber->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
    this->mainGrabber->display();

    if (this->karstPathInterface)
        this->karstPathInterface->display();
    if (this->spaceColonizationInterface)
        this->spaceColonizationInterface->display();

    if (this->isTakingScreenshots) {
#ifdef linux
        mode_t prevMode = umask(0011);
#endif
        if(this->screenshotIndex == 0 && voxelGrid)
        {
            std::ofstream outfile;
            outfile.open(this->screenshotFolder + "grid_data.json", std::ios_base::trunc);
            outfile << voxelGrid->toString();
            outfile.close();
        }
        this->window()->grab().save(QString::fromStdString(this->screenshotFolder + std::to_string(this->screenshotIndex++) + ".jpg"));
#ifdef linux
        chmod((this->screenshotFolder + std::to_string(this->screenshotIndex) + ".jpg").c_str(), 0666);
        umask(prevMode);
#endif
    }
}

void Viewer::mousePressEvent(QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
    if (curvesErosionConstructionMode && checkMouseOnVoxel()) {
        this->addCurvesControlPoint(this->mousePosWorld);
    }
    if (QApplication::keyboardModifiers().testFlag(Qt::AltModifier) == true)
    {
        this->throwRock();
    }
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
    // Defines the Alt+R shortcut.
    if (e->key() == Qt::Key_Z)
    {
        if (QApplication::keyboardModifiers().testFlag(Qt::ControlModifier) == true)
            this->voxelGrid->undo();
        else
            setViewerMode(ViewerMode::WIRE_MODE);
        update(); // Refresh display
    } else if (e->key() == Qt::Key_S)
    {
        setViewerMode(ViewerMode::FILL_MODE);
        update(); // Refresh display
    } else if (e->key() == Qt::Key_Q)
    {
        setMapMode(MapMode::VOXEL_MODE);
        update(); // Refresh display
    } else if (e->key() == Qt::Key_D)
    {
        setMapMode(MapMode::GRID_MODE);
        update(); // Refresh display
    } else if (e->key() == Qt::Key_F)
    {
        setMapMode(MapMode::LAYER_MODE);
        update(); // Refresh display
    } else if (e->key() == Qt::Key_R) {
        if (this->algorithm == NONE)
            setSmoothingAlgorithm(MARCHING_CUBES);
        else if (this->algorithm == MARCHING_CUBES)
            setSmoothingAlgorithm(NONE);
        // this->displayMessage(QString::fromStdString("Displaying using " + std::string(this->algorithm == MARCHING_CUBES ? " Marching cubes" : "no") + " algorithm") );
        update();
    } else if(e->key() == Qt::Key_V) {
        this->display_vertices = !this->display_vertices;
        update();
    } else if(e->key() == Qt::Key_P) {
        this->setAddingMatterMode(!this->addingMatterMode);
        // this->displayMessage( (addingMatterMode ? "Construction mode" : "Destruction mode") );
        update();
    } else if(e->key() == Qt::Key_Return) {
        erodeMap(e->modifiers() == Qt::ShiftModifier);
    } else if(e->key() == Qt::Key_Minus) {
        this->setManualErosionRocksSize(std::max(2, this->erosionSize - 2));
        // this->displayMessage(QString::fromStdString("Cursor size : " + std::to_string(this->manualErosionSize) ));
        update();
    } else if(e->key() == Qt::Key_Plus) {
        this->setManualErosionRocksSize(std::max(2, this->erosionSize + 2));
        // this->displayMessage(QString::fromStdString("Cursor size : " + std::to_string(this->manualErosionSize) ));
        update();
    } else if(e->key() == Qt::Key_Space) {
//        displayRockTrajectories = !displayRockTrajectories;
        // this->displayMessage(QString::fromStdString("Rock trajectories are : " + std::string(displayRockTrajectories ? "ON" : "OFF") ));
        update();
    } else if(e->key() == Qt::Key_0) {
        this->createGlobalGravity();
    } else if(e->key() == Qt::Key_Comma) {
        this->createSandGravity();
    } else if(e->key() == Qt::Key_1) {
        this->startStopRecording();
    } else if(e->key() == Qt::Key_2) {
        for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks) {
            vc->LoDIndex++;
            vc->needRemeshing = true;
        }
        this->voxelGrid->remeshAll();
        update();
    } else if(e->key() == Qt::Key_3) {
        // this->displayMessage( "Removing matter to create a tunnel" );
        createTunnel(true);
    } else if(e->key() == Qt::Key_4) {
        // this->displayMessage( "Adding matter to create a tunnel" );
        createTunnel(false);
    } else if(e->key() == Qt::Key_5) {
        this->setCamera(this->spaceColonizationInterface->visitingCamera);
    } else {
        QGLViewer::keyPressEvent(e);
    }
}
void Viewer::mouseMoveEvent(QMouseEvent* e)
{
    this->mousePos = e->pos();

    if (this->checkMouseOnVoxel())
    {
        this->mainGrabber->move(this->mousePosWorld);
        this->mainGrabber->setState(ACTIVE);
    } else {
        this->mainGrabber->setState(HIDDEN);
    }
    update();
    QGLViewer::mouseMoveEvent(e);
}

void Viewer::animate()
{
    /*if (voxelGrid) {
        voxelGrid->computeFlowfield();
        this->updateFlowfieldDebugMesh();
    }*/
    if (this->applyLetItFall)
        this->voxelGrid->makeItFall((this->applyLetSandFall ? -1.0 : 0.1));
    if (this->applyLetSandFall)
        this->voxelGrid->letGravityMakeSandFall();
}

Vector3 Viewer::minVoxelsShown()
{
    Vector3 minVec(minSliceMapX, minSliceMapY, minSliceMapZ);
    return /*Vector3(-this->voxelGrid->sizeX/2, -this->voxelGrid->sizeY/2, 0) + */Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ) * minVec;
}

Vector3 Viewer::maxVoxelsShown()
{
    Vector3 maxVec(maxSliceMapX, maxSliceMapY, maxSliceMapZ);
    return /*Vector3(-this->voxelGrid->sizeX/2, -this->voxelGrid->sizeY/2, 0) + */Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ) * maxVec;
}

void Viewer::erodeMap(bool sendFromCam)
{
//    this->voxelGrid->computeFlowfield();
    UnderwaterErosion erod(this->voxelGrid, this->erosionSize, this->erosionStrength, this->erosionQtt);
    std::shared_ptr<Vector3> pos = nullptr;
    std::shared_ptr<Vector3> dir = nullptr;
    if (sendFromCam)
    {
        dir = std::make_shared<Vector3>(1.f, .0f, .0f);/*
        Vec a;
        Vec b;
        camera()->convertClickToLine(QPoint(camera()->screenWidth()/2, camera()->screenHeight()/2), a, b);
//        pos = std::make_shared<Vector3>(a.x, a.y, a.z);
        dir = std::make_shared<Vector3>(b.x, b.y, b.z);
        // this->displayMessage( "Rocks launched from camera!" );*/
    } else {
        pos = nullptr;
        dir = std::make_shared<Vector3>(new Vector3(0.0, 0.0, 0.0));
        // this->displayMessage( "Rocks launched!" );
    }
    std::tie(this->lastRocksLaunched, this->lastFailedRocksLaunched) = erod.Apply(pos, dir, 10, this->erosionFlowfieldFactor, this->erosionFlowfieldRandomness, true);

    std::vector<Vector3> asOneVector;
    for(std::vector<Vector3>& coords : this->lastRocksLaunched) {
        asOneVector.insert(asOneVector.end(), coords.begin(), coords.end());
    }
    this->debugMeshes[ROCK_TRAILS].fromArray(asOneVector);
    this->debugMeshes[ROCK_TRAILS].update();
    asOneVector.clear();
    for(std::vector<Vector3>& coords : this->lastFailedRocksLaunched) {
        asOneVector.insert(asOneVector.end(), coords.begin(), coords.end());
    }
    this->debugMeshes[FAILED_ROCKS].fromArray(asOneVector);
    this->debugMeshes[FAILED_ROCKS].update();

    updateFlowfieldDebugMesh();
}

void Viewer::recomputeFlowfield()
{
    this->voxelGrid->computeFlowfield();
    updateFlowfieldDebugMesh();
}

void Viewer::setManualErosionRocksSize(int newSize)
{
    this->manualErosionSize = newSize;
    this->mainGrabber->radius = newSize / 2.f;
}

void Viewer::throwRock()
{
    if (this->mouseInWorld)
    {
        RockErosion rock(this->manualErosionSize, this->manualErosionStrength);
        rock.Apply(this->voxelGrid, this->mousePosWorld, addingMatterMode, true);
    }
    update();
}

void Viewer::computeLoD()
{
    for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks) {
        vc->LoDIndex = this->LoD;
        vc->needRemeshing = true;
    }
    this->voxelGrid->remeshAll();
    this->update();
}

void Viewer::swapCamera(qglviewer::Camera *altCamera, bool useAltCamera)
{
    if (altCamera != nullptr) {
        this->alternativeCamera = altCamera;
    }
    if (useAltCamera) {
        this->setCamera(this->alternativeCamera);
        /*
        if (this->spaceColonizationInterface->keyFramesPaths.size() > 0) {
//        if  (this->alternativeCamera->keyFrameInterpolator(0)) {
//            qglviewer::KeyFrameInterpolator *interp = this->alternativeCamera->keyFrameInterpolator(0);
            qglviewer::KeyFrameInterpolator *interp = this->spaceColonizationInterface->keyFramesPaths[0];
            interp->setInterpolationSpeed(1.0);
            interp->setLoopInterpolation(true);
            QObject::connect(interp, SIGNAL(interpolated()),
                             this, SLOT(frameInterpolated()));
            QObject::connect(interp, SIGNAL(endReached()),
                             this, SLOT(frameInterpolated()));
            for (int i = 0; i < interp->numberOfKeyFrames(); i++) {
                qglviewer::Vec v;
                interp->keyFrame(i).getPosition(v.x, v.y, v.z);
                std::cout << v.x << " " << v.y << " " << v.z << " -> t=" << interp->keyFrameTime(i) << std::endl;
//                interp->keyFrame(i).setPosition(v.x * 1000, v.y * 1000, v.z * 1000);
            }
//            this->camera()->setKeyFrameInterpolator(0, interp);
//            interp->setFrame(this->camera()->frame());
            interp->startInterpolation();
            std::cout << (interp->interpolationIsStarted() ? "Started" : "Not started") << std::endl;
            std::cout << "Should take " << interp->duration() << "s between " << interp->firstTime() << " and " << interp->lastTime() << std::endl;
//            this->camera()->playPath(0);
        }*/
    }
    else {
        this->setCamera(this->mainCamera);
    }
}

void Viewer::frameInterpolated()
{
    std::cout << "Interpolation" << std::endl;
}

void Viewer::addCurvesControlPoint(Vector3 pos, bool justUpdatePath)
{
    if (!justUpdatePath)
    {
        bool addTheNewPoint = true;
        for (auto& controls : this->debugControlPoints[TUNNEL_PATHS]) {
            if (controls->manipFrame.isManipulated()) {
                addTheNewPoint = false;
                break;
            }
        }
        if (addTheNewPoint) {
            this->debugControlPoints[TUNNEL_PATHS].push_back(new ControlPoint(pos, 5.f, INACTIVE));
            QObject::connect(this->debugControlPoints[TUNNEL_PATHS].back(), &ControlPoint::modified,
                             this, [&](){ this->addCurvesControlPoint(Vector3(), true); });
        }
    }
    this->currentTunnelPoints.clear();
    for (auto& controls : this->debugControlPoints[TUNNEL_PATHS]) {
        this->currentTunnelPoints.push_back(controls->position);
//        controls->onUpdate([=]{ this->addCurvesControlPoint(Vector3(), true); });
    }
    BSpline path(this->currentTunnelPoints);

    std::vector<Vector3> vertices = path.getPath(0.01);
    std::vector<Vector3> meshVertices;
    for (size_t i = 0; i < vertices.size() - 1; i++)
    {
        meshVertices.push_back(vertices[i]);
        meshVertices.push_back(vertices[i+1]);
    }
    this->debugMeshes[TUNNEL_PATHS].fromArray(meshVertices);
    this->debugMeshes[TUNNEL_PATHS].update();

    update();
}

void Viewer::clearTunnelPoints()
{
    this->curvesErosionConstructionMode = false;
    this->currentTunnelPoints.clear();
    this->debugControlPoints[TUNNEL_PATHS].clear();
    this->debugMeshes[TUNNEL_PATHS].clear();
    update();
}
void Viewer::createTunnel(bool removingMatter)
{
    this->curvesErosionConstructionMode = false;
    UnderwaterErosion erod(this->voxelGrid, this->curvesErosionSize, curvesErosionStrength, 10);
    if (this->currentTunnelPoints.empty())
        this->debugMeshes[TUNNEL_PATHS].fromArray(erod.CreateTunnel(3, !removingMatter));
    else
        this->debugMeshes[TUNNEL_PATHS].fromArray(erod.CreateTunnel(this->currentTunnelPoints, !removingMatter));
    this->currentTunnelPoints.clear();
    this->debugControlPoints[TUNNEL_PATHS].clear();
    this->debugMeshes[TUNNEL_PATHS].update();
    update();
}

bool Viewer::checkMouseOnVoxel()
{
    if (voxelGrid == nullptr)
        return false;
/*
    bool isFound = false;
    qglviewer::Vec pos = camera()->pointUnderPixel(mousePos, isFound);
    this->mousePosWorld = Vector3(pos.x, pos.y, pos.z);
    this->mouseInWorld = isFound;
    return isFound;


*/
    camera()->convertClickToLine(mousePos, orig, dir);
    float maxDist = std::max((int)camera()->distanceToSceneCenter(), std::max(voxelGrid->getSizeX(), std::max(voxelGrid->getSizeY(), voxelGrid->getSizeZ())));
    maxDist *= maxDist;

    Vector3 minPos = minVoxelsShown(), maxPos = maxVoxelsShown();

    bool found = false;
    Vector3 currPos(orig.x, orig.y, orig.z);
//    std::cout << "Click from " << currPos << " to ... ";
//    currPos += Vector3(voxelGrid->getSizeX()/2, voxelGrid->getSizeY()/2, 0.0);
    while((currPos / 2.f).norm2() < maxDist && !found)
    {
        currPos += Vector3(dir.x, dir.y, dir.z);
        if (minPos.x <= currPos.x && currPos.x <= maxPos.x && minPos.y <= currPos.y && currPos.y <= maxPos.y && minPos.z <= currPos.z && currPos.z <= maxPos.z) {
            float isoval = voxelGrid->getVoxelValue(currPos);
            if (isoval > 0.0)
                found = true;
        }
    }
    this->mouseInWorld = found;
    if (found) {
//        std::cout << "Click on " << currPos << std::endl;
        this->mousePosWorld = currPos;
        this->mainGrabber->position = currPos; // - Vector3(voxelGrid->getSizeX()/2, voxelGrid->getSizeY()/2, 0.0);
    }
//    std::cout << currPos << " (length=" << currPos.norm() << " over " << std::sqrt(maxDist) << ")" << std::endl;
    return found;
}

void Viewer::updateFlowfieldDebugMesh()
{
    std::vector<Vector3> normals;
    for (int x = this->voxelGrid->fluidSimRescale; x < this->voxelGrid->sizeX-1; x+= this->voxelGrid->fluidSimRescale) {
        for (int y = this->voxelGrid->fluidSimRescale; y < this->voxelGrid->sizeY-1; y+= this->voxelGrid->fluidSimRescale) {
            for (int z = this->voxelGrid->fluidSimRescale; z < this->voxelGrid->sizeZ - 1; z+= this->voxelGrid->fluidSimRescale) {
                normals.push_back(Vector3(x, y, z) + .5); // - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
                normals.push_back(Vector3(x, y, z) + (this->voxelGrid->getFlowfield(x, y, z) * (float)voxelGrid->fluidSimRescale) + .5); // - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
            }
        }
    }
    this->debugMeshes[FLOW_TRAILS].fromArray(normals);
    this->debugMeshes[FLOW_TRAILS].update();
    update();
}

void Viewer::closeEvent(QCloseEvent *e) {
    this->setCamera(this->mainCamera);
    if (this->isTakingScreenshots) this->startStopRecording();
    QGLViewer::closeEvent(e);
}

bool Viewer::createGlobalGravity()
{
    this->voxelGrid->makeItFall();
    update();
    return false;
    /*
    this->startAnimation();
    this->applyLetItFall = !this->applyLetItFall;
//    if (this->applyLetItFall)
        // this->displayMessage( "Gravity is making his job!" );
//    else
        // this->displayMessage( "Gravity stopped caring" );
    update();
    return this->applyLetItFall;*/
}

bool Viewer::createSandGravity()
{
    this->voxelGrid->letGravityMakeSandFall(true);
    update();
    return false;
    /*
    this->startAnimation();
    this->applyLetSandFall = !this->applyLetSandFall;
//    if (this->applyLetSandFall)
        // this->displayMessage( "Sand is falling!" );
//    else
        // this->displayMessage( "Sand stopped falling" );
    update();
    return this->applyLetSandFall;*/
}

bool Viewer::startStopRecording()
{
    if(!makedir(this->screenshotFolder)) {
        this->isTakingScreenshots = false;
        std::cerr << "Not possible to create folder " << this->screenshotFolder << std::endl;
        exit(-1);
    }
    this->isTakingScreenshots = !this->isTakingScreenshots;
    if (!isTakingScreenshots) {
        std::string command = "ffmpeg -f image2 -i ";
        command += this->screenshotFolder + "%d.jpg -framerate 10 " + this->screenshotFolder + "0.gif";
        if (this->screenshotIndex > 0) {
            int result = std::system(command.c_str());
            if (result != 0) {
                std::cerr << "Oups, the command `" << command << "` didn't finished as expected... maybe ffmpeg is not installed?" << std::endl;
            }
        }
        this->screenshotFolder += "__next-take";
        this->screenshotIndex = 0;
    }


    std::cout << (this->isTakingScreenshots ? "Smile, you're on camera" : "Ok, stop smiling, it's saved") << std::endl;
    update();
    return this->isTakingScreenshots;
}

void Viewer::loadMapUI()
{
    QString filename = QFileDialog::getOpenFileName(this, QString("Charger une carte"), QString::fromStdString(this->mapSavingFolder));

    const char* vShader_voxels = ":/src/Shaders/voxels_vertex_shader_blinn_phong.glsl";
    const char* fShader_voxels = ":/src/Shaders/voxels_fragment_shader_blinn_phong.glsl";
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    voxelGrid->retrieveMap(filename.toStdString());
    voxelGrid->fromIsoData();
    voxelGrid->displayWithMarchingCubes = (this->algorithm == MARCHING_CUBES);
    this->voxelGrid->createMesh();
    for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks)
        vc->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
}
void Viewer::saveMapUI()
{
    QString filename = QFileDialog::getSaveFileName(this, QString("Enregistrer la carte"), QString::fromStdString(this->mapSavingFolder));
    if (this->voxelGrid)
        voxelGrid->saveMap(filename.toStdString());
}
