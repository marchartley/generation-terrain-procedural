#include "Utils/Globals.h"
#include "Interface/Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>
#include <chrono>
#include "TerrainModification/UnderwaterErosion.h"
#include "DataStructure/Matrix.h"
#include "Utils/Utils.h"
#include <QTemporaryDir>
#ifdef linux
    #include "sys/stat.h"
#endif
Viewer::Viewer(QWidget *parent): Viewer(
        std::make_shared<Grid>(),
        std::make_shared<VoxelGrid>(),
        std::shared_ptr<LayerBasedGrid>(nullptr), // new LayerBasedGrid(10, 10, 50),
        VOXEL_MODE,
        FILL_MODE,
        parent
        )
{
    if (parent != nullptr)
        parent->installEventFilter(this);
    this->mainCamera = this->camera();
    this->flyingCamera = new qglviewer::Camera(*this->mainCamera);
}
Viewer::Viewer(std::shared_ptr<Grid> grid, std::shared_ptr<VoxelGrid> voxelGrid,
               std::shared_ptr<LayerBasedGrid> layerGrid, MapMode map,
               ViewerMode mode, QWidget *parent)
    : QGLViewer(parent), viewerMode(mode), mapMode(map), grid(grid), voxelGrid(voxelGrid), layerGrid(layerGrid)
{
    if (parent != nullptr)
        parent->installEventFilter(this);
    this->mainCamera = this->camera();
    this->flyingCamera = new qglviewer::Camera(*this->mainCamera);
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
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_3D);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    this->setBackgroundColor(QColor(127, 127, 127));

    this->camera()->setType(qglviewer::Camera::PERSPECTIVE);

    setTextIsEnabled(true);
    setMouseTracking(true);

    const char* vShader_voxels = ":/src/Shaders/voxels.vert";
    const char* fShader_voxels = ":/src/Shaders/voxels.frag";
    const char* vNoShader = ":/src/Shaders/no_shader.vert";
    const char* fNoShader = ":/src/Shaders/no_shader.frag";

    glEnable              ( GL_DEBUG_OUTPUT );
    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 );

    Shader::default_shader = std::make_shared<Shader>(vNoShader, fNoShader);
    ControlPoint::base_shader = std::make_shared<Shader>(vNoShader, fNoShader);

    ControlPoint::base_shader->setVector("color", std::vector<float>({160/255.f, 5/255.f, 0/255.f, 1.f}));
    this->mainGrabber = new ControlPoint(Vector3(), 1.f, ACTIVE, false);


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
        this->grid->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
    }
    if (layerGrid != nullptr) {
        this->layerGrid->createMesh();
        this->layerGrid->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
    }
    if (voxelGrid != nullptr) {

    }

    QObject::connect(this->spaceColonizationInterface.get(), &SpaceColonizationInterface::useAsMainCamera, this, &Viewer::swapCamera);
    QObject::connect(this->karstPathInterface.get(), &KarstPathGenerationInterface::useAsMainCamera, this, &Viewer::swapCamera);

    Mesh::setShaderToAllMeshesWithoutShader(*Shader::default_shader);

//    startAnimation();
    QGLViewer::init();
}

void Viewer::draw() {
    // Update the mouse position in the grid
    this->checkMouseOnVoxel();

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

    this->light.position = Vector3(voxelGrid->sizeX / 2.0, voxelGrid->sizeY, voxelGrid->sizeZ * 2.0);
    //this->light.position = Vector3(camera()->frame()->position()) + Vector3(0, 0, 100);

    float white[4] = {240/255.f, 240/255.f, 240/255.f, 1.f};
    Material ground_material(
                    new float[4] {220/255.f, 210/255.f, 110/255.f, 1.f}, // new float[4]{.48, .16, .04, 1.},
                    new float[4] { 70/255.f,  80/255.f,  70/255.f, 1.f}, // new float[4]{.60, .20, .08, 1.},
                    new float[4] {0/255.f, 0/255.f, 0/255.f, 1.f}, // new float[4]{.62, .56, .37, 1.},
                    1.f // 51.2f
                    );
    Material grass_material(
                    new float[4] { 70/255.f,  80/255.f,  70/255.f, 1.f}, // new float[4]{.28, .90, .00, 1.},
                    new float[4] {220/255.f, 210/255.f, 160/255.f, 1.f}, // new float[4]{.32, .80, .00, 1.},
                    new float[4] {0/255.f, 0/255.f, 0/255.f, 1.f}, // new float[4]{.62, .56, .37, 1.},
                    1.f // 51.2f
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
        shader->setBool("isSpotlight", this->usingSpotlight);
        if (this->usingSpotlight) {
            shader->setVector("light.position", this->camera()->position());
        } else {
            shader->setVector("light.position", (this->light.position + Vector3(this->camera()->position())) / 2.f);
        }
        shader->setBool("display_light_source", true);
        shader->setVector("min_vertice_positions", minVoxelsShown());
        shader->setVector("max_vertice_positions", maxVoxelsShown());
        shader->setFloat("fogNear", this->fogNear);
        shader->setFloat("fogFar", this->fogFar);
    });
    current_frame ++;
    if (this->terrainGenerationInterface)
        terrainGenerationInterface->display(this->mapMode, this->algorithm, this->displayParticles);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    this->mainGrabber->display();

    if (this->karstPathInterface)
        this->karstPathInterface->display();
    if (this->spaceColonizationInterface)
        this->spaceColonizationInterface->display();
    if (this->faultSlipInterface)
        this->faultSlipInterface->display();
    if (this->tunnelInterface)
        this->tunnelInterface->display();
    if (this->flowFieldInterface)
        this->flowFieldInterface->display();
    if (this->manualEditionInterface)
        this->manualEditionInterface->display();
    if (this->erosionInterface)
        this->erosionInterface->display();
    if (this->heightmapErosionInterface)
        this->heightmapErosionInterface->display();
    if (this->biomeInterface)
        this->biomeInterface->display();

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

bool Viewer::inFlyMode()
{
    Qt::Key keyBinded;
    Qt::KeyboardModifiers modifierBinded;
    Qt::MouseButton buttonBinded;
    this->getMouseActionBinding(CAMERA, ROTATE, false, keyBinded, modifierBinded, buttonBinded);
    return buttonBinded == Qt::NoButton;
}

void Viewer::mousePressEvent(QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
    checkMouseOnVoxel();
    if (curvesErosionConstructionMode && this->mouseInWorld) {
//        this->addCurvesControlPoint(this->mousePosWorld);
    }
    if (QApplication::keyboardModifiers().testFlag(Qt::AltModifier) == true)
    {
        this->throwRock();
    }
    if (this->mouseInWorld && e->button() == Qt::MouseButton::LeftButton) {
        std::cout << "Voxel (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << this->voxelGrid->getVoxelValue(this->mousePosWorld) << std::endl;

    }
    Q_EMIT this->mouseClickOnMap(this->mousePosWorld, this->mouseInWorld, e);
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
//        createTunnel(true);
    } else if(e->key() == Qt::Key_4) {
        // this->displayMessage( "Adding matter to create a tunnel" );
//        createTunnel(false);
    } else if(e->key() == Qt::Key_5) {
        this->setCamera(this->spaceColonizationInterface->visitingCamera);
    } else if(e->key() == Qt::Key_C) {
        this->toggleCameraMode(); //this->swapCamera(this->flyingCamera, !this->usingMainCamera);
    }
    QGLViewer::keyPressEvent(e);
}

void Viewer::keyReleaseEvent(QKeyEvent *e)
{
    QGLViewer::keyReleaseEvent(e);
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

    Q_EMIT this->mouseMovedOnMap((this->mouseInWorld ? this->mousePosWorld : Vector3(-10000, -10000, -10000)));
    update();
    QGLViewer::mouseMoveEvent(e);
}

void Viewer::animate()
{
    /*if (voxelGrid) {
        voxelGrid->computeFlowfield();
        this->updateFlowfieldDebugMesh();
    }*//*
    if (this->applyLetItFall)
        this->voxelGrid->makeItFall((this->applyLetSandFall ? -1.0 : 0.1));
    if (this->applyLetSandFall)
        this->voxelGrid->letGravityMakeSandFall();*/
    QGLViewer::animate();
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
{/*
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
        // this->displayMessage( "Rocks launched from camera!" );*//*
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

//    updateFlowfieldDebugMesh();*/
}
/*
void Viewer::recomputeFlowfield()
{
    this->voxelGrid->computeFlowfield();
    updateFlowfieldDebugMesh();
}
*/
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
        this->usingMainCamera = false;
        this->setCamera(this->alternativeCamera);
        this->displayParticles = true;
        this->fogNear = 5.f;
        this->fogFar = 30.f;
        this->usingSpotlight = true;

//        if (!this->inFlyMode())
//            this->toggleCameraMode();
    }
    else {
//        if (this->inFlyMode())
//            this->toggleCameraMode();
        this->usingMainCamera = true;
        this->setCamera(this->mainCamera);
        this->displayParticles = false;
        this->fogNear = 1000.f;
        this->fogFar = 5000.f;
        this->usingSpotlight = false;
    }

    if (alternativeCamera != nullptr && dynamic_cast<VisitingCamera*>(alternativeCamera) != nullptr) {
        dynamic_cast<VisitingCamera*>(alternativeCamera)->isVisiting = useAltCamera;
    }
}

void Viewer::frameInterpolated()
{
}
/*
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
        this->debugMeshes[TUNNEL_PATHS].fromArray(erod.CreateTunnel(this->currentTunnelPoints, !removingMatter, false));
    this->currentTunnelPoints.clear();
    this->debugControlPoints[TUNNEL_PATHS].clear();
    this->debugMeshes[TUNNEL_PATHS].update();
    update();
}
void Viewer::createCrack(bool removingMatter)
{
    if (this->currentTunnelPoints.size() < 2) return;

    this->curvesErosionConstructionMode = false;
    UnderwaterErosion erod(this->voxelGrid, this->curvesErosionSize, curvesErosionStrength, 10);
    this->debugMeshes[TUNNEL_PATHS].fromArray(erod.CreateCrack(this->currentTunnelPoints[0], this->currentTunnelPoints[1], true));
    this->currentTunnelPoints.clear();
    this->debugControlPoints[TUNNEL_PATHS].clear();
    this->debugMeshes[TUNNEL_PATHS].update();
    update();
}
*/
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
        this->mousePosWorld = currPos;
        this->mainGrabber->move(currPos); // - Vector3(voxelGrid->getSizeX()/2, voxelGrid->getSizeY()/2, 0.0);
    }
    return found;
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

    const char* vShader_voxels = ":/src/Shaders/voxels.vert";
    const char* fShader_voxels = ":/src/Shaders/voxels.frag";

    QString q_filename = QFileDialog::getOpenFileName(this, QString("Charger une carte"), QString::fromStdString(this->mapSavingFolder));
    if (q_filename.isEmpty()) return; // Cancel the action if no file has been selected
    std::string filename = q_filename.toStdString();
    std::string ext = toUpper(getExtention(filename));
    if (!this->grid)
        this->grid = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    if (ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        grid->loadFromHeightmap(filename, 127, 127, 255);
        voxelGrid->from2DGrid(*grid);
        voxelGrid->fromIsoData();
    } else if (ext == "JSON") {
        // The JSON file contains the list of actions made on a map
        std::ifstream file(filename);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json json_content = nlohmann::json::parse(content);
        if (!json_content.contains("actions"))
            return;
        for (auto action : json_content.at("actions")) {
            // Let all the interfaces try to replay their actions
//            karstPathInterface->replay(action);
            spaceColonizationInterface->replay(action);
            faultSlipInterface->replay(action);
            gravityInterface->replay(action);
            tunnelInterface->replay(action);
            flowFieldInterface->replay(action);
            manualEditionInterface->replay(action);
            undoRedoInterface->replay(action);
        }
    } else {
        // Then it's our custom voxel grid file
        voxelGrid->retrieveMap(filename);
        voxelGrid->fromIsoData();
        grid->fromVoxelGrid(*voxelGrid);
    }
    this->setSceneCenter(voxelGrid->getDimensions() * voxelGrid->getBlockSize() / 2.f);
    voxelGrid->displayWithMarchingCubes = (this->algorithm == MARCHING_CUBES);
    this->voxelGrid->createMesh();
    this->grid->createMesh();
    for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks)
        vc->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
    update();
}
void Viewer::saveMapUI()
{
    QString filename = QFileDialog::getSaveFileName(this, QString("Enregistrer la carte"), QString::fromStdString(this->mapSavingFolder));
    if (this->voxelGrid)
        voxelGrid->saveMap(filename.toStdString());
}

bool Viewer::eventFilter(QObject* obj, QEvent* event)
{
    if (event->type() == QEvent::KeyPress)
        this->keyPressEvent(static_cast<QKeyEvent *>(event));
    if (event->type() == QEvent::KeyRelease)
        this->keyReleaseEvent(static_cast<QKeyEvent *>(event));
    if (event->type() == QEvent::Shortcut)
        this->keyReleaseEvent(static_cast<QKeyEvent *>(event));
    if (event->type() == QEvent::MouseMove)
        this->mouseMoveEvent(static_cast<QMouseEvent *>(event));
    if (event->type() == QEvent::MouseButtonPress)
        this->mousePressEvent(static_cast<QMouseEvent *>(event));
    if (event->type() == QEvent::MouseButtonRelease)
        this->mouseReleaseEvent(static_cast<QMouseEvent *>(event));
    if (event->type() == QEvent::MouseButtonDblClick)
        this->mouseDoubleClickEvent(static_cast<QMouseEvent *>(event));
    if (event->type() == QEvent::Wheel)
        this->wheelEvent(static_cast<QWheelEvent *>(event));
    if (event->type() == QEvent::Timer)
        this->timerEvent(static_cast<QTimerEvent *>(event));

    // Don't block any event
    return false;
}

void Viewer::clipViewTemporarily(Vector3 direction, Vector3 center, bool active)
{
    if (direction.norm2() > 0) {
        Vector3 cameraDirection = this->camera()->viewDirection();
        if (cameraDirection.dot(direction) > 1e-2)
            direction *= -1.f; // The clipping will happen from the camera point of view

        this->clipPlaneDirection = direction;
        this->clipPlanePosition = center;
    }
    this->temporaryClipPlaneActivated = active;
    this->update();
}
