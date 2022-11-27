#include "Utils/Globals.h"
#include "Interface/Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>
#include <chrono>
#include "TerrainModification/UnderwaterErosion.h"
#include "DataStructure/Matrix.h"
#include "Utils/Utils.h"
#include "Interface/TerrainGenerationInterface.h"
#include "Interface/VisitingCamera.h"
#include <QTemporaryDir>
#include "Graphics/RayMarching.h"

#ifdef linux
    #include "sys/stat.h"
#endif
Viewer::Viewer(QWidget *parent): Viewer(
        std::make_shared<Grid>(),
        std::make_shared<VoxelGrid>(),
        std::make_shared<LayerBasedGrid>(),
        VOXEL_MODE,
        FILL_MODE,
        parent
        )
{
    if (parent != nullptr)
        parent->installEventFilter(this);
    this->mainCamera = this->camera();
}
Viewer::Viewer(std::shared_ptr<Grid> grid, std::shared_ptr<VoxelGrid> voxelGrid,
               std::shared_ptr<LayerBasedGrid> layerGrid, MapMode map,
               ViewerMode mode, QWidget *parent)
    : QGLViewer(parent), viewerMode(mode), mapMode(map), grid(grid), voxelGrid(voxelGrid), layerGrid(layerGrid)
{
    if (parent != nullptr)
        parent->installEventFilter(this);
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
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_3D);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    this->setBackgroundColor(QColor(127, 127, 127));

    this->camera()->setType(qglviewer::Camera::PERSPECTIVE);

    this->setShortcut(KeyboardAction::MOVE_CAMERA_DOWN, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_UP, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_LEFT, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_RIGHT, 0);

    setTextIsEnabled(true);
    setMouseTracking(true);

#ifdef linux
    std::string pathToShaders = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/src/Shaders/"; // ":/src/Shaders/"
#else
    std::string pathToShaders = "C:/codes/Qt/generation-terrain-procedural/src/Shaders/"; // ":/src/Shaders/"
#endif
    std::string vShader_voxels = pathToShaders + "voxels.vert";
    std::string fShader_voxels = pathToShaders + "voxels.frag";
    std::string vNoShader = pathToShaders + "no_shader.vert";
    std::string fNoShader = pathToShaders + "no_shader.frag";

    std::string vRayMarch = pathToShaders + "test_raymarching_voxels.vert";
    std::string fRayMarch = pathToShaders + "test_raymarching_voxels.frag";
    this->raymarchingShader = std::make_shared<Shader>(vRayMarch, fRayMarch);
    this->raymarchingShader->compileShadersFromSource({

                                                      });
    this->raymarchingQuad = Mesh({Vector3(-1.0f, -1.0f, 0.f),
                                  Vector3(1.0f, -1.0f, 0.f),
                                  Vector3(1.0f,  1.0f, 0.f),
                                 Vector3(-1.0f, -1.0f, 0.f),
                                 Vector3( 1.0f,  1.0f, 0.f),
                                 Vector3(-1.0f,  1.0f)}, raymarchingShader);

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

    Mesh::setShaderToAllMeshesWithoutShader(*Shader::default_shader);
    GlobalsGL::f()->glBindVertexArray(raymarchingQuad.vao);
//    sceneUBO = ShaderUBO("scene_buf", 1, sizeof(rt_scene));
//    sceneUBO.affectShader(raymarchingQuad.shader);
    QGLViewer::init();
}

void Viewer::draw() {
//    QGLViewer::draw();
    /*
    GlobalsGL::f()->glBindVertexArray(raymarchingQuad.vao);
    raymarchingQuad.shader->use();
    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);
    glm::mat4 mv_matrix(mvMatrix[0], mvMatrix[1], mvMatrix[2], mvMatrix[3], mvMatrix[4], mvMatrix[5], mvMatrix[6], mvMatrix[7], mvMatrix[8], mvMatrix[9], mvMatrix[10], mvMatrix[11], mvMatrix[12], mvMatrix[13], mvMatrix[14], mvMatrix[15]);
    glm::mat4 proj_matrix(pMatrix[0], pMatrix[1], pMatrix[2], pMatrix[3], pMatrix[4], pMatrix[5], pMatrix[6], pMatrix[7], pMatrix[8], pMatrix[9], pMatrix[10], pMatrix[11], pMatrix[12], pMatrix[13], pMatrix[14], pMatrix[15]);
    glm::vec3 cam_pos = Vector3(camera()->position());
    rt_scene scene("scene", raymarchingQuad.shader);
    scene.proj_matrix = proj_matrix;
    scene.mv_matrix = mv_matrix;
    scene.bg_color =  glm::vec4(1.f, 0.f, 0.f, 1.f);
    scene.canvas_width = this->width();
    scene.canvas_height = this->height();
    scene.camera_pos = glm::vec4(cam_pos, 0.0);
    scene.update();

    rt_material mat;
    mat.color = glm::vec4(1.f, 0.f, 0.f, 1.f);
    mat.absorb = glm::vec4(.5f, .5f, .5f, 0.f);
    mat.diffuse = 1.f;
    mat.refraction = 0.f;
    mat.reflection = 0.f;
    mat.specular = 16;
    mat.kd = .5f;
    mat.ks = .5f;

    rt_sphere sphere("spheres[0]", raymarchingQuad.shader);
    sphere.mat = mat;
    sphere.obj = glm::vec4(20.f, 0.f, 0.f, 20.f);
    sphere.quat_rotation = glm::vec4(0, 0, 0, 1);
    sphere.textureNum = 0;
    sphere.hollow = false;
    sphere.update();

    rt_sphere sphere2("spheres[1]", raymarchingQuad.shader);
    sphere2.mat = mat;
    sphere2.obj = glm::vec4(30.f, 0.f, 0.f, 20.f);
    sphere2.quat_rotation = glm::vec4(0, 0, 0, 1);
    sphere2.textureNum = 0;
    sphere2.hollow = false;
    sphere2.update();

    rt_box box("boxes[0]", raymarchingQuad.shader);
    box.mat = mat;
    box.pos = glm::vec4(0.f, 20.f, 0.f, 0.f);
    box.form = glm::vec4(10.f, 30.f, 100.f, 0.f);
    box.quat_rotation = glm::vec4(0, 0, 0, 1);
    box.textureNum = 0;
    box.update();

    rt_light_direct light("lights_direct[0]", raymarchingQuad.shader);
    light.direction = glm::vec4(.1f, .1f, .9f, 1.f);
    light.color =  glm::vec4(1.f, 1.f, 1.f, 1.f);
    light.intensity = 1.f;
    light.update();

    raymarchingQuad.shader->setTexture3D("dataFieldTex", 0, voxelGrid->getVoxelValues() / 6.f + .5f);

    Matrix3<int> materials; Matrix3<float> matHeights;
    std::tie(materials, matHeights) = voxelGrid->getLayersRepresentations();
    raymarchingQuad.shader->setTexture3D("matIndicesTex", 1, materials);
    raymarchingQuad.shader->setTexture3D("matHeightsTex", 2, matHeights);

    this->raymarchingQuad.display();*/
    this->drawingProcess();
}
void Viewer::drawingProcess() {
    // Update the mouse position in the grid
    this->checkMouseOnVoxel();

    this->frame_num ++;
    glClear(GL_DEPTH_BUFFER_BIT);
    bool displayFill = this->viewerMode != ViewerMode::WIRE_MODE;
    if (!displayFill) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);

    this->light.position = Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ);
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
        shader->setInt("voxels_displayed_on_borders", voxelsSmoothedOnBorders);
        shader->setFloat("fogNear", this->fogNear);
        shader->setFloat("fogFar", this->fogFar);
        shader->setBool("wireframeMode", !displayFill);
    });
    current_frame ++;
    if (this->interfaces.count("terrainGenerationInterface")) {
        static_cast<TerrainGenerationInterface*>(this->interfaces["terrainGenerationInterface"].get())->display(this->mapMode, this->algorithm, this->displayParticles);
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    this->mainGrabber->display();

    for (auto& actionInterface : this->interfaces)
        if (actionInterface.first != "terrainGenerationInterface")
            actionInterface.second->display();

    if (this->interfaces.count("terrainGenerationInterface")) {
        static_cast<TerrainGenerationInterface*>(this->interfaces["terrainGenerationInterface"].get())->displayWaterLevel();
    }
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

void Viewer::reloadAllShaders()
{
    Shader::applyToAllShaders([](std::shared_ptr<Shader> shader) {
        shader->compileShadersFromSource();
    });
    for (auto& [name, action] : this->interfaces) {
        action->reloadShaders();
    }
}


void Viewer::mousePressEvent(QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
    checkMouseOnVoxel();
    if (this->mouseInWorld && e->button() == Qt::MouseButton::LeftButton) {
        if (this->mapMode == MapMode::VOXEL_MODE) {
            std::cout << "Voxel (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << this->voxelGrid->getVoxelValue(this->mousePosWorld) << std::endl;
        } else if (this->mapMode == MapMode::LAYER_MODE) {
            auto [mat, height] = this->layerGrid->getMaterialAndHeight(mousePosWorld);
            std::cout << "Stack (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << LayerBasedGrid::densityFromMaterial(mat) << " for height " << height << std::endl;
        } else if (this->mapMode == MapMode::GRID_MODE) {
            std::cout << "Vertex (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has height " << this->grid->getHeight(mousePosWorld) << std::endl;
        }

    }
    Q_EMIT this->mouseClickOnMap(this->mousePosWorld, this->mouseInWorld, e);
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
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
    try {
        QGLViewer::mouseMoveEvent(e);
    }  catch (std::exception) {
        std::cout << "Catched this f***ing exception!" << std::endl;
    }

    update();
}

void Viewer::mouseDoubleClickEvent(QMouseEvent *e)
{
    QGLViewer::mouseDoubleClickEvent(e);
    checkMouseOnVoxel();
    if (this->mouseInWorld && e->button() == Qt::MouseButton::LeftButton) {
        std::cout << "Voxel (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << this->voxelGrid->getVoxelValue(this->mousePosWorld) << std::endl;

    }
    Q_EMIT this->mouseDoubleClickedOnMap(this->mousePosWorld, this->mouseInWorld, e);
}

void Viewer::animate()
{
    QGLViewer::animate();
}

Vector3 Viewer::minVoxelsShown()
{
    Vector3 minVec(minSliceMapX, minSliceMapY, minSliceMapZ);
    if (this->mapMode == MapMode::VOXEL_MODE)
        return voxelGrid->getDimensions() * minVec;
    else if (this->mapMode == MapMode::LAYER_MODE)
        return layerGrid->getDimensions() * minVec;
    else if (this->mapMode == MapMode::GRID_MODE)
        return grid->getDimensions() * minVec;
    return Vector3();
}

Vector3 Viewer::maxVoxelsShown()
{
    Vector3 maxVec(maxSliceMapX, maxSliceMapY, maxSliceMapZ);
    if (this->mapMode == MapMode::VOXEL_MODE)
        return voxelGrid->getDimensions() * maxVec;
    else if (this->mapMode == MapMode::LAYER_MODE)
        return layerGrid->getDimensions() * maxVec;
    else if (this->mapMode == MapMode::GRID_MODE)
        return grid->getDimensions() * maxVec;
    return Vector3();
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
    }
    else {
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
bool Viewer::checkMouseOnVoxel()
{
    if (voxelGrid == nullptr)
        return false;
    camera()->convertClickToLine(mousePos, orig, dir);
    /*
    float maxDist = std::max((int)camera()->distanceToSceneCenter(), std::max(voxelGrid->getSizeX(), std::max(voxelGrid->getSizeY(), voxelGrid->getSizeZ())));
    maxDist *= maxDist;

    Vector3 minPos = minVoxelsShown(), maxPos = maxVoxelsShown();

    bool found = false;
    Vector3 currPos(orig.x, orig.y, orig.z);
    auto values = voxelGrid->getVoxelValues();
    values.raiseErrorOnBadCoord = false;
    while((currPos / 2.f).norm2() < maxDist && !found)
    {
        currPos += Vector3(dir.x, dir.y, dir.z);
        if (minPos.x <= currPos.x && currPos.x <= maxPos.x && minPos.y <= currPos.y && currPos.y <= maxPos.y && minPos.z <= currPos.z && currPos.z <= maxPos.z) {
            float isoval = values.at(currPos);//voxelGrid->getVoxelValue(currPos);
            if (isoval > 0.0) // || (0.f <= currPos.z && currPos.z < 1.0) )
                found = true;
        }
    }*/
    Vector3 currPos(false);
    if (this->mapMode == MapMode::VOXEL_MODE) {
        currPos = voxelGrid->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    } else if (this->mapMode == MapMode::GRID_MODE) {
        currPos = grid->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    } else if (this->mapMode == MapMode::LAYER_MODE) {
        currPos = layerGrid->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    }
    bool found = currPos.isValid();
//    std::cout << found << std::endl;
    this->mouseInWorld = found;
    if (found) {
        this->mousePosWorld = currPos;
        this->mainGrabber->move(currPos);
    }
    return found;
}

void Viewer::closeEvent(QCloseEvent *e) {
    this->setCamera(this->mainCamera);
    if (this->isTakingScreenshots) this->startStopRecording();
    QGLViewer::closeEvent(e);
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
bool Viewer::eventFilter(QObject* obj, QEvent* event)
{/*
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
    return false;*/
    return QGLViewer::eventFilter(obj, event);
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
