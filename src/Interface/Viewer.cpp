#include "Interface/TerrainSavingInterface.h"
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

// ffmpeg -f image2 -i %d.png  -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -preset ultrafast -qp 0 out.mp4

#ifdef linux
    #include "sys/stat.h"
#endif
Viewer::Viewer(QWidget *parent): Viewer(
        std::make_shared<Heightmap>(),
        std::make_shared<VoxelGrid>(),
        std::make_shared<LayerBasedGrid>(),
        nullptr, // std::make_shared<ImplicitPatch>(new ImplicitPrimitive()),
        VOXEL_MODE,
        FILL_MODE,
        parent
        )
{
    if (parent != nullptr)
        parent->installEventFilter(this);
    this->mainCamera = this->camera();
}
Viewer::Viewer(std::shared_ptr<Heightmap> grid, std::shared_ptr<VoxelGrid> voxelGrid,
               std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch, MapMode map,
               ViewerMode mode, QWidget *parent)
    : QGLViewer(parent), viewerMode(mode), mapMode(map), heightmap(grid), voxelGrid(voxelGrid), layerGrid(layerGrid), implicitTerrain(implicitPatch)
{
    if (parent != nullptr)
        parent->installEventFilter(this);
    this->mainCamera = this->camera();
}
Viewer::Viewer(std::shared_ptr<Heightmap> g, QWidget *parent)
    : Viewer(g, nullptr, nullptr, nullptr, GRID_MODE, FILL_MODE, parent) {

}
Viewer::Viewer(std::shared_ptr<VoxelGrid> g, QWidget *parent)
    : Viewer(nullptr, g, nullptr, nullptr, LAYER_MODE, FILL_MODE, parent) {

}
Viewer::~Viewer()
{
}

void Viewer::init() {
    QGLViewer::init();
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
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    this->setBackgroundColor(QColor(127, 127, 127, 0));
//    this->setBackgroundColor(QColor(0, 0, 0));

    this->camera()->setType(qglviewer::Camera::PERSPECTIVE);

    this->setShortcut(KeyboardAction::MOVE_CAMERA_DOWN, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_UP, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_LEFT, 0);
    this->setShortcut(KeyboardAction::MOVE_CAMERA_RIGHT, 0);

    setTextIsEnabled(true);
    setMouseTracking(true);

    std::string pathToShaders = "src/Shaders/";

//    std::string vShader_voxels = pathToShaders + "voxels.vert";
//    std::string fShader_voxels = pathToShaders + "voxels.frag";
    std::string vNoShader = pathToShaders + "no_shader.vert";
    std::string fNoShader = pathToShaders + "no_shader.frag";

//    std::string vRayMarch = pathToShaders + "test_raymarching_voxels.vert";
//    std::string fRayMarch = pathToShaders + "test_raymarching_voxels.frag";
//    this->raymarchingShader = std::make_shared<Shader>(vRayMarch, fRayMarch);
//    this->raymarchingShader->compileShadersFromSource({
//                                                      });
//    this->raymarchingQuad = Mesh({Vector3(-1.0f, -1.0f, 0.f),
//                                  Vector3(1.0f, -1.0f, 0.f),
//                                  Vector3(1.0f,  1.0f, 0.f),
//                                 Vector3(-1.0f, -1.0f, 0.f),
//                                 Vector3( 1.0f,  1.0f, 0.f),
//                                 Vector3(-1.0f,  1.0f)}, raymarchingShader);

    glEnable              ( GL_DEBUG_OUTPUT );
//    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 ); // TODO : Add back

    Shader::default_shader = std::make_shared<Shader>(vNoShader, fNoShader);
    this->mainGrabber = std::make_unique<ControlPoint>(Vector3(), 1.f, ACTIVE, false);


    this->light = PositionalLight(
                {.5, .5, .5, 1.},
                {.2, .2, .2, 1.},
                {.5, .5, .5, 1.},
                Vector3(0.0, 0.0, -100.0)
                );

    this->setAnimationPeriod(0);

    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);


    this->main_screenshotFolder = "screenshots/";
    this->screenshotFolder = main_screenshotFolder;
    this->mapSavingFolder = "saved_maps/";
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
    }

    this->camera()->setViewDirection(qglviewer::Vec(-0.334813, -0.802757, -0.493438));
    this->camera()->setPosition(qglviewer::Vec(58.6367, 126.525002, 80.349899));
//    QGLViewer::init();
}

void Viewer::draw() {
    this->drawingProcess();
    this->window()->setWindowTitle("Simulation - " + QString::number(this->currentFPS()) + "FPS");
}

void Viewer::saveScreenshotPNG(std::string filename)
{
    int w = this->camera()->screenWidth();
    int h = this->camera()->screenHeight();
    int nbComp = 4;
    int size = w * h * nbComp;
    int newWidth = w * .5f;
    int newHeight = h * .5f;

    if (w != currentWidth || h != currentHeight) {
        if (currentWidth != -1 && currentHeight != -1){
            delete[] buffer;
            delete[] resized;
            delete[] flipped;
        }
        currentWidth = w;
        currentHeight = h;
        buffer = new GLubyte[size];
        resized = new GLubyte[newWidth * newHeight * nbComp];
        flipped = new GLubyte[newWidth * newHeight * nbComp];
    }

    GlobalsGL::f()->glPixelStorei(GL_PACK_ALIGNMENT, 1);
    GlobalsGL::f()->glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

    stbir_resize_uint8(buffer, w, h, 0, resized, newWidth, newHeight, 0, nbComp);

    // Flip the image vertically
    for (int y = 0; y < newHeight; y++) {
        for (int x = 0; x < newWidth; x++) {
            // Calculate the original and target positions
            int originalIndex = (y * newWidth + x) * nbComp;
            int targetIndex = ((newHeight - 1 - y) * newWidth + x) * nbComp;

            // Copy the pixel data
            for (int c = 0; c < nbComp; c++) {
                flipped[targetIndex + c] = resized[originalIndex + c];
            }
        }
    }
    stbi_write_png(filename.c_str(), newWidth, newHeight, nbComp, flipped, newWidth * nbComp);
}

void Viewer::copyLastScreenshotTo(std::string filename)
{
    std::ifstream  src(".tmp/screenshots/screen.png", std::ios::binary);
    std::ofstream  dst(filename,   std::ios::binary);

    dst << src.rdbuf();
    src.close();
    dst.close();
}

TerrainModel *Viewer::getCurrentTerrainModel()
{
    if (this->mapMode == VOXEL_MODE) {
        return this->voxelGrid.get();
    } else if (this->mapMode == GRID_MODE) {
        return this->heightmap.get();
    } else if (this->mapMode == LAYER_MODE) {
        return this->layerGrid.get();
    } else if (this->mapMode == IMPLICIT_MODE) {
        return this->implicitTerrain.get();
    }
    return nullptr;
}
void Viewer::drawingProcess() {
    /*GLuint depthBuffer; // Depth buffer for the FBO

    if (fbo == 0) {
        // Generate and bind the FBO
        GlobalsGL::f()->glGenFramebuffers(1, &fbo);
        GlobalsGL::f()->glBindFramebuffer(GL_FRAMEBUFFER, fbo);

        // Generate texture
        GlobalsGL::f()->glGenTextures(1, &texture);
        GlobalsGL::f()->glBindTexture(GL_TEXTURE_2D, texture);
        GlobalsGL::f()->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        GlobalsGL::f()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        GlobalsGL::f()->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // Attach the texture to the FBO
        GlobalsGL::f()->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

        // Create depth buffer and attach to FBO
        GlobalsGL::f()->glGenRenderbuffers(1, &depthBuffer);
        GlobalsGL::f()->glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
        GlobalsGL::f()->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
        GlobalsGL::f()->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);

        // Check if FBO is complete
        if (GlobalsGL::f()->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "Error setting up framebuffer!" << std::endl;
        }

        GlobalsGL::f()->glBindRenderbuffer(GL_RENDERBUFFER, 0); // Unbind the renderbuffer
    }

    GlobalsGL::f()->glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    GlobalsGL::f()->glViewport(0, 0, w, h); // Set the viewport
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
*/

    this->setSceneCenter(voxelGrid->getDimensions() / 2.f);
    auto allProcessStart = std::chrono::system_clock::now();
    std::map<std::shared_ptr<ActionInterface>, float> interfacesTimings;
//    std::chrono::milliseconds test;
    // Update the mouse position in the grid
    this->checkMouseOnVoxel();
    float mouseCheckTiming = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - allProcessStart).count();

    this->frame_num ++;
//    glClear(GL_DEPTH_BUFFER_BIT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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

    this->light.position = voxelGrid->getDimensions() * Vector3(.5f, .5f, 1.5f);
    Material ground_material(
                    {220/255.f, 210/255.f, 110/255.f, 1.f}, // new float[4]{.48, .16, .04, 1.},
                    { 70/255.f,  80/255.f,  70/255.f, 1.f}, // new float[4]{.60, .20, .08, 1.},
                    {0/255.f, 0/255.f, 0/255.f, 1.f}, // new float[4]{.62, .56, .37, 1.},
                    1.f // 51.2f
                    );
    Material grass_material(
                    { 70/255.f,  80/255.f,  70/255.f, 1.f}, // new float[4]{.28, .90, .00, 1.},
                    {220/255.f, 210/255.f, 160/255.f, 1.f}, // new float[4]{.32, .80, .00, 1.},
                    {0/255.f, 0/255.f, 0/255.f, 1.f}, // new float[4]{.62, .56, .37, 1.},
                    1.f // 51.2f
                    );
//    this->light.position = Vector3(100.0 * std::cos(this->frame_num / (float)10), 100.0 * std::sin(this->frame_num / (float)10), 0.0);
    Vector4 globalAmbiant = {.10, .10, .10, 1.0};

    Shader::applyToAllShaders([&](std::shared_ptr<Shader> shader) -> void {
        shader->setMatrix("proj_matrix", pMatrix);
        shader->setMatrix("mv_matrix", mvMatrix);
        shader->setPositionalLight("light", this->light);
        shader->setMaterial("ground_material", ground_material);
        shader->setMaterial("grass_material", grass_material);
        shader->setVector("globalAmbiant", globalAmbiant);
        shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
        shader->setBool("isSpotlight", this->usingSpotlight);
        if (this->usingSpotlight) {
            shader->setVector("light.position", this->camera()->position());
        } else {
            shader->setVector("light.position", this->light.position); //(this->light.position + Vector3(this->camera()->position())) / 2.f);
        }
        shader->setBool("display_light_source", true);
        shader->setVector("min_vertice_positions", minVoxelsShown());
        shader->setVector("max_vertice_positions", maxVoxelsShown());
        shader->setInt("voxels_displayed_on_borders", voxelsSmoothedOnBorders);
        shader->setFloat("fogNear", this->fogNear);
        shader->setFloat("fogFar", this->fogFar);
        shader->setBool("wireframeMode", !displayFill);


        Vector3 terrainMid = this->getCurrentTerrainModel()->getDimensions() * .5f;
        shader->setPositionalLight("lights[0]", this->light);
        shader->setVector("lights[0].position", terrainMid + Vector3(-100, 100, 200));
        shader->setPositionalLight("lights[1]", this->light);
        shader->setVector("lights[1].position", terrainMid + Vector3(100, 100, 0));
        shader->setPositionalLight("lights[2]", this->light);
        shader->setVector("lights[2].position", terrainMid + Vector3(0, 100, 200));
        shader->setPositionalLight("lights[3]", this->light);
        shader->setVector("lights[3].position", terrainMid + Vector3(0, -100, 100));
        shader->setPositionalLight("lights[4]", this->light);
        shader->setVector("lights[4].position", terrainMid + Vector3(100, -100, 200));
        shader->setPositionalLight("lights[5]", this->light);
        shader->setVector("lights[5].position", terrainMid + Vector3(-100, -100, -10));
    });
    current_frame ++;

    if (this->interfaces.count("terraingeneration")) {
        static_cast<TerrainGenerationInterface*>(this->interfaces["terraingeneration"].get())->setVisu(this->mapMode, this->algorithm, this->displayParticles);
        interfacesTimings[this->interfaces["terraingeneration"]] = timeIt([&]() { this->interfaces["terraingeneration"]->display();}); // std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    if (mouseDown) {
        this->mainGrabber->display();
    }

    for (auto& actionInterface : this->interfaces) {
        if (actionInterface.first != "terraingeneration") {
            interfacesTimings[actionInterface.second] = timeIt([&]() { actionInterface.second->display(this->camera()->position()); });
        }
    }

    if (this->interfaces.count("terraingeneration")) {
        interfacesTimings[this->interfaces["terraingeneration"]] += timeIt([&]() { static_cast<TerrainGenerationInterface*>(this->interfaces["terraingeneration"].get())->displayWaterLevel(); }); //std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    }

    auto startScreenSaving = std::chrono::system_clock::now();
    if (this->isTakingScreenshots) {
        if (makedir(".tmp/screenshots")) {
            this->saveScreenshotPNG(".tmp/screenshots/screen.png");
        }
        this->copyLastScreenshotTo(this->screenshotFolder + std::to_string(this->screenshotIndex++) + ".png");
    }
    float screenSavingTiming = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startScreenSaving).count();

    auto allProcessEnd = std::chrono::system_clock::now();

    bool displayTiming = false;
    if (displayTiming) {
        std::cout << "Total time/frame : " << showTime(std::chrono::duration_cast<std::chrono::milliseconds>(allProcessEnd - allProcessStart).count()) << std::endl;
        for (auto& [interf, time] : interfacesTimings) {
            std::cout << "\t" << interf->actionType << " : " << showTime(time) << std::endl;
        }
        std::cout << "\tMouse check : " << showTime(mouseCheckTiming) << std::endl;
        std::cout << "\tScreen saving : " << showTime(screenSavingTiming) << std::endl;
    }
//    GlobalsGL::f()->glBindFramebuffer(GL_FRAMEBUFFER, 0); // Unbind the FBO for subsequent rendering
}

void Viewer::reloadAllShaders()
{
    Shader::applyToAllShaders([](std::shared_ptr<Shader> shader) {
        shader->compileShadersFromSource();
    });
//    for (auto& [name, action] : this->interfaces) {
//        action->reloadShaders();
//    }
}

void Viewer::setupViewFromFile(std::string filename)
{
//    std::cout << "Using view from file " << filename << std::endl;
    this->setStateFileName(QString::fromStdString(filename));
    this->restoreStateFromFile();
    this->setStateFileName("");
//    this->saveStateToFile();
}

void Viewer::saveViewToFile(std::string filename)
{
    this->setStateFileName(QString::fromStdString(filename));
    this->saveStateToFile();
//    this->restoreStateFromFile();
    this->setStateFileName("");
}

void Viewer::screenshot()
{
    makedir(this->main_screenshotFolder + "shots/");
    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
    bool wasRecording = this->isTakingScreenshots;
    if (!wasRecording) {
        this->startRecording();
        this->draw();
    }
    this->copyLastScreenshotTo(this->main_screenshotFolder + "shots/" + s_time + ".png");
    if (!wasRecording)
        this->stopRecording();

    if (this->mapMode == MapMode::GRID_MODE) {
        this->heightmap->saveHeightmap(this->main_screenshotFolder + "shots/" + s_time + "-heightmap.png");
    } else if (this->mapMode == MapMode::VOXEL_MODE) {
        Heightmap tmp;
        tmp.fromVoxelGrid(*voxelGrid);
        tmp.saveHeightmap(this->main_screenshotFolder + "shots/" + s_time + "-heightmap_from_voxels.png");
//        this->voxelGrid->saveHeightmap(this->main_screenshotFolder + "shots/" + s_time + "-heightmap_from_voxels.png");
    }
    //    dynamic_cast<TerrainSavingInterface*>(this->interfaces["terrainSavingInterface"].get())->quickSaveAt(this->main_screenshotFolder + "shots", s_time, true, true, false);
}

void Viewer::resetScreenshotFolderName()
{
    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
    this->screenshotFolder = "screenshots/" + std::string(s_time) + "/";
    makedir(this->main_screenshotFolder);
}


void Viewer::mousePressEvent(QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
    this->mouseDown = true;
    checkMouseOnVoxel();
    Vector3 terrainScale = this->getCurrentTerrainModel()->scaling;
    Vector3 terrainTranslate = this->getCurrentTerrainModel()->translation;
    if (this->mouseInWorld && e->button() == Qt::MouseButton::LeftButton) {
        if (this->mapMode == MapMode::VOXEL_MODE) {
            std::cout << "Voxel (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << this->voxelGrid->getVoxelValue(this->mousePosWorld) << std::endl;
        } else if (this->mapMode == MapMode::LAYER_MODE) {
            auto [mat, height] = this->layerGrid->getMaterialAndHeight(mousePosWorld);
            std::cout << "Stack (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << LayerBasedGrid::densityFromMaterial(mat) << " for height " << height << std::endl;
        } else if (this->mapMode == MapMode::GRID_MODE) {
            std::cout << "Vertex (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has height " << this->heightmap->getHeight(mousePosWorld) << std::endl;
        } else if (this->mapMode == MapMode::IMPLICIT_MODE) {
            std::cout << "Implicit surface at " << mousePosWorld << std::endl;
        }
    }
    Q_EMIT this->mouseClickOnMap(this->mousePosWorld, this->mouseInWorld, e, this->getCurrentTerrainModel());
}

void Viewer::mouseReleaseEvent(QMouseEvent *e)
{
    QGLViewer::mouseReleaseEvent(e);
    this->mouseDown = false;
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
    auto start = std::chrono::system_clock::now();

    this->mousePos = e->pos();

    if (this->checkMouseOnVoxel())
    {
        this->mainGrabber->move(this->mousePosWorld);
        this->mainGrabber->setState(ACTIVE);
    } else {
        this->mainGrabber->setState(HIDDEN);
    }

    Vector3 terrainScale = this->getCurrentTerrainModel()->scaling;
    Vector3 terrainTranslate = this->getCurrentTerrainModel()->translation;
    Q_EMIT this->mouseMovedOnMap((this->mouseInWorld ? this->mousePosWorld : Vector3(-10000, -10000, -10000)), this->getCurrentTerrainModel());
    try {
        QGLViewer::mouseMoveEvent(e);
    }  catch (std::exception) {
        std::cout << "Catched this f***ing exception!" << std::endl;
    }

    update();

    auto end = std::chrono::system_clock::now();
//    std::cout << "Total time on MouseMoveEvent : " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}

void Viewer::mouseDoubleClickEvent(QMouseEvent *e)
{
    QGLViewer::mouseDoubleClickEvent(e);
    checkMouseOnVoxel();
    if (this->mouseInWorld && e->button() == Qt::MouseButton::LeftButton) {
        std::cout << "Voxel (" << int(mousePosWorld.x) << ", " << int(mousePosWorld.y) << ", " << int(mousePosWorld.z) << ") has value " << this->voxelGrid->getVoxelValue(this->mousePosWorld) << std::endl;

    }
    Q_EMIT this->mouseDoubleClickedOnMap(this->mousePosWorld, this->mouseInWorld, e, this->getCurrentTerrainModel());
}

void Viewer::animate()
{
    QGLViewer::animate();
    draw();
}

Vector3 Viewer::minVoxelsShown()
{
    Vector3 minVec(minSliceMapX, minSliceMapY, minSliceMapZ);
    if (this->mapMode == MapMode::VOXEL_MODE)
        return voxelGrid->getDimensions() * minVec;
    else if (this->mapMode == MapMode::LAYER_MODE)
        return layerGrid->getDimensions() * minVec;
    else if (this->mapMode == MapMode::GRID_MODE)
        return heightmap->getDimensions() * minVec;
    else if (this->mapMode == MapMode::IMPLICIT_MODE)
        return Vector3(); // implicitTerrain->getDimensions() * minVec;
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
        return heightmap->getDimensions() * maxVec;
    else if (this->mapMode == MapMode::IMPLICIT_MODE)
        return voxelGrid->getDimensions() * maxVec; //heightmap->getDimensions() * maxVec;
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

    Vector3 currPos(false);
    if (this->mapMode == MapMode::VOXEL_MODE) {
        currPos = voxelGrid->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    } else if (this->mapMode == MapMode::GRID_MODE) {
        currPos = heightmap->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    } else if (this->mapMode == MapMode::LAYER_MODE) {
        currPos = layerGrid->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
    } else if (this->mapMode == MapMode::IMPLICIT_MODE) {
        currPos = implicitTerrain->getIntersection(Vector3(orig.x, orig.y, orig.z), Vector3(dir.x, dir.y, dir.z), this->minVoxelsShown(), this->maxVoxelsShown());
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


bool Viewer::startRecording(std::string folderUsed)
{
    if (folderUsed != "")
        this->screenshotFolder = folderUsed;

    if(!makedir(this->screenshotFolder)) {
        this->isTakingScreenshots = false;
        std::cerr << "Not possible to create folder " << this->screenshotFolder << std::endl;
        exit(-1);
    }
    this->isTakingScreenshots = true;
    return this->isTakingScreenshots;
}

bool Viewer::stopRecording()
{
    std::string command = "ffmpeg -f image2 -i ";
//    command += this->screenshotFolder + "%d.png -framerate 10 " + this->screenshotFolder + "0.gif";
    command += this->screenshotFolder + "%d.png  -c:v libx264 -preset ultrafast -qp 0 " + this->screenshotFolder + "out.mp4";
    if (this->screenshotIndex > 0) {
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Oups, the command `" << command << "` didn't finished as expected... maybe ffmpeg is not installed?" << std::endl;
        }
    }
    this->screenshotFolder += "__next-take/";
    this->screenshotIndex = 0;

    this->isTakingScreenshots = false;
    return this->isTakingScreenshots;
}

bool Viewer::startStopRecording()
{
    if (!this->isTakingScreenshots) {
        std::cout << "Smile, you're on camera!" << std::endl;
        return this->startRecording();
    }
    else {
        std::cout << "Ok, you can stop smiling" << std::endl;
        return this->stopRecording();
    }
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

void Viewer::clipViewTemporarily(const Vector3& direction, const Vector3& center, bool active)
{
    if (direction.norm2() > 0) {
        Vector3 cameraDirection = this->camera()->viewDirection();
        if (cameraDirection.dot(direction) > 1e-2)
            this->clipPlaneDirection = -direction; // The clipping will happen from the camera point of view
        else
            this->clipPlaneDirection = direction;
        this->clipPlanePosition = center;
    }
    this->temporaryClipPlaneActivated = active;
    this->update();
}
