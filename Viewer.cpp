#include "Globals.h"
#include "Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <chrono>
#include "UnderwaterErosion.h"
#include "Matrix.h"

#include <sys/stat.h>


using namespace qglviewer;
using namespace std;

std::vector<std::string> split(std::string str, char c = ' ');
bool makedir(std::string path);

Viewer::Viewer(QWidget *parent):
    Viewer(
        nullptr, // new Grid(100, 100, 40, 1.0),
        new VoxelGrid(3, 3, 3, 1.0, .30),
        nullptr, // new LayerBasedGrid(10, 10, 50),
        VOXEL_MODE,
        FILL_MODE,
        parent
        )
{
}
Viewer::Viewer(Grid* grid, VoxelGrid* voxelGrid, LayerBasedGrid* layerGrid, MapMode map, ViewerMode mode, QWidget *parent)
    : QGLViewer(parent), viewerMode(mode), mapMode(map), grid(grid), voxelGrid(voxelGrid), layerGrid(layerGrid) {

}
Viewer::Viewer(Grid* g, QWidget *parent)
    : Viewer(g, nullptr, nullptr, GRID_MODE, FILL_MODE, parent) {

}
Viewer::Viewer(VoxelGrid* g, QWidget *parent)
    : Viewer(nullptr, g, nullptr, VOXEL_MODE, FILL_MODE, parent) {

}
Viewer::~Viewer()
{
}

void Viewer::init() {
    /*float totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for (VoxelChunk* vc : this->voxelGrid->chunks) {
            vc->applyToVoxels([](Voxel* v){});
        }
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Empty loop                     : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for (VoxelChunk* vc : this->voxelGrid->chunks) {
            vc->toFloat();
        }
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "To float                       : " << totalTime/10 << "s" << std::endl;
    std::vector<std::vector<std::vector<float>>> arr = voxelGrid->toFloat();
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for(int x = 0; x < voxelGrid->sizeX; x++)
            for(int y = 0; y < voxelGrid->sizeY; y++)
                for(int z = 0; z < voxelGrid->sizeZ; z++) {
                    Voxel* v = voxelGrid->getVoxel(x, y, z);
                    if(v) float a = v->getIsosurface();
                }
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "For on voxelGrid               : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for(int x = 0; x < voxelGrid->sizeX; x++)
            for(int y = 0; y < voxelGrid->sizeY; y++)
                for(int z = 0; z < voxelGrid->sizeZ; z++)
                    float a = arr[x][y][z];
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "For on float array             : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        voxelGrid->toFloat();
        for(int x = 0; x < voxelGrid->sizeX; x++)
            for(int y = 0; y < voxelGrid->sizeY; y++)
                for(int z = 0; z < voxelGrid->sizeZ; z++)
                    float a = arr[x][y][z];
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "For on float array + computing : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        voxelGrid->toVoxels();
        voxelGrid->toFloat();
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "To voxels + to array           : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        voxelGrid->toVoxels(arr);
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "To voxels with premade array   : " << totalTime/10 << "s" << std::endl;

    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for (VoxelChunk* vc : this->voxelGrid->chunks) {
            vc->applyToVoxels([](Voxel* v){
                v->parent->parent->getVoxel(v->globalPos());
            });
        }
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Fetching voxels from the grid  : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        for (VoxelChunk* vc : this->voxelGrid->chunks) {
            vc->applyMarchingCubes();
        }
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Marching cubes                 : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        this->voxelGrid->computeVoxelGroups();
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Computing groups               : " << totalTime/10 << "s" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        this->voxelGrid->displayWithMarchingCubes = true;
        this->voxelGrid->createMesh();
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Meshing with marching cubes    : " << totalTime/10 << "s" << std::endl;
    int nbTris = 0;
    for(VoxelChunk* vc : voxelGrid->chunks)
        nbTris += vc->mesh.vertexArray.size();
    std::cout << "(number of triangles = " << nbTris/3 << ")" << std::endl;
    totalTime = 0.0;
    for (int i = 0; i < 10; i ++) {
        auto startTime = std::chrono::system_clock::now();
        this->voxelGrid->displayWithMarchingCubes = false;
        this->voxelGrid->createMesh();
        auto endTime = std::chrono::system_clock::now();
        totalTime += std::chrono::duration<float>(endTime - startTime).count();
    }
    std::cout << "Meshing without marching cubes : " << totalTime/10 << "s" << std::endl;
    nbTris = 0;
        for(VoxelChunk* vc : voxelGrid->chunks)
            nbTris += vc->mesh.vertexArray.size();
        std::cout << "(number of triangles = " << nbTris/3 << ")" << std::endl;
    exit(0);*/

    restoreStateFromFile();
    setSceneRadius(500.0);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);

    this->camera()->setType(Camera::ORTHOGRAPHIC);

    setTextIsEnabled(true);
    setMouseTracking(true);

#ifdef _WIN32
    const char* vShader_grid = "C:/codes/Qt/generation-terrain-procedural/grid_vertex_shader_blinn_phong.glsl";
    const char* fShader_grid = "C:/codes/Qt/generation-terrain-procedural/grid_fragment_shader_blinn_phong.glsl";
    const char* vShader_voxels = "C:/codes/Qt/generation-terrain-procedural/voxels_vertex_shader_blinn_phong.glsl";
    const char* fShader_voxels = "C:/codes/Qt/generation-terrain-procedural/voxels_fragment_shader_blinn_phong.glsl";
    const char* vShader_layer = "C:/codes/Qt/generation-terrain-procedural/layer_based_vertex_shader.glsl";
    const char* fShader_layer = "C:/codes/Qt/generation-terrain-procedural/layer_based_fragment_shader.glsl";
    const char* vNoShader = "C:/codes/Qt/generation-terrain-procedural/no_vertex_shader.glsl";
    const char* fNoShader = "C:/codes/Qt/generation-terrain-procedural/no_fragment_shader.glsl";
#elif linux
    const char* vShader_grid = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/grid_vertex_shader_blinn_phong.glsl";
    const char* fShader_grid = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/grid_fragment_shader_blinn_phong.glsl";
    const char* vShader_voxels = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/voxels_vertex_shader_blinn_phong.glsl";
    const char* fShader_voxels = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/voxels_fragment_shader_blinn_phong.glsl";
    const char* vShader_layer = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/layer_based_vertex_shader.glsl";
    const char* fShader_layer = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/layer_based_fragment_shader.glsl";
    const char* vNoShader = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/no_vertex_shader.glsl";
    const char* fNoShader = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/no_fragment_shader.glsl";
#endif
    glEnable              ( GL_DEBUG_OUTPUT );
    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 );

    GlobalsGL::generateBuffers();
    this->rocksVBO = GlobalsGL::newBufferId();

    if (grid != nullptr) {
        this->grid->createMesh();
        this->grid->mesh.shader = new Shader(vShader_grid, fShader_grid);
    }
    if (layerGrid != nullptr) {
        this->layerGrid->createMesh();
        this->layerGrid->mesh.shader = new Shader(vShader_voxels, fShader_voxels);
    }
    if (voxelGrid != nullptr) {
        voxelGrid->displayWithMarchingCubes = (this->algorithm == MARCHING_CUBES);
        this->voxelGrid->createMesh();
        for(VoxelChunk* vc : this->voxelGrid->chunks)
            vc->mesh.shader = new Shader(vShader_voxels, fShader_voxels);
//        this->voxelGrid->mesh.shader = new Shader(vShader_voxels, fShader_voxels);
    }
    this->shader = Shader(vNoShader, fNoShader);
    this->rocksMeshes.shader = new Shader(vNoShader, fNoShader);
    this->flowDebugMeshes.shader = new Shader(vNoShader, fNoShader);

    this->matter_adder = RockErosion(this->erosionSize, 1.0);

    this->light = PositionalLight(
                new float[4]{.8, .8, .8, 1.},
                new float[4]{1., 1., 1., 1.},
                new float[4]{1., 1., 1., 1.},
                Vector3(100.0, 100.0, 0.0)
                );

    this->setAnimationPeriod(0);

    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);

#ifdef _WIN32
    this->screenshotFolder = "C:/codes/Qt/generation-terrain-procedural/screenshots/";
#elif linux
    this->screenshotFolder = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/screenshots/";
#endif
    if(!makedir(this->screenshotFolder)) {
        std::cerr << "Not possible to create folder " << this->screenshotFolder << std::endl;
        exit(-1);
    }
    if (this->voxelGrid != nullptr) {
        this->screenshotFolder += std::string(s_time) + "__" + voxelGrid->toShortString() + "/";
        this->displayMessage(QString::fromStdString(std::string("Screenshots will be saved in folder ") + std::string(this->screenshotFolder)));
    }

/*
    std::vector<Vector3> normals;
    for (int x = 0; x < this->voxelGrid->sizeX; x++) {
        for (int y = 0; y < this->voxelGrid->sizeY; y++) {
            for (int z = 0; z < this->voxelGrid->sizeZ; z++) {
                normals.push_back(Vector3(x, y, z) + .5 - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
                normals.push_back((Vector3(x, y, z) + this->voxelGrid->getFlowfield(x, y, z)) + .5 - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
            }
        }
    }
    this->flowDebugMeshes.fromArray(normals);
    this->flowDebugMeshes.update();*/
}

void Viewer::draw() {
    this->frame_num ++;
    glClear(GL_DEPTH_BUFFER_BIT);
    if (this->viewerMode == ViewerMode::WIRE_MODE)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//    this->shader.use(true);

    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);


    Material ground_material(
                    new float[4]{.24, .08, .02, 1.},
                    new float[4]{.30, .10, .04, 1.},
                    new float[4]{.62, .56, .37, 1.},
                    51.2f
                    );
    Material grass_material(
                    new float[4]{.14, .52, .00, 1.},
                    new float[4]{.16, .60, .00, 1.},
                    new float[4]{.62, .56, .37, 1.},
                    51.2f
                    );
//    this->light.position = Vector3(100.0 * std::cos(this->frame_num / (float)10), 100.0 * std::sin(this->frame_num / (float)10), 0.0);
    float globalAmbiant[4] = {.90, .90, .90, 1.0};

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
            for(VoxelChunk* vc : this->voxelGrid->chunks) {
                vc->mesh.shader->setMatrix("proj_matrix", pMatrix);
                vc->mesh.shader->setMatrix("mv_matrix", mvMatrix);
                vc->mesh.shader->setPositionalLight("light", this->light);
                vc->mesh.shader->setMaterial("ground_material", ground_material);
                vc->mesh.shader->setMaterial("grass_material", grass_material);
                vc->mesh.shader->setVector("globalAmbiant", globalAmbiant, 4);
                vc->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
                vc->mesh.shader->setBool("display_light_source", true);
                vc->mesh.shader->setFloat("offsetX", -this->voxelGrid->sizeX/2);
                vc->mesh.shader->setFloat("offsetY", -this->voxelGrid->sizeY/2);
            }
//            this->voxelGrid->mesh.shader->setMatrix("proj_matrix", pMatrix);
//            this->voxelGrid->mesh.shader->setMatrix("mv_matrix", mvMatrix);
//            this->voxelGrid->mesh.shader->setPositionalLight("light", this->light);
//            this->voxelGrid->mesh.shader->setMaterial("ground_material", ground_material);
//            this->voxelGrid->mesh.shader->setMaterial("grass_material", grass_material);
//            this->voxelGrid->mesh.shader->setVector("globalAmbiant", globalAmbiant, 4);
//            this->voxelGrid->mesh.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
//            this->voxelGrid->mesh.shader->setBool("display_light_source", true);
//            this->voxelGrid->mesh.shader->setFloat("offsetX", -this->voxelGrid->sizeX/2);
//            this->voxelGrid->mesh.shader->setFloat("offsetY", -this->voxelGrid->sizeY/2);
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
            this->layerGrid->mesh.shader->setFloat("offsetX", 0.0);
            this->layerGrid->mesh.shader->setFloat("offsetY", 0.0);
            this->layerGrid->display();
        }
    }

    this->rocksMeshes.shader->setMatrix("proj_matrix", pMatrix);
    this->rocksMeshes.shader->setMatrix("mv_matrix", mvMatrix);
    this->rocksMeshes.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
    this->flowDebugMeshes.shader->setMatrix("proj_matrix", pMatrix);
    this->flowDebugMeshes.shader->setMatrix("mv_matrix", mvMatrix);
    this->flowDebugMeshes.shader->setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
    if (displayRockTrajectories) {
        this->flowDebugMeshes.display(GL_LINES);
        this->rocksMeshes.display(GL_LINES);
    }


    if (this->isTakingScreenshots) {
#ifdef linux
        mode_t prevMode = umask(0011);
#endif
        if(this->screenshotIndex == 0)
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

    if (checkMouseOnVoxel())
    {
        Voxel* main_v = this->voxelGrid->getVoxel(this->mousePosWorld);
        this->displayMessage(QString::fromStdString(std::to_string( main_v->getIsosurface()) + " " + main_v->globalPos().toString() + " " + std::to_string( main_v->isOnGround) ));
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
        this->displayMessage(QString::fromStdString("Displaying using " + std::string(this->algorithm == MARCHING_CUBES ? " Marching cubes" : "no") + " algorithm") );
        update();
    } else if(e->key() == Qt::Key_V) {
        this->display_vertices = !this->display_vertices;
        update();
    } else if(e->key() == Qt::Key_P) {
        this->setAddingMatterMode(!this->addingMatterMode);
        this->displayMessage( (addingMatterMode ? "Construction mode" : "Destruction mode") );
        update();
    } else if(e->key() == Qt::Key_Return) {
        erodeMap(e->modifiers() == Qt::ShiftModifier);
    } else if(e->key() == Qt::Key_Minus) {
        this->setManualErosionRocksSize(max(2, this->erosionSize - 2));
        this->displayMessage(QString::fromStdString("Cursor size : " + std::to_string(this->manualErosionSize) ));
        update();
    } else if(e->key() == Qt::Key_Plus) {
        this->setManualErosionRocksSize(max(2, this->erosionSize + 2));
        this->displayMessage(QString::fromStdString("Cursor size : " + std::to_string(this->manualErosionSize) ));
        update();
    } else if(e->key() == Qt::Key_Space) {
        displayRockTrajectories = !displayRockTrajectories;
        this->displayMessage(QString::fromStdString("Rock trajectories are : " + std::string(displayRockTrajectories ? "ON" : "OFF") ));
        update();
    } else if(e->key() == Qt::Key_0) {
        this->createGlobalGravity();
    } else if(e->key() == Qt::Key_Comma) {
        this->createSandGravity();
    } else if(e->key() == Qt::Key_1) {
        this->startStopRecording();
    } else if(e->key() == Qt::Key_2) {
        for(VoxelChunk* vc : this->voxelGrid->chunks) {
            vc->LoDIndex++;
            vc->needRemeshing = true;
        }
        this->voxelGrid->remeshAll();
        update();
    } else if(e->key() == Qt::Key_3) {
        this->displayMessage( "Removing matter to create a tunnel" );
        createTunnel(true);
    } else if(e->key() == Qt::Key_4) {
        this->displayMessage( "Adding matter to create a tunnel" );
        createTunnel(false);
    } else {
        QGLViewer::keyPressEvent(e);
    }
}
void Viewer::mouseMoveEvent(QMouseEvent* e)
{
    this->mousePos = e->pos();

    QGLViewer::mouseMoveEvent(e);
}

void Viewer::animate()
{
    if (this->applyLetItFall)
        this->voxelGrid->makeItFall((this->applyLetSandFall ? -1.0 : 0.1));
    if (this->applyLetSandFall)
        this->voxelGrid->letGravityMakeSandFall();
}

void Viewer::erodeMap(bool sendFromCam)
{
    UnderwaterErosion erod(this->voxelGrid, this->erosionSize, this->erosionStrength, this->erosionQtt);
    Vector3 *pos = nullptr;
    Vector3* dir = nullptr;
    if (sendFromCam)
    {
        Vec a;
        Vec b;
        camera()->convertClickToLine(QPoint(camera()->screenWidth()/2, camera()->screenHeight()/2), a, b);
        pos = new Vector3(a.x, a.y, a.z);
        dir = new Vector3(b.x, b.y, b.z);
        this->displayMessage( "Rocks launched from camera!" );
    } else {
        pos = nullptr;
        dir = new Vector3(0.0, 0.0, 0.0);
        this->displayMessage( "Rocks launched!" );
    }
    this->lastRocksLaunched = erod.Apply(pos, dir, 10);
    this->rocksMeshes.vertexArrayFloat.clear();
    std::vector<float> oneThrow;
    for(std::vector<Vector3>& coords : this->lastRocksLaunched) {
        oneThrow = Vector3::toArray(coords);
        this->rocksMeshes.vertexArrayFloat.insert(this->rocksMeshes.vertexArrayFloat.end(), oneThrow.begin(), oneThrow.end());
    }
    this->rocksMeshes.update();
/*
    std::vector<Vector3> normals;
    for (int x = 0; x < this->voxelGrid->sizeX; x++) {
        for (int y = 0; y < this->voxelGrid->sizeY; y++) {
            for (int z = 0; z < this->voxelGrid->sizeZ; z++) {
                normals.push_back(Vector3(x, y, z) + .5 - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
                normals.push_back((Vector3(x, y, z) + this->voxelGrid->getFlowfield(x, y, z)) + .5 - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
            }
        }
    }
    this->flowDebugMeshes.fromArray(normals);
    this->flowDebugMeshes.update();
    update();*/
}

void Viewer::throwRock()
{
    if (this->mouseInWorld)
    {
        Voxel* main_v = this->voxelGrid->getVoxel(this->mousePosWorld);
        RockErosion rock(this->manualErosionSize, this->manualErosionStrength);
        rock.Apply(main_v, addingMatterMode);
    }
    update();
}
void Viewer::createTunnel(bool removingMatter)
{
    UnderwaterErosion erod(this->voxelGrid, this->curvesErosionSize, curvesErosionStrength, 10);
    this->rocksMeshes.fromArray(erod.CreateTunnel(nullptr, nullptr, 3, !removingMatter));
    this->rocksMeshes.update();
    update();
}

bool Viewer::checkMouseOnVoxel()
{
    if (voxelGrid == nullptr)
        return false;
    camera()->convertClickToLine(mousePos, orig, dir);
    float maxDist = max((int)camera()->distanceToSceneCenter(), max(voxelGrid->getSizeX(), max(voxelGrid->getSizeY(), voxelGrid->getSizeZ())));
    maxDist *= maxDist;

    bool found = false;
    Vector3 currPos(orig.x, orig.y, orig.z);
    currPos += Vector3(voxelGrid->getSizeX()/2, voxelGrid->getSizeY()/2, 0.0);
    while(currPos.norm() < maxDist)
    {
        currPos += Vector3(dir.x, dir.y, dir.z);
        Voxel* v = voxelGrid->getVoxel(currPos.x, currPos.y, currPos.z);
        if (v != nullptr && (bool)*v) {
            found = true;
            break;
        }
    }
    this->mouseInWorld = found;
    if (found) {
        this->mousePosWorld = currPos;
    }
    return found;
}

void Viewer::closeEvent(QCloseEvent *e) {
    QGLViewer::closeEvent(e);

    std::string command = "ffmpeg -f image2 -i ";
    command += this->screenshotFolder + "%d.jpg -framerate 10 " + this->screenshotFolder + "0.gif";
    if (this->screenshotIndex > 0) {
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Oups, the command `" << command << "` didn't finished as expected... maybe ffmpeg is not installed?" << std::endl;
        }
    }
    delete this->grid;
    delete this->voxelGrid;
    delete this->layerGrid;
}

bool Viewer::createGlobalGravity()
{
    this->startAnimation();
    this->applyLetItFall = !this->applyLetItFall;
    if (this->applyLetItFall)
        this->displayMessage( "Gravity is making his job!" );
    else
        this->displayMessage( "Gravity stopped caring" );
    update();
    return this->applyLetItFall;
}

bool Viewer::createSandGravity()
{
    this->startAnimation();
    this->applyLetSandFall = !this->applyLetSandFall;
    if (this->applyLetSandFall)
        this->displayMessage( "Sand is falling!" );
    else
        this->displayMessage( "Sand stopped falling" );
    update();
    return this->applyLetSandFall;
}

bool Viewer::startStopRecording()
{
    this->isTakingScreenshots = !this->isTakingScreenshots;
    if(!makedir(this->screenshotFolder)) {
        this->isTakingScreenshots = false;
        this->displayMessage(QString::fromStdString("Not possible to create folder " + this->screenshotFolder ));
        exit(-1);
    }
    this->displayMessage(this->isTakingScreenshots ? "Smile, you're on camera" : "Ok, stop smiling");
    update();
    return this->isTakingScreenshots;
}

std::vector<std::string> split(std::string str, char c)
{
    std::vector<std::string> result;
    size_t pos = str.npos;
    do {
        pos = str.rfind(c, pos);
        if (pos != str.npos) {
            std::string sub = str.substr(pos + 1);
            if(sub != "")
                result.insert(result.begin(), sub);
            str = str.substr(0, pos);
        }
    } while(pos != str.npos);
    result.insert(result.begin(), str);
    return result;
}

bool makedir(std::string path)
{
    std::vector<std::string> splitted = split(path, '/');
    int result = 0;
    std::string currentPath = "";
    for (size_t i = 0; i < splitted.size(); i++) {
        currentPath += splitted[i] + '/';
        struct stat info;
        if(stat(currentPath.c_str(), &info) != 0) { // Folder doesn't exist
#ifdef linux
            mode_t prevMode = umask(0011);
            result = mkdir(currentPath.c_str(), 0666); // Create with full permission for linux
            chmod(currentPath.c_str(), 0666);
            umask(prevMode);
#elif _WIN32
            result = mkdir(currentPath.c_str()); // Create for windows
#endif
            if (result != 0)
                return false;
        }
    }
    return true;
}
