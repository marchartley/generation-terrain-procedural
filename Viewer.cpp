#include "Globals.h"
#include "Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <chrono>
#include "UnderwaterErosion.h"
#include "Matrix.h"

#include <sys/stat.h>


using namespace qglviewer;
using namespace std;

std::vector<std::string> split(std::string str, char c = ' ')
{
    std::vector<std::string> result;
    bool insertHeadingSlash = str.size() > 0 && str[0] == '/';
    int pos = str.npos;
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
//    if(insertHeadingSlash)
//        result.insert(result.begin(), "/");// Case of the heading "/" for Linux
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

Viewer::Viewer(Grid* grid, VoxelGrid* voxelGrid, MapMode map, ViewerMode mode)
    : QGLViewer(), viewerMode(mode), mapMode(map), grid(grid), voxelGrid(voxelGrid) {

}
Viewer::Viewer(Grid* g)
    : Viewer(g, NULL, GRID_MODE, FILL_MODE) {

}
Viewer::Viewer(VoxelGrid* g)
    : Viewer(NULL, g, VOXEL_MODE, FILL_MODE) {

}

void Viewer::init() {
    restoreStateFromFile();
    setSceneRadius(200.0);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);

//    this->camera()->setType(Camera::ORTHOGRAPHIC);

    setMouseTracking(true);

#ifdef _WIN32
    const char* vShader = "C:/codes/Qt/generation-terrain-procedural/vertex_shader_blinn_phong.glsl";
    const char* fShader = "C:/codes/Qt/generation-terrain-procedural/fragment_shader_blinn_phong.glsl";
#elif linux
    const char* vShader = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/vertex_shader_blinn_phong.glsl";
    const char* fShader = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/fragment_shader_blinn_phong.glsl";
#endif
    this->shader = Shader(vShader, fShader);
    glEnable              ( GL_DEBUG_OUTPUT );
    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 );

    this->rocksVBO = GlobalsGL::newBufferId();

    this->grid->createMesh();
    voxelGrid->displayWithMarchingCubes = this->algorithm == MARCHING_CUBES;
    this->voxelGrid->createMesh();

    this->matter_adder = RockErosion(this->selectionSize, 1.0);

    this->light = PositionalLight(
                new float[4]{.8, .8, .8, 1.},
                new float[4]{1., 1., 1., 1.},
                new float[4]{1., 1., 1., 1.},
                Vector3(100.0, 100.0, -100.0)
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
    this->screenshotFolder += std::string(s_time) + "__" + voxelGrid->toShortString() + "/";
    std::cout << "Screenshots will be saved in folder " << this->screenshotFolder << std::endl;
}

void Viewer::draw() {
    glClear(GL_DEPTH_BUFFER_BIT);
    if (this->viewerMode == ViewerMode::WIRE_MODE)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    this->shader.use(true);

    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);

    this->shader.setMatrix("proj_matrix", pMatrix);
    this->shader.setMatrix("mv_matrix", mvMatrix);

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
    float globalAmbiant[4] = {.90, .90, .90, 1.0};

    this->shader.setPositionalLight("light", this->light);
    this->shader.setMaterial("ground_material", ground_material);
    this->shader.setMaterial("grass_material", grass_material);
    this->shader.setVector("globalAmbiant", globalAmbiant, 4);
    this->shader.setMatrix("norm_matrix", Matrix(4, 4, mvMatrix).transpose().inverse());
    this->shader.setBool("display_light_source", true);

    if (this->mapMode == GRID_MODE) {
        this->grid->display(true);
    }
    else if (this->mapMode == VOXEL_MODE) {
        this->voxelGrid->display();
    }

    if (displayRockTrajectories) {
        glPushMatrix();
        GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->rocksVBO]);
        this->shader.setFloat("offsetX", 0.f);
        this->shader.setFloat("offsetY", 0.f);
//        glTranslatef(this->voxelGrid->getSizeX()/2.0, this->voxelGrid->getSizeY()/2.0, 0); //, -this->voxelGrid->getSizeZ()/2.0);
        for(std::vector<Vector3>& coords : this->lastRocksLaunched) {
            glBegin(GL_LINE_STRIP);
            glColor4f(0.0, 0.0, 0.0, .5);
            for(Vector3& pos : coords) {
                glVertex3f(pos.x, pos.y, pos.z);
            }
            glEnd();
        }
        glPopMatrix();
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
        std::cout << main_v->getIsosurface() << " " << main_v->globalPos() << std::endl;
    }
    if (QApplication::keyboardModifiers().testFlag(Qt::AltModifier) == true)
    {
        if (this->mouseInWorld)
        {
            Voxel* main_v = this->voxelGrid->getVoxel(this->mousePosWorld);
            RockErosion rock(this->selectionSize, .50);
            rock.Apply(main_v, addingMatterMode);
        }
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
    } else if (e->key() == Qt::Key_R) {
        if (this->algorithm == NONE)
            setSmoothingAlgorithm(MARCHING_CUBES);
        else if (this->algorithm == MARCHING_CUBES)
            setSmoothingAlgorithm(NONE);
//        else if (this->algorithm == DUAL_CONTOURING)
//            setSmoothingAlgorithm(NONE);
        std::cout << "Displaying using " << (this->algorithm == MARCHING_CUBES ? " Marching cubes" : "no") << " algorithm" << std::endl;
        voxelGrid->displayWithMarchingCubes = this->algorithm == MARCHING_CUBES;
        voxelGrid->createMesh();
        update();
    } else if(e->key() == Qt::Key_V) {
        this->display_vertices = !this->display_vertices;
        update();
    } else if(e->key() == Qt::Key_P) {
        this->addingMatterMode = !this->addingMatterMode;
        std::cout << (addingMatterMode ? "Construction mode" : "Destruction mode") << std::endl;
        update();
    } else if(e->key() == Qt::Key_Return) {
        UnderwaterErosion erod(this->voxelGrid, 10, 2.0, 10);
        if (e->modifiers() == Qt::ShiftModifier)
        {
            Vector3 pos(camera()->position().x, camera()->position().y, camera()->position().z);
            this->lastRocksLaunched = erod.Apply(pos);
            std::cout << "Rocks launched from camera!" << std::endl;
        } else {
            this->lastRocksLaunched = erod.Apply();
            std::cout << "Rocks launched!" << std::endl;
        }
        update();
    } else if(e->key() == Qt::Key_Minus) {
        this->selectionSize = max(2, this->selectionSize - 2);
        std::cout << "Cursor size : " << this->selectionSize << std::endl;
        update();
    } else if(e->key() == Qt::Key_Plus) {
        this->selectionSize = max(2, this->selectionSize + 2);
        std::cout << "Cursor size : " << this->selectionSize << std::endl;
        update();
    } else if(e->key() == Qt::Key_Space) {
        displayRockTrajectories = !displayRockTrajectories;
        std::cout << "Rock trajectories are : " << (displayRockTrajectories ? "ON" : "OFF") << std::endl;
        update();
    } else if(e->key() == Qt::Key_0) {
        if (this->animationIsStarted())
            this->stopAnimation();
        else
            this->startAnimation();
//        this->voxelGrid->makeItFall();
        std::cout << "It's falling!" << std::endl;
        update();
    } else if(e->key() == Qt::Key_1) {
        this->isTakingScreenshots = !this->isTakingScreenshots;
        if(!makedir(this->screenshotFolder)) {
            std::cerr << "Not possible to create folder " << this->screenshotFolder << std::endl;
            exit(-1);
        }
        std::cout << (this->isTakingScreenshots ? "Smile, you're on camera" : "Ok, stop smiling") << std::endl;
        update();
    } else {
        QGLViewer::keyPressEvent(e);
    }
}
std::vector<std::vector<Vector3>> Viewer::getSphereVertices(int rings, int halves) {
    std::vector<std::vector<Vector3>> coords;
    std::vector<Vector3> tris;
    for(int i = 0; i < halves - 1; i++)
    {
        for (int j = 0; j < rings -1 ; j++)
        {
            Vector3 v1(i, j, 1.0);
            Vector3 v2(i+1, j, 1.0);
            Vector3 v3(i+1, j+1, 1.0);
            Vector3 v4(i, j+1, 1.0);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();

            v1 = Vector3(i, j, -1.0);
            v2 = Vector3(i+1, j, -1.0);
            v3 = Vector3(i+1, j+1, -1.0);
            v4 = Vector3(i, j+1, -1.0);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();

            v1 = Vector3(1.0, i, j);
            v2 = Vector3(1.0, i+1, j);
            v3 = Vector3(1.0, i+1, j+1);
            v4 = Vector3(1.0, i, j+1);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();

            v1 = Vector3(-1.0, i, j);
            v2 = Vector3(-1.0, i+1, j);
            v3 = Vector3(-1.0, i+1, j+1);
            v4 = Vector3(-1.0, i, j+1);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();

            v1 = Vector3(i, 1.0, j);
            v2 = Vector3(i+1, 1.0, j);
            v3 = Vector3(i+1, 1.0, j+1);
            v4 = Vector3(i, 1.0,  j+1);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();

            v1 = Vector3(i, -1.0, j);
            v2 = Vector3(i+1, -1.0, j);
            v3 = Vector3(i+1, -1.0, j+1);
            v4 = Vector3(i, -1.0,  j+1);
            v1.normalize(); v2.normalize(); v3.normalize(); v4.normalize();

            tris.push_back(v1); tris.push_back(v2); tris.push_back(v3);
            coords.push_back(tris);
            tris.clear();
        }
    }
    return coords;
}
void Viewer::drawSphere(float radius, int rings, int halves)
{
    std::vector<std::vector<Vector3>> coords = getSphereVertices(rings, halves);
    glBegin(GL_TRIANGLES);
    for (std::vector<Vector3> triangle : coords)
    {
        for(Vector3 vert : triangle)
            glVertex3f(vert.x * radius, vert.y * radius, vert.z * radius);
    }
    glEnd();
}
void Viewer::mouseMoveEvent(QMouseEvent* e)
{
    this->mousePos = e->pos();

    QGLViewer::mouseMoveEvent(e);
}

void Viewer::animate()
{
    this->voxelGrid->makeItFall();
}

bool Viewer::checkMouseOnVoxel()
{
    camera()->convertClickToLine(mousePos, orig, dir);
    float maxDist = max((int)camera()->distanceToSceneCenter(), max(voxelGrid->getSizeX(), max(voxelGrid->getSizeY(), voxelGrid->getSizeZ())));
    maxDist *= maxDist;

    bool found = false;
    Vector3 currPos(orig.x, orig.y, orig.z);
//    currPos += Vector3(voxelGrid->getSizeX()/2 - 1, voxelGrid->getSizeY()/2 - 1, 0.0); //, voxelGrid->getSizeZ()/2 - 1);
    currPos += Vector3(voxelGrid->getSizeX()/2, voxelGrid->getSizeY()/2, 0.0); //, voxelGrid->getSizeZ()/2 - 1);
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
        this->mousePosWorld = currPos; // + Vector3(voxelGrid->getSizeX()/2.0, voxelGrid->getSizeY()/2.0, voxelGrid->getSizeZ()/2.0);
    }
    return found;
}

void Viewer::closeEvent(QCloseEvent *e) {
    QGLViewer::closeEvent(e);

    std::string command = "ffmpeg -f image2 -i ";
    command += this->screenshotFolder + "%d.jpg -framerate 10 0.gif";
    if (this->screenshotIndex > 0) {
        int result = std::system(command.c_str());
        if (result != 0) {
            std::cerr << "Oups, the command `" << command << "` didn't finished as expected... maybe ffmpeg is not installed?" << std::endl;
        }
    }
}
