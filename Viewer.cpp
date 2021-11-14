#include "Globals.h"
#include "Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <chrono>
#include "UnderwaterErosion.h"



using namespace qglviewer;
using namespace std;

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
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable( GL_BLEND );
//    glEnable(GL_AUTO_NORMAL);
//    this->camera()->setType(Camera::ORTHOGRAPHIC);

    setMouseTracking(true);

//    GlobalsGL::createShaderProgram("vertex_shader.glsl", "fragment_shader.glsl");
    this->shader = Shader("vertex_shader.glsl", "fragment_shader.glsl");
    glEnable              ( GL_DEBUG_OUTPUT );
    GlobalsGL::f()->glDebugMessageCallback( GlobalsGL::MessageCallback, 0 );

    this->grid->createMesh();
    this->voxelGrid->createMesh();

    this->matter_adder = RockErosion(this->selectionSize, 1.0);
}

void Viewer::draw() {
    glClear(GL_DEPTH_BUFFER_BIT);
//    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    this->shader.use();

    float pMatrix[16];
    float mvMatrix[16];
    camera()->getProjectionMatrix(pMatrix);
    camera()->getModelViewMatrix(mvMatrix);

    this->shader.setMatrix("proj_matrix", pMatrix);
    this->shader.setMatrix("mv_matrix", mvMatrix);

    voxelGrid->display();
//    drawAxis();
    //grid->floatArrayMesh.size()/3);

/*

    if (this->viewerMode == ViewerMode::WIRE_MODE)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    // Place light at camera position
    const Vec cameraPos = camera()->position();
//    const GLfloat pos[4] = {(float)cameraPos[0], (float)cameraPos[1],
//                            (float)cameraPos[2] + 10.f, 1.0f};
    const GLfloat pos[4] = {0, 0.0, 50.0, 1.0};
    glLightfv(GL_LIGHT1, GL_POSITION, pos);

    // Orientate light along view direction
//    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, camera()->viewDirection());
    const GLfloat orient[4] = {0.1, 0.1, 1.0};
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, orient);


    // Light default parameters
    const GLfloat light_ambient[4] = {.92, .54, .20, .01};
    const GLfloat light_specular[4] = {.8, 1.0, .8, .001};
    const GLfloat light_diffuse[4] = {1.0, 1.0, 1.0, .01};

//    glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 30.0);
//    glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 100.0);
//    glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, .000001f);
//    glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, .003f);
//    glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.003f);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
//    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);

//    glPushMatrix();
//    glRotatef(180.0, 1.0, 0.0, 0.0);
    drawAxis();
//    glPopMatrix();
//    drawLight(GL_LIGHT1);

    if (displayRockTrajectories) {
        glPushMatrix();
        glTranslatef(-this->voxelGrid->getSizeX()/2.0, -this->voxelGrid->getSizeY()/2.0, -this->voxelGrid->getSizeZ()/2.0);
        for(std::vector<Vector3> coords : this->lastRocksLaunched) {
            glBegin(GL_LINE_STRIP);
            glColor4f(1.0, 1.0, 1.0, .5);
            for(Vector3 pos : coords) {
                glVertex3f(pos.x, pos.y, pos.z);
            }
            glEnd();
        }
        glPopMatrix();
    }

    if (this->mapMode == GRID_MODE) {
        this->grid->display(true);
    }
    else if (this->mapMode == VOXEL_MODE) {
        this->voxelGrid->display((this->algorithm & MARCHING_CUBES), this->display_vertices, 0.0);
    }
*/
    /*if(this->mouseInWorld)
    {
        glPushMatrix();
        glTranslatef(-this->voxelGrid->getSizeX()/2.0, -this->voxelGrid->getSizeY()/2.0, -this->voxelGrid->getSizeZ()/2.0);
        glTranslatef(mousePosWorld.x, mousePosWorld.y, mousePosWorld.z);
        glColor4f(1., 1., 1., .5);
        this->drawSphere(this->selectionSize, 8, 8);
        glColor4f(1., 1., 1., 1.);
        glPopMatrix();
    }*/
}


void Viewer::drawWithNames()
{
    this->draw();
}
void Viewer::postSelection(const QPoint &point)
{
    auto full_start = std::chrono::system_clock::now();
    camera()->convertClickToLine(point, orig, dir);

    // Find the selectedPoint coordinates, using camera()->pointUnderPixel().
    bool found;
    selectedPoint = camera()->pointUnderPixel(point, found);
    selectedPoint -= 0.01f * dir; // Small offset to make point clearly visible.
    // Note that "found" is different from (selectedObjectId()>=0) because of the
    // size of the select region.

    if (selectedName() == -1) {
        std::cout << "Nope." << std::endl;
    }
    else {
        intptr_t ptr_int = selectedName();
        if (ptr_int > 0) {
            Voxel* main_v = reinterpret_cast<Voxel*>(ptr_int);
            RockErosion rock(this->selectionSize, 1.0);
            rock.Apply(main_v, addingMatterMode);

        }
    }
    std::cout << "Total duration : " << (std::chrono::duration<double>(std::chrono::system_clock::now() - full_start)).count() << "s" << std::endl;
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
        update();
    } else if(e->key() == Qt::Key_V) {
        this->display_vertices = !this->display_vertices;
        update();
    } else if(e->key() == Qt::Key_P) {
        this->addingMatterMode = !this->addingMatterMode;
        std::cout << (addingMatterMode ? "Construction mode" : "Destruction mode") << std::endl;
        update();
    } else if(e->key() == Qt::Key_Return) {
        UnderwaterErosion erod(this->voxelGrid, 10, 0.05, 2000);
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
    /*qglviewer::Camera *c = camera();
    this->mousePos = e->pos();
    camera()->convertClickToLine(mousePos, orig, dir);
    float maxDist = max((int)camera()->distanceToSceneCenter(), max(voxelGrid->getSizeX(), max(voxelGrid->getSizeY(), voxelGrid->getSizeZ())));
    maxDist *= maxDist;

    bool found = false;
    Vector3 currPos(orig.x, orig.y, orig.z);
    currPos += Vector3(voxelGrid->getSizeX()/2.0, voxelGrid->getSizeY()/2.0, voxelGrid->getSizeZ()/2.0);
    while(currPos.norm() < maxDist)
    {
        currPos += Vector3(dir.x, dir.y, dir.z);
        Voxel* v = voxelGrid->getVoxel(currPos.x, currPos.y, currPos.z);
        if (v != nullptr && v->type != TerrainTypes::AIR) {
            found = true;
            break;
        }
    }
    this->mouseInWorld = found;
    if (found) {
        this->mousePosWorld = currPos; // + Vector3(voxelGrid->getSizeX()/2.0, voxelGrid->getSizeY()/2.0, voxelGrid->getSizeZ()/2.0);
//        std::cout << mousePosWorld << std::endl;
    }*/

/*
    bool found;
    this->mousePosWorld = camera()->pointUnderPixel(mousePos, found);
    this->mouseInWorld = found;
    if (found)
        std::cout << "Found" << std::endl;

    if (selectedName() == -1) {
        std::cout << "Nope." << std::endl;
    }
    else {
        std::cout << "Here" << std::endl;
    }

    float depth;
    // Qt uses upper corner for its origin while GL uses the lower corner.
    std::cout << camera()->screenHeight() <<" " << camera()->screenWidth() << std::endl;
    glReadPixels(mousePos.x(), camera()->screenHeight() - 1 - mousePos.y(), 1, 1,
                 GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    found = static_cast<double>(depth)< 1.0;
    Vec point(mousePos.x(), mousePos.y(), static_cast<double>(depth));
    point = camera()->unprojectedCoordinatesOf(point);
*/

    QGLViewer::mouseMoveEvent(e);
}
