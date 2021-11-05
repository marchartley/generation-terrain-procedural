#include "Viewer.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <chrono>
#include "UnderwaterErosion.h"

using namespace qglviewer;
using namespace std;

Viewer::Viewer(Grid* grid, VoxelGrid* voxelGrid, MapMode map, ViewerMode mode)
    : QGLViewer(), mapMode(map), viewerMode(mode), grid(grid), voxelGrid(voxelGrid) {

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

    this->matter_adder = RockErosion(this->selectionSize, 1.0);
    this->selectionShape = this->matter_adder.attackMask;
//    std::cout.precision(2);
//    selectionShape = new float**[selectionSize];
//    for (int i = 0; i < selectionSize; i++) {
//        selectionShape[i] = new float*[selectionSize];
//        for (int j = 0; j < selectionSize; j++) {
//            selectionShape[i][j] = new float[selectionSize];
//            for (int k = 0; k < selectionSize; k++) {
//                float i_i = (i - (selectionSize-1)/2.0) / ((selectionSize-1)/2.0);
//                float j_i = (j - (selectionSize-1)/2.0) / ((selectionSize-1)/2.0);
//                float k_i = (k - (selectionSize-1)/2.0) / ((selectionSize-1)/2.0);
//                selectionShape[i][j][k] = (sqrt(3) - sqrt(i_i*i_i + j_i*j_i + k_i*k_i))/(sqrt(3));
////                std::cout << selectionShape[i][j][k] << "\t" << std::flush;
////                continue;
//            }
////            std::cout << std::endl;
//        }
////        std::cout << std::endl;
//    }

}

void Viewer::draw() {
    if (this->viewerMode == ViewerMode::WIRE_MODE)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    // Place light at camera position
    const Vec cameraPos = camera()->position();
    const GLfloat pos[4] = {(float)cameraPos[0], (float)cameraPos[1],
                            (float)cameraPos[2], 1.0f};
//    const GLfloat pos[4] = {0, -5.0, -5.0, 1.0};
    glLightfv(GL_LIGHT1, GL_POSITION, pos);

    // Orientate light along view direction
//    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, camera()->viewDirection());
    const GLfloat orient[4] = {0.1, 0.1, 1.0};
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, orient);


    // Light default parameters
    const GLfloat light_ambient[4] = {1.0, .9, .9, 1.0};
    const GLfloat light_specular[4] = {1.0, 1.0, 1.0, 1.0};
    const GLfloat light_diffuse[4] = {1.0, 1.0, 1.0, 100.0};

    glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 30.0);
    glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 100.0);
//    glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, .1f);
    glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, .00003f);
    glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.0003f);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
//    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
//    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);

//    glPushMatrix();
//    glRotatef(180.0, 1.0, 0.0, 0.0);
    drawAxis();
//    glPopMatrix();
//    drawLight(GL_LIGHT1);

    if (this->mapMode == GRID_MODE) {
        this->grid->display(true);
    }
    else if (this->mapMode == VOXEL_MODE) {
        this->voxelGrid->display((this->algorithm & MARCHING_CUBES), this->display_vertices);
    }
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
            auto total_time = std::chrono::duration<double>();
            Voxel* main_v = reinterpret_cast<Voxel*>(ptr_int);
            this->matter_adder.Apply(main_v, addingMatterMode);
            /*
            for (int x = -(int(selectionSize/2)); x <= int(selectionSize/2); x++) {
                for (int y = -(int(selectionSize/2)); y <= int(selectionSize/2); y++) {
                    for (int z = -(int(selectionSize/2)); z <= int(selectionSize/2); z++){
                        auto start = std::chrono::system_clock::now();
                        int v_x = main_v->x + x, v_y = main_v->y + y, v_z = main_v->z + z;
                        VoxelChunk* current_parent = main_v->parent;
                        if ((main_v->x + x < 0 && current_parent->x == 0) || (main_v->x + x >= current_parent->sizeX && current_parent->lastChunkOnX))
                            continue;
                        if ((main_v->y + y < 0 && current_parent->y == 0) || (main_v->y + y >= current_parent->sizeY && current_parent->lastChunkOnY))
                            continue;
                        if (main_v->z + z < 0 || main_v->z + z >= current_parent->height)
                            continue;

                        if (main_v->x + x < 0) {
                            current_parent = current_parent->neighboring_chunks[LEFT];
                            v_x += current_parent->sizeX;
                        }
                        if (main_v->y + y < 0) {
                            current_parent = current_parent->neighboring_chunks[FRONT];
                            v_y += current_parent->sizeY;
                        }
                        if (main_v->x + x >= current_parent->sizeX) {
                            current_parent = current_parent->neighboring_chunks[RIGHT];
                            v_x -= current_parent->sizeX;
                        }
                        if (main_v->y + y >= current_parent->sizeY) {
                            current_parent = current_parent->neighboring_chunks[BACK];
                            v_y -= current_parent->sizeY;
                        }
                        Voxel* v = current_parent->voxels[v_x][v_y][v_z];
//                        v->type = TerrainTypes::AIR;
                        v->manual_isosurface += this->selectionShape[x+int(selectionSize/2)][y+int(selectionSize/2)][z+int(selectionSize/2)] * (this->addingMatterMode ? 1 : -1);
                        if (v->getIsosurface() < 0.0)
                            v->type = TerrainTypes::AIR;
                        else
                            v->type = TerrainTypes::DIRT;
                        v->parent->data[v->x][v->y][v->z] = v->type;
                        auto end = std::chrono::system_clock::now();
//                        std::cout << (end - start).count() << std::endl;
                        total_time += (end - start);
                    }
                }
            }
            std::cout << "Time in loop : " << total_time.count() << std::endl;
            auto start = std::chrono::system_clock::now();
            main_v->parent->createMesh();
            auto end = std::chrono::system_clock::now();
            total_time = (end - start);
            std::cout << "Remeshing : " << total_time.count() << std::endl;*/

        }
    }
    std::cout << "Total duration : " << (std::chrono::duration<double>(std::chrono::system_clock::now() - full_start)).count() << std::endl;
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
        update();
    } else if(e->key() == Qt::Key_V) {
        this->display_vertices = !this->display_vertices;
        update();
    } else if(e->key() == Qt::Key_P) {
        this->addingMatterMode = !this->addingMatterMode;
        update();
    } else if(e->key() == Qt::Key_Return) {
        UnderwaterErosion erod(this->voxelGrid, 10, 2.0, 100);
        erod.Apply();
        update();
    } else {
        QGLViewer::keyPressEvent(e);
    }
}
