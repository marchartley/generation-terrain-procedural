#ifndef GRID_H
#define GRID_H

#include <QGLViewer/qglviewer.h>
#include "vector3.h"

class Grid;

#include "VoxelGrid.h"

class Grid {
public:
    Grid();
    Grid(int nx, int ny, float tileSize = 0.1);

    void display(bool displayAsWires = false, bool displayNormals = false);

    int getSizeX() {return this->sizeX;}
    int getSizeY() {return this->sizeY;}

    float getHeight(int x, int y) { return this->vertices[x][y].z; }
    Vector3 getNormal(int x, int y) { return this->normals[x][y]; }

    void fromVoxelGrid(VoxelGrid& voxelGrid);

protected:
    Vector3** vertices;
    Vector3** normals;
    int sizeX, sizeY;
    float tileSize;

    void computeNormals();
};

#endif // GRID_H
