#ifndef GRID_H
#define GRID_H

#include "DataStructure/Vector3.h"

class Grid;

#include "TerrainGen/VoxelGrid.h"
#include "Graphics/Mesh.h"

class Grid {
public:
    Grid();
    Grid(int nx, int ny, float max_height, float tileSize = 0.1);

    void display(bool displayNormals = false);

    int getSizeX() {return this->sizeX;}
    int getSizeY() {return this->sizeY;}

    float getHeight(int x, int y) { return this->vertices[x][y].z; }
    float getMaxHeight() { return maxHeight; }
    float getTileSize() { return this->tileSize; }
    Vector3 getNormal(int x, int y) { return this->normals[x][y]; }

    void createMesh();

    void fromVoxelGrid(VoxelGrid& voxelGrid);

//protected:
    void computeNormals();
    Vector3** vertices;
    Vector3** normals;
    int sizeX, sizeY;
    float maxHeight;
    float tileSize;
    std::vector<float> vertexArrayFloat;
    Mesh mesh;
};

#endif // GRID_H
