#ifndef GRID_H
#define GRID_H

#include "DataStructure/Vector3.h"

class Grid;

#include "TerrainGen/VoxelGrid.h"
#include "Graphics/Mesh.h"
#include "DataStructure/Matrix3.h"

class Grid {
public:
    Grid();
    Grid(int nx, int ny, float max_height, float tileSize = 0.1);
    Grid(std::string heightmap_filename, int nx = -1, int ny = -1, float max_height = -1, float tileSize = 0.1);

    void display(bool displayNormals = false);

    int getSizeX() {return vertices.sizeX;}
    int getSizeY() {return vertices.sizeY;}

    float getHeight(int x, int y) { return this->vertices.at(x, y).z; }
    float getMaxHeight();
    float getTileSize() { return this->tileSize; }
    Vector3 getNormal(int x, int y) { return this->normals.at(x, y); }

    void createMesh();

    void fromVoxelGrid(VoxelGrid& voxelGrid);

    void loadFromHeightmap(std::string heightmap_filename, int nx = -1, int ny = -1, float max_height = -1, float tileSize = 0.1);
    void saveHeightmap(std::string heightmap_filename);

//protected:
    void computeNormals();
    Matrix3<Vector3> vertices;
    Matrix3<Vector3> normals;
    float maxHeight;
    float tileSize;
    std::vector<float> vertexArrayFloat;
    Mesh mesh;
};

#endif // GRID_H
