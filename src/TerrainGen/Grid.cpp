#include "Utils/Globals.h"
#include "TerrainGen/Grid.h"

#include "Utils/FastNoiseLit.h"
//#include <QtOpenGLExtensions/QtOpenGLExtensions>


Grid::Grid(int nx, int ny, float maxHeight, float tileSize)
    : sizeX(nx), sizeY(ny), maxHeight(maxHeight), tileSize(tileSize) {
    // Create and configure FastNoise object
    FastNoiseLite noise;
    noise.SetFrequency(0.01);
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(10);

    // Make a first pass, to get some noise
    float min = 1000, max = -1000;
    this->vertices = new Vector3*[this->sizeX];
    this->normals = new Vector3*[this->sizeX];
    for (int x = 0; x < this->sizeX; x++) {
        this->vertices[x] = new Vector3[this->sizeY];
        this->normals[x] = new Vector3[this->sizeY];
        for (int y = 0; y < this->sizeY; y++) {
            float z = noise.GetNoise((float)x, (float)y);
            this->vertices[x][y] = Vector3(x, y, z);
            min = min < z ? min : z;
            max = max > z ? max : z;
        }
    }
    // Sadely, we need a second pass because the noise function is a little f*cked up... we don't always have values in [-1;1]
    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            this->vertices[x][y].z -= min;
            this->vertices[x][y].z /= (max - min);
            this->vertices[x][y].z *= maxHeight;
//            this->vertices[x][y].z = maxHeight;
        }
    }
    this->computeNormals();
}
Grid::Grid() : Grid(10, 10, 5.0) {

}

void Grid::createMesh()
{
    std::vector<Vector3> vecs;

    for(int x = 0; x < this->sizeX - 1; x++) {
        for (int y = 0; y < this->sizeY - 1; y++) {
            vecs.push_back(Vector3(this->vertices[x][y].x, this->vertices[x][y].y, this->vertices[x][y].z));
            vecs.push_back(Vector3(this->vertices[x][y+1].x, this->vertices[x][y+1].y, this->vertices[x][y+1].z));
            vecs.push_back(Vector3(this->vertices[x+1][y+1].x, this->vertices[x+1][y+1].y, this->vertices[x+1][y+1].z));

            vecs.push_back(Vector3(this->vertices[x][y].x, this->vertices[x][y].y, this->vertices[x][y].z));
            vecs.push_back(Vector3(this->vertices[x+1][y].x, this->vertices[x+1][y].y, this->vertices[x+1][y].z));
            vecs.push_back(Vector3(this->vertices[x+1][y+1].x, this->vertices[x+1][y+1].y, this->vertices[x+1][y+1].z));
        }
    }
    this->mesh.fromArray(vecs);
//    return this->vertexArrayFloat;
}

void Grid::computeNormals() {
//    Vector3 upVector(0.0, 0.0, 1.0);
//    for (int x = 0; x < this->sizeX; x++)
//        for (int y = 0; y < this->sizeY; y++)
//            this->normals[x][y] = upVector;


    for (int x = 0; x < this->sizeX; x++) {
        for (int y = 0; y < this->sizeY; y++) {
            this->normals[x][y] = Vector3();

            Vector3 v1 = x > 0 ? this->vertices[x - 1][y    ] - this->vertices[x][y] : Vector3();
            Vector3 v2 = x > 0 && y > 0 ? this->vertices[x - 1][y - 1] - this->vertices[x][y] : Vector3();
            Vector3 v3 = y > 0 ? this->vertices[x    ][y - 1] - this->vertices[x][y] : Vector3();
            Vector3 v4 = x < sizeX - 1 && y > 0 ? this->vertices[x + 1][y - 1] - this->vertices[x][y] : Vector3();
            Vector3 v5 = x < sizeX - 1 ? this->vertices[x + 1][y    ] - this->vertices[x][y] : Vector3();
            Vector3 v6 = x < sizeX - 1 && y < sizeY - 1 ?this->vertices[x + 1][y + 1] - this->vertices[x][y] : Vector3();
            Vector3 v7 = y < sizeY - 1 ? this->vertices[x    ][y + 1] - this->vertices[x][y] : Vector3();
            Vector3 v8 = x > 0 && y < sizeY - 1 ? this->vertices[x - 1][y + 1] - this->vertices[x][y] : Vector3();

            this->normals[x][y] += v1.cross(v2);
            this->normals[x][y] += v2.cross(v3);
            this->normals[x][y] += v3.cross(v4);
            this->normals[x][y] += v4.cross(v5);
            this->normals[x][y] += v5.cross(v6);
            this->normals[x][y] += v6.cross(v7);
            this->normals[x][y] += v7.cross(v8);
            this->normals[x][y] += v8.cross(v1);
            this->normals[x][y].normalize();
        }
    }
}

void Grid::display(bool displayNormals) {
    this->mesh.display();
    /*
    glPushMatrix();

    float maxi = -9999999, mini = 9999999;
    for (int x = 0; x < this->sizeX; x++)
        for(int y = 0; y < this->sizeY; y++){
            if (this->vertices[x][y].z > maxi)
                maxi = this->vertices[x][y].z;
            if (this->vertices[x][y].z < mini)
                mini = this->vertices[x][y].z;
        }
    glScalef(1/this->tileSize, 1/this->tileSize, 1/this->tileSize);
    glTranslatef(-(sizeX - 1) / 2.0, - (sizeY - 1) / 2.0, - (maxi + mini) / 2.0);

    glBegin(GL_TRIANGLES);
//    glColor4f(1.0, 0.5, 0.5, 0.0);
    for(int x = 0; x < this->sizeX - 1; x++) {
        for (int y = 0; y < this->sizeY - 1; y++) {
            glVertex3f(this->vertices[x][y].x, this->vertices[x][y].y, this->vertices[x][y].z);
            glVertex3f(this->vertices[x][y+1].x, this->vertices[x][y+1].y, this->vertices[x][y+1].z);
            glVertex3f(this->vertices[x+1][y+1].x, this->vertices[x+1][y+1].y, this->vertices[x+1][y+1].z);

            glVertex3f(this->vertices[x][y].x, this->vertices[x][y].y, this->vertices[x][y].z);
            glVertex3f(this->vertices[x+1][y].x, this->vertices[x+1][y].y, this->vertices[x+1][y].z);
            glVertex3f(this->vertices[x+1][y+1].x, this->vertices[x+1][y+1].y, this->vertices[x+1][y+1].z);
        }
    }
    glEnd();

    if (displayNormals) {
        this->computeNormals();
        for(int x = 0; x < this->sizeX - 1; x++) {
            for (int y = 0; y < this->sizeY - 1; y++) {
                glBegin(GL_LINES);
                glVertex3f(this->vertices[x][y].x, this->vertices[x][y].y, this->vertices[x][y].z);
                glVertex3f(this->vertices[x][y].x + this->normals[x][y].x, this->vertices[x][y].y + this->normals[x][y].y, this->vertices[x][y].z + this->normals[x][y].z);
                glEnd();
            }
        }
    }
    glPopMatrix();*/
}

void Grid::fromVoxelGrid(VoxelGrid &voxelGrid) {
    for (int x = 0; x < voxelGrid.getSizeX(); x++) {
        for (int y = 0; y < voxelGrid.getSizeY(); y++) {
            this->vertices[x][y] = Vector3(x, y, voxelGrid.getHeight(x, y)/ (voxelGrid.getSizeZ() / 2.0));
        }
    }
    this->computeNormals();
}
