#include "Voxel.h"

#include "MarchingCubes.h"


Voxel::Voxel(int x, int y, int z, TerrainTypes type, float blockSize)
    : x(x), y(y), z(z), type(type), blockSize(blockSize) {
    this->has_neighbors[TOP] = true;
    this->has_neighbors[BOTTOM] = true;
    this->has_neighbors[FRONT] = true;
    this->has_neighbors[BACK] = true;
    this->has_neighbors[RIGHT] = true;
    this->has_neighbors[LEFT] = true;

    this->vertices[0] = Vertex(0.1, 0.9, 0.1);
    this->vertices[1] = Vertex(0.9, 0.9, 0.1);
    this->vertices[2] = Vertex(0.9, 0.1, 0.1);
    this->vertices[3] = Vertex(0.1, 0.1, 0.1);
    this->vertices[4] = Vertex(0.1, 0.9, 0.9);
    this->vertices[5] = Vertex(0.9, 0.9, 0.9);
    this->vertices[6] = Vertex(0.9, 0.1, 0.9);
    this->vertices[7] = Vertex(0.1, 0.1, 0.9);

//    this->vertices[0] = Vertex(0.0, 1.0, 0.0);
//    this->vertices[1] = Vertex(1.0, 1.0, 0.0);
//    this->vertices[2] = Vertex(1.0, 0.0, 0.0);
//    this->vertices[3] = Vertex(0.0, 0.0, 0.0);
//    this->vertices[4] = Vertex(0.0, 1.0, 1.0);
//    this->vertices[5] = Vertex(1.0, 1.0, 1.0);
//    this->vertices[6] = Vertex(1.0, 0.0, 1.0);
//    this->vertices[7] = Vertex(0.0, 0.0, 1.0);

    for (int i = 0; i < 8; i++)
        this->isosurfaces[i] = &this->vertices[i].isosurface;
}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0) {

}

void Voxel::display(bool apply_marching_cubes, bool display_vertices) {
    glPushMatrix();
    glTranslatef(this->x, this->y, this->z);

    if (apply_marching_cubes) {
        int cube_index = 0;
        for (int i = 0; i < 8; i++){
            if (this->vertices[i].isosurface < 0)
                cube_index ^= 1 << i;
            if (display_vertices)
                this->vertices[i].display();
        }
        int* edgesForTriangles = MarchingCubes::triangleTable[cube_index];
        glBegin(GL_TRIANGLES);
        Vertex originalVertex;
        Vertex firstVertex;
        Vertex secondVertex;
        for (int i = 0; i < 16; i++) {
            if (edgesForTriangles[i] == -1)
                continue;
            Vertex& v1 = this->vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][0]];
            Vertex& v2 = this->vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][1]];
            Vertex midpoint = (v1 + v2);
            midpoint /= 2.0;
            if (i % 3 == 0)
                originalVertex = midpoint;
            else if (i % 3 == 1)
                firstVertex = midpoint;
            else {
                secondVertex = midpoint;
                Vector3 normal = (firstVertex - originalVertex).cross((secondVertex - originalVertex)).normalize();
                glColor3f(abs(normal.x), abs(normal.y), abs(normal.z));/*
                std::cout << normal << " = " << firstVertex << " - " << originalVertex << std::endl;
                std::cout << normal.x << ", " << 1-normal.x << ", 0.0" << std::endl;*/
            }
            glVertex3f(midpoint.x, midpoint.y, midpoint.z);
        }
        glEnd();
    } else {
        if (this->type != TerrainTypes::AIR) {
            // BOTTOM
            if (!this->has_neighbors[BOTTOM]) {
                glColor3f(1.0, 1.0, 1.0);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 0);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(0, 1, 0);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // RIGHT
            if (!this->has_neighbors[RIGHT]) {
                glColor3f(0.8, 0.8, 0.8);
                glBegin(GL_QUADS);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(1, 1, 1);
                glVertex3f(1, 0, 1);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // TOP
            if (!this->has_neighbors[TOP]) {
                glColor3f(0.3, 0.3, 0.3);
                glBegin(GL_QUADS);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 1, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 0, 1);
                glEnd();
            }

            // LEFT
            if (!this->has_neighbors[LEFT]) {
                glColor3f(0.6, 0.6, 0.6);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 1, 0);
                glVertex3f(0, 0, 0);
                glEnd();
            }

            // FRONT
            if (!this->has_neighbors[FRONT]) {
                glColor3f(0.6, 0.6, 0.6);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 0);
                glVertex3f(0, 0, 1);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 0, 0);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // BACK
            if (!this->has_neighbors[BACK]) {
                glColor3f(0.8, 0.8, 0.8);
                glBegin(GL_QUADS);
                glVertex3f(0, 1, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(1, 1, 1);
                glVertex3f(0, 1, 1);
                glEnd();
            }
        }
    }
    glPopMatrix();
}
