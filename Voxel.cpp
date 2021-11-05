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
    this->isosurface = 0;
    for (int i = 0; i < 8; i++) {
        this->isosurfaces[i] = &this->vertices[i].isosurface;
        this->isosurface += *this->isosurfaces[i];
    }
    this->isosurface /= 8;
}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0) {

}

float Voxel::getIsosurface() {
    this->isosurface = this->type == TerrainTypes::AIR ? -1.0 : 1.0;
//    this->isosurface = 0;
//    for (int i = 0; i < 8; i++) {
//        this->isosurface += *this->isosurfaces[i];
//    }
//    this->isosurface /= 8;
    return this->isosurface + manual_isosurface;
//    return this->isosurface * this->isosurface * (this->isosurface > 0 ? 1 : -1) + this->manual_isosurface;
}
void Voxel::display(bool apply_marching_cubes, bool display_vertices) {
    glPushMatrix();
    glTranslatef(this->x, this->y, this->z);

    if (apply_marching_cubes) {
        glPushName(reinterpret_cast<intptr_t>(this));
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
            glPushName(reinterpret_cast<intptr_t>(this));
            glBegin(GL_QUADS);
            // BOTTOM
            if (!this->has_neighbors[BOTTOM]) {
                glColor3f(1.0, 1.0, 1.0);
                glVertex3f(0, 0, 0);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(0, 1, 0);
                glColor3f(1.0, 1.0, 1.0);
            }

            // RIGHT
            if (!this->has_neighbors[RIGHT]) {
                glColor3f(0.8, 0.8, 0.8);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(1, 1, 1);
                glVertex3f(1, 0, 1);
                glColor3f(1.0, 1.0, 1.0);
            }

            // TOP
            if (!this->has_neighbors[TOP]) {
                glColor3f(0.3, 0.3, 0.3);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 1, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 0, 1);
            }

            // LEFT
            if (!this->has_neighbors[LEFT]) {
                glColor3f(0.6, 0.6, 0.6);
                glVertex3f(0, 0, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 1, 0);
                glVertex3f(0, 0, 0);
            }

            // FRONT
            if (!this->has_neighbors[FRONT]) {
                glColor3f(0.6, 0.6, 0.6);
                glVertex3f(0, 0, 0);
                glVertex3f(0, 0, 1);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 0, 0);
                glColor3f(1.0, 1.0, 1.0);
            }

            // BACK
            if (!this->has_neighbors[BACK]) {
                glColor3f(0.8, 0.8, 0.8);
                glVertex3f(0, 1, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(1, 1, 1);
                glVertex3f(0, 1, 1);
            }
            glEnd();
        }
    }
    glPopName();
    glPopMatrix();
}


std::ostream& operator<<(std::ostream& io, const Voxel& v)
{
    io << "Voxel (" << v.x << ", " << v.y << ", " << v.z << ")";
    return io;
}
std::ostream& operator<<(std::ostream& io, Voxel* v)
{
    io << "Voxel (" << v->x << ", " << v->y << ", " << v->z << ")";
    return io;
}

bool Voxel::contains(Vector3 v) {
    return this->contains(v.x, v.y, v.z);
}

bool Voxel::contains(float x, float y, float z) {
    return (this->globalX() <= x && x < this->globalX() + 1 && this->globalY() <= y && y < this->globalY() + 1 && this->globalZ() <= z && z < this->globalZ() + 1);
}

float Voxel::globalX()  {
    return this->x + this->parent->x;
}
float Voxel::globalY()  {
    return this->y + this->parent->y;
}
float Voxel::globalZ()  {
    return this->z;
}
