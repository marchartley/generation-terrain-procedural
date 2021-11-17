#include "Voxel.h"

#include "MarchingCubes.h"

std::vector<std::set<int>> Voxel::voxelGroups;
int Voxel::currentLabelIndex = 0;

Voxel::Voxel(int x, int y, int z, TerrainTypes type, float blockSize, float isosurface)
    : x(x), y(y), z(z), type(type), blockSize(blockSize), isosurface(isosurface) {
    this->has_neighbors[TOP] = false;
    this->has_neighbors[BOTTOM] = false;
    this->has_neighbors[FRONT] = false;
    this->has_neighbors[BACK] = false;
    this->has_neighbors[RIGHT] = false;
    this->has_neighbors[LEFT] = false;

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
//    this->vertices[2] = Vestd::vector<rtex(1.0, 0.0, 0.0);
//    this->vertices[3] = Vertex(0.0, 0.0, 0.0);
//    this->vertices[4] = Vertex(0.0, 1.0, 1.0);
//    this->vertices[5] = Vertex(1.0, 1.0, 1.0);
//    this->vertices[6] = Vertex(1.0, 0.0, 1.0);
//    this->vertices[7] = Vertex(0.0, 0.0, 1.0);
//    this->isosurface = 0;
//    for (int i = 0; i < 8; i++) {
//        this->isosurfaces[i] = &this->vertices[i].isosurface;
//        this->isosurface += *this->isosurfaces[i];
//    }
//    this->isosurface /= 8;
}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0, 0.0) {

}
void Voxel::addNeighbor(Voxel* neighbor)
{
    if(neighbor == nullptr || neighbor->type == TerrainTypes::AIR || this->type == TerrainTypes::AIR)
        return;
    if (this->globalX() < neighbor->globalX()) {
        this->neighbors[RIGHT] = neighbor;
        this->has_neighbors[RIGHT] = true;
        neighbor->neighbors[LEFT] = this;
        neighbor->has_neighbors[LEFT] = true;
    }
    if (this->globalX() > neighbor->globalX()) {
        this->neighbors[LEFT] = neighbor;
        this->has_neighbors[LEFT] = true;
        neighbor->neighbors[RIGHT] = this;
        neighbor->has_neighbors[RIGHT] = true;
    }
    if (this->globalY() < neighbor->globalY()) {
        this->neighbors[BACK] = neighbor;
        this->has_neighbors[BACK] = true;
        neighbor->neighbors[FRONT] = this;
        neighbor->has_neighbors[FRONT] = true;
    }
    if (this->globalY() > neighbor->globalY()) {
        this->neighbors[FRONT] = neighbor;
        this->has_neighbors[FRONT] = true;
        neighbor->neighbors[BACK] = this;
        neighbor->has_neighbors[BACK] = true;
    }
    if (this->globalZ() < neighbor->globalZ()) {
        this->neighbors[TOP] = neighbor;
        this->has_neighbors[TOP] = true;
        neighbor->neighbors[BOTTOM] = this;
        neighbor->has_neighbors[BOTTOM] = true;
    }
    if (this->globalZ() > neighbor->globalZ()) {
        this->neighbors[BOTTOM] = neighbor;
        this->has_neighbors[BOTTOM] = true;
        neighbor->neighbors[TOP] = this;
        neighbor->has_neighbors[TOP] = true;
    }
    if (this->group == -1)
        this->group = neighbor->group;
    this->group = std::min(this->group, neighbor->group);
//    if(this->group == 0)
//        this->group = neighbor->group;
//    if(this->group != neighbor->group)
//    {
//        Voxel::voxelGroups[std::max(this->group, neighbor->group)] = std::min(this->group, neighbor->group);
//    }
}
void Voxel::removeNeighbor(Voxel* neighbor)
{
    if (neighbor == nullptr)
        return;
    if (this->globalX() < neighbor->globalX()) {
        neighbor->neighbors[RIGHT] = nullptr;
        neighbor->has_neighbors[RIGHT] = false;
        this->neighbors[LEFT] = nullptr;
        this->has_neighbors[LEFT] = false;
    }
    if (this->globalX() > neighbor->globalX()) {
        neighbor->neighbors[LEFT] = nullptr;
        neighbor->has_neighbors[LEFT] = false;
        this->neighbors[RIGHT] = nullptr;
        this->has_neighbors[RIGHT] = false;
    }
    if (this->globalY() < neighbor->globalY()) {
        neighbor->neighbors[FRONT] = nullptr;
        neighbor->has_neighbors[FRONT] = false;
        this->neighbors[BACK] = nullptr;
        this->has_neighbors[BACK] = false;
    }
    if (this->globalY() > neighbor->globalY()) {
        neighbor->neighbors[BACK] = nullptr;
        neighbor->has_neighbors[BACK] = false;
        this->neighbors[FRONT] = nullptr;
        this->has_neighbors[FRONT] = false;
    }
    if (this->globalZ() < neighbor->globalZ()) {
        neighbor->neighbors[BOTTOM] = nullptr;
        neighbor->has_neighbors[BOTTOM] = false;
        this->neighbors[TOP] = nullptr;
        this->has_neighbors[TOP] = false;
    }
    if (this->globalZ() > neighbor->globalZ()) {
        neighbor->neighbors[TOP] = nullptr;
        neighbor->has_neighbors[TOP] = false;
        this->neighbors[BOTTOM] = nullptr;
        this->has_neighbors[BOTTOM] = false;
    }
}
void Voxel::resetNeighbors() {
    for(std::map<VOXEL_NEIGHBOR, Voxel*>::iterator it = this->neighbors.begin(); it != this->neighbors.end(); it++)
        removeNeighbor(it->second);
}

float Voxel::getIsosurface() {
    return this->isosurface + manual_isosurface;
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
                glColor3f(.9, .9, .9);
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
                glColor3f(.2, .5, .2);
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


std::vector<Vector3> Voxel::getMeshVertices()
{
    std::vector<Vector3> vertex;
    if (this->type != TerrainTypes::AIR) {
        // BOTTOM
        if (!this->has_neighbors[BOTTOM]) {
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 0 + globalZ()));
        }

        // RIGHT
        if (!this->has_neighbors[RIGHT]) {
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 1 + globalZ()));
        }

        // TOP
        if (!this->has_neighbors[TOP]) {
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 1 + globalZ()));
        }

        // LEFT
        if (!this->has_neighbors[LEFT]) {
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 0 + globalZ()));
        }

        // FRONT
        if (!this->has_neighbors[FRONT]) {
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 0 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 0 + globalY(), 0 + globalZ()));
        }

        // BACK
        if (!this->has_neighbors[BACK]) {
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 0 + globalZ()));
            vertex.push_back(Vector3(1 + globalX(), 1 + globalY(), 1 + globalZ()));
            vertex.push_back(Vector3(0 + globalX(), 1 + globalY(), 1 + globalZ()));
        }
    }
    std::vector<Vector3> returnArray;
    for (size_t i = 0; i < vertex.size(); i+=4)
    {
        returnArray.push_back(vertex[i + 0]);
        returnArray.push_back(vertex[i + 1]);
        returnArray.push_back(vertex[i + 2]);
        returnArray.push_back(vertex[i + 2]);
        returnArray.push_back(vertex[i + 3]);
        returnArray.push_back(vertex[i + 0]);
    }
    return returnArray;
}
