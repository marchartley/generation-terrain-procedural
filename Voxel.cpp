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

    this->vertices[0] = Vertex(x, y+1, z);
    this->vertices[1] = Vertex(x+1, y+1, z);
    this->vertices[2] = Vertex(x+1, y, z);
    this->vertices[3] = Vertex(x, y, z);
    this->vertices[4] = Vertex(x, y+1, z+1);
    this->vertices[5] = Vertex(x+1, y+1, z+1);
    this->vertices[6] = Vertex(x+1, y, z+1);
    this->vertices[7] = Vertex(x, y, z+1);
}

Voxel::Voxel() : Voxel(0, 0, 0, TerrainTypes::AIR, 1.0) {

}
void Voxel::addNeighbor(Voxel &neighbor) {
    if (abs(neighbor.x - this->x) == 1) {
        if (neighbor.x > this->x) {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(LEFT, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(RIGHT, *this));
        } else {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(RIGHT, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(LEFT, *this));
        }
    }
    if (abs(neighbor.y - this->y) == 1) {
        if (neighbor.y > this->y) {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(FRONT, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(BACK, *this));
        } else {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(BACK, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(FRONT, *this));
        }
    }
    if (abs(neighbor.z - this->z) == 1) {
        if (neighbor.z > this->z) {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(TOP, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(BOTTOM, *this));
        } else {
            this->neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(BOTTOM, neighbor));
            neighbor.neighbors.insert(std::pair<VOXEL_NEIGHBOR, Voxel&>(TOP, *this));
        }
    }
}
void Voxel::removeNeighbor(Voxel &neighbor) {
    if (abs(neighbor.x - this->x) == 1) {
        if (neighbor.x > this->x) {
            this->neighbors.erase(LEFT);
            neighbor.neighbors.erase(RIGHT);
        } else {
            this->neighbors.erase(RIGHT);
            neighbor.neighbors.erase(LEFT);
        }
    }
    if (abs(neighbor.y - this->y) == 1) {
        if (neighbor.y > this->y) {
            this->neighbors.erase(FRONT);
            neighbor.neighbors.erase(BACK);
        } else {
            this->neighbors.erase(BACK);
            neighbor.neighbors.erase(FRONT);
        }
    }
    if (abs(neighbor.z - this->z) == 1) {
        if (neighbor.z > this->z) {
            this->neighbors.erase(TOP);
            neighbor.neighbors.erase(BOTTOM);
        } else {
            this->neighbors.erase(BOTTOM);
            neighbor.neighbors.erase(TOP);
        }
    }
}


void Voxel::display() {
    bool marching_cubes = false;
    glPushMatrix();
    glTranslatef(-this->x, -this->y, -this->z);

    if (marching_cubes) {
        int cube_index = 0;
        for (int i = 0; i < 8; i++){
            if (this->vertices[i].isosurface > 0)
                cube_index ^= 1 << i;
//            glPushMatrix();
//            glRotatef(180, 0.0, 0.0, 1.0);
//            glTranslatef(this->x, this->y, this->z);
            this->vertices[i].display();
//            glPopMatrix();
        }
//        std::cout << cube_index << std::endl;
//        if(rand() % 200 == 0)
//            exit(0);
        int* edgesForTriangles = MarchingCubes::triangleTable[cube_index];
        glBegin(GL_TRIANGLES);
        glPushMatrix();
//        glRotatef(90.0, 0.0, 0.0, 1.0);
//        glRotatef(90.0, 0.0, 1.0, 0.0);
//        glRotatef(90.0, 1.0, 0.0, 0.0);
//        glRotatef(-90.0, 0.0, 0.0, 1.0);
//        glRotatef(-90.0, 0.0, 1.0, 0.0);
//        glRotatef(-90.0, 1.0, 0.0, 0.0);
//        glRotatef(180.0, 0.0, 0.0, 1.0);
//        glRotatef(180.0, 0.0, 1.0, 0.0);
//        glRotatef(180.0, 1.0, 0.0, 0.0);
        for (int i = 0; i < 16; i++) {
            if (edgesForTriangles[i] == -1)
                continue;
            Vertex& v1 = this->vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][0]];
            Vertex& v2 = this->vertices[MarchingCubes::edgeToCorner[edgesForTriangles[i]][1]];
            Vertex midpoint = (v1 + v2);
            midpoint /= 2.0;
            midpoint -= Vector3(this->getX(), this->getY(), this->getZ());
            glVertex3f(midpoint.x, midpoint.y, midpoint.z);
        }
        glEnd();
        glPopMatrix();
    } else {
        if (this->type != TerrainTypes::AIR) {
            // TOP
            if (!this->has_neighbors[TOP]) {
                glColor3f(1.0, 1.0, 1.0);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 0);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(0, 1, 0);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // LEFT
            if (!this->has_neighbors[LEFT]) {
                glColor3f(0.8, 0.8, 0.8);
                glBegin(GL_QUADS);
                glVertex3f(1, 0, 0);
                glVertex3f(1, 1, 0);
                glVertex3f(1, 1, 1);
                glVertex3f(1, 0, 1);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // BOTTOM
            if (!this->has_neighbors[BOTTOM]) {
                glColor3f(0.3, 0.3, 0.3);
                glBegin(GL_QUADS);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 1, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 0, 1);
                glEnd();
            }

            // RIGHT
            if (!this->has_neighbors[RIGHT]) {
                glColor3f(0.6, 0.6, 0.6);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 1);
                glVertex3f(0, 1, 1);
                glVertex3f(0, 1, 0);
                glVertex3f(0, 0, 0);
                glEnd();
            }

            // BACK
            if (!this->has_neighbors[BACK]) {
                glColor3f(0.6, 0.6, 0.6);
                glBegin(GL_QUADS);
                glVertex3f(0, 0, 0);
                glVertex3f(0, 0, 1);
                glVertex3f(1, 0, 1);
                glVertex3f(1, 0, 0);
                glEnd();
                glColor3f(1.0, 1.0, 1.0);
            }

            // FRONT
            if (!this->has_neighbors[FRONT]) {
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
