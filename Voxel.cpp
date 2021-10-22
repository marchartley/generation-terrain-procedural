#include "Voxel.h"


Voxel::Voxel(int x, int y, int z, TerrainTypes type, float blockSize)
    : x(x), y(y), z(z), type(type), blockSize(blockSize) {

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
    if (this->type == TerrainTypes::AIR)
        return;
    glPushMatrix();
    glTranslatef(-this->x, -this->y, -this->z);
    glScalef(1/this->blockSize, 1/this->blockSize, 1/this->blockSize);

    // BOTTOM
    if (!this->has_neighbors[BOTTOM]) {
        glColor3f(0.3, 0.2, 1.0);
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
        glColor3f(1.0, 0.2, 0.3);
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
        glBegin(GL_QUADS);
        glVertex3f(1, 0, 1);
        glVertex3f(1, 1, 1);
        glVertex3f(0, 1, 1);
        glVertex3f(0, 0, 1);
        glEnd();
    }

    // LEFT
    if (!this->has_neighbors[LEFT]) {
        glBegin(GL_QUADS);
        glVertex3f(0, 0, 1);
        glVertex3f(0, 1, 1);
        glVertex3f(0, 1, 0);
        glVertex3f(0, 0, 0);
        glEnd();
    }

    // FRONT
    if (!this->has_neighbors[FRONT]) {
        glColor3f(0.3, 1.0, 0.2);
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
        glBegin(GL_QUADS);
        glVertex3f(0, 1, 0);
        glVertex3f(1, 1, 0);
        glVertex3f(1, 1, 1);
        glVertex3f(0, 1, 1);
        glEnd();
    }
    glPopMatrix();
}
