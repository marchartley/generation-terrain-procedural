#ifndef VERTEX_H
#define VERTEX_H

#include "DataStructure/Vector3.h"

class Vertex : public Vector3
{
public:
    Vertex();
    Vertex(float x, float y, float z, float isosurface = -1.0);
    Vertex(const Vector3& v, float isosurface = -1.0);

    void display();

    float isosurface = -1.0;
};

#endif // VERTEX_H
