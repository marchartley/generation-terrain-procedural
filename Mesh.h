#ifndef MESH_H
#define MESH_H

#include "Globals.h"
#include "Vector3.h"
#include <vector>

class Mesh
{
public:
    Mesh();
    Mesh(std::vector<Vector3> _vertexArray);
    Mesh(std::vector<float> _vertexArrayFloat);

    void update();
    void pushToBuffer();
    void display();

    unsigned int bufferID;
    bool bufferReady = false;

//protected:
    std::vector<Vector3> vertexArray;
    std::vector<float> vertexArrayFloat;
};

#endif // MESH_H
