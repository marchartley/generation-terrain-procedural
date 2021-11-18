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
    Mesh fromArray(std::vector<Vector3> vertices);

    void update();
    void pushToBuffer();
    void display();

    unsigned int bufferID;
    bool bufferReady = false;


//protected:
    void computeNormals();
    std::vector<Vector3> vertexArray;
    std::vector<float> vertexArrayFloat;
    std::vector<Vector3> normalsArray;
    std::vector<float> normalsArrayFloat;
    std::vector<float> colorArrayFloat;
};

#endif // MESH_H
