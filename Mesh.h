#ifndef MESH_H
#define MESH_H

#include "Globals.h"
#include "Vector3.h"
#include "Shader.h"
#include <vector>
#include <map>

class Mesh
{
public:
    Mesh();
    Mesh(std::vector<Vector3> _vertexArray);
    Mesh(std::vector<float> _vertexArrayFloat);
    Mesh fromArray(std::vector<Vector3> vertices, std::vector<int> indices = std::vector<int>());
    Mesh fromArray(std::vector<float> vertices, std::vector<int> indices = std::vector<int>());

    void clear();
    Mesh merge(std::shared_ptr<Mesh> other, bool recomputeIndices = true);
    Mesh merge(std::vector<std::shared_ptr<Mesh>> others);

    void update();
    void pushToBuffer();
    void display(GLenum shape = GL_TRIANGLES);


    unsigned int bufferID;
    bool bufferReady = false;

    bool useIndices = true;

//protected:
    void computeIndices(std::vector<int> indices = std::vector<int>());
    void computeNormals();
    void computeColors();
    std::vector<Vector3> vertexArray;
    std::vector<float> vertexArrayFloat;
    std::vector<Vector3> normalsArray;
    std::vector<Vector3> normalsArrayIndex;
    std::vector<float> normalsArrayFloat;
    std::vector<Vector3> colorsArray;
    std::vector<Vector3> colorsArrayIndex;
    std::vector<float> colorsArrayFloat;
    std::map<int, int> indices;
    std::map<std::tuple<int, int, int>, int> vectorToIndex;
    std::shared_ptr<Shader> shader = nullptr;

};

#endif // MESH_H
