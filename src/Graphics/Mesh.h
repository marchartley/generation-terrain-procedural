#ifndef MESH_H
#define MESH_H

#include "Utils/Globals.h"
#include "DataStructure/Vector3.h"
#include "Graphics/Shader.h"
#include <vector>
#include <map>

class Mesh
{
public:
    Mesh(std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<Vector3> _vertexArray, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<float> _vertexArrayFloat, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh fromArray(std::vector<Vector3> vertices, std::vector<int> indices = std::vector<int>());
    Mesh fromArray(std::vector<float> vertices, std::vector<int> indices = std::vector<int>());

    void clear();
    Mesh merge(std::shared_ptr<Mesh> other, bool recomputeIndices = true);
    Mesh merge(std::vector<std::shared_ptr<Mesh>> others);

    void update();
    void pushToBuffer();
    void display(GLenum shape = 0);


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
    bool isDisplayed;
    GLenum displayShape;

};

#endif // MESH_H
