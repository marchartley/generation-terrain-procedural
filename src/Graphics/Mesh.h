#ifndef MESH_H
#define MESH_H

#include "Utils/Globals.h"
#include "DataStructure/Vector3.h"
#include "Graphics/Shader.h"
#include <vector>
#include <map>

class Mesh //: public std::enable_shared_from_this<Mesh>
{
public:
    Mesh(std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
//    Mesh(Shader shader, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<Vector3> _vertexArray, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<float> _vertexArrayFloat, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh& fromArray(std::vector<Vector3> vertices, std::vector<int> indices = std::vector<int>());
    Mesh& fromArray(std::vector<float> vertices, std::vector<int> indices = std::vector<int>());

    Mesh& fromStl(std::string filename);

    Mesh& normalize();
    Mesh& scale(float factor);
    Mesh& scale(float factor_x, float factor_y, float factor_z);
    Mesh& scale(Vector3 factor);
    Mesh& translate(Vector3 translation);
    Mesh& translate(float translation_x, float translation_y, float translation_z);

    void clear();
    Mesh merge(std::shared_ptr<Mesh> other, bool recomputeIndices = true);
    Mesh merge(std::vector<std::shared_ptr<Mesh>> others);

    void update();
    void pushToBuffer();
    void display(GLenum shape = -1, float lineWeight = 1);
    void displayNormals();

    void shareShader(std::shared_ptr<Shader> sharedShader) { this->shader = sharedShader; }
    void shareShader(const Mesh& otherMesh) { this->shader = otherMesh.shader; }

    bool isHidden() { return !this->isDisplayed; }
    void hide() { this->isDisplayed = false; }
    void show() { this->isDisplayed = true; }

    void setShader(std::shared_ptr<Shader> shader);

    std::vector<std::vector<Vector3>> getTriangles(std::vector<int> indices = std::vector<int>());

    std::string toOFF();


    unsigned int bufferID;
    bool bufferReady = false;

    bool useIndices = true;

    bool needToUpdatePositions = true;
    bool needToUpdateTextures = true;
    bool needToUpdateNormals = true;
    bool needToUpdateColors = true;

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
    bool isDisplayed = true;
    GLenum displayShape;

    static void setShaderToAllMeshesWithoutShader(Shader newShader);

    static std::vector<Mesh*> all_meshes;

    GLuint vao;
    GLuint vbo[4];

//    bool useModernOpenGL = true;
};

#endif // MESH_H
