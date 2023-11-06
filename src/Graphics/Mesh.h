#ifndef MESH_H
#define MESH_H

#include "Utils/Globals.h"
#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "Graphics/Shader.h"
#include <vector>
#include <map>

class Mesh
{
public:
    Mesh(std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<Vector3> _vertexArray, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
    Mesh(std::vector<float> _vertexArrayFloat, std::shared_ptr<Shader> shader = nullptr, bool isDisplayed = true, GLenum displayShape = GL_TRIANGLES);
//    Mesh& fromArray(std::vector<std::vector<Vector3>> triangles, std::vector<int> indices = std::vector<int>());
    Mesh& fromArray(std::vector<Vector3> vertices, std::vector<int> indices = std::vector<int>());
    Mesh& fromArray(std::vector<float> vertices, std::vector<int> indices = std::vector<int>());

    Mesh& fromStl(std::string filename);
    Mesh& fromFBX(std::string filename);

    Mesh& normalize();
    Mesh& scale(float factor);
    Mesh& scale(float factor_x, float factor_y, float factor_z);
    Mesh& scale(const Vector3& factor);
    Mesh& translate(const Vector3& translation);
    Mesh& translate(float translation_x, float translation_y, float translation_z);
    Mesh& rotate(const Vector3& rotation);
    Mesh& rotate(float rotation_x, float rotation_y, float rotation_Z);

    void clear();
    Mesh merge(const Mesh &other, bool recomputeIndices = true);
    Mesh merge(std::vector<Mesh> others);

    void update();
    void pushToBuffer();
    void display(GLenum shape = -1, float lineWeight = 1);
    void displayWithOutlines(std::vector<float> faceColor, GLenum shape = -1, float lineWeight = 1);
    void displayNormals();

    void displayAsScalarField(GridF field, const Vector3& cameraPosition, std::vector<float> isoValues = {0.5f});
    void displayAsVectorField(GridV3 field, const Vector3& finalDimensions = Vector3(false), float maxMaginitude = -1, bool normalize = false);

    void shareShader(std::shared_ptr<Shader> sharedShader) { this->shader = sharedShader; }
    void shareShader(const Mesh& otherMesh) { this->shader = otherMesh.shader; }

    bool isHidden() { return !this->isDisplayed; }
    void hide() { this->isDisplayed = false; }
    void show() { this->isDisplayed = true; }

    void setShader(std::shared_ptr<Shader> shader);

    void reorderVertices(const Vector3& camPos);
    void reorderLines(const Vector3& camPos);
    void reorderTriangles(const Vector3& camPos);
    void reorderAny(const Vector3& camPos, int nbVertexToUse);

    std::vector<std::vector<Vector3>> getTriangles(std::vector<int> indices = std::vector<int>()) const;

    std::string toOBJ();
    std::string toOFF();
    std::string toSTL();

    static Mesh applyMarchingCubes(const GridF &values);

    GridI voxelize(const Vector3& dimensions) const;
    GridI voxelizeSurface(const Vector3& dimensions) const;

    bool isWatertight();


    static std::vector<Vector3> getPointsForArrow(const Vector3& from, const Vector3& to);
    static Mesh createVectorField(GridV3 field, const Vector3& finalDimensions = Vector3(false), Mesh *mesh = nullptr, float maxMaginitude = -1, bool normalize = false, bool displayArrow = false);

    static void displayScalarField(GridF field, Mesh& mesh, const Vector3& cameraPosition, std::vector<float> isoValues = {0.5f});


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

    bool cullFace = true;

    Mesh extractGeometryFromShaders(GridF &values);

    static void setShaderToAllMeshesWithoutShader(Shader newShader);

    static std::vector<Mesh*> all_meshes;

    GLuint vao;
    GLuint vbo[4];
};

#endif // MESH_H
