#include "Graphics/Mesh.h"

Mesh::Mesh(std::shared_ptr<Shader> shader, bool isDisplayed, GLenum displayShape)
    : Mesh(std::vector<Vector3>(), shader, isDisplayed, displayShape)
{

}
Mesh::Mesh(std::vector<Vector3> _vertexArray, std::shared_ptr<Shader> shader, bool isDisplayed, GLenum displayShape)
    : vertexArray(_vertexArray), shader(shader), isDisplayed(isDisplayed), displayShape(displayShape)
{
    this->fromArray(_vertexArray);
}
Mesh::Mesh(std::vector<float> _vertexArrayFloat, std::shared_ptr<Shader> shader, bool isDisplayed, GLenum displayShape)
    : shader(shader), isDisplayed(isDisplayed), displayShape(displayShape)
{/*
    for(size_t i = 0; i < _vertexArrayFloat.size(); i+=3)
    {
        this->vertexArray.push_back(Vector3(_vertexArrayFloat[i],
                                            _vertexArrayFloat[i + 1],
                                            _vertexArrayFloat[i + 2]));
    }*/
    this->fromArray(_vertexArrayFloat); // This is time consuming, but less code to do
}

void Mesh::clear()
{
    this->vertexArray.clear();
    this->colorsArray.clear();
}
Mesh Mesh::merge(std::shared_ptr<Mesh> other, bool recomputeIndices)
{
    this->vertexArray.insert(this->vertexArray.end(), other->vertexArray.begin(), other->vertexArray.end());
    this->colorsArray.insert(this->colorsArray.end(), other->colorsArray.begin(), other->colorsArray.end());

    if (recomputeIndices) {
        this->fromArray(this->vertexArray);
    }
    return *this;
}
Mesh Mesh::merge(std::vector<std::shared_ptr<Mesh>> others)
{
    for (std::vector<std::shared_ptr<Mesh>>::iterator it = others.begin(); it != others.end(); it++)
        merge(*it, false);

    this->fromArray(this->vertexArray);
    return *this;
}

Mesh Mesh::fromArray(std::vector<Vector3> vertices, std::vector<int> indices)
{
    this->vertexArray = vertices;
    computeIndices(indices);
    this->vertexArrayFloat = Vector3::toArray(vertices);

    this->computeNormals();
    this->computeColors();
    return *this;
}
Mesh Mesh::fromArray(std::vector<float> vertices, std::vector<int> indices)
{
    std::vector<Vector3> vecVertices;
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        vecVertices.push_back(Vector3(vertices[i], vertices[i+1], vertices[i+2]));
    }
    return this->fromArray(vecVertices, indices);
}

void Mesh::computeIndices(std::vector<int> indices)
{
    if (!useIndices)
        return;
    if (indices.size() > 0) {
        for (size_t i = 0; i < this->vertexArray.size(); i++) {
            this->indices[i] = indices[i];
        }
    } else {
        this->indices.clear();
        this->vectorToIndex.clear();
        int index = 0;
        int mainIndex = 0;
        for (Vector3& pos : this->vertexArray) {
            std::tuple<int, int, int> vecIndex = (pos*100).toIntTuple();
            std::map<std::tuple<int, int, int>, int>::iterator it = this->vectorToIndex.find(vecIndex);
            if (it == this->vectorToIndex.end()) {
                // This index will be the main index for all points at this position
                this->vectorToIndex[vecIndex] = index;
                mainIndex = index;
            } else {
                mainIndex = it->second;

            }
            this->indices[index] = mainIndex;
            index++;
        }
    }
}

void Mesh::computeNormals()
{
    this->normalsArray.clear();
    this->normalsArrayIndex.clear();
    this->normalsArrayFloat.clear();
    this->normalsArrayIndex.resize(this->indices.size());
    for (size_t i = 0; i < this->vertexArray.size(); i+=3)
    {
        Vector3 normal0 = (this->vertexArray[i+2] - this->vertexArray[i]).cross(
                    (this->vertexArray[i+1] - this->vertexArray[i]));
        Vector3 normal1 = (this->vertexArray[i] - this->vertexArray[i+1]).cross(
                    (this->vertexArray[i+2] - this->vertexArray[i+1]));
        Vector3 normal2 = (this->vertexArray[i+1] - this->vertexArray[i+2]).cross(
                    (this->vertexArray[i] - this->vertexArray[i+2]));
        if(useIndices) {
            this->normalsArrayIndex[this->indices[i  ]] += normal0;
            this->normalsArrayIndex[this->indices[i+1]] += normal1;
            this->normalsArrayIndex[this->indices[i+2]] += normal2;
        } else {
            std::vector<float> n0 = Vector3::toArray(normal0);
            std::vector<float> n1 = Vector3::toArray(normal1);
            std::vector<float> n2 = Vector3::toArray(normal2);
            this->normalsArrayFloat.insert(this->normalsArrayFloat.end(), n0.begin(), n0.end());
            this->normalsArrayFloat.insert(this->normalsArrayFloat.end(), n1.begin(), n1.end());
            this->normalsArrayFloat.insert(this->normalsArrayFloat.end(), n2.begin(), n2.end());
        }
    }
    if(useIndices) {
        for (size_t i = 0; i < this->vertexArray.size(); i++)
        {
            std::vector<float> meanNormal = Vector3::toArray(this->normalsArrayIndex[this->indices[i]].normalize());
            this->normalsArrayFloat.insert(this->normalsArrayFloat.end(), meanNormal.begin(), meanNormal.end());
        }
    }
}
void Mesh::computeColors()
{
    this->colorsArrayIndex.clear();
    this->colorsArrayFloat.clear();
    if (!useIndices) {
        for (size_t i = 0; i < this->colorsArray.size(); i++)
        {
            std::vector<float> color = Vector3::toArray(this->colorsArray[i]);
            this->colorsArrayFloat.insert(this->colorsArrayFloat.end(), color.begin(), color.end());
        }
    } else {
        this->colorsArrayIndex.resize(this->indices.size());
        if (this->indices.size() == 0)
            return;
        std::vector<int> numVertexPerIndex(this->indices.size(), 0);
        for (size_t i = 0; i < this->colorsArray.size(); i++)
        {
            this->colorsArrayIndex[this->indices[i]] += this->colorsArray[i];
            numVertexPerIndex[this->indices[i]] ++;
        }
        for (size_t i = 0; i < this->colorsArray.size(); i++)
        {
            std::vector<float> meanNormal = Vector3::toArray(this->colorsArrayIndex[this->indices[i]] / (float)numVertexPerIndex[this->indices[i]]);
            this->colorsArrayFloat.insert(this->colorsArrayFloat.end(), meanNormal.begin(), meanNormal.end());
        }
    }
}

void Mesh::update()
{
    if (!bufferReady)
    {
        this->bufferID = GlobalsGL::newBufferId();
        GlobalsGL::generateBuffers();
    }
    pushToBuffer();

}

void Mesh::pushToBuffer()
{
    // Vertex
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->vertexArrayFloat.size() * sizeof(float), &this->vertexArrayFloat.front(), GL_STATIC_DRAW);

    // Textures

    // Normals
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 2]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->normalsArrayFloat.size() * sizeof(float), &this->normalsArrayFloat.front(), GL_STATIC_DRAW);

    // Random stuff
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 3]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->colorsArrayFloat.size() * sizeof(float), &this->colorsArrayFloat.front(), GL_STATIC_DRAW);
    this->bufferReady = true;
}
void Mesh::display(GLenum shape)
{
    if (!isDisplayed) return;
    if (shape != 0) this->displayShape = shape;
    this->update();
    if(this->shader != nullptr)
        this->shader->use();
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID]);
    GlobalsGL::f()->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    GlobalsGL::f()->glEnableVertexAttribArray(0);

    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 2]);
    GlobalsGL::f()->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
    GlobalsGL::f()->glEnableVertexAttribArray(2);

    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 3]);
    GlobalsGL::f()->glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
    GlobalsGL::f()->glEnableVertexAttribArray(3);
    GlobalsGL::checkOpenGLError();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDrawArrays(this->displayShape, 0, this->vertexArrayFloat.size()/3);

    if (this->shader != nullptr) {
        GlobalsGL::printShaderErrors(this->shader->vShader);
        GlobalsGL::printShaderErrors(this->shader->fShader);
//        GlobalsGL::printShaderErrors(this->shader->gShader);
        GlobalsGL::printProgramErrors(this->shader->programID);
        GlobalsGL::checkOpenGLError();
    }
}
