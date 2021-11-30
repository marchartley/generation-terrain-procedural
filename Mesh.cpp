#include "Mesh.h"

Mesh::Mesh()
{

}
Mesh::Mesh(std::vector<Vector3> _vertexArray)
    : vertexArray(_vertexArray)
{
    this->fromArray(_vertexArray);
}
Mesh::Mesh(std::vector<float> _vertexArrayFloat)
//    : vertexArrayFloat(_vertexArrayFloat)
{
    for(size_t i = 0; i < _vertexArrayFloat.size(); i+=3)
    {
        this->vertexArray.push_back(Vector3(_vertexArrayFloat[i],
                                            _vertexArrayFloat[i + 1],
                                            _vertexArrayFloat[i + 2]));
    }
    this->fromArray(this->vertexArray); // This is time consuming, but less code to do
}
Mesh Mesh::fromArray(std::vector<Vector3> vertices)
{
    this->vertexArray = vertices;
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
    this->vertexArrayFloat = Vector3::toArray(vertices);

    this->computeNormals();
    this->computeColors();
    return *this;
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
        this->normalsArrayIndex[this->indices[i  ]] += normal0;
        this->normalsArrayIndex[this->indices[i+1]] += normal1;
        this->normalsArrayIndex[this->indices[i+2]] += normal2;
    }
    for (size_t i = 0; i < this->vertexArray.size(); i++)
    {
        std::vector<float> meanNormal = Vector3::toArray(this->normalsArrayIndex[this->indices[i]].normalize());
        this->normalsArrayFloat.insert(this->normalsArrayFloat.end(), meanNormal.begin(), meanNormal.end());
    }
}
void Mesh::computeColors()
{
    this->colorsArrayIndex.clear();
    this->colorsArrayFloat.clear();
    this->colorsArrayIndex.resize(this->indices.size());
    if (this->indices.size() == 0)
        return;
    int* numVertexPerIndex = new int[this->indices.size()]{0};
    for (size_t i = 0; i < this->vertexArray.size(); i++)
    {
        this->colorsArrayIndex[this->indices[i]] += this->colorsArray[i];
        numVertexPerIndex[this->indices[i]] ++;
    }
    for (size_t i = 0; i < this->vertexArray.size(); i++)
    {
        std::vector<float> meanNormal = Vector3::toArray(this->colorsArrayIndex[this->indices[i]] / (float)numVertexPerIndex[this->indices[i]]);
        this->colorsArrayFloat.insert(this->colorsArrayFloat.end(), meanNormal.begin(), meanNormal.end());
    }
    delete[] numVertexPerIndex;
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
    glDrawArrays(shape, 0, this->vertexArrayFloat.size()/3);

    if (this->shader != nullptr) {
        GlobalsGL::printShaderErrors(this->shader->vShader);
        GlobalsGL::printShaderErrors(this->shader->fShader);
//        GlobalsGL::printShaderErrors(this->shader->gShader);
        GlobalsGL::printProgramErrors(this->shader->programID);
        GlobalsGL::checkOpenGLError();
    }
}
