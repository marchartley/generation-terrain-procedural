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
    : vertexArrayFloat(_vertexArrayFloat)
{
    for(size_t i = 0; i < _vertexArrayFloat.size(); i+=3)
    {
        this->vertexArray.push_back(Vector3(_vertexArrayFloat[i],
                                            _vertexArrayFloat[i + 1],
                                            _vertexArrayFloat[i + 2]));
    }
}
Mesh Mesh::fromArray(std::vector<Vector3> vertices)
{
    this->vertexArray = vertices;
    this->vertexArrayFloat = Vector3::toArray(vertices);

    this->computeNormals();
    return *this;
}

void Mesh::computeNormals()
{
    this->normalsArray.clear();
    for (size_t i = 0; i < this->vertexArray.size(); i+=3)
    {
        Vector3 normal = (this->vertexArray[i+1] - this->vertexArray[i]).cross(
                    (this->vertexArray[i+2] - this->vertexArray[i])).normalize();
//        for (int r = 0; r < 3; r++)
        this->normalsArray.push_back((this->vertexArray[i+1] - this->vertexArray[i]).cross((this->vertexArray[i+2] - this->vertexArray[i])).normalize() * -1.0);
        // Small fix, the third one is counted twice...
//        this->normalsArray.push_back((this->vertexArray[i] - this->vertexArray[i+1]).cross((this->vertexArray[i+2] - this->vertexArray[i+1])).normalize());
        this->normalsArray.push_back((this->vertexArray[i] - this->vertexArray[i+2]).cross((this->vertexArray[i+1] - this->vertexArray[i+2])).normalize() * -1.0);
        this->normalsArray.push_back((this->vertexArray[i] - this->vertexArray[i+2]).cross((this->vertexArray[i+1] - this->vertexArray[i+2])).normalize() * -1.0);
    }
    this->normalsArrayFloat = Vector3::toArray(this->normalsArray);
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
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->colorArrayFloat.size() * sizeof(float), &this->colorArrayFloat.front(), GL_STATIC_DRAW);
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
