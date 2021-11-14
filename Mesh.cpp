#include "Mesh.h"

Mesh::Mesh()
{

}
Mesh::Mesh(std::vector<Vector3> _vertexArray)
    : vertexArray(_vertexArray)
{
    vertexArrayFloat = Vector3::toArray(_vertexArray);
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

void Mesh::update()
{
    if (!bufferReady)
    {
        this->pushToBuffer();
    } else {
//        GlobalsGL::f()->glBufferData(this->)
    }
}

void Mesh::pushToBuffer()
{
    this->bufferID = GlobalsGL::newBufferId();
    GlobalsGL::f()->glGenVertexArrays(1, GlobalsGL::vao);
    GlobalsGL::f()->glBindVertexArray(GlobalsGL::vao[this->bufferID]);
    GlobalsGL::f()->glGenBuffers(numVBOs, GlobalsGL::vbo);

    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->vertexArrayFloat.size() * sizeof(float), &this->vertexArrayFloat.front(), GL_STATIC_DRAW);
    this->bufferReady = true;
}
void Mesh::display()
{
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID]);
    GlobalsGL::f()->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    GlobalsGL::f()->glEnableVertexAttribArray(0);
    GlobalsGL::checkOpenGLError();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDrawArrays(GL_TRIANGLES, 0, this->vertexArrayFloat.size()/3);
}
