#include "CubeMesh.h"

std::vector<Vector3> CubeMesh::cubesVertices = {
    {0, 1, 0}, {1, 1, 0}, {1, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0},
    {1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 0}, {1, 1, 1}, {1, 0, 1},
    {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {1, 0, 1}, {0, 1, 1}, {0, 0, 1},
    {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0},
    {0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 0}, {1, 0, 1}, {0, 0, 1},
    {0, 1, 0}, {0, 1, 1}, {1, 1, 1}, {0, 1, 0}, {1, 1, 1}, {1, 1, 0}
};
std::vector<Vector3> CubeMesh::cubesEdgesVertices = {
    {0, 0, 0}, {0, 1, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 0}, {1, 0, 0}, {1, 0, 0}, {0, 0, 0},
    {0, 0, 1}, {0, 1, 1}, {0, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 0, 1}, {1, 0, 1}, {0, 0, 1},
    {0, 0, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0, 1, 0}, {0, 1, 1}
};
CubeMesh::CubeMesh()
{
    this->fromArray(CubeMesh::cubesVertices);
}

void CubeMesh::update(std::vector<Vector3> positions)
{
    if (positions.empty()) positions.push_back(Vector3(0, 0, 0));
    if (!bufferReady)
    {
        this->bufferID = GlobalsGL::newBufferId();
//        GlobalsGL::generateBuffers();
    }
    pushToBuffer(positions);

}

void CubeMesh::pushToBuffer(std::vector<Vector3> positions)
{
    this->positions = positions;
    this->positionsFloat = Vector3::toArray(positions);

    GlobalsGL::f()->glBindVertexArray(GlobalsGL::vao[this->bufferID]);

    // Vertex
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID * 10 + 0]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->vertexArrayFloat.size() * sizeof(float), &this->vertexArrayFloat.front(), GL_STATIC_DRAW);
    GlobalsGL::f()->glEnableVertexAttribArray(0);
    GlobalsGL::f()->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    // Textures

    // Normals
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID * 10 + 2]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->normalsArrayFloat.size() * sizeof(float), &this->normalsArrayFloat.front(), GL_STATIC_DRAW);
    GlobalsGL::f()->glEnableVertexAttribArray(2);
    GlobalsGL::f()->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

    // Colors
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID * 10 + 3]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->colorsArrayFloat.size() * sizeof(float), &this->colorsArrayFloat.front(), GL_STATIC_DRAW);
    GlobalsGL::f()->glEnableVertexAttribArray(3);
    GlobalsGL::f()->glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);

    // Offsets
    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID * 10 + 4]);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->positionsFloat.size() * sizeof(float), &this->positionsFloat.front(), GL_STATIC_DRAW);
    GlobalsGL::f()->glEnableVertexAttribArray(4);
    GlobalsGL::f()->glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 0, 0);
    GlobalsGL::f()->glVertexAttribDivisor(4, 1);

    this->bufferReady = true;
}
void CubeMesh::display(GLenum shape)
{
    this->shader = nullptr;
    if (!isDisplayed) return;
    if (shape != 0) this->displayShape = shape;
    this->update(this->positions);
    if(this->shader != nullptr)
        this->shader->use();
//    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 0]);
//    GlobalsGL::f()->glEnableVertexAttribArray(0);
//    GlobalsGL::f()->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

//    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 2]);
//    GlobalsGL::f()->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
//    GlobalsGL::f()->glEnableVertexAttribArray(2);

//    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 3]);
//    GlobalsGL::f()->glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
//    GlobalsGL::f()->glEnableVertexAttribArray(3);

//    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, GlobalsGL::vbo[this->bufferID + 4]);
//    GlobalsGL::f()->glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 0, 0);
//    GlobalsGL::f()->glEnableVertexAttribArray(4);
//    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, 0);
//    GlobalsGL::f()->glVertexAttribDivisor(4, 1);
    GlobalsGL::f()->glBindVertexArray(GlobalsGL::vao[this->bufferID]);
    GlobalsGL::checkOpenGLError();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    if (this->positions.empty()) return; // To delete
    GlobalsGL::f()->glDrawArraysInstanced(GL_TRIANGLES, 0, 1, 1); //this->vertexArrayFloat.size()/3, this->positions.size());

    if (this->shader != nullptr) {
        GlobalsGL::printShaderErrors(this->shader->vShader);
        GlobalsGL::printShaderErrors(this->shader->fShader);
//        GlobalsGL::printShaderErrors(this->shader->gShader);
        GlobalsGL::printProgramErrors(this->shader->programID);
        GlobalsGL::checkOpenGLError();
    }
}

std::vector<Vector3> CubeMesh::createTriangles(const Vector3& minPos, const Vector3& maxPos)
{
    std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
    Vector3 dim = maxPos - minPos;
    for (auto& p : box) {
        p = p * dim + minPos;
    }
    return box;
}

std::vector<Vector3> CubeMesh::createTriangles(AABBox box)
{
    return CubeMesh::createTriangles(box.min(), box.max());
}
