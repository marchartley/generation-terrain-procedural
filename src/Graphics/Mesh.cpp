#include "Graphics/Mesh.h"

#include "Utils/stl_reader.h"

std::vector<Mesh*> Mesh::all_meshes;

Mesh::Mesh(std::shared_ptr<Shader> shader, bool isDisplayed, GLenum displayShape)
    : Mesh(std::vector<Vector3>(), shader, isDisplayed, displayShape)
{

}
Mesh::Mesh(std::vector<Vector3> _vertexArray, std::shared_ptr<Shader> shader, bool isDisplayed, GLenum displayShape)
    : vertexArray(_vertexArray), shader(shader), isDisplayed(isDisplayed), displayShape(displayShape)
{
    this->fromArray(_vertexArray);
    Mesh::all_meshes.push_back(this);
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
    Mesh::all_meshes.push_back(this);
}

void Mesh::clear()
{
    vertexArray.clear();
    vertexArrayFloat.clear();
    normalsArray.clear();
    normalsArrayIndex.clear();
    normalsArrayFloat.clear();
    colorsArray.clear();
    colorsArrayIndex.clear();
    colorsArrayFloat.clear();
    indices.clear();
    vectorToIndex.clear();
    this->needToUpdatePositions = true;
    this->needToUpdateColors = true;
    this->needToUpdateNormals = true;
    this->needToUpdateTextures = true;
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

Mesh& Mesh::fromArray(std::vector<Vector3> vertices, std::vector<int> indices)
{
    if (indices.empty()) {
        this->vertexArray = vertices;
        computeIndices(indices);
        this->vertexArrayFloat = Vector3::toArray(this->vertexArray);
    }
    this->needToUpdatePositions = true;

    this->computeNormals();
    this->computeColors();
    return *this;
}
Mesh& Mesh::fromArray(std::vector<float> vertices, std::vector<int> indices)
{
    std::vector<Vector3> vecVertices;
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        vecVertices.push_back(Vector3(vertices[i], vertices[i+1], vertices[i+2]));
    }
    return this->fromArray(vecVertices, indices);
}

Mesh &Mesh::fromStl(std::string filename)
{
    try {
        stl_reader::StlMesh <float, unsigned int> stl_model (filename);
        std::vector<Vector3> vertices;
        for (size_t ivert = 0; ivert < stl_model.num_vrts(); ivert++) {
            Vector3 v(stl_model.vrt_coords(ivert));
            vertices.push_back(v);
        }
        std::vector<Vector3> trianglesVec3;
        std::vector<int> indices;
        for (size_t itri = 0; itri < stl_model.num_tris(); itri++) {
            // Todo : Put that back in the right order... But this means I have to do that for ALL triangles
            int index = stl_model.tri_corner_ind(itri, 0);
            trianglesVec3.push_back(Vector3(vertices[index]));
            indices.push_back(index);

            index = stl_model.tri_corner_ind(itri, 1);
            trianglesVec3.push_back(Vector3(vertices[index]));
            indices.push_back(index);

            index = stl_model.tri_corner_ind(itri, 2);
            trianglesVec3.push_back(Vector3(vertices[index]));
            indices.push_back(index);
        }
        this->fromArray(trianglesVec3); //vertices, indices);
    }
    catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    return *this;
}

Mesh& Mesh::normalize()
{
    Vector3 minAABBox(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vector3 maxAABBox(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    Vector3 center;

    for (const auto& p : this->vertexArray) {
        center += p;
        minAABBox = Vector3(std::min(minAABBox.x, p.x), std::min(minAABBox.y, p.y), std::min(minAABBox.z, p.z));
        maxAABBox = Vector3(std::max(maxAABBox.x, p.x), std::max(maxAABBox.y, p.y), std::max(maxAABBox.z, p.z));
    }
    center /= (float)vertexArray.size();

    this->translate(center * -1.f);
    this->scale(1.f / (maxAABBox - minAABBox).norm());

    return *this;
}

Mesh& Mesh::scale(float factor)
{
    return scale(factor, factor, factor);
}

Mesh& Mesh::scale(float factor_x, float factor_y, float factor_z)
{
    return scale(Vector3(factor_x, factor_y, factor_z));
}

Mesh& Mesh::scale(Vector3 factor)
{
    for (auto& p : this->vertexArray)
        p *= factor;

    this->vertexArrayFloat = Vector3::toArray(this->vertexArray);
    this->needToUpdatePositions = true;

    return *this;
}

Mesh& Mesh::translate(Vector3 translation)
{
    for (auto& p : this->vertexArray)
        p += translation;

    this->vertexArrayFloat = Vector3::toArray(this->vertexArray);
    this->needToUpdatePositions = true;

    return *this;
}

Mesh& Mesh::translate(float translation_x, float translation_y, float translation_z)
{
    return translate(Vector3(translation_x, translation_y, translation_z));
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

    for (size_t i = 0; i < normalsArrayFloat.size(); i += 3) {
        this->normalsArray.push_back(Vector3(normalsArrayFloat[i], normalsArrayFloat[i+1], normalsArrayFloat[i+2]));
    }
    this->needToUpdateNormals = true;
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
    this->needToUpdateColors = true;
}

void Mesh::setShaderToAllMeshesWithoutShader(Shader newShader)
{
    /*for (Mesh*& mesh : Mesh::all_meshes) {
        if (!mesh->shader || mesh->shader->vShader < 0)
            mesh->shader = std::make_shared<Shader>(newShader);
    }*/
}

void Mesh::update()
{
    if (!bufferReady)
    {
        GlobalsGL::f()->glGenVertexArrays(1, &this->vao);
        GlobalsGL::f()->glGenBuffers(4, this->vbo);
        GlobalsGL::checkOpenGLError();
    }
    if (!this->shader && Shader::default_shader) {
        this->shader = std::make_shared<Shader>(*Shader::default_shader);
    }
    if (this->shader && Shader::allShaders.find(this->shader) == Shader::allShaders.end()) {
        Shader::allShaders.insert(this->shader);
    }
    if (this->shader) {
        pushToBuffer();
    }

}

void Mesh::pushToBuffer()
{
#if useModernOpenGL
        this->shader->use();
//        if(!this->shader->use()) return;
        GlobalsGL::f()->glBindVertexArray(this->vao); //GlobalsGL::vao[this->bufferID]);
        // Vertex
        GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, this->vbo[0]); // GlobalsGL::vbo[this->bufferID * 10 + 0]);
        if (needToUpdatePositions)
            GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->vertexArrayFloat.size() * sizeof(float), &this->vertexArrayFloat.front(), GL_STATIC_DRAW);
        GlobalsGL::f()->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        GlobalsGL::f()->glEnableVertexAttribArray(0);

        // Textures

        // Normals
        GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, this->vbo[2]); // GlobalsGL::vbo[this->bufferID * 10 + 2]);
        if (needToUpdateNormals)
            GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->normalsArrayFloat.size() * sizeof(float), &this->normalsArrayFloat.front(), GL_STATIC_DRAW);
        GlobalsGL::f()->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
        GlobalsGL::f()->glEnableVertexAttribArray(2);

        // Colors
        GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, this->vbo[3]); // GlobalsGL::vbo[this->bufferID * 10 + 3]);
        if (needToUpdateColors)
            GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, this->colorsArrayFloat.size() * sizeof(float), &this->colorsArrayFloat.front(), GL_STATIC_DRAW);
        GlobalsGL::f()->glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
        GlobalsGL::f()->glEnableVertexAttribArray(3);
//        if(this->shader->use()) return;
#endif
    this->bufferReady = true;
}
void Mesh::display(GLenum shape, float lineWeight)
{
    if (!isDisplayed)
        return;
    if (shape != -1) this->displayShape = shape;
    this->update();
    if(this->shader != nullptr) {
        this->shader->use();
        this->shader->setBool("cullFace", this->cullFace);
    }
    if(!this->shader->use()) return;
    //    glEnable(GL_DEPTH_TEST);
    //    glDepthFunc(GL_LEQUAL);
    GLfloat previousLineWidth[1];
    glGetFloatv(GL_LINE_WIDTH, previousLineWidth);
    glLineWidth(lineWeight);
    if (this->shader != nullptr) {
        if (this->shader->vShader != -1)
            GlobalsGL::printShaderErrors(this->shader->vShader);
        if (this->shader->fShader != -1)
            GlobalsGL::printShaderErrors(this->shader->fShader);
        if (this->shader->gShader != -1)
            GlobalsGL::printShaderErrors(this->shader->gShader);
        GlobalsGL::printProgramErrors(this->shader->programID);
        GlobalsGL::checkOpenGLError();
    }

#if useModernOpenGL
    GlobalsGL::f()->glBindVertexArray(this->vao); // GlobalsGL::vao[this->bufferID]);
    GlobalsGL::checkOpenGLError();

    if (this->displayShape == GL_TRIANGLES)
        glDrawArrays(this->displayShape, 0, this->vertexArrayFloat.size()/3);
    else if (this->displayShape == GL_LINE_STRIP)
        glDrawArrays(this->displayShape, 0, this->vertexArrayFloat.size() - 1);
    else
        glDrawArrays(this->displayShape, 0, this->vertexArrayFloat.size());
#else
        glEnable(GL_RESCALE_NORMAL);
        glBegin(this->displayShape);
        for (size_t i = 0; i < this->vertexArray.size(); i++) {
            if (this->colorsArray.size() > i)
                glColor3f(colorsArray[i].x, colorsArray[i].y, colorsArray[i].z);
            if (this->normalsArray.size() > i)
                glNormal3f(normalsArray[i].x, normalsArray[i].y, normalsArray[i].z);
            glVertex3f(vertexArray[i].x, vertexArray[i].y, vertexArray[i].z);
        }
        glEnd();
#endif
    glLineWidth(previousLineWidth[0]);
}

void Mesh::displayNormals()
{
    glBegin(GL_LINES);
    for (size_t i = 0; i < this->vertexArray.size(); i++) {
        glVertex3f(this->vertexArray[i].x, this->vertexArray[i].y, this->vertexArray[i].z);
        glVertex3f(this->vertexArray[i].x + this->normalsArray[i].x, this->vertexArray[i].y + this->normalsArray[i].y, this->vertexArray[i].z + this->normalsArray[i].z);
    }
    glEnd();
}

void Mesh::setShader(std::shared_ptr<Shader> shader)
{

}

std::vector<std::vector<Vector3> > Mesh::getTriangles(std::vector<int> indices)
{
    if (indices.empty()) {
        indices = std::vector<int>(this->vertexArray.size());
        for (size_t i = 0; i < indices.size(); i++)
            indices[i] = i;
    }

    std::vector<std::vector<Vector3>> triangles;
    for (size_t i = 0; i < indices.size(); i += 3) {
        triangles.push_back(std::vector<Vector3>({
                                                     this->vertexArray[indices[i+0]],
                                                     this->vertexArray[indices[i+1]],
                                                     this->vertexArray[indices[i+2]],
                                                 }));
    }
    return triangles;
}


std::string Mesh::toOFF()
{
    std::string out = "OFF\n# ICAR Team - Terrain generation\n";

    out += std::to_string(vertexArray.size()) + " " + std::to_string(vertexArray.size() / 3) + " 0\n# Vertices\n";
    for (size_t i = 0; i < vertexArray.size(); i++)
        out += std::to_string(vertexArray[i].x) + " " + std::to_string(vertexArray[i].y) + " " + std::to_string(vertexArray[i].z) + "\n";

    out += "# Faces\n";
    for (size_t i = 0; i < vertexArray.size(); i += 3)
        out += "3 " + std::to_string(i + 0) + " " + std::to_string(i + 1) + " " + std::to_string(i + 2) + "\n";
    return out;
}
