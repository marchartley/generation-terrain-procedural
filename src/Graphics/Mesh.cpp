#include "Graphics/Mesh.h"

#include "Graphics/MarchingCubes.h"
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

Mesh& Mesh::rotate(Vector3 rotation)
{
    for (auto& p : this->vertexArray)
        p.rotate(rotation);


    this->vertexArrayFloat = Vector3::toArray(this->vertexArray);
    this->needToUpdatePositions = true;

    return *this;
}

Mesh& Mesh::rotate(float rotation_x, float rotation_y, float rotation_Z)
{
    return this->rotate(Vector3(rotation_x, rotation_y, rotation_Z));
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
        if (vertexArray.size() < i + 3) continue;
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

Mesh Mesh::extractGeometryFromShaders(Matrix3<float>& values)
{
    /*
    GLuint vertices_positions_buffer;
    GlobalsGL::f()->glGenBuffers(1, &vertices_positions_buffer);
    GlobalsGL::f()->glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, vertices_positions_buffer);

    GLuint feedback_buffer;
    GlobalsGL::f()->glGenTransformFeedbacks(1, &feedback_buffer);

    GLuint myTBO;               // Transform Buffer Object (TBO)
    GlobalsGL::f()->glGenBuffers(1, &myTBO);          // Generate an OpenGL buffer

    GlobalsGL::f()->glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, feedback_buffer);
    GlobalsGL::f()->glBeginTransformFeedback(GL_TRIANGLES);

    this->display();

    GlobalsGL::f()->glEndTransformFeedback();

    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, vertices_positions_buffer);

    GLint dataSize_feedback = 0;
    GLint dataSize_vertices = 0;
    GlobalsGL::f()->glGetNamedBufferParameteriv(feedback_buffer, GL_BUFFER_SIZE, &dataSize_feedback);
    GlobalsGL::f()->glGetNamedBufferParameteriv(vertices_positions_buffer, GL_BUFFER_SIZE, &dataSize_vertices);
    std::cout << "Feedback has " << dataSize_feedback << "b data \t" << "Vertices has " << dataSize_feedback << "b data" << std::endl;
    *//*
    int maxNumEdgesFeedback = this->vertexArray.size() * 5; // 3 * (unitSphere.GetNumTriangles() + unitSphere.GetNumTrianglesInSlice());
    int sizeOfEdgeFeedback = 3 * sizeof(float); //3 * sizeof(int) + 4 * sizeof(float);
    int maxSizeOfFeedback = sizeOfEdgeFeedback* maxNumEdgesFeedback;

    GlobalsGL::f()->glBindBuffer(GL_ARRAY_BUFFER, myTBO);
    GlobalsGL::f()->glBufferData(GL_ARRAY_BUFFER, maxSizeOfFeedback, (void*)0, GL_STATIC_READ);
    GlobalsGL::f()->glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, myTBO);
    unsigned int myTBOquery;          // Will hold an OpenGL query object
    GlobalsGL::f()->glGenQueries(1, &myTBOquery);
    Mesh copy;

    GlobalsGL::f()->glEnable(GL_RASTERIZER_DISCARD);
//    GlobalsGL::f()->glUseProgram(copy.shader->programID);
//    this->shader->use();
//    GlobalsGL::f()->glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, feedback_buffer);
    GlobalsGL::f()->glBeginQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, myTBOquery);
    GlobalsGL::f()->glBeginTransformFeedback(GL_TRIANGLES);
    this->display();
//    GlobalsGL::f()->glFlush();
    GlobalsGL::f()->glEndTransformFeedback();
    GlobalsGL::f()->glEndQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN);
    unsigned int dataSize;
    GlobalsGL::f()->glGetQueryObjectuiv(myTBOquery, GL_QUERY_RESULT, &dataSize);
    GlobalsGL::f()->glDisable(GL_RASTERIZER_DISCARD);

//    GLint dataSize = 0;
//    GlobalsGL::f()->glGetNamedBufferParameteriv(feedback_buffer, GL_BUFFER_SIZE, &dataSize);
//    dataSize /= sizeof(float);
//    float *data = new float[dataSize / sizeof(float)];
    float *data = new float[dataSize * 1];

    GlobalsGL::f()->glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeOfEdgeFeedback*dataSize, data);
//    GlobalsGL::f()->glGetNamedBufferSubData(feedback_buffer, 0, dataSize, data);
    */
    Mesh copy = this->applyMarchingCubes(values);/*
    copy.vertexArray.clear();
    for (size_t i = 0; i < dataSize / (3 * sizeof(float)); i++) {
        Vector3 pos = Vector3(
                    data[3 * i + 0],
                    data[3 * i + 1],
                    data[3 * i + 2]);
        copy.vertexArray.push_back(pos);
    }
    std::cout << copy.vertexArray.size() / 3 << " triangles to store" << std::endl;*/
    std::ofstream off;
    off.open("test.off");
    off << copy.toOFF();
    off.close();

    std::ofstream stl;
    stl.open("test.stl");
    stl << copy.toSTL();
    stl.close();

//    GlobalsGL::f()->glDeleteTransformFeedbacks(1, &feedback_buffer);
//    GlobalsGL::f()->glDeleteBuffers(1, &myTBO);
//    GlobalsGL::f()->glDeleteBuffers(1, &vertices_positions_buffer);

//    GlobalsGL::f()->glDeleteQueries(1, &myTBOquery);
//    return copy;
    return *this;
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

std::string Mesh::toOBJ()
{
    std::string out = "# ICAR Team - Terrain generation\n# mesh.obj\n\no terrain\n\n";

    for (size_t i = 0; i < vertexArray.size(); i++)
        out += "v " + std::to_string(vertexArray[i].x) + " " + std::to_string(vertexArray[i].z) + " " + std::to_string(vertexArray[i].y) + "\n";

//    out += "\n";
//    for (size_t i = 0; i < vertexArray.size(); i++)
//        out += "vt 0.0 0.0\n";

//    out += "\n";
//    for (size_t i = 0; i < normalsArray.size(); i++)
//        out += "vn " + std::to_string(normalsArray[i].x) + " " + std::to_string(normalsArray[i].y) + " " + std::to_string(normalsArray[i].z) + "\n";

//    out += "\n";
//    for (size_t i = 0; i < vertexArray.size(); i += 3)
//        out += "f " + std::to_string(i + 0) + "/" + std::to_string(i + 0) + "/" + std::to_string(i + 0) + " " + std::to_string(i + 1) + "/" + std::to_string(i + 1) + "/" + std::to_string(i + 1) + " " + std::to_string(i + 2) + "/" + std::to_string(i + 2) + "/" + std::to_string(i + 2) + "\n";
    for (size_t i = 0; i < vertexArray.size(); i += 3)
        out += "f " + std::to_string(i + 1) + " " + std::to_string(i + 2) + " " + std::to_string(i + 3) + "\n";
    return out;
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

std::string Mesh::toSTL()
{
    std::string objName = "mesh-output";
    //binary file
    std::string header_info = "solid " + objName;

    std::ostringstream file;

    file << header_info << std::endl;
    bool wasUsingIndices = this->useIndices;
    this->useIndices = false;
    this->computeNormals();
    for (size_t _i = 0; _i < this->vertexArray.size() / 3; _i++) { // (const auto& t : Triangles) {
        size_t i = 3 * _i;
        Vector3 n = (this->normalsArray[i] + this->normalsArray[i+1] + this->normalsArray[i+2])/3.f;
        Vector3 v1 = this->vertexArray[i + 0];
        Vector3 v2 = this->vertexArray[i + 1];
        Vector3 v3 = this->vertexArray[i + 2];
        file << "\t" << "facet normal" << " " << n.x << " " << n.y << " " << n.z << "\n";
        file << "\t\t" << "outer loop" << "\n";
        file << "\t\t\t" << " " << "vertex" << " " << v1.x << " " << v1.y << " " << v1.z << "\n";
        file << "\t\t\t" << " " << "vertex" << " " << v2.x << " " << v2.y << " " << v2.z << "\n";
        file << "\t\t\t" << " " << "vertex" << " " << v3.x << " " << v3.y << " " << v3.z << "\n";
        file << "\t\t" << "endloop" << "\n";
        file << "\t" << "endfacet" << "\n";
    }
    file << "endsolid" << " " << objName << "\n";
    this->useIndices = wasUsingIndices;
    return file.str();
    /*char head[80];
    std::strncpy(head,header_info.c_str(),sizeof(head)-1);
    char attribute[2] = "0";
    unsigned long nTriLong = triangles.size() ;

//    std::ofstream myfile;
    char* myText = new char[80 + 4 + 50 * nTriLong]; // Head = 80 + number facets = 4 + 50 per triangle


    myfile.open((Filename + "-out.stl").c_str(),  std::ios::out | std::ios::binary);
    myfile.write(head,sizeof(head));
    myfile.write((char*)&nTriLong,4);

    //write down every triangle
    for (std::vector<tri>::iterator it = triangles.begin(); it!=triangles.end(); it++){
        //normal vector coordinates

        myfile.write((char*)&it->m_n.m_x, 4);
        myfile.write((char*)&it->m_n.m_y, 4);
        myfile.write((char*)&it->m_n.m_z, 4);

        //p1 coordinates
        myfile.write((char*)&it->m_p1.m_x, 4);
        myfile.write((char*)&it->m_p1.m_y, 4);
        myfile.write((char*)&it->m_p1.m_z, 4);

        //p2 coordinates
        myfile.write((char*)&it->m_p2.m_x, 4);
        myfile.write((char*)&it->m_p2.m_y, 4);
        myfile.write((char*)&it->m_p2.m_z, 4);

        //p3 coordinates
        myfile.write((char*)&it->m_p3.m_x, 4);
        myfile.write((char*)&it->m_p3.m_y, 4);
        myfile.write((char*)&it->m_p3.m_z, 4);

        myfile.write(attribute,2);
    }

    myfile.close();*/
}

Mesh Mesh::applyMarchingCubes(Matrix3<float>& values)
{
    auto triTable = MarchingCubes::triangleTable;
    auto edges = MarchingCubes::cubeEdges;
    float offsetX = 0.f;
    float offsetY = 0.f;
    float offsetZ = 0.f;
    Vector3 scale(1.f, 1.f, 1.f);
    float isolevel = 0.f;

    Vector3 vertDecals[8] = {
        Vector3(0.0, 0.0, 0.0),
        Vector3(1.0, 0.0, 0.0),
        Vector3(1.0, 1.0, 0.0),
        Vector3(0.0, 1.0, 0.0),
        Vector3(0.0, 0.0, 1.0),
        Vector3(1.0, 0.0, 1.0),
        Vector3(1.0, 1.0, 1.0),
        Vector3(0.0, 1.0, 1.0)
    };
    std::function cubePos = [&](Vector3 voxelPos, int i) -> Vector3 {
        return voxelPos + vertDecals[i];
    };

    //Get vertex i value within current marching cube
    std::function cubeVal = [&](Vector3 pos) -> float {
        if (!values.checkCoord(pos)) return -1.f;
        return values.at(pos);
    };
    //Get vertex i value within current marching cube
    std::function cubeVali = [&](Vector3 voxelPos, int i) -> float {
        return cubeVal(cubePos(voxelPos, i));
    };

    //Get triangle table value
    std::function triTableValue = [&](int i, int j) -> int{
        return triTable[i][j];
    };

    //Compute interpolated vertex along an edge
    std::function vertexInterp = [&](float isolevel, Vector3 v0, float l0, Vector3 v1, float l1) -> Vector3 {
        float iso = std::clamp((isolevel-l0)/(l1-l0), 0.f, 1.f);
        return v0 * iso + v1 * (1.f - iso);
//        return mix(v0, v1, clamp() * scale + vec3(offsetX, offsetY, offsetZ);
    };

    std::function getPosition = [&](Vector3 position, Vector3 _offset) -> Vector3 {
    //    return position + vec4(_offset, 0.0);
        _offset += Vector3(offsetX, offsetY, offsetZ);
        position *= scale;

//        float distToLimits = (voxels_displayed_on_borders > 1 ? min(mincomp(abs(position.xyz - min_vertice_positions)), mincomp(abs(position.xyz + vec3(1.0) - max_vertice_positions))) : 1.0);
        Vector3 off = _offset * 1.f; // (clamp(distToLimits / float(voxels_displayed_on_borders), 0.0, 1.0));
        return position + off; //clamp(position + vec4(off, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
    //    return clamp (position + vec4(_offset, 0.0), vec4(min_vertice_positions, 1.0), vec4(max_vertice_positions, 1.0));
    };

    std::function getCubeIndex = [&](Vector3 voxPos, Vector3 normal) -> int {
        int cubeindex = 0;
        float cubeVal0 = cubeVali(voxPos, 0);
        float cubeVal1 = cubeVali(voxPos, 1);
        float cubeVal2 = cubeVali(voxPos, 2);
        float cubeVal3 = cubeVali(voxPos, 3);
        float cubeVal4 = cubeVali(voxPos, 4);
        float cubeVal5 = cubeVali(voxPos, 5);
        float cubeVal6 = cubeVali(voxPos, 6);
        float cubeVal7 = cubeVali(voxPos, 7);
        float refined_isolevel = isolevel + 0.0001;
        //Determine the index into the edge table which
        //tells us which vertices are inside of the surface
        cubeindex  = int(cubeVal0 < refined_isolevel);
        cubeindex += int(cubeVal1 < refined_isolevel)*2;
        cubeindex += int(cubeVal2 < refined_isolevel)*4;
        cubeindex += int(cubeVal3 < refined_isolevel)*8;
        cubeindex += int(cubeVal4 < refined_isolevel)*16;
        cubeindex += int(cubeVal5 < refined_isolevel)*32;
        cubeindex += int(cubeVal6 < refined_isolevel)*64;
        cubeindex += int(cubeVal7 < refined_isolevel)*128;

        normal = Vector3(0, 0, 0);

        if (cubeindex != 0 && cubeindex != 255) {
            Vector3 vertlist[12];

            //Find the vertices where the surface intersects the cube
            vertlist[0] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 1), cubeVal1);
            vertlist[1] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 2), cubeVal2);
            vertlist[2] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 3), cubeVal3);
            vertlist[3] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 0), cubeVal0);
            vertlist[4] = vertexInterp(refined_isolevel, cubePos(voxPos, 4), cubeVal4, cubePos(voxPos, 5), cubeVal5);
            vertlist[5] = vertexInterp(refined_isolevel, cubePos(voxPos, 5), cubeVal5, cubePos(voxPos, 6), cubeVal6);
            vertlist[6] = vertexInterp(refined_isolevel, cubePos(voxPos, 6), cubeVal6, cubePos(voxPos, 7), cubeVal7);
            vertlist[7] = vertexInterp(refined_isolevel, cubePos(voxPos, 7), cubeVal7, cubePos(voxPos, 4), cubeVal4);
            vertlist[8] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 4), cubeVal4);
            vertlist[9] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 5), cubeVal5);
            vertlist[10] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 6), cubeVal6);
            vertlist[11] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 7), cubeVal7);


//            vec3 edge1 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 1)];
//            vec3 edge2 = vertlist[triTableValue(cubeindex, 0)] - vertlist[triTableValue(cubeindex, 2)];
//            normal = normalize(cross(edge1, edge2));
        }
        return cubeindex;
    };

    Mesh marched;
    float refined_isolevel = isolevel + 0.0001;
    Vector3 vertlist[12];
    for (int x = 0; x < values.sizeX; x++) {
        for (int y = 0; y < values.sizeY; y++) {
            for (int z = 0; z < values.sizeZ; z++) {
                Vector3 position = Vector3(x, y, z);
                Vector3 voxPos = position;

                float cubeVal0 = cubeVali(voxPos, 0);
                float cubeVal1 = cubeVali(voxPos, 1);
                float cubeVal2 = cubeVali(voxPos, 2);
                float cubeVal3 = cubeVali(voxPos, 3);
                float cubeVal4 = cubeVali(voxPos, 4);
                float cubeVal5 = cubeVali(voxPos, 5);
                float cubeVal6 = cubeVali(voxPos, 6);
                float cubeVal7 = cubeVali(voxPos, 7);

                Vector3 normal;
                int cubeindex = getCubeIndex(voxPos, normal);

                //Cube is entirely in/out of the surface
                if (cubeindex == 0 || cubeindex == 255)
                    continue;

                //Find the vertices where the surface intersects the cube
                vertlist[0] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 1), cubeVal1);
                vertlist[1] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 2), cubeVal2);
                vertlist[2] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 3), cubeVal3);
                vertlist[3] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 0), cubeVal0);
                vertlist[4] = vertexInterp(refined_isolevel, cubePos(voxPos, 4), cubeVal4, cubePos(voxPos, 5), cubeVal5);
                vertlist[5] = vertexInterp(refined_isolevel, cubePos(voxPos, 5), cubeVal5, cubePos(voxPos, 6), cubeVal6);
                vertlist[6] = vertexInterp(refined_isolevel, cubePos(voxPos, 6), cubeVal6, cubePos(voxPos, 7), cubeVal7);
                vertlist[7] = vertexInterp(refined_isolevel, cubePos(voxPos, 7), cubeVal7, cubePos(voxPos, 4), cubeVal4);
                vertlist[8] = vertexInterp(refined_isolevel, cubePos(voxPos, 0), cubeVal0, cubePos(voxPos, 4), cubeVal4);
                vertlist[9] = vertexInterp(refined_isolevel, cubePos(voxPos, 1), cubeVal1, cubePos(voxPos, 5), cubeVal5);
                vertlist[10] = vertexInterp(refined_isolevel, cubePos(voxPos, 2), cubeVal2, cubePos(voxPos, 6), cubeVal6);
                vertlist[11] = vertexInterp(refined_isolevel, cubePos(voxPos, 3), cubeVal3, cubePos(voxPos, 7), cubeVal7);

                int i = 0;
                while(true){
                    if(triTableValue(cubeindex, i) != -1){
                        position = vertlist[triTableValue(cubeindex, i+0)];
                        marched.vertexArray.push_back(position);

                        position = vertlist[triTableValue(cubeindex, i+1)];
                        marched.vertexArray.push_back(position);

                        position = vertlist[triTableValue(cubeindex, i+2)];
                        marched.vertexArray.push_back(position);
                    }else{
                        break;
                    }

                    i = i + 3;
                }
            }
        }
    }
    return marched;
}
