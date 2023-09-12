#include "Graphics/Mesh.h"

#include "DataStructure/BVH.h"
#include "Graphics/MarchingCubes.h"
#include "Utils/stl_reader.h"
#include "Graphics/ofbx.h"

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


//Mesh &Mesh::fromArray(std::vector<std::vector<Vector3> > triangles, std::vector<int> indices)
//{
//    std::vector<Vector3> flattenVertices(triangles.size() * 3);
//    for (size_t i = 0; i < triangles.size(); i++) {
//        flattenVertices[3 * i + 0] = triangles[i][0];
//        flattenVertices[3 * i + 1] = triangles[i][1];
//        flattenVertices[3 * i + 2] = triangles[i][2];
//    }
//    return this->fromArray(flattenVertices, indices);
//}
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
        this->normalize();
    }
    catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    return *this;
}

Mesh &Mesh::fromFBX(std::string filename)
{
    this->clear();

    FILE* fp = fopen(filename.c_str(), "rb");

    if (!fp) return *this;

    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    auto* content = new ofbx::u8[file_size];
    fread(content, 1, file_size, fp);

    ofbx::LoadFlags flags =
            ofbx::LoadFlags::TRIANGULATE |
    //		ofbx::LoadFlags::IGNORE_MODELS |
            ofbx::LoadFlags::IGNORE_BLEND_SHAPES |
            ofbx::LoadFlags::IGNORE_CAMERAS |
            ofbx::LoadFlags::IGNORE_LIGHTS |
            ofbx::LoadFlags::IGNORE_TEXTURES |
            ofbx::LoadFlags::IGNORE_SKIN |
            ofbx::LoadFlags::IGNORE_BONES |
            ofbx::LoadFlags::IGNORE_PIVOTS |
            ofbx::LoadFlags::IGNORE_MATERIALS |
            ofbx::LoadFlags::IGNORE_POSES |
            ofbx::LoadFlags::IGNORE_VIDEOS |
            ofbx::LoadFlags::IGNORE_LIMBS |
    //		ofbx::LoadFlags::IGNORE_MESHES |
            ofbx::LoadFlags::IGNORE_ANIMATIONS;

    auto scene = ofbx::load((ofbx::u8*)content, file_size, (ofbx::u16)flags);
    if (scene->getGeometryCount() == 0)
        return *this;
//    for (int i = 0; i < scene->getGeometryCount(); i++) {
    int i = 0;
        auto geom = scene->getGeometry(i);
        auto verts = geom->getVertices();
        auto norms = geom->getNormals();
        int nbVertices = geom->getVertexCount();
        std::vector<Vector3> vertices(nbVertices);
        std::vector<Vector3> normals(nbVertices);
        for (int j = 0; j < nbVertices; j++) {
            vertices[j] = Vector3(verts[j].x, verts[j].y, verts[j].z);
            normals[j] = Vector3(norms[j].x, norms[j].y, norms[j].z);
        }
        this->normalsArray.insert(normalsArray.end(), normals.begin(), normals.end());
        this->vertexArray.insert(vertexArray.end(), vertices.begin(), vertices.end());
//    }
    this->fromArray(vertexArray);
    this->normalize();
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

Mesh& Mesh::scale(const Vector3& factor)
{
    for (auto& p : this->vertexArray)
        p *= factor;

    this->vertexArrayFloat = Vector3::toArray(this->vertexArray);
    this->needToUpdatePositions = true;

    return *this;
}

Mesh& Mesh::translate(const Vector3& translation)
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

Mesh& Mesh::rotate(const Vector3& rotation)
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

Mesh Mesh::extractGeometryFromShaders(GridF& values)
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
    Mesh copy = Mesh::applyMarchingCubes(values);
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

void Mesh::displayWithOutlines(std::vector<float> faceColor, GLenum shape, float lineWeight)
{
    this->shader->setVector("color", faceColor);
    this->display(shape, lineWeight);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    this->shader->setVector("color", std::vector{0.f, 0.f, 0.f, faceColor[3]});
    this->display(shape, lineWeight);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    this->shader->setVector("color", std::vector{0.f, 0.f, 0.f, faceColor[3]});
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

void Mesh::displayAsScalarField(GridF field, const Vector3& cameraPosition, std::vector<float> isoValues)
{
    std::vector<Vector3> positions(field.size());
    for (size_t i = 0; i < positions.size(); i++) {
        positions[i] = field.getCoordAsVector3(i);
    }
    this->fromArray(positions);
    this->update();
    GlobalsGL::f()->glBindVertexArray(this->vao);
    this->shader->setTexture3D("dataFieldTex", 0, field + .5f);
    this->shader->setInt("dataFieldTex", 0);
    this->shader->setInt("edgeTableTex", 1);
    this->shader->setInt("triTableTex", 2);
    this->shader->setFloat("isolevel", 0.f);
    this->shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    this->shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    this->shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    this->shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    this->shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    this->shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    this->shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    this->shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
    this->shader->setBool("useMarchingCubes", true);
    //Edge Table texture//
    //This texture store the 256 different configurations of a marching cube.
    //This is a table accessed with a bitfield of the 8 cube edges states
    //(edge cut by isosurface or totally in or out).
    //(cf. MarchingCubes.cpp)

    GLuint edgeTableTex, triTableTex;
    GlobalsGL::f()->glGenTextures(1, &edgeTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, edgeTableTex);
    //Integer textures must use nearest filtering mode

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 256, 1, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::cubeEdges));

    //Triangle Table texture//
    //This texture store the vertex index list for
    //generating the triangles of each configurations.
    //(cf. MarchingCubes.cpp)

    glGenTextures(1, &triTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE2);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, triTableTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 16, 256, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::triangleTable));


    GlobalsGL::f()->glActiveTexture(GL_TEXTURE0);

    // Ignore parameters to hide some voxels
//    this->shader->setVector("min_vertice_positions", Vector3::min());
//    this->shader->setVector("max_vertice_positions", Vector3::max());
//    this->shader->setFloat("min_isolevel", -1000.f); // 3.5f);
//    this->shader->setFloat("max_isolevel", 1000.f);

    for (size_t i = 0; i < isoValues.size(); i++) {
        float iso = isoValues[i];
        Vector3 color = HSVtoRGB(i / float(isoValues.size()), 1.f, 1.f);
        this->shader->setVector("color", std::vector<float> {color.x, color.y, color.z, .3f});
        this->shader->setFloat("isolevel", iso);

        // display the mesh
        this->reorderVertices(cameraPosition);
        this->display(GL_POINTS);
    }
}

void Mesh::displayAsVectorField(GridV3 field, const Vector3& finalDimensions, float maxMaginitude, bool normalize)
{
    if (maxMaginitude > 0.f) {
        for (auto& v : field)
            v.maxMagnitude(maxMaginitude);
    }
    if (normalize) {
        field.normalize();
    }
    std::vector<Vector3> normals;
    for (int x = 0; x < field.sizeX - 1; x++) {
        for (int y = 0; y < field.sizeY - 1; y++) {
            for (int z = 0; z < field.sizeZ - 1; z++) {
                normals.push_back(Vector3(x, y, z) + Vector3(.5f, .5f, .5f));
                normals.push_back(Vector3(x, y, z) + field.at(x, y, z) + Vector3(.5f, .5f, .5f));
            }
        }
    }
    if (finalDimensions.isValid()) {
        Vector3 ratio = finalDimensions / field.getDimensions();
        for (auto& n : normals) {
            n *= ratio;
        }
    }
    this->fromArray(normals);
    this->display(GL_LINES);
}

void Mesh::setShader(std::shared_ptr<Shader> shader)
{

}

void Mesh::reorderVertices(const Vector3& camPos)
{
    this->reorderAny(camPos, 1);
}

void Mesh::reorderLines(const Vector3& camPos)
{
    this->reorderAny(camPos, 2);
}

void Mesh::reorderTriangles(const Vector3& camPos)
{
    this->reorderAny(camPos, 3);
}

void Mesh::reorderAny(const Vector3& camPos, int nbVertexToUse)
{
    std::vector<Vector3> lineCenters(vertexArray.size() / nbVertexToUse);
    std::vector<int> newOrder(lineCenters.size());
    std::vector<float> distances(lineCenters.size());
    for (size_t i = 0; i < lineCenters.size(); i++) {
        for (int n = 0; n < nbVertexToUse; n++)
            lineCenters[i] += vertexArray[nbVertexToUse * i + n];
        lineCenters[i] = (lineCenters[i] / (float)(nbVertexToUse)) - camPos;
        distances[i] = lineCenters[i].norm2();
        newOrder[i] = i;
    }
    std::sort(newOrder.begin(), newOrder.end(), [&](const int& idA, const int& idB) { return distances[idA] > distances[idB]; });

    for (size_t i = 0; i < lineCenters.size(); i++) {
        for (int iVertex = 0; iVertex < nbVertexToUse; iVertex++) {
            size_t vectorIndex = nbVertexToUse * newOrder[i] + iVertex;
            size_t flatIndex = 3 * (nbVertexToUse * i + iVertex);
            if (vertexArray.size() > vectorIndex) {
                vertexArrayFloat[flatIndex + 0] = vertexArray[vectorIndex].x;
                vertexArrayFloat[flatIndex + 1] = vertexArray[vectorIndex].y;
                vertexArrayFloat[flatIndex + 2] = vertexArray[vectorIndex].z;
            }
            if (normalsArray.size() > vectorIndex) {
                normalsArrayFloat[flatIndex + 0] = normalsArray[vectorIndex].x;
                normalsArrayFloat[flatIndex + 1] = normalsArray[vectorIndex].y;
                normalsArrayFloat[flatIndex + 2] = normalsArray[vectorIndex].z;
            }
            if (colorsArray.size() > vectorIndex) {
                colorsArrayFloat[flatIndex + 0] = colorsArray[vectorIndex].x;
                colorsArrayFloat[flatIndex + 1] = colorsArray[vectorIndex].y;
                colorsArrayFloat[flatIndex + 2] = colorsArray[vectorIndex].z;
            }
        }
    }
    // Retrieve the new order as Vector3's
    for (size_t i = 0; i < lineCenters.size(); i++) {
        for (int iVertex = 0; iVertex < nbVertexToUse; iVertex++) {
            size_t vectorIndex = nbVertexToUse * i + iVertex;
            size_t flatIndex = 3 * (nbVertexToUse * i + iVertex);
            if (vertexArray.size() > vectorIndex) {
                vertexArray[vectorIndex].x = vertexArrayFloat[flatIndex + 0];
                vertexArray[vectorIndex].y = vertexArrayFloat[flatIndex + 1];
                vertexArray[vectorIndex].z = vertexArrayFloat[flatIndex + 2];
            }
            if (normalsArray.size() > vectorIndex) {
                normalsArray[vectorIndex].x = normalsArrayFloat[flatIndex + 0];
                normalsArray[vectorIndex].y = normalsArrayFloat[flatIndex + 1];
                normalsArray[vectorIndex].z = normalsArrayFloat[flatIndex + 2];
            }
            if (colorsArray.size() > vectorIndex) {
                colorsArray[vectorIndex].x = colorsArrayFloat[flatIndex + 0];
                colorsArray[vectorIndex].y = colorsArrayFloat[flatIndex + 1];
                colorsArray[vectorIndex].z = colorsArrayFloat[flatIndex + 2];
            }
        }
        /*for (int j = 0; j < 3 * nbVertexToUse; j++) {
            if (vertexArray.size() > )
            vertexArray[nbVertexToUse * i + j / 3][j % 3] = vertexArrayFloat[nbVertexToUse * 3 * i + j];

            normalsArray[nbVertexToUse * i + j / 3][j % 3] = normalsArrayFloat[nbVertexToUse * 3 * i + j];

            colorsArray[nbVertexToUse * i + j / 3][j % 3] = colorsArrayFloat[nbVertexToUse * 3 * i + j];
        }*/
    }
    this->needToUpdatePositions = true;
    this->needToUpdateNormals = true;
    this->needToUpdateColors = true;
}

std::vector<std::vector<Vector3> > Mesh::getTriangles(std::vector<int> indices) const
{
    if (indices.empty()) {
        indices = std::vector<int>(this->vertexArray.size());
        for (size_t i = 0; i < indices.size(); i++)
            indices[i] = i;
    }

    std::vector<std::vector<Vector3>> triangles(indices.size() / 3);
    for (size_t i = 0; i < indices.size(); i += 3) {
        triangles[i / 3] = std::vector<Vector3>({
                                                     this->vertexArray[indices[i+0]],
                                                     this->vertexArray[indices[i+1]],
                                                     this->vertexArray[indices[i+2]],
                                                 });
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

Mesh Mesh::applyMarchingCubes(const GridF& values)
{
    auto triTable = MarchingCubes::triangleTable;
    auto edges = MarchingCubes::edgeToCorner;
    Vector3 scale(1.f, 1.f, 1.f);
    float isolevel = 0.f;
    float refined_isolevel = isolevel + 0.0001;

    Vector3 vertDecals[8] = {
/*        Vector3(0.f, 0.f, 0.f),
        Vector3(1.f, 0.f, 0.f),
        Vector3(1.f, 1.f, 0.f),
        Vector3(0.f, 1.f, 0.f),
        Vector3(0.f, 0.f, 1.f),
        Vector3(1.f, 0.f, 1.f),
        Vector3(1.f, 1.f, 1.f),
        Vector3(0.f, 1.f, 1.f)
*/

        Vector3(-.5f, -.5f, -.5f),
        Vector3(0.5f, -.5f, -.5f),
        Vector3(0.5f, 0.5f, -.5f),
        Vector3(-.5f, 0.5f, -.5f),
        Vector3(-.5f, -.5f, 0.5f),
        Vector3(0.5f, -.5f, 0.5f),
        Vector3(0.5f, 0.5f, 0.5f),
        Vector3(-.5f, 0.5f, 0.5f)
    };

    Mesh marched;
    marched.vertexArray.reserve(100000);

//    #pragma omp parallel for
    for (int z = 0; z < values.sizeZ; z++) {
        for (int y = 0; y < values.sizeY; y++) {
            for (int x = 0; x < values.sizeX; x++) {
                Vector3 position = Vector3(x, y, z);
                Vector3 voxPos = position;

                std::vector<float> cubeVals(8);
                for (int i = 0; i < 8; i++)
                    cubeVals[i] = values.checkCoord(voxPos + vertDecals[i]) && (voxPos + vertDecals[i]).minComp() > 0 ? values.at(voxPos + vertDecals[i]) : -100000.f;

                int cubeindex = 0;
                for (int i = 0; i < 8; i++)
                    cubeindex += int(cubeVals[i] < refined_isolevel) << i;

                if (cubeindex == 0 || cubeindex == 255)
                    continue;

                Vector3 vertlist[12];
                for (int iEdge = 0; iEdge < 12; iEdge++)  {
                    int p0 = edges[iEdge][0];
                    int p1 = edges[iEdge][1];
                    float l0 = cubeVals[p0];
                    float l1 = cubeVals[p1];
                    float epsilon = 1e-3;
                    if (std::abs(l0 - l1) < epsilon) {
                        vertlist[iEdge] = (voxPos + vertDecals[p0] + voxPos + vertDecals[p1]) * .5f;
                    } else {
                        float iso = 1.f - std::clamp((isolevel-l0)/(l1-l0), 0.f, 1.f);
                        vertlist[iEdge] = (voxPos + vertDecals[p0]) * iso + (voxPos + vertDecals[p1]) * (1.f - iso);
                    }
                }

                int i = 0;
                while(triTable[cubeindex][i] != -1){
//                    #pragma omp critical
                    {
                        marched.vertexArray.push_back(vertlist[triTable[cubeindex][i]]);
                        marched.vertexArray.push_back(vertlist[triTable[cubeindex][i+1]]);
                        marched.vertexArray.push_back(vertlist[triTable[cubeindex][i+2]]);
                    }
                    i += 3;
                }
            }
        }
    }

    marched.scale(values.getDimensions() / values.getDimensions()); //.translate(-Vector3(.5, .5, .5));
//    marched.scale((values.getDimensions() + Vector3(1, 1, 1)) / values.getDimensions()).translate(-Vector3(.5, .5, .5));
    return marched;
}


GridI Mesh::voxelize(const Vector3& dimensions) const
{
    AABBox myDims(this->vertexArray);
    GridI res(dimensions);

    BVHTree tree;
    auto triangles = this->getTriangles();
    std::cout << "Mesh dimensions: " << myDims << std::endl;
    for (auto& tri : triangles) {
        for (auto& p : tri) {
            p = (p - myDims.mini) / myDims.dimensions() /** Vector3(1.5, 1.5, 4.5)*/ * dimensions;
        }
    }
    tree.build(Triangle::vectorsToTriangles(triangles));

    int dimX = dimensions.x;
    int dimY = dimensions.y;
    int dimZ = dimensions.z;

#pragma omp parallel for collapse(3)
    for (int x = 0; x < dimX; x++) {
        for (int y = 0; y < dimY; y++) {
            for (int z = 0; z < dimZ; z++) {
                Vector3 pos = Vector3(x, y, z); // /(dimensions / myDims.dimensions()) + myDims.min();
                Vector3 ray = Vector3(pos.x + 1, pos.y + 3, pos.z + 500); // Vector3(pos.x, pos.y, pos.z + myDims.dimensions().z + 1);
//                res.at(x, y, z) = (tree.checkIntersection(pos, ray) ? 1 : 0);
//                res.at(x, y, z) = (tree.getIntersection(pos, ray).isValid() ? 1 : 0);
                res.at(pos) = (tree.getAllIntersections(pos, ray).size() % 2 == 0 ? 0.f : 1.f);
            }
        }
    }
    return res;
}

GridI Mesh::voxelizeSurface(const Vector3& dimensions) const
{
    AABBox myDims(this->vertexArray);
    GridI res(dimensions, -1.f);
    auto triangles = this->getTriangles();

    for (auto& t : triangles)
        for (auto& p : t)
            p = (p + Vector3(.5f, .5f, .5f)) * dimensions;

    BVHTree tree;
    tree.build(Triangle::vectorsToTriangles(triangles));

    int dimX = dimensions.x;
    int dimY = dimensions.y;
    int dimZ = dimensions.z;

#pragma omp parallel for collapse(3)
    for (int x = 0; x < dimX; x++) {
        for (int y = 0; y < dimY; y++) {
            for (int z = 0; z < dimZ; z++) {
                Vector3 pos = Vector3(x, y, z);
                if (tree.getIntersection(pos, pos + Vector3(1, 1, 1)).isValid())
                    res.at(pos) = 1.f;
            }
        }
    }
    return res;
}

bool Mesh::isWatertight()
{
    bool previousVal = this->useIndices;
    this->useIndices = true;
    this->computeIndices();
    this->useIndices = previousVal;

//    auto triangles = this->getTriangles();
    std::map<std::pair<int, int>, int> edgesCount;
    for (size_t i = 0; i < this->vertexArray.size(); i += 3) {
        std::vector<int> ids({indices[i], indices[i + 1], indices[i + 2]});
        std::sort(ids.begin(), ids.end());
        if (edgesCount.count({ids[0], ids[1]}) == 0)
            edgesCount[{ids[0], ids[1]}] = 0;
        edgesCount[{ids[0], ids[1]}] ++;
        if (edgesCount.count({ids[0], ids[2]}) == 0)
            edgesCount[{ids[0], ids[2]}] = 0;
        edgesCount[{ids[0], ids[2]}] ++;
        if (edgesCount.count({ids[1], ids[2]}) == 0)
            edgesCount[{ids[1], ids[2]}] = 0;
        edgesCount[{ids[1], ids[2]}] ++;
    }

    for (auto& edge : edgesCount) {
        int count = edge.second;
        if (count != 2)
            return false;
    }
    return true;


}

Mesh Mesh::createVectorField(GridV3 field, const Vector3& finalDimensions, Mesh* mesh, float maxMaginitude, bool normalize, bool displayArrow)
{
    if (maxMaginitude > 0.f) {
        for (auto& v : field)
            v.maxMagnitude(maxMaginitude);
    }
    if (normalize) {
        field.normalize();
    }

    Vector3 offsetToCenter(.5f, .5f, .5f);
    std::vector<Vector3> normals;
    for (int x = 0; x < field.sizeX; x++) {
        for (int y = 0; y < field.sizeY; y++) {
            for (int z = 0; z < field.sizeZ; z++) {
                Vector3 pos(x, y, z);
                Vector3 value = field.at(x, y, z);
                normals.push_back(pos + offsetToCenter);
                normals.push_back(pos + offsetToCenter + value);

                if (displayArrow) {
                    // Arrow part
                    Vector3 cross = value.cross(Vector3(0, 0, 1)).setMag(value.norm() * .2f);
                    normals.push_back(pos + offsetToCenter + value * .8f - cross);
                    normals.push_back(pos + offsetToCenter + value);
                    normals.push_back(pos + offsetToCenter + value);
                    normals.push_back(pos + offsetToCenter + value * .8f + cross);
                }
            }
        }
    }
    if (finalDimensions.isValid()) {
        Vector3 ratio = finalDimensions / field.getDimensions();
        for (auto& n : normals) {
            n *= ratio;
        }
    }
    if (mesh == nullptr)
        return Mesh(normals, nullptr, true, GL_LINES);
    else {
        return mesh->fromArray(normals);
    }
}

void Mesh::displayScalarField(GridF field, Mesh &mesh, const Vector3& cameraPosition, std::vector<float> isoValues)
{
    std::vector<Vector3> positions(field.size());
    for (size_t i = 0; i < positions.size(); i++) {
        positions[i] = field.getCoordAsVector3(i);
    }
    mesh.fromArray(positions);
    mesh.update();
    GlobalsGL::f()->glBindVertexArray(mesh.vao);
    mesh.shader->setTexture3D("dataFieldTex", 0, field + .5f);
    mesh.shader->setInt("dataFieldTex", 0);
    mesh.shader->setInt("edgeTableTex", 1);
    mesh.shader->setInt("triTableTex", 2);
    mesh.shader->setFloat("isolevel", 0.f);
    mesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    mesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    mesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    mesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    mesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    mesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    mesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    mesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
    mesh.shader->setBool("useMarchingCubes", true);
    //Edge Table texture//
    //This texture store the 256 different configurations of a marching cube.
    //This is a table accessed with a bitfield of the 8 cube edges states
    //(edge cut by isosurface or totally in or out).
    //(cf. MarchingCubes.cpp)

    GLuint edgeTableTex, triTableTex;
    GlobalsGL::f()->glGenTextures(1, &edgeTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, edgeTableTex);
    //Integer textures must use nearest filtering mode

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 256, 1, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::cubeEdges));

    //Triangle Table texture//
    //This texture store the vertex index list for
    //generating the triangles of each configurations.
    //(cf. MarchingCubes.cpp)

    glGenTextures(1, &triTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE2);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, triTableTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 16, 256, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::triangleTable));


    GlobalsGL::f()->glActiveTexture(GL_TEXTURE0);

    // Ignore parameters to hide some voxels
//    mesh.shader->setVector("min_vertice_positions", Vector3::min());
//    mesh.shader->setVector("max_vertice_positions", Vector3::max());
//    mesh.shader->setFloat("min_isolevel", -1000.f); // 3.5f);
//    mesh.shader->setFloat("max_isolevel", 1000.f);

    std::sort(isoValues.begin(), isoValues.end(), [&](float A, float B) { return A > B; });
    for (size_t i = 0; i < isoValues.size(); i++) {
        float iso = isoValues[i];
        Vector3 color = HSVtoRGB(i / float(isoValues.size()), 1.f, 1.f);
        mesh.shader->setVector("color", std::vector<float> {color.x, color.y, color.z, .3f});
        mesh.shader->setFloat("isolevel", iso);

        // display the mesh
        mesh.reorderVertices(cameraPosition);
        mesh.display(GL_POINTS);
    }
}
