#include "Graphics/Shader.h"

#include <QFile>
#include <QTextStream>

std::shared_ptr<Shader> Shader::default_shader = nullptr;
std::set<std::shared_ptr<Shader>> Shader::allShaders;

Shader::Shader()
    : Shader("", "", "")
{

}
Shader::Shader(const std::string& vertexShaderFilename)
    : Shader(vertexShaderFilename, "", "")
{

}
Shader::Shader(const std::string& vertexShaderFilename, std::string fragmentShaderFilename)
    : Shader(vertexShaderFilename, fragmentShaderFilename, "")
{

}
Shader::Shader(const std::string& vertexShaderFilename, std::string fragmentShaderFilename,
       std::string geometryShaderFilename)
    : vertexShaderFilename(vertexShaderFilename), fragmentShaderFilename(fragmentShaderFilename),
      geometryShaderFilename(geometryShaderFilename)
{
    // GPU compilation is deferred until initializeGL()/use().
}
Shader::Shader(const Shader& copy) {
    this->vertexShaderFilename = copy.vertexShaderFilename;
    this->fragmentShaderFilename = copy.fragmentShaderFilename;
    this->geometryShaderFilename = copy.geometryShaderFilename;
    this->definitions = copy.definitions;
    // GPU compilation is deferred until initializeGL()/use().
}
/*
Shader::Shader()
    : Shader(nullptr, nullptr, nullptr)
{

}
Shader::Shader(const char* vertexShaderFilename)
    : Shader(vertexShaderFilename, nullptr, nullptr)
{

}
Shader::Shader(const char* vertexShaderFilename, const char* fragmentShaderFilename)
    : Shader(vertexShaderFilename, fragmentShaderFilename, nullptr)
{

}
Shader::Shader(const char* vertexShaderFilename, const char* fragmentShaderFilename,
       const char* geometryShaderFilename)
    : vertexShaderFilename(vertexShaderFilename), fragmentShaderFilename(fragmentShaderFilename),
      geometryShaderFilename(geometryShaderFilename)
{
    // GPU compilation is deferred until initializeGL()/use().
}

Shader::Shader(const std::string& vertexShaderFilename)
    : Shader(vertexShaderFilename.c_str())
{

}

Shader::Shader(const std::string& vertexShaderFilename, std::string fragmentShaderFilename)
    : Shader(vertexShaderFilename.c_str(), fragmentShaderFilename.c_str())
{

}

Shader::Shader(const std::string& vertexShaderFilename, std::string fragmentShaderFilename, std::string geometryShaderFilename)
    : Shader(vertexShaderFilename.c_str(), fragmentShaderFilename.c_str(), geometryShaderFilename.c_str())
{

}*/
Shader& Shader::operator=(const Shader& copy) {
    if (this == &copy) return *this;
    gpuResources.reset();
    programID = 0;
    vShader = fShader = gShader = -1;
    textureSlotIndices.clear();
    volumeTextures.clear();
    invalidateProgramCache();
    vertexShaderFilename = copy.vertexShaderFilename;
    fragmentShaderFilename = copy.fragmentShaderFilename;
    geometryShaderFilename = copy.geometryShaderFilename;
    definitions = copy.definitions;
    return *this;
}

void Shader::ensureResources() {
    auto* ctx = GlobalsGL::context();
    if (!gpuResources || gpuResources->released || gpuResources->owner != ctx) {
        gpuResources = GlobalsGL::Resources::create();
        programID = 0;
        vShader = fShader = gShader = -1;
        textureSlotIndices.clear();
        volumeTextures.clear();
        invalidateProgramCache();
    }
}

void Shader::compileShadersFromSource(std::map<std::string, std::string> addedDefinitions, std::vector<std::string> feedbackValues)
{
#if useModernOpenGL || !useModernOpenGL
    GlobalsGL::ContextScope contextScope;
    if (addedDefinitions.empty()) addedDefinitions = definitions;
    definitions = addedDefinitions;
    ensureResources();
    vShader = fShader = gShader = -1;
    if (this->programID != 0) {
        GlobalsGL::f()->glDeleteProgram(this->programID);
        gpuResources->programs.clear();
    }
    this->programID = GlobalsGL::f()->glCreateProgram();
    gpuResources->programs.push_back(programID);
    if (vertexShaderFilename != "")
    {
        std::string content = Shader::addDefinitionsToSource(Shader::readShaderSource(vertexShaderFilename), addedDefinitions);
        if (!content.empty()) {
            this->vShader = GlobalsGL::f()->glCreateShader(GL_VERTEX_SHADER);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(this->vShader, 1, &src, NULL);
            GlobalsGL::f()->glCompileShader(this->vShader);
            GlobalsGL::f()->glAttachShader(this->programID, this->vShader);
            GlobalsGL::printShaderErrors(this->vShader);
            GlobalsGL::checkOpenGLError();
        } else {
            vertexShaderFilename = "";
        }
    }
    if (fragmentShaderFilename != "")
    {
        std::string content = Shader::addDefinitionsToSource(Shader::readShaderSource(fragmentShaderFilename), addedDefinitions);
        if (!content.empty()) {
            this->fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(this->fShader, 1, &src, NULL);
            GlobalsGL::f()->glCompileShader(this->fShader);
            GlobalsGL::f()->glAttachShader(this->programID, this->fShader);
            GlobalsGL::printShaderErrors(this->fShader);
            GlobalsGL::checkOpenGLError();
        } else {
            fragmentShaderFilename = "";
        }
    }
    if (geometryShaderFilename != "")
    {
        std::string content = Shader::addDefinitionsToSource(Shader::readShaderSource(geometryShaderFilename), addedDefinitions);
        if (!content.empty()) {
            this->gShader = GlobalsGL::f()->glCreateShader(GL_GEOMETRY_SHADER);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(this->gShader, 1, &src, NULL);
            GlobalsGL::f()->glCompileShader(this->gShader);
            GlobalsGL::f()->glAttachShader(this->programID, this->gShader);
            GlobalsGL::printShaderErrors(this->gShader);
            GlobalsGL::checkOpenGLError();
        } else {
            geometryShaderFilename = "";
        }
    }

    if (!feedbackValues.empty())
    {
        std::vector<const GLchar*> variables(feedbackValues.size());
        for (size_t i = 0; i < feedbackValues.size(); i++)
            variables[i] = (GLchar*)feedbackValues[i].c_str();
        GlobalsGL::f()->glTransformFeedbackVaryings(this->programID, feedbackValues.size(), variables.data(), GL_INTERLEAVED_ATTRIBS);
    }

    GlobalsGL::checkOpenGLError();
    GlobalsGL::f()->glLinkProgram(this->programID);
    invalidateProgramCache();
    GLint linked = GL_FALSE;
    GlobalsGL::f()->glGetProgramiv(programID, GL_LINK_STATUS, &linked);
    programLinked = linked == GL_TRUE;
    if (this->vShader > -1)
        GlobalsGL::f()->glDeleteShader(this->vShader);
    if (this->fShader > -1)
        GlobalsGL::f()->glDeleteShader(this->fShader);
    if (this->gShader > -1)
        GlobalsGL::f()->glDeleteShader(this->gShader);
    vShader = fShader = gShader = -1;
    GlobalsGL::printProgramErrors(programID);
    GlobalsGL::checkOpenGLError();

#endif
}

void Shader::invalidateProgramCache()
{
    programLinked = false;
    uniformLocations.clear();
    uniformValues.clear();
}

bool Shader::ensureLinked()
{
    ensureResources();
    if (programID == 0) initializeGL();
    return programLinked;
}

GLint Shader::uniformLocation(const std::string& name)
{
    const auto found = uniformLocations.find(name);
    if (found != uniformLocations.end()) return found->second;
    const GLint location = GlobalsGL::f()->glGetUniformLocation(programID, name.c_str());
    uniformLocations.emplace(name, location);
    return location;
}

bool Shader::uniformChanged(GLint location, const void* data, size_t bytes)
{
    if (location < 0) return false;
    auto& previous = uniformValues[location];
    if (previous.size() == bytes && std::memcmp(previous.data(), data, bytes) == 0) return false;
    const auto* begin = static_cast<const unsigned char*>(data);
    previous.assign(begin, begin + bytes);
    return true;
}

bool Shader::use(bool update_source_file)
{
    GlobalsGL::ContextScope contextScope;
    if (update_source_file) initializeGL();
    if (!ensureLinked()) return false;
    GlobalsGL::f()->glUseProgram(programID);
    return true;
}

void Shader::setBool(const std::string& pname, bool value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const int integer = int(value);
    const GLint location = uniformLocation(pname);
    if (uniformChanged(location, &integer, sizeof(integer)))
        GlobalsGL::f()->glProgramUniform1i(programID, location, integer);
}
void Shader::setInt(const std::string& pname, int value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const GLint location = uniformLocation(pname);
    if (uniformChanged(location, &value, sizeof(value)))
        GlobalsGL::f()->glProgramUniform1i(programID, location, value);
}

void Shader::setFloat(const std::string& pname, float value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const GLint location = uniformLocation(pname);
    if (uniformChanged(location, &value, sizeof(value)))
        GlobalsGL::f()->glProgramUniform1f(programID, location, value);
}

void Shader::setVector(const std::string& pname, const Vector3& value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    this->setVector(pname, value.data(), 3);
}

void Shader::setVector(const std::string& pname, const Vector4 &value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const float components[] = {value.x(), value.y(), value.z(), value.w};
    this->setVector(pname, components, 4);
}

void Shader::setVector(const std::string& pname, glm::vec2 value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const GLint location = uniformLocation(pname);
    if (uniformChanged(location, &value[0], sizeof(float) * 2))
        GlobalsGL::f()->glProgramUniform2fv(programID, location, 1, &value[0]);
}

//void Shader::setVector(const std::string& pname, glm::vec3 value)
//{
//    if (!this->use()) return;
//    GlobalsGL::f()->glUniform3fv(uniformLocation(pname), 1, &value[0]);
//}

/*void Shader::setVector(const std::string& pname, glm::vec4 value)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    GlobalsGL::f()->glUniform4fv(uniformLocation(pname), 1, &value[0]);
}*/
void Shader::setVector(const std::string& pname, const float* value, int n)
{
    GlobalsGL::ContextScope contextScope;
    if (!this->use()) return;
    const GLint loc = uniformLocation(pname);
    if (!uniformChanged(loc, value, sizeof(float) * n)) return;
    switch (n) {
    case 1:
        GlobalsGL::f()->glProgramUniform1fv(programID, loc, 1, value);
        break;
    case 2:
        GlobalsGL::f()->glProgramUniform2fv(programID, loc, 1, value);
        break;
    case 3:
        GlobalsGL::f()->glProgramUniform3fv(programID, loc, 1, value);
        break;
    case 4:
        GlobalsGL::f()->glProgramUniform4fv(programID, loc, 1, value);
        break;
    }
}
void Shader::setLightSource(const std::string& pname, LightSource &value)
{
    this->setVector((pname + ".ambiant").c_str(), value.ambiant);
    this->setVector((pname + ".diffuse").c_str(), value.diffuse);
    this->setVector((pname + ".specular").c_str(), value.specular);
}

void Shader::addLightSource(const std::string& pname, LightSource &value)
{
    this->setLightSource(pname + "[" + std::to_string(lightCount) +"]", value);
    lightCount ++;
    this->setInt(pname + "_count", this->lightCount);
}

void Shader::clearLightSources(const std::string& pname)
{
    this->lightCount = 0;
    this->setInt(pname + "_count", 0);
}

void Shader::setPositionalLight(const std::string& pname, PositionalLight &value)
{
    this->setLightSource(pname, value);
    this->setVector((pname + ".position").c_str(), value.position);
    this->setBool(pname + ".castsShadows", value.castsShadows);
}

void Shader::setMaterial(const std::string& pname, Material &value)
{
    this->setVector((pname + ".ambiant").c_str(), value.ambiant);
    this->setVector((pname + ".diffuse").c_str(), value.diffuse);
    this->setVector((pname + ".specular").c_str(), value.specular);
    this->setFloat((pname + ".shininness").c_str(), value.shininess);
}

// Todo : change the integer type to a template (or select int or float)
void Shader::setTexture2D(const std::string& pname, int index, GridI texture)
{
    setTexture2D(pname, index, texture.sizeX, texture.sizeY, texture.data.data());
}

void Shader::setTexture2D(const std::string& pname, int index, int width, int height, int *data)
{
    GlobalsGL::ContextScope contextScope;
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    bool justUpdateTexture = true;
    if (!this->use()) return;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        gpuResources->textures.push_back(texIndex);
        textureSlotIndices[textureSlot] = texIndex;
        justUpdateTexture = false;
    }
    texIndex = textureSlotIndices[textureSlot];
    GLint previousTextureUnit = GL_TEXTURE0;
    glGetIntegerv(GL_ACTIVE_TEXTURE, &previousTextureUnit);
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glBindTexture(GL_TEXTURE_2D, texIndex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    if (justUpdateTexture) {
        GlobalsGL::f()->glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height,
        GL_ALPHA_INTEGER_EXT, GL_INT, data);
    } else {
        glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, width, height, 0,
        GL_ALPHA_INTEGER_EXT, GL_INT, data);
    }
    this->setInt(pname, index);
    GlobalsGL::f()->glActiveTexture(previousTextureUnit);
}

void Shader::setTexture2D(const std::string& pname, int index, int width, int height, int **data)
{
    // OpenGL requires contiguous pixels, not an array of pointers to rows.
    std::vector<int> pixels(size_t(width) * size_t(height));
    for (int y = 0; y < height; ++y)
        std::copy(data[y], data[y] + width, pixels.begin() + size_t(y) * width);
    setTexture2D(pname, index, width, height, pixels.data());
}

// Todo : change the integer type to a template (or select int or float)
void Shader::setTexture3D(const std::string& pname, int index, const GridF& texture)
{
    GlobalsGL::ContextScope contextScope;
    if (!use() || texture.empty()) return;
    const int textureSlot = GL_TEXTURE0 + index;
    auto& cached = volumeTextures[textureSlot];
    const bool allocate = cached.dimensions != texture.getDimensions();
    const bool upload = allocate || std::memcmp(cached.samples.data.data(), texture.data.data(), texture.size() * sizeof(float)) != 0;
    GLuint texIndex;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        gpuResources->textures.push_back(texIndex);
        textureSlotIndices[textureSlot] = texIndex;
    }
    texIndex = textureSlotIndices.at(textureSlot);
    GLint previousTextureUnit = GL_TEXTURE0;
    glGetIntegerv(GL_ACTIVE_TEXTURE, &previousTextureUnit);
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glBindTexture(GL_TEXTURE_3D, texIndex);
    if (allocate) {
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        GlobalsGL::f()->glTexImage3D(GL_TEXTURE_3D, 0, GL_ALPHA32F_ARB,
            texture.sizeX, texture.sizeY, texture.sizeZ, 0, GL_ALPHA, GL_FLOAT, texture.data.data());
        ++volumeAllocationCount;
    } else if (upload) {
        GlobalsGL::f()->glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0,
            texture.sizeX, texture.sizeY, texture.sizeZ, GL_ALPHA, GL_FLOAT, texture.data.data());
    }
    if (upload) {
        cached.dimensions = texture.getDimensions();
        cached.samples = texture;
        ++volumeUploadCount;
    }
    setInt(pname, index);
    GlobalsGL::f()->glActiveTexture(previousTextureUnit);
}

//void Shader::setMatrix(const std::string& pname, Matrix value)
//{
//    this->setMatrix(pname, value.)
//}


void Shader::applyToAllShaders(std::function<void (std::shared_ptr<Shader>)> func)
{
    GlobalsGL::ContextScope contextScope;
    for (auto it = allShaders.begin(); it != allShaders.end();) {
        // The registry must not keep shaders from discarded meshes alive forever.
        if (!*it || it->use_count() == 1) {
            it = allShaders.erase(it);
        } else {
            auto shader = *it++;
            func(shader);
        }
    }
}


std::string Shader::readShaderSource(const std::string& filename)
{
    std::string content = "";
    QString qFilename = QString::fromStdString(filename);
    if (!QFile::exists(qFilename))
        qFilename = ":" + qFilename;
    if (!QFile::exists(qFilename)) {
        std::cerr << "The shader " << filename << " doesn't exist!" << std::endl;
        return "";
    }
    QFile file(qFilename);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    std::string line;
    QTextStream in(&file);
    while (!in.atEnd()) {
        line = in.readLine().toStdString();
        content += line + " \n";
    }
    file.close();
    return content;
}

std::string Shader::addDefinitionsToSource(std::string src, std::map<std::string, std::string> newDefinitions)
{
    auto start = src.find("#version");
    if (start == std::string::npos)
        return src;

    start = src.find("\n", start) + 1;
    for (auto& [name, val] : newDefinitions) {
        src.insert(start, "#define " + name + " " + val + "\n");
    }
    return src;
}
