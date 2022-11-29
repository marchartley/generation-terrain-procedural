#include "Graphics/Shader.h"

#include <QFile>
#include <QTextStream>

std::shared_ptr<Shader> Shader::default_shader = nullptr;
std::set<std::shared_ptr<Shader>> Shader::allShaders;

Shader::Shader()
    : Shader("", "", "")
{

}
Shader::Shader(std::string vertexShaderFilename)
    : Shader(vertexShaderFilename, "", "")
{

}
Shader::Shader(std::string vertexShaderFilename, std::string fragmentShaderFilename)
    : Shader(vertexShaderFilename, fragmentShaderFilename, "")
{

}
Shader::Shader(std::string vertexShaderFilename, std::string fragmentShaderFilename,
       std::string geometryShaderFilename)
    : vertexShaderFilename(vertexShaderFilename), fragmentShaderFilename(fragmentShaderFilename),
      geometryShaderFilename(geometryShaderFilename)
{
    this->compileShadersFromSource();
}
Shader::Shader(Shader& copy) {
    this->vertexShaderFilename = copy.vertexShaderFilename;
    this->fragmentShaderFilename = copy.fragmentShaderFilename;
    this->geometryShaderFilename = copy.geometryShaderFilename;
    this->compileShadersFromSource();
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
    this->compileShadersFromSource();
}

Shader::Shader(std::string vertexShaderFilename)
    : Shader(vertexShaderFilename.c_str())
{

}

Shader::Shader(std::string vertexShaderFilename, std::string fragmentShaderFilename)
    : Shader(vertexShaderFilename.c_str(), fragmentShaderFilename.c_str())
{

}

Shader::Shader(std::string vertexShaderFilename, std::string fragmentShaderFilename, std::string geometryShaderFilename)
    : Shader(vertexShaderFilename.c_str(), fragmentShaderFilename.c_str(), geometryShaderFilename.c_str())
{

}*/
void Shader::compileShadersFromSource(std::map<std::string, std::string> addedDefinitions)
{
#if useModernOpenGL || !useModernOpenGL
    this->programID = GlobalsGL::f()->glCreateProgram();
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
        } else {
            geometryShaderFilename = "";
        }
    }

    GlobalsGL::f()->glLinkProgram(this->programID);

#endif
}

bool Shader::use(bool update_source_file)
{
    if (update_source_file)
    {
        this->compileShadersFromSource();
    }

    GlobalsGL::checkOpenGLError();
    if (programID > 0) {
        GlobalsGL::f()->glUseProgram(this->programID);
        GlobalsGL::printProgramErrors(this->programID);
        GlobalsGL::printShaderErrors(this->vShader);
        GlobalsGL::printShaderErrors(this->fShader);
//        return true;
    }
    GlobalsGL::checkOpenGLError();
    return programID > 0;
}

void Shader::setBool(std::string pname, bool value)
{
    if (!this->use()) return;
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), (int)value);
}
void Shader::setInt(std::string pname, int value)
{
    if (!this->use()) return;
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}

void Shader::setFloat(std::string pname, float value)
{
    if (!this->use()) return;
    GlobalsGL::f()->glUniform1f(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}

void Shader::setVector(std::string pname, Vector3 value)
{
    if (!this->use()) return;
    this->setVector(pname, value, 3);
}

void Shader::setVector(std::string pname, glm::vec2 value)
{
    if (!this->use()) return;
    GlobalsGL::f()->glUniform2fv(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), 1, &value[0]);
}

//void Shader::setVector(std::string pname, glm::vec3 value)
//{
//    if (!this->use()) return;
//    GlobalsGL::f()->glUniform3fv(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), 1, &value[0]);
//}

void Shader::setVector(std::string pname, glm::vec4 value)
{
    if (!this->use()) return;
    GlobalsGL::f()->glUniform4fv(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), 1, &value[0]);
}
void Shader::setVector(std::string pname, float value[], int n)
{
    if (!this->use()) return;
    GLuint loc = GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str());
    switch (n) {
    case 1:
        GlobalsGL::f()->glUniform1fv(loc, 1, value);
        break;
    case 2:
        GlobalsGL::f()->glUniform2fv(loc, 1, value);
        break;
    case 3:
        GlobalsGL::f()->glUniform3fv(loc, 1, value);
        break;
    case 4:
        GlobalsGL::f()->glUniform4fv(loc, 1, value);
        break;
    }
}
void Shader::setLightSource(std::string pname, LightSource &value)
{
    this->setVector((pname + ".ambiant").c_str(), value.ambiant, 4);
    this->setVector((pname + ".diffuse").c_str(), value.diffuse, 4);
    this->setVector((pname + ".specular").c_str(), value.specular, 4);
}

void Shader::addLightSource(std::string pname, LightSource &value)
{
    this->setLightSource(pname + "[" + std::to_string(lightCount) +"]", value);
    lightCount ++;
    this->setInt(pname + "_count", this->lightCount);
}

void Shader::clearLightSources(std::string pname)
{
    this->lightCount = 0;
    this->setInt(pname + "_count", 0);
}

void Shader::setPositionalLight(std::string pname, PositionalLight &value)
{
    this->setLightSource(pname, value);
    this->setVector((pname + ".position").c_str(), value.position);
}

void Shader::setMaterial(std::string pname, Material &value)
{
    this->setVector((pname + ".ambiant").c_str(), value.ambiant, 4);
    this->setVector((pname + ".diffuse").c_str(), value.diffuse, 4);
    this->setVector((pname + ".specular").c_str(), value.specular, 4);
    this->setFloat((pname + ".shininness").c_str(), value.shininess);
}

// Todo : change the integer type to a template (or select int or float)
void Shader::setTexture2D(std::string pname, int index, Matrix3<int> texture)
{
    int **data = new int*[texture.sizeY];
    for (int i = 0; i < texture.sizeY; i++) {
        data[i] = new int[texture.sizeX];
        for (int j = 0; j < texture.sizeX; j++) {
            data[i][j] = texture.at(j, i);
        }
    }

    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    bool justUpdateTexture = true;
    if (!this->use()) return;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        textureSlotIndices[textureSlot] = texIndex;
        justUpdateTexture = false;
    }
    texIndex = textureSlotIndices[textureSlot];
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texIndex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    if (justUpdateTexture) {
        GlobalsGL::f()->glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, texture.sizeX, texture.sizeY,
        GL_ALPHA_INTEGER_EXT, GL_INT, data);
    } else {
        glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, texture.sizeX, texture.sizeY, 0,
        GL_ALPHA_INTEGER_EXT, GL_INT, data);
    }
    this->setInt(pname, index);

    for (int i = 0; i < texture.sizeX; i++)
        delete[] data[i];
    delete[] data;
}

void Shader::setTexture2D(std::string pname, int index, int width, int height, int *data)
{
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    bool justUpdateTexture = true;
    if (!this->use()) return;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        textureSlotIndices[textureSlot] = texIndex;
        justUpdateTexture = false;
    }
    texIndex = textureSlotIndices[textureSlot];
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glEnable(GL_TEXTURE_2D);
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
}

void Shader::setTexture2D(std::string pname, int index, int width, int height, int **data)
{
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    bool justUpdateTexture = true;
    if (!this->use()) return;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        textureSlotIndices[textureSlot] = texIndex;
        justUpdateTexture = false;
    }
    texIndex = textureSlotIndices[textureSlot];
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glEnable(GL_TEXTURE_2D);
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
}

// Todo : change the integer type to a template (or select int or float)
void Shader::setTexture3D(std::string pname, int index, Matrix3<float> texture)
{
    for (auto& val : texture)
        val = std::max(val, 0.f);
    glEnable(GL_TEXTURE_3D);
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    bool justUpdateTexture = false; // true;
    if (!this->use()) return;
    if (textureSlotIndices.count(textureSlot) == 0) {
        glGenTextures(1, &texIndex);
        textureSlotIndices[textureSlot] = texIndex;
        justUpdateTexture = false;
    }
    texIndex = textureSlotIndices[textureSlot];

    GlobalsGL::f()->glActiveTexture(textureSlot);
    glBindTexture(GL_TEXTURE_3D, texIndex);
//    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    if (justUpdateTexture) {
        GlobalsGL::f()->glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, texture.sizeX, texture.sizeY, texture.sizeZ,
        GL_ALPHA, GL_FLOAT, texture.data.data());
    } else {
        GlobalsGL::f()->glTexImage3D( GL_TEXTURE_3D, 0, GL_ALPHA32F_ARB, texture.sizeX, texture.sizeY, texture.sizeZ, 0,
        GL_ALPHA, GL_FLOAT, texture.data.data());
    }
    this->setInt(pname, index);
}

//void Shader::setMatrix(std::string pname, Matrix value)
//{
//    this->setMatrix(pname, value.)
//}


void Shader::applyToAllShaders(std::function<void (std::shared_ptr<Shader>)> func)
{
    std::set<std::shared_ptr<Shader>>::iterator it;
    for (it = Shader::allShaders.begin(); it != Shader::allShaders.end(); it++) {
        if (bool(*it)) {
            func(*it);
        }
    }
}


std::string Shader::readShaderSource(std::string filename)
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
