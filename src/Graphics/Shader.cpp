#include "Graphics/Shader.h"

#include <QFile>
#include <QTextStream>

std::shared_ptr<Shader> Shader::default_shader = nullptr;
std::set<std::shared_ptr<Shader>> Shader::allShaders;

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
Shader::Shader(Shader& copy) {
    this->vertexShaderFilename = copy.vertexShaderFilename;
    this->fragmentShaderFilename = copy.fragmentShaderFilename;
    this->geometryShaderFilename = copy.geometryShaderFilename;
    this->compileShadersFromSource();
}
void Shader::compileShadersFromSource()
{
#if useModernOpenGL || !useModernOpenGL
    this->programID = GlobalsGL::f()->glCreateProgram();
    if (vertexShaderFilename != nullptr)
    {
        this->vShader = GlobalsGL::f()->glCreateShader(GL_VERTEX_SHADER);
        std::string content = Shader::readShaderSource(vertexShaderFilename);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(this->vShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(this->vShader);
        GlobalsGL::f()->glAttachShader(this->programID, this->vShader);
        GlobalsGL::printShaderErrors(this->vShader);
    }
    if (fragmentShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(fragmentShaderFilename);
        this->fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(this->fShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(this->fShader);
        GlobalsGL::f()->glAttachShader(this->programID, this->fShader);
        GlobalsGL::printShaderErrors(this->fShader);
    }
    if (geometryShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(geometryShaderFilename);
        this->gShader = GlobalsGL::f()->glCreateShader(GL_GEOMETRY_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(this->gShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(this->gShader);
        GlobalsGL::f()->glAttachShader(this->programID, this->gShader);
        GlobalsGL::printShaderErrors(this->gShader);
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

    if (programID > 0) {
        GlobalsGL::checkOpenGLError();
        GlobalsGL::f()->glUseProgram(this->programID);
        GlobalsGL::printProgramErrors(this->programID);
        GlobalsGL::printShaderErrors(this->vShader);
        GlobalsGL::printShaderErrors(this->fShader);
        GlobalsGL::checkOpenGLError();
        return true;
    }
    return false;
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
    /// WARNING: It should work as GL_TEXTURE1 is defined just after GL_TEXTURE0, but it may not be stable...
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    if (!this->use()) return;
    GlobalsGL::f()->glGenTextures(1, &texIndex);
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glEnable(GL_TEXTURE_2D);
    GlobalsGL::f()->glBindTexture(GL_TEXTURE_2D, texIndex);
    //Integer textures must use nearest filtering mode
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, texture.sizeX, texture.sizeY, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, texture.data.data());
    this->setInt(pname, index);

}

// Todo : change the integer type to a template (or select int or float)
void Shader::setTexture3D(std::string pname, int index, Matrix3<float> texture)
{
    /// WARNING: It should work as GL_TEXTURE1 is defined just after GL_TEXTURE0, but it may not be stable...
    int textureSlot = GL_TEXTURE0 + index;

    GLuint texIndex;
    if (!this->use()) return;
    glGenTextures(1, &texIndex);
    GlobalsGL::f()->glActiveTexture(textureSlot);
    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, texIndex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    GlobalsGL::f()->glTexImage3D( GL_TEXTURE_3D, 0, GL_ALPHA32F_ARB, texture.sizeX, texture.sizeY, texture.sizeZ, 0,
    GL_ALPHA, GL_FLOAT, texture.data.data());
    this->setInt(pname, index);
}


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
    if (!QFile::exists(qFilename))
        std::cerr << "The shader " << filename << " doesn't exist!" << std::endl;
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
