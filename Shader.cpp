#include "Shader.h"

#include <QFile>
#include <QTextStream>

std::vector<Shader*> Shader::shaders;

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
    Shader::shaders.push_back(this);
}
void Shader::compileShadersFromSource()
{
    this->programID = GlobalsGL::f()->glCreateProgram();
    if (vertexShaderFilename != nullptr)
    {
        GLuint vShader = GlobalsGL::f()->glCreateShader(GL_VERTEX_SHADER);
        std::string content = Shader::readShaderSource(vertexShaderFilename);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(vShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(vShader);
        GlobalsGL::f()->glAttachShader(this->programID, vShader);
        GlobalsGL::printShaderErrors(vShader);
    }
    if (fragmentShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(fragmentShaderFilename);
        GLuint fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(fShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(fShader);
        GlobalsGL::f()->glAttachShader(this->programID, fShader);
        GlobalsGL::printShaderErrors(fShader);
    }
    if (geometryShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(geometryShaderFilename);
        GLuint gShader = GlobalsGL::f()->glCreateShader(GL_GEOMETRY_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(gShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(gShader);
        GlobalsGL::f()->glAttachShader(this->programID, gShader);
        GlobalsGL::printShaderErrors(gShader);
    }

    GlobalsGL::f()->glLinkProgram(this->programID);
}

void Shader::use(bool update_source_file)
{
    if (update_source_file)
    {
        this->compileShadersFromSource();
    }
    GlobalsGL::checkOpenGLError();
    GlobalsGL::f()->glUseProgram(this->programID);
    GlobalsGL::checkOpenGLError();
}

void Shader::setBool(std::string pname, bool value)
{
    this->use();
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), (int)value);
}
void Shader::setInt(std::string pname, int value)
{
    this->use();
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}

void Shader::setFloat(std::string pname, float value)
{
    this->use();
    GlobalsGL::f()->glUniform1f(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}

void Shader::setVector(std::string pname, Vector3 value)
{
    this->use();
    this->setVector(pname, value, 3);
}
void Shader::setVector(std::string pname, float value[], int n)
{
    this->use();
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


std::string Shader::readShaderSource(std::string filename)
{
    std::string content = "";
    QString qFilename = QString::fromStdString(filename);
    if (!QFile::exists(qFilename))
        qFilename = ":" + qFilename;
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
