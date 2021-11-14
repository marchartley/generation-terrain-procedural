#include "Shader.h"

#include <QFile>
#include <QTextStream>

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
    }
    if (fragmentShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(fragmentShaderFilename);
        GLuint fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(fShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(fShader);
        GlobalsGL::f()->glAttachShader(this->programID, fShader);
    }
    if (geometryShaderFilename != nullptr)
    {
        std::string content = Shader::readShaderSource(geometryShaderFilename);
        GLuint gShader = GlobalsGL::f()->glCreateShader(GL_GEOMETRY_SHADER);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(gShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(gShader);
        GlobalsGL::f()->glAttachShader(this->programID, gShader);
    }

    GlobalsGL::f()->glLinkProgram(this->programID);
}

void Shader::use()
{
    GlobalsGL::f()->glUseProgram(this->programID);
}

void Shader::setBool(std::string pname, bool value)
{
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), (int)value);
}
void Shader::setInt(std::string pname, int value)
{
    GlobalsGL::f()->glUniform1i(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}

void Shader::setFloat(std::string pname, float value)
{
    GlobalsGL::f()->glUniform1f(GlobalsGL::f()->glGetUniformLocation(this->programID, pname.c_str()), value);
}


std::string Shader::readShaderSource(std::string filename)
{
    std::string content = "";
    QString qFilename = ":" + QString::fromStdString(filename);
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
