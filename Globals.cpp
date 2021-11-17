#include "Globals.h"

#include <string>
#include <fstream>

std::default_random_engine random_gen::random_generator;


QOpenGLContext* GlobalsGL::_context;
QOpenGLFunctions* GlobalsGL::_f;
QOpenGLExtraFunctions* GlobalsGL::_ef;
GLuint GlobalsGL::_renderingProgram;
GLuint GlobalsGL::vao[numVAOs];
GLuint GlobalsGL::vbo[numVBOs];
GLuint GlobalsGL::currentBufferId = 0;
bool GlobalsGL::buffersGenerated = false;

QOpenGLContext* GlobalsGL::context() {
    if (GlobalsGL::_context == nullptr)
        GlobalsGL::_context = QOpenGLContext::currentContext();
    return GlobalsGL::_context;
}
QOpenGLExtraFunctions* GlobalsGL::f() {
    return GlobalsGL::ef();
//    if (GlobalsGL::_f == nullptr)
//        GlobalsGL::_f = GlobalsGL::context()->functions();
//    return GlobalsGL::_f;
}
QOpenGLExtraFunctions* GlobalsGL::ef() {
    if (GlobalsGL::_ef == nullptr)
        GlobalsGL::_ef = GlobalsGL::context()->extraFunctions();
    return GlobalsGL::_ef;
}/*
GLuint GlobalsGL::renderingProgram() {
    return GlobalsGL::_renderingProgram;
}

std::string GlobalsGL::readShaderSource(std::string filename)
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

GLuint GlobalsGL::createShaderProgram(std::string vertexShaderFile, std::string fragmentShaderFile)
{
    GLuint vShader = GlobalsGL::f()->glCreateShader(GL_VERTEX_SHADER);
    GLuint fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
    if (vertexShaderFile != "")
    {
        std::string content = GlobalsGL::readShaderSource(vertexShaderFile);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(vShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(vShader);
    }
    if (fragmentShaderFile != "")
    {
        std::string content = GlobalsGL::readShaderSource(fragmentShaderFile);
        const char* src = content.c_str();
        GlobalsGL::f()->glShaderSource(fShader, 1, &src, NULL);
        GlobalsGL::f()->glCompileShader(fShader);
    }
    GLuint vProgram = GlobalsGL::f()->glCreateProgram();

    if(vertexShaderFile != "")
        GlobalsGL::f()->glAttachShader(vProgram, vShader);
    if(fragmentShaderFile != "")
        GlobalsGL::f()->glAttachShader(vProgram, fShader);
    GlobalsGL::f()->glLinkProgram(vProgram);

    GlobalsGL::_renderingProgram = vProgram;

    GlobalsGL::checkOpenGLError();
    GlobalsGL::printShaderErrors(vShader);
    GlobalsGL::printShaderErrors(fShader);
    GlobalsGL::printProgramErrors(vProgram);

    return GlobalsGL::_renderingProgram;
}
*/
void GlobalsGL::generateBuffers()
{
    if(GlobalsGL::buffersGenerated)
        return;
    GlobalsGL::f()->glGenVertexArrays(numVAOs, GlobalsGL::vao);
    GlobalsGL::f()->glBindVertexArray(GlobalsGL::vao[0]);
    GlobalsGL::f()->glGenBuffers(numVBOs, GlobalsGL::vbo);
    GlobalsGL::buffersGenerated = true;
}
GLuint GlobalsGL::newBufferId()
{
    return GlobalsGL::currentBufferId += 4; // Gives space for vertex, texture and normals
}
bool GlobalsGL::checkOpenGLError()
{
    bool error = false;
    int glErr = glGetError();
    while(glErr != GL_NO_ERROR)
    {
        std::cout << "[OpenGL] Error: " << glErr << std::endl;
        error = true;
        glErr = glGetError();
    }
    return !error;
}

bool GlobalsGL::printShaderErrors(GLuint shader)
{
    int state = 0;
    GlobalsGL::f()->glGetShaderiv(shader, GL_COMPILE_STATUS, &state);
    if (state == 1)
        return true;
    int len = 0;
    int chWritten = 0;
    char* log;
    GlobalsGL::f()->glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
    if (len > 0)
    {
        log = (char*)malloc(len);
        GlobalsGL::f()->glGetShaderInfoLog(shader, len, &chWritten, log);
        std::cout << "[OpenGL] Shader error: " << log << std::endl;
        free(log);
    }
    return false;
}
bool GlobalsGL::printProgramErrors(int program)
{
    int state = 0;
    GlobalsGL::f()->glGetProgramiv(program, GL_LINK_STATUS, &state);
    if (state == 1)
        return true;
    int len = 0;
    int chWritten = 0;
    char* log;
    GlobalsGL::f()->glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
    if (len > 0)
    {
        log = (char*)malloc(len);
        GlobalsGL::f()->glGetProgramInfoLog(program, len, &chWritten, log);
        std::cout << "[OpenGL] Program error: " << log << std::endl;
        free(log);
    }
    return false;
}
