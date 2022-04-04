#include "Utils/Globals.h"

#include <string>
#include <fstream>

#define UNUSED(expr) (void)(expr)

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
}
QOpenGLExtraFunctions* GlobalsGL::ef() {
    if (GlobalsGL::_ef == nullptr)
        GlobalsGL::_ef = GlobalsGL::context()->extraFunctions();
    return GlobalsGL::_ef;
}
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
    return GlobalsGL::currentBufferId += 1; // Gives space for vertex, texture and normals
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
void GLAPIENTRY GlobalsGL::MessageCallback( GLenum source, GLenum type,
                                            GLuint id, GLenum severity,
                                            GLsizei length, const GLchar* message,
                                            const void* userParam )
{
    UNUSED(source);
    UNUSED(type);
    UNUSED(length);
    UNUSED(userParam);
    if (id == 131154) return; // Ignore "Pixel-path performance warning: Pixel transfer is synchronized with 3D rendering." due to screenshots
    if (severity == GL_DEBUG_SEVERITY_HIGH || severity == GL_DEBUG_SEVERITY_MEDIUM || severity == GL_DEBUG_SEVERITY_LOW) {
        std::string s_severity = (severity == GL_DEBUG_SEVERITY_HIGH ? "High" : severity == GL_DEBUG_SEVERITY_MEDIUM ? "Medium" : "Low");
        std::cout << "Error " << id << " [severity=" << s_severity << "]: " << message << std::endl;
    }
}
