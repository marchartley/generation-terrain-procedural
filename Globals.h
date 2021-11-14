#ifndef GLOBALS_H
#define GLOBALS_H

//class random_gen;
//class GlobalsGL;

#include <GL/glew.h>
#include <random>
//#include <QOpenGLContext>
#include <QtWidgets>
//#include <QOpenGLFunctions>
#include <fstream>
#include <iostream>
//#include "Shader.h"

#define numVAOs 4000
#define numVBOs 4000

/*
GLuint renderingProgram;
GLuint vao[numVAOs];
GLuint vbo[numVBOs];*/

class random_gen {
public:
    static std::default_random_engine random_generator;
    static float generate(float min = 0.0, float max = 1.0) {
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(random_gen::random_generator);
    }
};

class GlobalsGL {
public:
    static QOpenGLContext *_context;
    static QOpenGLContext *context();
    static QOpenGLFunctions* _f;
    static QOpenGLExtraFunctions *f(); // Alias for ef()
    static QOpenGLExtraFunctions* _ef;
    static QOpenGLExtraFunctions *ef();
    static GLuint _renderingProgram;
    static GLuint renderingProgram();

    static GLuint vao[numVAOs];
    static GLuint vbo[numVBOs];
    static bool buffersGenerated;

    static GLuint currentBufferId;

    static std::string readShaderSource(std::string filename);


    static GLuint createShaderProgram(std::string vertexShaderFile = "", std::string fragmentShaderFile = "");

    static GLuint newBufferId();

    static void generateBuffers();

    static bool checkOpenGLError();
    static bool printShaderErrors(GLuint shader);
    static bool printProgramErrors(int program);

    static void GLAPIENTRY
    MessageCallback( GLenum source,
                     GLenum type,
                     GLuint id,
                     GLenum severity,
                     GLsizei length,
                     const GLchar* message,
                     const void* userParam )
    {
        std::cout << "Error [severity=" << severity << "]: " << message << std::endl;
//      fprintf( stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
//               ( type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : "" ),
//                type, severity, message );
    }
};
/*
GLuint CreateShaderProgram(std::string vertexShaderFile, std::string fragmentShaderFile)
{
    GLuint vShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fShader = glCreateShader(GL_FRAGMENT_SHADER);
    if (vertexShaderFile != "")
    {
        std::string content = GlobalsGL::readShaderSource(vertexShaderFile);
        const char* src = content.c_str();
        glShaderSource(vShader, 1, &src, NULL);
        glCompileShader(vShader);
    }
    if (fragmentShaderFile != "")
    {
        std::string content = GlobalsGL::readShaderSource(fragmentShaderFile);
        const char* src = content.c_str();
        glShaderSource(fShader, 1, &src, NULL);
        glCompileShader(fShader);
    }
    GLuint vProgram = glCreateProgram();

    if(vertexShaderFile != "")
        glAttachShader(vProgram, vShader);
    if(fragmentShaderFile != "")
        glAttachShader(vProgram, fShader);
    glLinkProgram(vProgram);

    renderingProgram = vProgram;

    GlobalsGL::checkOpenGLError();
    GlobalsGL::printShaderErrors(vShader);
    GlobalsGL::printShaderErrors(fShader);
    GlobalsGL::printProgramErrors(vProgram);

    return renderingProgram;
}
*/
#endif // GLOBALS_H
