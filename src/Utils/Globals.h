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
//#include "Graphics/Shader.h"

#define numVAOs 1000
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

    static void GLAPIENTRY MessageCallback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam );
};
#endif // GLOBALS_H
