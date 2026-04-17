#ifndef GLOBALS_H
#define GLOBALS_H

#include "Utils/stb_image.h"
#include "Utils/stb_image_write.h"
#include "Utils/stb_image_resize.h"
#define useModernOpenGL 1

#include <QOpenGLFunctions_4_5_Compatibility>
#include <QtWidgets>
#include <fstream>
#include <iostream>
#include "Utils/Random.h"

#define numVAOs 1000
#define numVBOs 4000

class GlobalsGL {
public:
    static QOpenGLContext *_context;
    static QOpenGLContext *context();
    static QOpenGLFunctions* _f;
//    static QOpenGLFunctions_4_6_Compatibility *f(); // Alias for ef()
    static QOpenGLFunctions_4_5_Compatibility *f(); // Alias for ef()
    static QOpenGLExtraFunctions* _ef;
    static QOpenGLFunctions_4_5_Compatibility* f45();
    static QOpenGLFunctions_4_5_Compatibility* _ef45;
//    static QOpenGLFunctions_4_6_Compatibility* f46();
//    static QOpenGLFunctions_4_6_Compatibility* _ef46;
    static QOpenGLExtraFunctions *ef();
    static GLuint _renderingProgram;
    static GLuint renderingProgram();

    static GLuint vao[numVAOs];
    static GLuint vbo[numVBOs];
    static bool buffersGenerated;

    static GLuint currentBufferId;

    static std::string readShaderSource(const std::string& filename);


    static GLuint createShaderProgram(const std::string& vertexShaderFile = "", std::string fragmentShaderFile = "");

    static GLuint newBufferId();

    static void generateBuffers();

    static bool checkOpenGLError();
    static bool printShaderErrors(GLuint shader);
    static bool printProgramErrors(int program);

//    static void GLAPIENTRY MessageCallback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam );
};
#endif // GLOBALS_H
