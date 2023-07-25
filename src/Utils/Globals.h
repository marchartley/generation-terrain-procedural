#ifndef GLOBALS_H
#define GLOBALS_H

#include "Utils/stb_image.h"
#include "Utils/stb_image_write.h"
#define useModernOpenGL 1

//class random_gen;
//class GlobalsGL;

//#include <GL/glew.h>
#include <random>
//#include <QOpenGLContext>
// #include <glad/glad.h>
#include <QOpenGLFunctions_4_5_Core>
#include <QtWidgets>
//#include <QOpenGLFunctions>
#include <fstream>
#include <iostream>
//#include "Graphics/Shader.h"
//#include "Utils/Utils.h"
#include "Utils/FastNoiseLit.h"

#define numVAOs 1000
#define numVBOs 4000

/*
GLuint renderingProgram;
GLuint vao[numVAOs];
GLuint vbo[numVBOs];*/

class random_gen {
public:
    static std::default_random_engine random_generator;
    static FastNoiseLite perlinNoise;
    static float generate(float min, float max) {
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(random_gen::random_generator);
    }
    static float generate(float max) {
        std::uniform_real_distribution<float> distribution(0.f, max);
        return distribution(random_gen::random_generator);
    }
    static float generate() {
        std::uniform_real_distribution<float> distribution(0.f, 1.f);
        return distribution(random_gen::random_generator);
    }

    static float generate_normal(float mu, float sigma) {
        std::normal_distribution<float> distribution(mu, sigma);
        return distribution(random_gen::random_generator);
    }
    static float generate_normal() {
        std::uniform_real_distribution<float> distribution(0.f, 1.f);
        return distribution(random_gen::random_generator);
    }

    static float generate_perlin(float x, float y, float z = 0) {
        return perlinNoise.GetNoise(x, y, z);
    }
};

class GlobalsGL {
public:
    static QOpenGLContext *_context;
    static QOpenGLContext *context();
    static QOpenGLFunctions* _f;
    static QOpenGLFunctions_4_5_Core *f(); // Alias for ef()
    static QOpenGLExtraFunctions* _ef;
    static QOpenGLFunctions_4_5_Core* f45();
    static QOpenGLFunctions_4_5_Core* _ef45;
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

//    static void GLAPIENTRY MessageCallback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam );
};
#endif // GLOBALS_H
