#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>
#include <QOpenGLContext>
#include <QtWidgets>
#include <QOpenGLFunctions>
#include <fstream>
#include <iostream>

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
    static QOpenGLFunctions* _f;
    static QOpenGLContext *context() {
        if (GlobalsGL::_context == nullptr)
            GlobalsGL::_context = QOpenGLContext::currentContext();
        return GlobalsGL::_context;
    }
    static QOpenGLFunctions *f() {
        if (GlobalsGL::_f == nullptr)
            GlobalsGL::_f = GlobalsGL::context()->functions();
        return GlobalsGL::_f;
    }

    static std::string readShaderSource(std::string filename)
    {
        std::string content = "";
        QString qFilename = QCoreApplication::applicationDirPath() + ".." + QDir::separator() + QString::fromStdString(filename);
        QFile file(qFilename);
        file.open(QIODevice::ReadOnly | QIODevice::Text);
        std::string line;
        QTextStream in(&file);
        while (!file.atEnd()) {
            line = in.readLine().toStdString();
            content += line + " \n";
        }
        file.close();
        return content;
    }

    static GLuint createShaderProgram(std::string vertexShaderFile = "", std::string fragmentShaderFile = "")
    {
        if (vertexShaderFile != "")
        {
            GLuint vShader = GlobalsGL::f()->glCreateShader(GL_VERTEX_SHADER);
            std::string content = GlobalsGL::readShaderSource(vertexShaderFile);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(vShader, 1, &src, NULL);
        }
        if (fragmentShaderFile != "")
        {
            GLuint fShader = GlobalsGL::f()->glCreateShader(GL_FRAGMENT_SHADER);
            std::string content = GlobalsGL::readShaderSource(fragmentShaderFile);
            const char* src = content.c_str();
            GlobalsGL::f()->glShaderSource(fShader, 1, &src, NULL);
        }
    }
};
#endif // GLOBALS_H
