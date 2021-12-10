#include "Globals.h"
#include "Viewer.h"

#include "Vector3.h"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>
#include "Matrix.h"
#include "Interface.h"
//#define ViewerInterface Viewer

int main(int argc, char *argv[])
{
//    Vector3 a (0, 0, 0);
//    Vector3 b (1, 1, 1);
//    Vector3 c (0, 0, 0);
//    std::hash<Vector3> x;
//    std::cout << x(a) << std::endl;
//    std::cout << x(b) << std::endl;
//    std::cout << x(c) << std::endl;
//    exit(0);

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
//    QGuiApplication app(argc, argv);
    QApplication app(argc, argv);

    QGLFormat glFormat;
    glFormat.setVersion(4, 1);
    glFormat.setProfile(QGLFormat::CompatibilityProfile);
    glFormat.setSampleBuffers(true);
    glFormat.setDefaultFormat(glFormat);
    glFormat.setSwapInterval(1);
    QGLWidget widget(glFormat);
    widget.makeCurrent();

    const QOpenGLContext *context = GlobalsGL::context();

    qDebug() << "Context valid: " << context->isValid();
    qDebug() << "Really used OpenGl: " << context->format().majorVersion() << "." << context->format().minorVersion();
    qDebug() << "OpenGl information: VENDOR:       " << (const char*)glGetString(GL_VENDOR);
    qDebug() << "                    RENDERDER:    " << (const char*)glGetString(GL_RENDERER);
    qDebug() << "                    VERSION:      " << (const char*)glGetString(GL_VERSION);
    qDebug() << "                    GLSL VERSION: " << (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);

    ViewerInterface vi;
    vi.show();

    return app.exec();
}
