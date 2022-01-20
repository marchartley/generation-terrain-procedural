//#define OPENVDB_DLL

#include "Utils/Globals.h"
#include "Interface/Viewer.h"

#include "DataStructure/Vector3.h"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>
#include "DataStructure/Matrix.h"
#include "Interface/Interface.h"
//#define ViewerInterface Viewer

//#include <openvdb/openvdb.h>
#include "DataStructure/Matrix3.h"
#include "Graph/FastPoissonGraph.h"
#include "Karts/KarstPathsGeneration.h"

int main(int argc, char *argv[])
{
//    BSpline b({Vector3(0, 0), Vector3(0, 1), Vector3(1, 1), Vector3(1, 0)});

//    Vector3 point;
//    for (float t = 0; t <= 1.0; t += 0.2) {
//        point = b.getPoint(t);
//        std::cout << point << std::endl;
//    }
//    openvdb::initialize();

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
