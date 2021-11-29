#include "Globals.h"
#include "Viewer.h"

#include "Vector3.h"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
//#include <qapplication.h>
#include <iostream>
#include "Matrix.h"



int main(int argc, char *argv[])
{
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

    Grid* grid = nullptr; // new Grid(100, 100, 40, 1.0);
    LayerBasedGrid* lGrid = new LayerBasedGrid(10, 10, 50);
    VoxelGrid* vGrid = new VoxelGrid(10, 10, 50, 1.0, -0.0);
//    VoxelGrid* vGrid = new VoxelGrid(*grid);
//    grid->fromVoxelGrid(*vGrid);

    Viewer view(grid, vGrid, lGrid, MapMode::LAYER_MODE, ViewerMode::FILL_MODE);
    view.show();

    return app.exec();
}
