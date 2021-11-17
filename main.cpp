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
    QGuiApplication app(argc, argv);
    QApplication a(argc, argv);

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

    Grid* grid = new Grid(10, 10, 40, 1.0);
    VoxelGrid* vGrid = new VoxelGrid(20, 20, 160, 1.0);
//    VoxelGrid* vGrid = new VoxelGrid(*grid);
//    grid->fromVoxelGrid(*vGrid);

    Viewer view(grid, vGrid, MapMode::VOXEL_MODE, ViewerMode::FILL_MODE);
    view.show();

    return app.exec();
}
