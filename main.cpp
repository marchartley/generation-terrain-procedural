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
//    float data[3][3] = {{0, 0, 1}, {2, -1, 3}, {1, 1, 4}};
//    std::vector<std::vector<float>> data = {{0, 0, 1}, {2, -1, 3}, {1, 1, 4}};
//    Matrix test(4, 4);
//    test[0][0] = 0;
//    test[0][1] = 2;
//    test[0][2] = 1;
//    test[0][3] = 4;
//    test[1][0] = 5;
//    test[1][1] = 6;
//    test[1][2] = 7;
//    test[1][3] = 8;
//    test[2][0] = 9;
//    test[2][1] = 9;
//    test[2][2] = 8;
//    test[2][3] = 7;
//    test[3][0] = 6;
//    test[3][1] = 5;
//    test[3][2] = 4;
//    test[3][3] = 3;
//    Matrix test(3, 3);
//    test[0][0] = 0;
//    test[0][1] = 2;
//    test[0][2] = 3;
//    test[1][0] = 4;
//    test[1][1] = 5;
//    test[1][2] = 6;
//    test[2][0] = 7;
//    test[2][1] = 8;
//    test[2][2] = 9;
    Matrix test(3, 3);
    test[0][0] = 0;
    test[0][1] = 0;
    test[0][2] = 1;
    test[1][0] = 2;
    test[1][1] = -1;
    test[1][2] = 3;
    test[2][0] = 1;
    test[2][1] = 1;
    test[2][2] = 4;

    std::cout << test << std::endl;
//    std::cout << test.det() << std::endl;
//    std::cout << test.cofactors() << std::endl;
    std::cout << test.inverse() << std::endl;

    exit(0);


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
    qDebug() << "endstuff\n";

    Grid* grid = new Grid(80, 80, 40, 1.0);
    VoxelGrid* vGrid = new VoxelGrid(100, 100, 40, 1.0);
//    VoxelGrid* vGrid = new VoxelGrid(*grid);
//    grid->fromVoxelGrid(*vGrid);

    Viewer view(grid, vGrid, MapMode::GRID_MODE, ViewerMode::FILL_MODE);
    view.show();

    return app.exec();
}
