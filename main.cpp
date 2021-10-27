#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>

#include "Viewer.h"

#include "Vector3.h"


int main(int argc, char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

    QGuiApplication app(argc, argv);
    QApplication a(argc, argv);

    Grid* grid = new Grid(200, 200, 20, 1.0);
    VoxelGrid* vGrid = new VoxelGrid(40, 40, 20, 1.0);
//    VoxelGrid* vGrid = new VoxelGrid(*grid);
//    grid->fromVoxelGrid(*vGrid);
//    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Viewer view(grid, vGrid, ViewerMode::GRID_MODE | ViewerMode::WIRE_MODE);

    view.show();

    return app.exec();
}
