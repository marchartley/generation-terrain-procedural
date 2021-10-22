#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>

#include "Viewer.h"

#include "Vector3.h"


int main(int argc, char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

    QGuiApplication app(argc, argv);
    QApplication a(argc, argv);

    Grid* grid = new Grid(200, 200, 100, 10.0 / 100);
    VoxelGrid* vGrid = new VoxelGrid(*grid);
//    grid->fromVoxelGrid(*vGrid);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Viewer view(grid, vGrid, ViewerMode::VOXEL_MODE);

    view.show();

    return app.exec();
}
