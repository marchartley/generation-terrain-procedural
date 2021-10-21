#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>

#include "Viewer.h"

#include "vector3.h"


int main(int argc, char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

    QGuiApplication app(argc, argv);
    QApplication a(argc, argv);

    Grid grid(100, 100, 10.0 / 100);
    VoxelGrid vGrid(grid);
//    grid.fromVoxelGrid(vGrid);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    Viewer view(vGrid);

    view.show();

    return app.exec();
}
