#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

//#define OPENVDB_DLL
#include "Utils/Globals.h"
//#include "sim-fluid-loganzartman/Game.hpp"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>
#include "Interface/Interface.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
    QApplication app(argc, argv);


    QGLFormat glFormat;
    glFormat.setVersion(4, 5);
    glFormat.setProfile(QGLFormat::CoreProfile);
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

/*
    int isize = 64;
    int jsize = 64;
    int ksize = 64;
    float dx = 0.125;
    FluidSimulation fluidsim(isize, jsize, ksize, dx);

    fluidsim.setSurfaceSubdivisionLevel(2);

    float x, y, z;
    fluidsim.getSimulationDimensions(&x, &y, &z);
    fluidsim.addImplicitFluidPoint(x/2, y/2, z/2, 7.0);

    fluidsim.addBodyForce(0.0, -25.0, 0.0);
    fluidsim.initialize();

    float timestep = 1.0 / 30.0;
    for (;;) {
        fluidsim.update(timestep);
    }

    return 0;*/



    /*
    std::cout << timeIt([]() {
        Fluid fluid;
        gfx::Program texture_copy_program;

        fluid.init();
        texture_copy_program.vertex({"screen_quad.vs.glsl"}).fragment({"texture_copy.fs.glsl"}).compile();

    //    fluid.resize(40, 40, 40);

        for (int i = 0; i < 100; i++) {
            fluid.step();
            fluid.ssbo_barrier();

//            if (i % 100 != 0) continue;
//            const auto particles = fluid.particle_ssbo.map_buffer_readonly<Particle>();
            auto grid = fluid.grid_ssbo.map_buffer<GridCell>();

            int count = 0;
            for (int y = fluid.grid_dimensions.y-1; y >= 0; --y) {
                for (int x = 0; x < fluid.grid_dimensions.x; ++x) {
                    for (int z = 0; z < fluid.grid_dimensions.z; ++z) {
                        const int i = fluid.get_grid_index({x, y, z});
                        if (z == 5) std::cout << (grid[i].type == GRID_FLUID ? "O" : (grid[i].type == GRID_SOLID ? "#" : "."));
                        count += (grid[i].type == 2 ? 1 : 0);
                    }
                }
                std::cout << "\n";
            }
            std::cout << "\n" << count << "\n" << std::endl;
        }
    }) << "ms" << std::endl;

    return 0;
    */
    /*Vector3 A = Vector3(1, 0, 0);
    Vector3 B = Vector3(0.5, 1, 1).normalize();
    Vector3 UP = Vector3(0, 0, 1); // A.cross(B);
    Vector3 LEFT = UP.cross(A);

    Vector3 C = Vector3(
                    std::acos(B.y),
                    std::acos(B.z),
                    std::acos(B.x)
                ) * (180.f / 3.141592f);

    std::cout << C << std::endl;
    return 0;*/
/*
    AABBox bbox(Vector3(0, 0, 0), Vector3(2, 2, 1));

    Matrix3<float> M = Matrix3<float>({
                                          {0, 1},
                                          {1, 0}
                                          });
//    std::cout << M.displayValues() << std::endl;

    Matrix3<float> m(10, 10);
    Vector3 ratio = (M.getDimensions() - Vector3(1, 1, 0)) / (bbox.dimensions() - Vector3(1, 1, 0));
    for (int _x = 0; _x < m.sizeX; _x++) {
        for (int _y = 0; _y < m.sizeY; _y++) {
            for (int _z = 0; _z < m.sizeZ; _z++) {
                float x = _x, y = _y, z = _z;
                Vector3 pos(x, y, z);
                Vector3 query = bbox.normalize(pos); // * M.getDimensions(); //pos * ratio;
                query.z = 0;
                auto val = M.interpolate(query);
                m.at(pos) = val;
            }
        }
    }
    std::cout << M.displayValues() << std::endl;
    std::cout << m.displayValues() << std::endl;
    return 0;*/

/*
    ShapeCurve A = ShapeCurve({
                                  Vector3(0, 0.5, 0),
                                  Vector3(0.5, 0, 0),
                                  Vector3(1, 0.1, 0),
                                  Vector3(0.9, 0.2, 0),
                                  Vector3(0.5, 0.5, 0),
                                  Vector3(0.9, 0.8, 0),
                                  Vector3(1, 0.9, 0),
                                  Vector3(0.5, 1, 0)
                              });
    ShapeCurve B = ShapeCurve({
                                  Vector3(2, 0.5, 0),
                                  Vector3(1.5, 0, 0),
                                  Vector3(1, 0.1, 0),
                                  Vector3(1.1, 0.2, 0),
                                  Vector3(1.5, 0.5, 0),
                                  Vector3(1.1, 0.8, 0),
                                  Vector3(1, 0.9, 0),
                                  Vector3(1.5, 1, 0)
                              });
    ShapeCurve AB = merge(A, B);

    Plotter::init();
    Plotter::getInstance()->addPlot(A.closedPath(), "A", Qt::blue);
    Plotter::getInstance()->addPlot(B.closedPath(), "B", Qt::green);
    Plotter::getInstance()->addPlot(AB.closedPath(), "AB", Qt::red);
    Plotter::getInstance()->addScatter(AB.points, "");
    std::cout << A << std::endl;
    std::cout << B << std::endl;
    std::cout << AB << std::endl;
    return Plotter::getInstance()->exec();*/

//    RegularSimplicialComplex grid(10, 10);
//    grid.getNode(2, 3)->value = 0;
//    grid.getNode(1, 3)->value = 0;
//    grid.getNode(1, 4)->value = 0;
//    grid.removeUnavailableLinks();
//    grid.display();

//    return 0;

//    auto g = CombinMap();
//    g.addFace(5, {}, {100, 101, 102, 103, 104});
//    g.addFace(4, {g.root->beta2, g.root->beta2->beta1}, {0});
//    g.addFace(4, {g.root->beta2, g.root->beta2->beta1}, {1});
//    g.debug();
//    Graph<int> G = g.toGraph().forceDrivenPositioning();
//    auto dual = g.dual(g.root->beta1->beta2);
//    dual.debug();
//    G = dual.toGraph().forceDrivenPositioning();
//    return 0;

    ViewerInterface vi;
    vi.show();

    return app.exec();
}
