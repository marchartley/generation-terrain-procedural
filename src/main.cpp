//#define OPENVDB_DLL
#include "Utils/Globals.h"
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

#include "third-party/glm/glm.hpp"
#include "third-party/glm/ext.hpp"

//#include "Graphics/DisplayGraphics.h"
//#include "Graph/TopoMap.h"
//#include "Graph/RegularSimplicialComplex.h"
/*
float eval (Vector3 pos, float width, float depth, float height) {
    Vector3 minPos = Vector3();
    Vector3 maxPos = Vector3(width, depth, height);
    float distanceToBBox = Vector3::signedDistanceToBoundaries((pos / maxPos), minPos, Vector3(1, 1, 1)); //maxPos);
//        float distanceFactor = std::max(1.f, maxPos.maxComp() * .5f); //(supportWidth - width).maxComp() * .5f); // Don't get a 0 here
//        distanceToBBox /= distanceFactor; // Make it depending on the support area
//    float distanceFalloff = interpolation::wyvill(std::clamp(distanceToBBox * 2.f, 0.f, 1.f));
    float distanceFalloff = 1.f - (distanceToBBox + 0.5f); // std::clamp(, 0.f, 1.f);
    float evaluation = distanceFalloff; // std::sin(distanceFalloff * 10.f); //std::clamp(distanceFalloff, 0.f, 1.f); // Maybe...
    return evaluation; // - 0.5f ?
}*/
int main(int argc, char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
    QApplication app(argc, argv);

/*
    float width = 5.f, depth = 5.f, height = 5.f;
    Matrix3<float> m(3 * width, 3 * depth);
    m.raiseErrorOnBadCoord = false;
    for (int x = -width; x < 2 * width; x++) {
        for (int y = -width; y < 2 * depth; y++) {
            m.at(width + x, height + y) = eval(Vector3(x, y, height * .5f), width, depth, height);
        }
    }
    std::cout << std::setw(2) << m.displayValues() << std::endl;
    Plotter plt;
    plt.addImage(m);
    plt.exec();*/

    /*
    Matrix3<Vector3> img(10, 10);
    for (int i = 0; i < img.sizeX; i++) {
        for (int j = 0; j < img.sizeY; j++) {
            img.at(i, j) = Vector3(float(i) / float(img.sizeX), float(j) / float(img.sizeY), 0.5f);//Vector3::random().abs();
        }
    }
    Plotter plt;
    plt.addImage(img);
    return plt.exec();
    */
    /*
    std::function blend = [](float valA, float valB) -> float {
        float n = 2.f;
        return std::pow(std::pow(valA, n) + std::pow(valB, n), 1/n);
    };

    int X = 20, Y = 20;
    Matrix3<float> matA(X, Y);
    Matrix3<float> matB(X, Y);
    Matrix3<float> matAB(X, Y);
    float delta = 0.001;
    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            float _x = x / float(X - 1), _y = 1.f - (y / float(Y - 1));
//            float blendVal = blend(_x, _y);
            Vector3 derivee = Vector3(
                        blend(_x + delta, _y) - blend(_x - delta, _y),
                        blend(_x, _y + delta) - blend(_x, _y - delta)
                        ).normalize();

            float contribA = derivee.x / (derivee.x + derivee.y);
            float contribB = derivee.y / (derivee.x + derivee.y);
            matA.at(x, y) = contribA;
            matB.at(x, y) = contribB;
        }
    }

    std::cout << matA.displayAsPlot(0.f, 1.f, {}, {{.5f, "/"}}) << std::endl;

    return 0;
*/
    /*BSpline opCurve = BSpline({
                                  Vector3(.0f, .5f),
                                  Vector3(.3f, .5f),
                                  Vector3(1.f, 1.f),
//                                  Vector3(2.f, 2.f),
//                                  Vector3(1.f, 1.f),
                                  Vector3(.5f, .3f),
                                  Vector3(.5f, .0f)
                              });
    Plotter* plot = new Plotter;
    plot->addPlot(opCurve.getPath(50));
    plot->exec();
    delete plot;
    return 0;*/
    /*
    Vector3 gridRes(151, 151, 1);
    Matrix3<float> grid(gridRes);
    for (int x = 0; x < gridRes.x; x++) {
        for (int y = 0; y < gridRes.y; y++) {
            Vector3 pos = Vector3(x, y) / gridRes;
            float dist = opCurve.estimateSignedDistanceFrom(pos);
            grid.at(x, y) = 1.f - std::clamp(dist + .5f, 0.f, 1.f);
        }
    }
    std::cout << grid.displayAsPlot(0, 0, {"-", "+"}, {{.5f, "0"}}, 0.005f) << std::endl;
    return 0;*/

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

    ViewerInterface vi;
    vi.show();

    return app.exec();
}
