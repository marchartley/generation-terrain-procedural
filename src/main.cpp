#include "Utils/Table.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION

//#define OPENVDB_DLL
#include "Utils/Globals.h"
//#include "sim-fluid-loganzartman/Game.hpp"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>
#include "Interface/Interface.h"

#include "EnvObject/EnvObject.h"


#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <set>

using namespace std;

int main(int argc, char *argv[])
{
    /*
    GridF values = GridF::random(100, 100, 50) - .8f;
    for (int x = 0; x < values.sizeX; x++) {
        for (int y =0; y < values.sizeY; y++) {
            for (int z = 0; z < values.sizeZ; z++) {
                float h = float(z) / float(values.sizeZ);
                values.at(x, y, z) -= h * .5f;
            }
        }
    }
//    BVHTree tree;
    Mesh m = Mesh::applyMarchingCubes(values);
    auto triangles = Triangle::vectorsToTriangles(m.getTriangles());

    std::ofstream file;
    file.open("TEST.stl");
    file << m.toSTL();
    file.close();

    int intersectionCount;
    int nbRays = 500000;
    std::vector<Vector3> starts(nbRays), ends(nbRays);
    for (int i = 0; i < nbRays; i++) {
        starts[i] = Vector3::random(Vector3(), values.getDimensions());
        ends[i] = Vector3::random(Vector3(), values.getDimensions());
    }

    std::vector<std::string> names;
    std::vector<std::string> columns = {"nbTriangles", "parallel", "build", "eval", "total"};
    std::vector<std::vector<float>> timings;
    std::vector<float> timing;

    float timeBuild, timeEval;

    for (auto tris : {2, 4, 16, 64, 256, 1024}) {
//    for (auto tris : {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048}) {
        for (auto parallel : {false, true}) {
            timing = {float(tris), float(parallel)};
            names.push_back("SAH");
            std::cout << "Time for " << tris << " triangles per leaf, parallel = " << parallel << ", SAH     ... " << std::flush;
            BVHTree tree1;
            tree1.maxTrianglesPerLeaves = tris;
            tree1.useParallel = parallel;
            tree1.useSAH = true;
            tree1.useQuickSelect = false;
            timeBuild = timeIt([&]() {
                tree1.build(triangles);
            });
            timeEval = timeIt([&]() {
                intersectionCount = 0;
                #pragma omp parallel for
                for (int i = 0; i < nbRays; i++) {
                    tree1.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]);
//                    intersectionCount += int(tree1.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]).size());
                }
            });
            timings.push_back(vectorMerge(timing, {timeBuild, timeEval, timeBuild + timeEval}));
            std::cout << intersectionCount << " intersections in " << showTime(timeBuild + timeEval) << std::endl;




            names.push_back("Quick");
            timing = {float(tris), float(parallel)};
            std::cout << "Time for " << tris << " triangles per leaf, parallel = " << parallel << ", Quick   ... " << std::flush;
            BVHTree tree2;
            tree2.maxTrianglesPerLeaves = tris;
            tree2.useParallel = parallel;
            tree2.useSAH = false;
            tree2.useQuickSelect = true;
            timeBuild = timeIt([&]() {
                tree2.build(triangles);
            });
            timeEval = timeIt([&]() {
                intersectionCount = 0;
                #pragma omp parallel for
                for (int i = 0; i < nbRays; i++) {
                    tree2.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]);
                //                    intersectionCount += int(tree2.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]).size());
                }
            });
            timings.push_back(vectorMerge(timing, {timeBuild, timeEval, timeBuild + timeEval}));
            std::cout << intersectionCount << " intersections in " << showTime(timeBuild + timeEval) << std::endl;




            names.push_back("Midpoint");
            timing = {float(tris), float(parallel)};
            std::cout << "Time for " << tris << " triangles per leaf, parallel = " << parallel << ", Midpoint... " << std::flush;
            BVHTree tree3;
            tree3.maxTrianglesPerLeaves = tris;
            tree3.useParallel = parallel;
            tree3.useSAH = false;
            tree3.useQuickSelect = false;
            timeBuild = timeIt([&]() {
                tree3.build(triangles);
            });
            timeEval = timeIt([&]() {
                intersectionCount = 0;
                #pragma omp parallel for
                for (int i = 0; i < nbRays; i++) {
                    tree3.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]);
                //                    intersectionCount += int(tree3.getAllIntersectionsAndTrianglesIndices(starts[i], ends[i]).size());
                }
            });
            timings.push_back(vectorMerge(timing, {timeBuild, timeEval, timeBuild + timeEval}));
            std::cout << intersectionCount << " intersections in " << showTime(timeBuild + timeEval) << std::endl;

        }
    }
    Table results(timings, columns, names);
    std::cout << results.sortBy("total").displayTable() << std::endl;
    return 0;*/
/*
    ImplicitPrimitive* primA = new ImplicitPrimitive;
    ImplicitPrimitive* primB = new ImplicitPrimitive;
    ImplicitNaryOperator* nary = new ImplicitNaryOperator;
    ImplicitBinaryOperator* binary = new ImplicitBinaryOperator;
    ImplicitUnaryOperator* unary = new ImplicitUnaryOperator;

    nary->addChild(binary);
    binary->addChild(primA, 0);
    binary->addChild(unary, 1);
    unary->addChild(primA);

    primA->dimensions = Vector3(2, 2, 0);
    primA->position = Vector3(5, 5, 0) - primA->dimensions*.5f;

    primB->dimensions = Vector3(1, 1, 0);
    primB->position = Vector3(10, 10, 0) - primB->dimensions*.5f;

//    unary->scale(Vector3(3, 3, 1));
//    unary->translate(Vector3(1, 0, 0));
//    unary->rotate(0, 0, M_PI / 2.f);

//    std::cout << nary->getBBox() << std::endl;

    Vector3 pos, newPos;

    pos = Vector3(0, 0, 0);
    newPos = primA->getGlobalPositionOf(pos);
    std::cout << "Position " << pos << " becomes " << newPos << std::endl;

    pos = Vector3(3, 0, 0);
    newPos = primA->getGlobalPositionOf(pos);
    std::cout << "Position " << pos << " becomes " << newPos << std::endl;

    pos = Vector3(0, 3, 0);
    newPos = primA->getGlobalPositionOf(pos);
    std::cout << "Position " << pos << " becomes " << newPos << std::endl;

    pos = Vector3(3, 3, 0);
    newPos = primA->getGlobalPositionOf(pos);
    std::cout << "Position " << pos << " becomes " << newPos << std::endl;

    return 0;*/


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
    glFormat.setAlpha(true);
    QGLWidget widget(glFormat);
    widget.makeCurrent();

    const QOpenGLContext *context = GlobalsGL::context();

    qDebug() << "Context valid: " << context->isValid();
    qDebug() << "Really used OpenGl: " << context->format().majorVersion() << "." << context->format().minorVersion();
    qDebug() << "OpenGl information: VENDOR:       " << (const char*)glGetString(GL_VENDOR);
    qDebug() << "                    RENDERDER:    " << (const char*)glGetString(GL_RENDERER);
    qDebug() << "                    VERSION:      " << (const char*)glGetString(GL_VERSION);
    qDebug() << "                    GLSL VERSION: " << (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);

    EnvObject::readFile("saved_maps/primitives.json");
    EnvObject::sandDeposit = GridF(100, 100);


    ViewerInterface vi;
    vi.show();

    return app.exec();

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
}
