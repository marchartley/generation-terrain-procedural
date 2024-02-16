#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION

#undef interface

//#define OPENVDB_DLL
#include "Utils/Globals.h"
//#include "sim-fluid-loganzartman/Game.hpp"
#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <qapplication.h>
#include <iostream>
#include "Interface/Interface.h"
#include "EnvObject/EnvObject.h"
#include "FluidSimulation/OpenFoamParser.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>

#include "TerrainModification/particleErosion.h"
using namespace std;

std::map<std::string, std::string> getAllEnvironmentVariables() {
    std::string cmd = "/bin/bash -i -c 'source ~/.bashrc && env'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Unable to run command." << std::endl;
        return {};
    }

    std::string currentKey = "";

    char buffer[1024];
    std::map<std::string, std::string> results;
    std::string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result = std::string(buffer);
        if (!result.empty() && result[result.length() - 1] == '\n') {
            result.erase(result.length() - 1);
        }
        size_t pos = result.find("=");
        if (pos != result.npos) {
            currentKey = result.substr(0, pos);
            results[currentKey] = result.substr(pos + 1);

        } else {
            results[currentKey] += result;
        }
    }
    pclose(pipe);
    return results;
}

int main(int argc, char *argv[])
{
    auto allVars = getAllEnvironmentVariables();
    for (auto& [key, val] : allVars) {
        auto lowerKey = toLower(key);
//        if (lowerKey == "path" || lowerKey == "ld_library_path" || lowerKey.find("foam") != lowerKey.npos) {
            std::string s_cmd = key + "=" + val;
            auto cmd = const_cast<char*>(s_cmd.c_str());
            setenv(key.c_str(), val.c_str(), 1);
//        }
    }
    //OpenFoamParser::createSimulationFile("OpenFOAM/simple", GridF());
    //return 0;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
    QApplication app(argc, argv);


    QGLFormat glFormat;
    glFormat.setVersion(4, 5);
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


    /*
     * Unit test: BSpline with duplicate values such that it doesn't return any NaN
    GridV3 vals(200, 200);
    BSpline spline;
    // Without duplicate positions
    spline = BSpline({Vector3(0, 0, 0), Vector3(100, 50, 0), Vector3(100, 25, 0), Vector3(150, 100, 0), Vector3(200, 50, 0)});
    for (int i = 0; i <= 1000; i++) {
        float t = float(i) / 1000.f;
        auto p = spline.getPoint(t);
        vals(p).x = 1;
        if (p.x != p.x) {
            std::cerr << "Without duplicates, NaN found for t = " << t << std::endl;
        }
    }
    // With duplicate positions
    spline = BSpline({Vector3(0, 0, 0), Vector3(100, 50, 0), Vector3(100, 25, 0), Vector3(100, 25, 0), Vector3(150, 100, 0), Vector3(200, 50, 0)});
    for (int i = 0; i <= 1000; i++) {
        float t = float(i) / 1000.f;
        auto p = spline.getPoint(t);
        vals(p).y = 1;
        if (p.x != p.x) {
            std::cerr << "With duplicates, NaN found for t = " << t << std::endl;
        }
    }
    Plotter::get()->addImage(vals);
    return Plotter::get()->exec();
    */

    /*
     * Unit test: testing the color palette
    GridV3 colors(1000, 1);
    colors.iterate([&] (const Vector3& p) {
        colors(p) = colorPalette(p.x / colors.sizeX, {Vector3(1, 0, 0), Vector3(0, 0, 1), Vector3(0, 0, 1), Vector3(1, 1, 0)});
    });
    Plotter::get()->addImage(colors);
    return Plotter::get()->exec();
    */
    /*
     * Unit test: divergence of a vector field
     * There should be the value "2" everywhere, except on the right and bottom borders (and 1.5 on top and left borders)
    GridV3 vels(50, 50);
    vels.iterateParallel([&](const Vector3& p) {
        vels(p) = p;
    });
    GridF divergence = vels.divergence();

    Plotter::get()->addImage(divergence);
    Plotter::get()->exec();
    return 0;
    */
    /*
     * Unit test: adding value in cells (with interpolation)
    GridF values(10, 10);
    values.addValueAt(1.f, Vector3(1.5f, 1.5f));
    std::cout << values.displayValues() << std::endl;
    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, Plotter::get(), [&](const Vector3& clickPos, Vector3 value) {
        std::cout << "Adding at " << clickPos << std::endl;
        values.reset();
        values.addValueAt(1.f, clickPos);
        std::cout << values.displayValues() << std::endl;
        Plotter::get()->addImage(values);
        Plotter::get()->draw();
    });

    Plotter::get()->addImage(values);
    Plotter::get()->exec();
    return 0;
    */
    /*
     * Unit test: Save and load vector fields
    GridV3 input(100, 100, 100);
    Vector3 center(50, 50, 10);
    input.iterateParallel([&](const Vector3& p) {
        input(p).setValid((p - center).norm2() < 30*30);
    });

    VectorFieldDataFile data(input);
    data.write("test_vectors.raw");

    Plotter::get()->addImage(input);
    Plotter::get()->exec();

    VectorFieldDataFile outData;
    outData.load("test_vectors.raw");
    GridV3 output = outData.data;

    output.iterateParallel([&](size_t i) {
        output[i].x = (output[i].isValid() ? 1.f : 0.f);
    });
    Plotter::get()->addImage(output);
    Plotter::get()->exec();
    return 0;
    */
    /*
     * Unit tests: FFT and iFFT on arbitrary sizes (limited to 3 dimensions max)
//    GridF vals(128, 1);
//    vals.iterateParallel([&](const Vector3& p) {
//        vals(p) = std::sin(p.x * .2f) + std::sin(p.x * .1f);
//    });
//    GridF vals(1, 128);
//    vals.iterateParallel([&](const Vector3& p) {
//        vals(p) = std::sin(p.y * .2f) + std::sin(p.y * .1f);
//    });
//    GridF vals(1, 1, 128);
//    vals.iterateParallel([&](const Vector3& p) {
//        vals(p) = std::sin(p.z * .2f) + std::sin(p.z * .1f);
//    });
    GridF vals(100, 200);
    vals.iterateParallel([&](const Vector3& p) {
        vals(p) = std::sin(p.x * .2f) + std::sin(p.y * .1f);
    });
    for (size_t i = 0; i < vals.size(); i++)
        std::cout << vals[i] << "\n";
    Plotter::get()->addImage(vals);
    Plotter::get()->exec();
    auto FFT = vals.FFT();
    auto iFFT = FFT.iFFT();

    std::cout << " ----------------- FFT -------------------- " << std::endl;
    for (size_t i = 0; i < FFT.size(); i++)
        std::cout << FFT[i] << "\n";

    std::cout << " ---------------- iFFT -------------------- " << std::endl;
    for (size_t i = 0; i < iFFT.size(); i++)
        std::cout << iFFT[i] << "\n";
    std::cout << std::endl;

    GridF absolutes(FFT.getDimensions());
    absolutes.iterateParallel([&](size_t i) {
        absolutes[i] = std::abs(FFT[i]);
    });
    Plotter::get()->addImage(absolutes);
    Plotter::get()->exec();

    absolutes.iterateParallel([&](size_t i) {
        absolutes[i] = iFFT[i].real();
    });
    Plotter::get()->addImage(absolutes);
    Plotter::get()->exec();

    GridF diff = absolutes - GridF(absolutes.getDimensions()).paste(vals);
    Plotter::get()->addImage(diff);
    Plotter::get()->exec();
    return 0;*/


    /*
    // Unit test : Slerp between vectors
    Vector3 A(.5f, .5f, 0);
    Vector3 B(-.5f, .5f, 0);

    GridV3 check(100, 100, 1);
    int nbSamples = 500;
    float timeForLerp = 0.f;
    float timeForSlerp = 0.f;
    for (int i = 0; i < nbSamples; i++) {
        float t = float(i) / float(nbSamples);
        Vector3 C_slerp, C_lerp;
        timeForSlerp += timeIt([&]() { for (int _ = 0; _ < 100000; _++) C_slerp = Vector3::slerp(t, A, B); });
        timeForLerp += timeIt([&]() { for (int _ = 0; _ < 100000; _++) C_lerp = Vector3::lerp(t, A, B); });

        check(C_slerp * 50 + Vector3(50, 50)) = Vector3(t, 0, 0);
        check(C_lerp * 50 + Vector3(50, 50)) = Vector3(0, t, 0);
    }
    std::cout << "Lerp : " << showTime(timeForLerp) << "\nSlerp: " << showTime(timeForSlerp) << std::endl;
    Plotter::get()->addImage(check);
    Plotter::get()->exec();
    return 0;
    */
    /*
    // Unit test : get angle between vectors
    Vector3 A(5, 0, 0);

    Vector3 B;

    GridF check(100, 100, 1);
    int nbSamples = 500;
    for (int i = 0; i < nbSamples; i++) {
        float t = float(i) / float(nbSamples);
        float angle = t * 2.f * M_PI;
        B = Vector3(1, 0, 0).rotate(0, 0, angle);
        std::cout << rad2deg(A.getAngleWith(B) * sign(A.cross(B).z)) << " deg from " << A << " to " << B << std::endl;
        check(B * 50 + Vector3(50, 50)) = rad2deg(A.getAngleWith(B) * sign(A.cross(B).z));
    }
    Plotter::get()->addImage(check);
    Plotter::get()->exec();
    return 0;
    */
    /*
    GridI grid({
                   {0, 1, 1, 0, 0, 0, 0, 0, 0},
                   {0, 1, 1, 0, 0, 0, 0, 0, 0},
                   {0, 1, 1, 0, 0, 1, 0, 0, 0},
                   {0, 0, 1, 1, 1, 1, 0, 0, 0},
                   {0, 0, 1, 1, 1, 1, 1, 0, 0},
                   {0, 0, 1, 1, 1, 1, 1, 1, 0},
                   {0, 0, 1, 1, 1, 1, 1, 1, 0},
                   {0, 0, 0, 1, 1, 1, 1, 1, 0},
                   {0, 0, 0, 0, 1, 1, 1, 0, 0},
                   {0, 0, 0, 0, 0, 0, 0, 0, 0}
               });

    std::vector<ShapeCurve> contours = grid.findContoursAsCurves();
    grid = GridI(grid.getDimensions() * Vector3(2.f, 2.f, 1.f));
    for (size_t iCurve = 0; iCurve < contours.size(); iCurve++) {
        auto& contour = contours[iCurve];
        contour = contour.simplifyByRamerDouglasPeucker(.2f);
        for (size_t i = 0; i < contour.size(); i++) {
            grid(contour[i] * 2.f) = iCurve + 1;
        }
    }
    std::cout << grid.displayValues() << std::endl;

    for (size_t iCurve = 0; iCurve < contours.size(); iCurve++) {
        auto& contour = contours[iCurve].scale(2.f);
        std::cout << contour << std::endl;
    }
    return 0;
    */
    /*
    Vector3 dims(1000, 1000, 1);
    GridF mask = GridF::random(dims);
    GridF grid = GridF::random(dims);
    GridF res(dims);
    GridF res2(dims);

    float timeParallel = timeIt([&]() {
        for (int iter = 0; iter < 1000; iter++) {
            FastNoiseLite noise(iter);
            res = GridF::fbmNoise1D(noise, dims.x, dims.y, dims.z);
//#pragma omp parallel for
//            for (size_t i = 0; i < grid.size(); i++) {
//                res[i] = mask[i] + grid[i];
//            }
//            grid.convolution(mask);
//            grid.iterateParallel([&](float x, float y, float z) {
//                grid(x, y, z) = x + y + z;
//            });
        }
    });
    float timeSeq = timeIt([&]() {
        for (int iter = 0; iter < 1000; iter++) {
//            res2 = mask + grid;
//            for (size_t i = 0; i < grid.size(); i++) {
//                res[i] = mask[i] + grid[i];
//            }
//            grid.convolutionSlow(mask);
//            grid.iterate([&](float x, float y, float z) {
//                grid(x, y, z) = x + y + z;
//            });
        }
    });
//    std::cout << grid.displayAsPlot() << std::endl;
    std::cout << "Parallel: " << showTime(timeParallel) << "\nSequential: " << showTime(timeSeq) << std::endl;
    std::cout << "Check: " << (res == res2 ? "OK" : "Not OK!!!") << std::endl;
    return 0;
    */
    /*std::string initialPath = "erosionsTests";
    std::string destFolder = "asOBJ";
    std::string finalPath = initialPath + "/" + destFolder;
    std::vector<std::string> filenames;
    QDirIterator it(QString::fromStdString(initialPath), QDir::Files, QDirIterator::Subdirectories);
    while (it.hasNext()) {
        QString dir = it.next();
        filenames.push_back(dir.toStdString());
    }
    size_t nbFiles = filenames.size();

    if (!checkPathExists(finalPath)) {
        makedir(finalPath);
    }
//    std::vector<Mesh> meshes(nbFiles);
    std::cout << showTime(timeIt([&]() {
//    #pragma omp parallel for
        for (size_t i = 0; i < nbFiles; i++) {
            std::string dir = filenames[i];
            auto path = split(dir, "/");
            std::string previousName = path.back();
            std::string basename = split(previousName, ".")[0] + ".obj";
            std::string newFilename = finalPath + "/" + basename;

    //        std::cout << dir << " --> " << newFilename << std::endl;
            std::ofstream file(newFilename);
            file << Mesh().fromStl(dir).toOBJ();
            file.close();
        }
    })) << std::endl;
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
/*
    BSpline spline({
                       Vector3(0, 10, 0),
                       Vector3(50, 1, 0),
                       Vector3(30, 25, 0),
                       Vector3(50, 50, 0),
                       Vector3(0, 50, 0),
                       Vector3(0, 0, 0),
                   });
//    BSpline spline;
//    float r = 20;
//    for (int i = 0; i < 50; i++) {
//        float t = float(i) / 48.f;
//        float a = t * M_PI * 2.f;
//        spline.points.push_back(Vector3(
//                                    std::cos(a) * r + 50,
//                                    std::sin(a) * r + 50,
//                                    0
//                                    ));
//    }
//    spline.resamplePoints();
    GridF grid(105, 105, 1, 0);
    grid.iterateParallel([&](const Vector3& pos) {
        if (spline.estimateDistanceFrom(pos) < 10.f){
            float x = spline.estimateClosestTime(pos);
            float curve = spline.getCurvature(x);
            grid(pos) = curve;
        }
    });

    for (int i = 0; i < 200; i++) {
        float t = float(i) / 200.f;

        std::cout << spline.getCurvature(t) << " " << spline.getDerivative(t) << " " << spline.getSecondDerivative(t) << std::endl;
    }
    std::cout << grid.displayAsPlot() << "\n" << spline << std::endl;
    Plotter::get()->addImage(grid);
    Plotter::get()->exec();
    return 0;
    */
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
    std::cout << showTime(timeIt([]() {
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
                std::cout << "/n";
            }
            std::cout << "/n" << count << "/n" << std::endl;
        }
    })) << std::endl;

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

    /*
     * Unit test: smoothing 1D data
    GridF data(20, 1);
    for (int i = 0; i < data.size(); i++) {
        data[i] = i;
    }
    data[19] = 0;
    for (int _ = 0; _ < 20; _++) {
        data = data.meanSmooth(3, 3, 3, true);
        Plotter::get()->addImage(data);
        Plotter::get()->exec();
    }
    std::cout << data.displayValues() << std::endl;
    return 0;
    */
    /*
     * Unit test: median blur on 1D data
     *
    GridF data = GridF::fromImageBW("saved_maps/lena.png"); // = GridF(100, 1)
//    for (int i = 0; i < data.size(); i++) {
//        data[i] = std::sin(i / 10.f);
//    }
//    data[19] = 0;
    for (int _ = 0; _ < 20; _++) {
        Plotter::get()->addImage(data);
        Plotter::get()->exec();
        data = data.medianBlur(3, 3, 3, true);
    }
    return 0;
    */

    EnvObject::readFile("saved_maps/primitives.json");

    ViewerInterface vi;
    vi.show();

    return app.exec();
}
