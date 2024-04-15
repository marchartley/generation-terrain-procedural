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

#include "TerrainModification/ParticleErosion.h"
#include "Utils/Delaunay.h"
#include "Graph/RegularSimplicialComplex.h"
#include "Graph/TopoMap.h"
#include "Utils/Table.h"
#include "Graph/WaveFunctionCollapse.h"
#include "DataStructure/Kelvinlet.h"
#include "DataStructure/Image.h"
#include "EnvObject/ExpressionParser.h"
#include "Utils/HotreloadFile.h"

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


    //((d((island - pos)?{0 0 0}) min 20 max 10) - 10
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
     * Unit test: extracting the contours of a binary map
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
     * Unit test : checking that the naive "convolution" method is slower than the parallel one
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
    /*
     * Unit test : Generate OBJ files from a folder filled with STL files
    std::string initialPath = "erosionsTests";
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
     * Unit test : getting global and local positions of implicit patches relative to their parents
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
     * Unit test : extract curvature from a curve (?)
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
     * Unit test : tried to run a fluid simulation on GPU (?)
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
    /*
     * Unit test : extract euler angles between 2 vectors
    Vector3 A = Vector3(1, 0, 0);
    Vector3 B = Vector3(0.5, 1, 1).normalize();
    Vector3 UP = Vector3(0, 0, 1); // A.cross(B);
    Vector3 LEFT = UP.cross(A);

    Vector3 C = Vector3(
                    std::acos(B.y),
                    std::acos(B.z),
                    std::acos(B.x)
                ) * (180.f / M_PIf);

    std::cout << C << std::endl;
    return 0;*/
    /*
     * No clue what I'm trying to test here...
    AABBox bbox(Vector3(0, 0, 0), Vector3(2, 2, 1));

    Matrix3<float> M = Matrix3<float>({
                                          {0, 1},
                                          {1, 0}
                                          });
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
    return 0;
    */

    /*
     * Unit test : Merging 2 shapes (boolean OR)
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
    ShapeCurve AB = A.merge(B);
    Plotter::getInstance()->addPlot(A.closedPath(), "A", Qt::blue);
    Plotter::getInstance()->addPlot(B.closedPath(), "B", Qt::green);
    Plotter::getInstance()->addPlot(AB.closedPath(), "AB", Qt::red);
    Plotter::getInstance()->addScatter(AB.points, "");
    std::cout << A << std::endl;
    std::cout << B << std::endl;
    std::cout << AB << std::endl;
    return Plotter::getInstance()->exec();
    */

    /*
     * Unit test : Create a regular simplicial and fill some areas
    RegularSimplicialComplex grid(10, 10);
    grid.getNode(1, 3)->value = 0;
    grid.getNode(2, 3)->value = 0;
    grid.getNode(1, 4)->value = 0;
    grid.getNode(2, 4)->value = 0;
    grid.getNode(3, 4)->value = 0;
    grid.getNode(1, 5)->value = 0;
    grid.getNode(2, 5)->value = 0;
    grid.getNode(3, 5)->value = 0;

    grid.getNode(5, 3)->value = 1;
    grid.getNode(6, 3)->value = 1;
    grid.getNode(5, 4)->value = 1;
    grid.getNode(6, 4)->value = 1;
    grid.getNode(7, 4)->value = 1;
    grid.getNode(5, 5)->value = 1;
    grid.getNode(6, 5)->value = 1;
    grid.getNode(7, 5)->value = 1;
    grid.removeUnavailableLinks();
    grid.display();

    return Plotter::get()->exec();
    */

    /*
     * Unit test : create a topological map, extract a dual, and provide a geometry using forces
     * Soiler: the display looks all broken...
    auto g = CombinMap();
    g.addFace(5, {}, {100, 101, 102, 103, 104});
    g.addFace(4, {g.root->beta2, g.root->beta2->beta1}, {0});
    g.addFace(4, {g.root->beta2, g.root->beta2->beta1}, {1});
    g.debug();
    Graph<int> G = g.toGraph().forceDrivenPositioning();
    G.draw();
    Plotter::get()->exec();
    auto dual = g.dual(g.root->beta1->beta2);
    dual.debug();
    G = dual.toGraph().forceDrivenPositioning();
    G.draw();
    return Plotter::get()->exec();
    */
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
    /*
     * Unit test : Looking at Voronoi and Delaunay timings
     * Spoiler: Voronoi is in O(n) but since I encode the graphs as matrices (non-sparse), the complexity of Delaunay gets in O(n^2)...
    for (auto& nbPoints : {20, 100, 500, 2000, 10000}) {
        displayProcessTime("Voronoi + Delaunay on " + std::to_string(nbPoints) + " with 10 relaxations ...", [&]() {
            auto points = ShapeCurve({Vector3(0, 0, 0), Vector3(100, 0, 0), Vector3(100, 100, 0), Vector3(0, 100, 0)}).randomPointsInside(nbPoints);
            Voronoi voro(points, Vector3(100, 100));
            voro.solve(true, 10);
            Delaunay delaunay(voro);
        });
    }
    return 0;
    */

    /*
     * Unit test : Compute and display Voronoi and Delaunay graphs after multiple relaxations
     *
    Vector3 size(100, 100, 1);
    Voronoi voro(100, size.xy());

    for (int relax = 0; relax < 10; relax++) {
        voro.solve(true, relax);
        Delaunay delaunay(voro);

        GridV3 screen(size);

        screen.iterateParallel([&](const Vector3& pos) {
            for (auto& p : voro.pointset) {
                if ((pos - p).norm2() < 4.f) {
                    screen(pos).y = 1.f;
                }
            }

            for (auto& region : voro.areas) {
                float distance = region.estimateSqrDistanceFrom(pos);
                if (0 < distance && distance < 1.f) {
                    screen(pos).x = 1.f;
                }
            }

            for (auto& nodeA : delaunay.graph.nodes) {
                for (auto& [nodeB, weight] : nodeA->neighbors) {
                    BSpline path({nodeA->pos, nodeB->pos});
                    float distance = path.estimateSqrDistanceFrom(pos);
                    if (0 < distance && distance < 1.f) {
                        screen(pos).z = 1.f;
                    }
                }
            }
        });

        std::cout << "Relaxation " << relax << std::endl;
        Plotter::get()->addImage(screen);
        Plotter::get()->exec();
    }
    return Plotter::get()->exec();
    */

    /*
     * Unit test : High number of relaxation on Voronoi must look like a Poisson-Disk sampling
    Vector3 size(1000, 1000, 1);
    Voronoi voro(1000, size.xy());
    voro.solve(true, 100);
    Delaunay delaunay(voro);

    GridV3 screen(size);
    displayProcessTime("Comuting some distances for map of size " + std::to_string(size.x * size.y) + "... ", [&]() {
        screen.iterateParallel([&](const Vector3& pos) {
            for (auto& p : voro.pointset) {
                if ((pos - p).norm2() < 2000.f)
                    screen(pos).x += .1f;
            }
        });
    });
    Plotter::get()->addImage(screen);
    Plotter::get()->setNormalizedModeImage(true);
    return Plotter::get()->exec();*/

    /*
     * Unit test : Checking that sparse graphs are faster than matrix graphs
    std::vector<int> nodesID(500);
    for (int i = 0; i < nodesID.size(); i++) {
        nodesID[i] = i;
    }


    Graph g1(false);
    Graph g2(true);

    displayProcessTime("Generate List graph... ", [&]() {
        g1.addNodes(nodesID);
    });
    displayProcessTime("Random neighbors... ", [&]() {
        for (int i = 0; i < g1.nodes.size(); i++) {
            for (int j = 0; j < g1.nodes.size(); j++) {
                if (random_gen::generate() < 2.f) {
                    g1.addConnection(i, j);
                }
            }
        }
    });
    displayProcessTime("Generate Matrix graph... ", [&]() {
        g2.addNodes(nodesID);
    });
    displayProcessTime("Random neighbors... ", [&]() {
        for (int i = 0; i < g2.nodes.size(); i++) {
            for (int j = 0; j < g2.nodes.size(); j++) {
                if (random_gen::generate() < 0.2f) {
                    g2.addConnection(i, j);
                }
            }
        }
    });

    g1.circularLayout();
    return 0;
    */

    /*
     * Unit test : Shortest path algorithms
     * Graph source : https://upload.wikimedia.org/wikipedia/commons/thumb/2/29/DijkstraBis01.svg/330px-DijkstraBis01.svg.png
    Graph G(false);
    G.addNodes({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    G.addConnection(0, 1, 85);
    G.addConnection(0, 2, 217);
    G.addConnection(0, 4, 173);
    G.addConnection(1, 5, 80);
    G.addConnection(2, 6, 186);
    G.addConnection(2, 7, 103);
    G.addConnection(3, 7, 183);
    G.addConnection(5, 8, 250);
    G.addConnection(7, 9, 167);
    G.addConnection(4, 9, 502);
    G.addConnection(8, 9, 84);

    std::vector<float> expectedValues = {0, 85, 217, 503, 173, 165, 403, 320, 415, 487};

    float distAStar;
    std::vector<float> distancesDjikstra, distancesBellman;
    std::vector<int> predecessorAStar, predecessorDjikstra, predecessorBellman;
    GridF distancesFloyd, distancesFloyd2, distancesJohnson;
    GridI predecessorFloyd, predecessorFloyd2, predecessorJohnson;

    std::tie(distancesDjikstra, predecessorDjikstra) = Pathfinding::Djikstra(G, 0);
    std::tie(distancesBellman, predecessorBellman) = Pathfinding::BellmanFord(G, 0);
    std::tie(distancesFloyd, predecessorFloyd) = Pathfinding::FloydWarshall(G);
    std::tie(distancesFloyd2, predecessorFloyd2) = Pathfinding::FloydWarshallImproved(G);
    std::tie(distancesJohnson, predecessorJohnson) = Pathfinding::Johnson(G);


    std::vector<std::vector<float>> results(G.size(), std::vector<float>(7));
    std::vector<std::string> names(G.size());
    for (int iNode = 0; iNode < distancesDjikstra.size(); iNode++) {
        auto distDjikstra = distancesDjikstra[iNode];
        auto distBellman = distancesBellman[iNode];
        auto distFloyd = distancesFloyd(0, iNode);
        auto distFloyd2 = distancesFloyd2(0, iNode);
        auto distJohnson = distancesJohnson(0, iNode);
        auto node = iNode;
        std::tie(distAStar, predecessorAStar) = Pathfinding::AStar(G, 0, iNode);
//        std::cout << "A to " << char('A' + node) << " -> A* : " << distAStar << " -- Djikstra : " << distDjikstra << " -- Bellman : " << distBellman << " -- Floyd : " << distFloyd << " -- Floyd2 : " << distFloyd2 << " -- Johnson : " << distJohnson <<  " -- GT : " << expectedValues[node] << "\n";
        results[iNode] = {distAStar, distDjikstra, distBellman, distFloyd, distFloyd2, distJohnson, expectedValues[iNode]};
        names[iNode] = "A to " + std::string(1, 'A' + node);
    }

    Table tab(results, {"A*", "Djikstra", "Bellman", "Floyd", "Floyd (2)", "Johnson", "Ground truth"}, names);
    std::cout << tab.displayTable() << std::endl;
    std::cout << "End." << std::endl;
    return 0;*/
    /*
     * Unit test : Shortest path algorithms should give the same output.
    std::vector<Vector3> randomPositions(100);
    for (auto& p : randomPositions)
        p = Vector3::random(Vector3(), Vector3(100, 100));

    Voronoi voro(randomPositions);
    voro.solve(false, 1);
    Delaunay delaunay;
    Graph& G = delaunay.fromVoronoi(voro).graph;

    float distAStar;
    std::vector<float> distancesDjikstra, distancesBellman;
    std::vector<int> predecessorAStar, predecessorDjikstra, predecessorBellman;
    GridF distancesFloyd, distancesFloyd2, distancesJohnson;
    GridI predecessorFloyd, predecessorFloyd2, predecessorJohnson;

    std::tie(distancesDjikstra, predecessorDjikstra) = Pathfinding::Djikstra(G, 0);
    std::tie(distancesBellman, predecessorBellman) = Pathfinding::BellmanFord(G, 0);
    std::tie(distancesFloyd, predecessorFloyd) = Pathfinding::FloydWarshall(G);
    std::tie(distancesFloyd2, predecessorFloyd2) = Pathfinding::FloydWarshallImproved(G);
    std::tie(distancesJohnson, predecessorJohnson) = Pathfinding::Johnson(G);

    float longest = 0;
    int bestTarget;
    for (int iNode = 0; iNode < distancesDjikstra.size(); iNode++) {
        std::tie(distAStar, predecessorAStar) = Pathfinding::AStar(G, 0, iNode);
        if (longest < distAStar) {
            longest = distAStar;
            bestTarget = iNode;
        }
    }

    int last = bestTarget;
    std::vector<Vector3> pathAStar, pathDjikstra, pathBellman, pathFloyd1, pathFloyd2, pathJohnson;
//    for (int idx : Pathfinding::getPath(last, predecessorDjikstra)) {
//        pathDjikstra.push_back(G.nodes[idx]->pos);
//    }
    std::cout << "\nDjikstra: ";
    for (int idx : Pathfinding::getPath(last, predecessorDjikstra)) {
        std::cout << idx << " ";
        pathDjikstra.push_back(G.nodes[idx]->pos +Vector3(0, 0, 0));
    }
    std::cout << "\nBellman: ";
    for (int idx : Pathfinding::getPath(last, predecessorBellman)) {
        std::cout << idx << " ";
        pathBellman.push_back(G.nodes[idx]->pos + Vector3(1, 1));
    }
    std::cout << "\nFloyd  : ";
    for (int idx : Pathfinding::getPath(0, last, predecessorFloyd)) {
        std::cout << idx << " ";
        pathFloyd1.push_back(G.nodes[idx]->pos + Vector3(2, 2));
    }
    std::cout << "\nFloyd2 : ";
    for (int idx : Pathfinding::getPath(0, last, predecessorFloyd2)) {
        std::cout << idx << " ";
        pathFloyd2.push_back(G.nodes[idx]->pos + Vector3(3, 3));
    }
    std::cout << "\nJohnson : ";
    for (int idx : Pathfinding::getPath(0, last, predecessorJohnson)) {
        std::cout << idx << " ";
        pathJohnson.push_back(G.nodes[idx]->pos + Vector3(4, 4));
    }
    std::cout << std::endl;
    G.draw();
    Plotter::get()->addPlot(pathDjikstra, "Djikstra", Qt::red);
    Plotter::get()->addPlot(pathBellman, "Bellman", Qt::black);
    Plotter::get()->addPlot(pathFloyd1, "Floyd", Qt::blue);
    Plotter::get()->addPlot(pathFloyd2, "Floyd improved", Qt::green);
    Plotter::get()->addPlot(pathJohnson, "Johnson", Qt::cyan);
    Plotter::get()->exec();
    std::cout << "End." << std::endl;
    return 0;
    */

    /*
     * Unit test : Check speed of different shortest path algorithms
    std::vector<int> sizes = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120};
    std::vector<std::string> names(sizes.size());
    std::vector<std::vector<std::string>> results(sizes.size(), std::vector<std::string>(6));
    std::vector<bool> toIgnore(6, false);
    for (int _i = 0; _i < sizes.size(); _i++) {
        int graphSize = sizes[_i];
        names[_i] = std::to_string(graphSize);
        Graph G(false);

        for (int iNode = 0; iNode < graphSize; iNode++) {
            G.addNode(iNode);
        }
        for (int iNodeA = 0; iNodeA < graphSize; iNodeA++) {
            for (int iNodeB = 0; iNodeB < graphSize; iNodeB++) {
                if (random_gen::generate() < .1f) {
                    G.addConnection(iNodeA, iNodeB, random_gen::generate(2.f, 10.f));
                }
            }
        }

        float distAStar;
        std::vector<float> distancesDjikstra, distancesBellman;
        std::vector<int> predecessorAStar, predecessorDjikstra, predecessorBellman;
        GridF distancesFloyd, distancesFloyd2, distancesJohnson;
        GridI predecessorFloyd, predecessorFloyd2, predecessorJohnson;

        double AStarTime, djikstraTime, bellmanTime, floyd1Time, floyd2Time, johnsonTime;
        djikstraTime = timeIt([&]() {
            std::tie(distancesDjikstra, predecessorDjikstra) = Pathfinding::Djikstra(G, 0);
        }, (toIgnore[1] ? 0 : 5));
        bellmanTime = timeIt([&]() {
            std::tie(distancesBellman, predecessorBellman) = Pathfinding::BellmanFord(G, 0);
        }, (toIgnore[2] ? 0 : 5));
        floyd1Time = timeIt([&]() {
            std::tie(distancesFloyd, predecessorFloyd) = Pathfinding::FloydWarshall(G);
        }, (toIgnore[3] ? 0 : 2));
        floyd2Time = timeIt([&]() {
            std::tie(distancesFloyd2, predecessorFloyd2) = Pathfinding::FloydWarshallImproved(G);
        }, (toIgnore[4] ? 0 : 2));
        johnsonTime = timeIt([&]() {
            std::tie(distancesJohnson, predecessorJohnson) = Pathfinding::Johnson(G);
        }, (toIgnore[5] ? 0 : 1));

        for (int iNode = 0; iNode < graphSize; iNode++) {
            AStarTime += timeIt([&]() {
                std::tie(distAStar, predecessorAStar) = Pathfinding::AStar(G, 0, iNode);
            }, (toIgnore[0] ? 0 : 1));
        }
        toIgnore[0] = AStarTime > 1*1e9;
        toIgnore[1] = djikstraTime > 1*1e9;
        toIgnore[2] = bellmanTime > 1*1e9;
        toIgnore[3] = floyd1Time > 1*1e9;
        toIgnore[4] = floyd2Time > 1*1e9;
        toIgnore[5] = johnsonTime > 1*1e9;
        results[_i] = {showTime(AStarTime), showTime(djikstraTime), showTime(bellmanTime), showTime(floyd1Time/float(graphSize)), showTime(floyd2Time/float(graphSize)), showTime(johnsonTime/float(graphSize))};
        std::cout << "Iteration " << _i+1 << "/" << sizes.size() << std::endl;
    }

    Table tab(results, {"A*", "Djikstra", "Bellman", "Floyd", "Floyd (2)", "Johnson"}, names);
    std::cout << tab.displayTable() << std::endl;
    std::cout << "NB: time is the avg. computation time for shortest path from 1 node to all others" << std::endl;
    return 0;
    */

    /*
     * Unit test : COmputation of the curl on a 2D grid
     *
    GridV3 cyclone(200, 200, 1);
    cyclone.iterateParallel([&](const Vector3& p) {
        Vector3 pos = p - cyclone.getDimensions().xy() * .5f;
        float r = pos.norm();
        Vector3 dir = Vector3(-pos.y, pos.x).normalized(); // * 20.f * r;
        cyclone(p) = dir;
    });
    GridV3 curl;
    for (float radius : {1.f, 5.f, 10.f, 20.f}) {
        displayProcessTime("Curl with radius " + std::to_string(radius) + " : ", [&]() {
            curl = cyclone.curl(radius).meanSmooth(30, 30, 1, true);
        });
        Plotter::get()->addImage(curl);
        Plotter::get()->exec();
    }
    return 0;
    */

    /*
     * Unit test: Fast Gaussian filter on a Matrix3
     * Notes: Did not try using FFT as We should have to cast things and it's tiring
    GridF grid(100, 100, 1);
    GridF ones(100, 100, 1, 1.f);
    grid.paste(ones, grid.getDimensions() * .5f - ones.getDimensions() * .5f);

    for (auto sigma : {1.f, 5.f, 10.f, 50.f, 100.f}) {
        displayProcessTime("Gaussian for sigma=" + std::to_string(sigma), [&]() {
            grid.gaussianSmooth(sigma, false);
        });
    }
    Plotter::get()->addImage(grid.gaussianSmooth(50.f, true));
    Plotter::get()->exec();

    return 0;
    */

    /*
     * Not a unit test, just being bored...
    GridV3 res(1000, 1000, 1);

    auto recompute = [&](std::complex<float> c) {
        res.iterateParallel([&](const Vector3& p) {
            std::complex<float> z1 = std::complex<float>(p.x - res.sizeX * .5f, p.y - res.sizeY * .5f) / std::complex<float>(res.sizeX * .25f, res.sizeY * .25f);
            auto z = z1;
            int iter = 0;
            int maxIter = 200;
            for (iter = 0; iter < maxIter; iter++) {
                z = z * z + c; //z1;
                auto diff = z - z1;
                if (diff.real() * diff.real() + diff.imag() * diff.imag() > 5)
                    break;
            }
            res(p) = colorPalette(iter / float(maxIter), {Vector3(0, 0, 0), Vector3(1, 1, 1), Vector3(0, 0, 1)});
        });
        Plotter::get()->addImage(res);
        Plotter::get()->show();
    };


    QObject::connect(Plotter::get()->chartView, &ChartView::mouseMoved, [&](Vector3 p) {
        recompute({(p.x - .5f) * 2.f, (p.y - .5f) * 2.f});
    });
    recompute({0.f, 0.f});
    Plotter::get()->exec();
    return 0;
    */

    /*
     * Unit test: Resampling the paths on BSplines
    BSpline curve({
        Vector3(0, 0, 0),
        Vector3(1, 0, 0),
        Vector3(2, 0, 0),
        Vector3(10, 1, 0),
        Vector3(11, 1, 0),
        Vector3(12, 1, 0)
    });

    curve = curve.resamplePoints(200).resamplePoints(10);
    auto path = curve.points; // curve.getPath(100, true);
    Plotter::get()->addScatter(path);
    Plotter::get()->addPlot(path);
    return Plotter::get()->exec();
    */

    /*
    int nbCandidates = 100;
    Vector3 startPosition(0, 0, 0);
    Vector3 direction(1, -1, 0);
    float radius = 100;
    float angle = deg2rad(30);

    std::vector<Vector3> points(nbCandidates);
    float initialAngle = direction.getSignedAngleWith(Vector3(1, 0, 0));
    for (int i = 0; i < nbCandidates; i++) {
        float phi = interpolation::inv_linear(random_gen::generate(), -angle, angle);
        float r = std::sqrt(random_gen::generate()) * radius; // Use square root to svoid bias towards center of the disk

        points[i] = Vector3(r, std::sin(phi) * r).rotate(0, 0, initialAngle) + startPosition;
    }

    Plotter::get()->addScatter(points);
    return Plotter::get()->exec();
    */

    /*
     * Unit test: Matrix products
    {
        std::cout << "Test 1" << std::endl;
        Matrix A(std::vector<std::vector<float>>{{1, 2, 3}});
        Matrix B(std::vector<std::vector<float>>{
                     {1, 0, 0}, {0, 2, 0}, {0, 0, 3}
                 });
        std::cout << A << "\n*\n" << B << "\n=\n" << std::flush;
        std::cout << A.product(B) << std::endl;
    }
    {
        std::cout << "Test 2" << std::endl;
        Matrix A(std::vector<std::vector<float>>{{1}, {2}, {3}});
        Matrix B = A.transpose();
        std::cout << A << "\n*\n" << B << "\n=\n" << std::flush;
        std::cout << A.product(B) << std::endl;
    }
    {
        std::cout << "Test 3" << std::endl;
        Matrix A(std::vector<std::vector<float>>{
                     {1, 0, 0}, {0, 2, 0}, {0, 0, 3}
                 });
        Matrix B(std::vector<std::vector<float>>{{1}, {2}, {3}});
        std::cout << A << "\n*\n" << B << "\n=\n" << std::flush;
        std::cout << A.product(B) << std::endl;
    }

    return 0;
    */


    /*
     * Unit test: Kelvinlets on an image
    Vector3 size = Vector3(200, 200, 1);
    Vector3 center = size.xy() * .5f;

    GridV3 initialImage = Image::readFromFile("poster/profile.png").colorImage.resize(size) / 255.f;

    std::vector<Kelvinlet*> operations;

    auto f = [&](float r, float s, const Vector3& mousePos) {
//        TwistKelvinlet* t = new TwistKelvinlet;
//        TranslateKelvinlet* t = new TranslateKelvinlet;
//        ScaleKelvinlet* t = new ScaleKelvinlet;
        PinchKelvinlet* t = new PinchKelvinlet;
        Vector3 offset = (mousePos - center);
        t->pos = center;
//        center = mousePos;
//        t->force = Vector3(0, 0, -offset.getSignedAngleWith(Vector3(1, 0, 0)));
        t->force = offset * 10.f; //Vector3(0, s, 0);
        // t->force = Vector3(20, 0, 0);
//        t->scale = offset.norm();
        t->radialScale = r;
        std::cout << t->force << std::endl;

        for (auto& op : operations)
            delete op;
        operations.clear();
        operations.push_back(t);

        GridV3 img(size);
        displayProcessTime("Computing " + std::to_string(operations.size()) + " operations... ", [&]() {
            img.iterateParallel([&](const Vector3& p) {
                if (p == center) return;
//                float rEps = t->densityFunction((t->pos - p).norm());
//                img(p) = Vector3(1, 1, 1) * rEps;
                Vector3 pos = p;
                for (int i = operations.size() - 1; i >= 0; i--) {
                    const auto& t = operations[i];
                    Vector3 warp = t->evaluate(pos);
                    pos -= warp;
                }
//                pos = p - warp;
//                img(p) = warp;
    //            img(p) = pos - mousePos;
//                img(p) = pos.xy() / size;
//                img(p) = pos * Vector3(1, 1, 1) * std::sin(pos.y * .1f) * std::sin(pos.x * .1f);
                img(p) = initialImage.interpolate(pos);

            });
        });
        Plotter::get()->addImage(img);
        Plotter::get()->show();
    };
    QObject::connect(Plotter::get()->chartView, &ChartView::mouseMoved, [&](Vector3 p) {
        float r = p.x * 1.f * size.x;
        float s = p.y * 1.f * size.y;
        f(5.f, s, p * size);
    });
    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, [&](const Vector3& clickPos, Vector3 value) {
        center = clickPos;
    });

    f(10, 0, center + Vector3(0, 10, 0));

    return Plotter::get()->exec();
    */

    /*
     * Unit test: Kelvinlets on curves
     * Notes: I guess only "TranslateKelvinletCurve" has a meaning...
    Vector3 size = Vector3(200, 200, 1);
    Vector3 center = size.xy() * .5f;

    GridV3 initialImage = Image::readFromFile("poster/profile.png").colorImage.resize(size) / 255.f;

    std::vector<Kelvinlet*> operations;

    auto f = [&](float r, float s, const Vector3& mousePos) {
//        TwistKelvinletCurve* t = new TwistKelvinletCurve;
//        TranslateKelvinletCurve* t = new TranslateKelvinletCurve;
//        ScaleKelvinletCurve* t = new ScaleKelvinletCurve;
        PinchKelvinletCurve* t = new PinchKelvinletCurve;
        t->curve = BSpline({center, center + Vector3(0, size.y * .25f), mousePos});
        t->force = 3000.f;
        t->radialScale = r;

        for (auto& op : operations)
            delete op;
        operations.clear();
        operations.push_back(t);

        GridV3 img(size);
        displayProcessTime("Computing " + std::to_string(operations.size()) + " operations... ", [&]() {
            img.iterateParallel([&](const Vector3& p) {
                if (p == center) return;
                Vector3 pos = p;
                for (int i = operations.size() - 1; i >= 0; i--) {
                    const auto& t = operations[i];
                    Vector3 warp = t->evaluate(pos);
                    pos -= warp;
                }
//                img(p) = pos - p;
    //            img(p) = pos - mousePos;
//                img(p) = pos.xy() / size;
//                img(p) = pos * Vector3(1, 1, 1) * std::sin(pos.y * .1f) * std::sin(pos.x * .1f);
                img(p) = initialImage.interpolate(pos);

            });
        });
        Plotter::get()->addImage(img);
        Plotter::get()->show();
    };
    QObject::connect(Plotter::get()->chartView, &ChartView::mouseMoved, [&](Vector3 p) {
        float r = p.x * 1.f * size.x;
        float s = p.y * 1.f * size.y;
        f(5.f, s, p * size);
    });
    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, [&](const Vector3& clickPos, Vector3 value) {
        center = clickPos;
    });

    f(10, 0, center + Vector3(0, 10, 0));

    return Plotter::get()->exec();
    */

    /*
     * Not unit test: Displaying the divergence of the resulting vector fields after Kelvinlets
     *
    Vector3 size = Vector3(200, 200, 1);
    Vector3 center = size.xy() * .5f;

    GridV3 initialImage = Image::readFromFile("poster/profile.png").colorImage.resize(size) / 255.f;

    std::vector<Kelvinlet*> operations;

    auto f = [&](float r, float s, const Vector3& mousePos) {
//        TwistKelvinletCurve* t = new TwistKelvinletCurve;
//        TranslateKelvinletCurve* t = new TranslateKelvinletCurve;
//        ScaleKelvinletCurve* t = new ScaleKelvinletCurve;
        PinchKelvinletCurve* t = new PinchKelvinletCurve;
        t->curve = BSpline({center, center + Vector3(0, size.y * .25f), mousePos});
        t->force = -3000.f;
        t->radialScale = r;

        for (auto& op : operations)
            delete op;
        operations.clear();
        operations.push_back(t);

        GridV3 img(size);
        GridV3 distortion(size);
        displayProcessTime("Computing " + std::to_string(operations.size()) + " operations... ", [&]() {
            img.iterateParallel([&](const Vector3& p) {
                if (p == center) return;
                Vector3 pos = p;
                for (int i = operations.size() - 1; i >= 0; i--) {
                    const auto& t = operations[i];
                    Vector3 warp = t->evaluate(pos);
                    pos -= warp;
                }
                distortion(p) = pos - p;
//                img(p) = pos - p;
    //            img(p) = pos - mousePos;
//                img(p) = pos.xy() / size;
                img(p) = Vector3(1, 1, 1) * (std::sin(pos.y * .1f) * std::sin(pos.x * .1f) > 0.f ? 1.f : 0.f);
//                img(p) = initialImage.interpolate(pos);

            });
            GridF divergence = distortion.divergence().meanSmooth(3, 3, 1);
//            divergence.normalize();
            float divMin = divergence.min(), divMax = divergence.max();
            divergence.iterateParallel([&](size_t i){
                divergence[i] = ((divergence[i] / (divergence[i] < 0 ? std::abs(divMin) : divMax)) + 1.f) * .5f;
            });
//            std::cout << divMin << " " << divMax << " " << divergence.min() << " " << divergence.max() << std::endl;
            img.iterateParallel([&](size_t i) {
                img[i] *= colorPalette(divergence[i], {Vector3(1, 0, 0), Vector3(1, 1, 1), Vector3(0, 0, 1)});
            });
        });
        Plotter::get()->addImage(img);
        Plotter::get()->show();
    };
    QObject::connect(Plotter::get()->chartView, &ChartView::mouseMoved, [&](Vector3 p) {
        float r = p.x * 1.f * size.x;
        float s = p.y * 1.f * size.y;
        f(5.f, s, p * size);
    });
    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, [&](const Vector3& clickPos, Vector3 value) {
        center = clickPos;
    });

    f(10, 0, center + Vector3(0, 10, 0));

    return Plotter::get()->exec();
    */

    /*
    GridF img(100, 100, 1);

    Vector3 center = img.getDimensions().xy() * .5f;
    float radius = 10.f;
    ShapeCurve shape = ShapeCurve({
                                      Vector3(10, 10),
                                      Vector3(20, 50),
                                      Vector3(50, 70),
                                      Vector3(70, 70),
                                      Vector3(50, 10)
                                  });
//    ShapeCurve shape = ShapeCurve::circle(30, center, 10);
    img.iterateParallel([&](const Vector3& p) {
        Vector3 pos = p;
        pos += Vector3(random_gen::generate_perlin(p.x * 2.f, p.y * 2.f, 0), random_gen::generate_perlin(p.x * 2.f, p.y * 2.f, 100), 0) * 100.f;
        float val = (random_gen::generate_perlin(pos.x, pos.y, pos.z) + 1.f) * .5f;

        float signDistance = BSpline(shape.closedPath()).estimateDistanceFrom(p, true); // * (shape.containsXY(p, false) ? 1.f : -1.f);
        float falloff = std::clamp((radius - signDistance) / radius, 0.f, 1.f); // shape.containsXY(p, false) ? 1.f : std::max(0.f, (radius - shape.estimateDistanceFrom(p)) / radius);
//        float falloff = interpolation::smooth((shape.containsXY(p) ? 1.f : std::max(0.f, radius - shape.estimateDistanceFrom(p) / radius)));
//        float falloff = interpolation::smooth(std::max(0.f, (radius - (p - center).norm()) / radius));
//        img(p) = val * falloff;
        img(p) = falloff;
    });
    Heightmap heightmap;
    heightmap.heights = img * 100.f;
    std::ofstream off("test.stl");
    off << heightmap.getGeometry().toSTL();
    off.close();
//    return 0;
    Plotter::get()->addImage(img);
    return Plotter::get()->exec();
    */


    /*Matrix mat({
                   {1, 2, 3, 11},
                   {4, 5, 6, 22},
                   {7, 8, 9, 33}
               });
    GridF mat2(mat);

    std::cout << mat << "\n\n" << mat2.displayValues() << std::endl;

    std::cout << "As an array: " << std::flush;
    float* values = mat.toArray();
    for (int i = 0; i < mat.rows() * mat.cols(); i++) {
        std::cout << values[i] << " ";
    }
    delete[] values;
    std::cout << std::endl;
    return 0;*/

    /*
    HotreloadFile file("/home/marc/generation-terrain-procedural/saved_maps/envMaterials.json");

    file.onChange([&](std::string content) {
        std::cout << content << std::endl;
    });
    while (true) {
        file.check();
        sleep(1000);
    }

    return 0;
    */
    /*
    Vector3 size(300, 300, 1);
    GridV3 img(size, Vector3(0, 1, 0));
    GridV3 initialImage = Image::readFromFile("poster/profile.png").colorImage.resize(size) / 255.f;


    float radialScale = 10.f;
    float force = 100.f;
    Vector3 center = size.xy() * .5f;

    TwistKelvinlet k1;
    k1.pos = center + Vector3(-radialScale, 0, 0);
    k1.force = Vector3(0, 0, force);
    k1.radialScale = radialScale;

    TwistKelvinlet k2;
    k2.pos = center + Vector3(radialScale, 0, 0);
    k2.force = Vector3(0, 0, -force);
    k2.radialScale = radialScale;

    img.iterateParallel([&](const Vector3& p) {
        img(p) = initialImage(p + k1.evaluate(p + k2.evaluate(p))); // initialImage(p + k1.evaluate(p) + k2.evaluate(p));
    });

    Plotter::get()->addImage(img);
    return Plotter::get()->exec();
    */

    /*
     * Unit test: Stringify and loading Matrix3 (float and Vector3 only)
     *
    GridF testF(4, 4, 4);
    testF.iterateParallel([&](size_t i) {
        testF[i] = i;
    });
//    auto resultF = loadGridF(stringifyGridF(testF));
//    std::cout << "Initial values:\n" << testF.displayValues() << "\nResulting values:\n" << resultF.displayValues() << "\nDifference:\n" << (resultF - testF).displayValues() << std::endl;
    std::cout << "Check floats binary: " << ((testF - loadGridF(stringifyGridF(testF, true), true)).sum() == 0 ? "OK" : "Error") << "\nCheck float plain: " << ((testF - loadGridF(stringifyGridF(testF, false), false)).sum() == 0 ? "OK" : "Error") << std::endl;

    GridV3 testV3(4, 4, 4);
    testV3.iterateParallel([&](size_t i) {
        testV3[i] = Vector3(i, i, i);
    });
//    auto resultV3 = loadGridV3(stringifyGridV3(testV3));
//    std::cout << "Initial values:\n" << testV3.displayValues() << "\nResulting values:\n" << resultV3.displayValues() << "\nDifference:\n" << (resultV3 - testV3).displayValues() << std::endl;
    std::cout << "Check vec3 binary: " << ((testV3 - loadGridV3(stringifyGridV3(testV3, true), true)).sum() == Vector3() ? "OK" : "Error") << "\nCheck vec3 plain: " << ((testV3 - loadGridV3(stringifyGridV3(testV3, false), false)).sum() == Vector3() ? "OK" : "Error") << std::endl;
    return 0;*/

    /*
    GridF testF(4, 4, 1);
    testF.iterateParallel([&](size_t i) {
        testF[i] = 0;
    });
    testF(Vector3(2, 1, 0)) = 3.f;
    auto resultF = loadGridF(stringifyGridF(testF));

    std::cout << testF.displayValues() << "\n\n" << resultF.displayValues() << std::endl;
    return 0;*/


    EnvObject::readEnvMaterialsFile("saved_maps/envMaterials.json");
    EnvObject::readEnvObjectsFile("saved_maps/primitives.json");

    ViewerInterface vi;
    vi.show();

    return app.exec();
}
