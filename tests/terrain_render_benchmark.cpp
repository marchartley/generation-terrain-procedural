#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "Interface/Interface.h"
#include "Interface/TerrainGenerationInterface.h"
#include <QApplication>
#include <QTemporaryDir>
#include <QTimer>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

struct RenderCase {
    const char* name;
    MapMode mode;
    SmoothingAlgorithm algorithm;
    float occlusion;
    bool shadows;
};

static void benchmark(ViewerInterface& window, const QString& imagePrefix)
{
    auto terrain = std::dynamic_pointer_cast<TerrainGenerationInterface>(
        window.viewer->interfaces.at("terraingeneration"));
    terrain->heightmap->heights = GridF(64, 64, 1);
    terrain->heightmap->heights.iterate([&](int x, int y, int) {
        const float dx = x - 31.5f, dy = y - 31.5f;
        terrain->heightmap->heights(x, y) = 4.f + 20.f * std::exp(-(dx * dx + dy * dy) / 500.f)
            + .6f * std::sin(x * .4f) * std::cos(y * .3f);
    });
    terrain->voxelGrid->setVoxelValues(GridF(64, 64, 32, -1.f));
    terrain->heightmapToAll();
    terrain->setWaterLevel(.25f);
    window.viewer->setSceneCenter(Vector3(32, 32, 12));
    window.viewer->setSceneRadius(64);
    window.viewer->showEntireScene();
    const std::vector<RenderCase> cases = {
        {"heightmap", GRID_MODE, MARCHING_CUBES, 0, false},
        {"voxels", VOXEL_MODE, MARCHING_CUBES, 0, false},
        {"layers", LAYER_MODE, MARCHING_CUBES, 0, false},
        {"implicit", IMPLICIT_MODE, MARCHING_CUBES, 0, false},
        {"heightmap_ao", GRID_MODE, MARCHING_CUBES, .5f, false},
        {"voxels_ao", VOXEL_MODE, MARCHING_CUBES, .5f, false},
        {"layers_raw", LAYER_MODE, NONE, 0, false},
        {"voxels_raw", VOXEL_MODE, NONE, 0, false},
        {"voxels_shadows", VOXEL_MODE, MARCHING_CUBES, 0, true}
    };
    for (const auto& test : cases) {
        window.viewer->setMapMode(test.mode);
        window.viewer->setSmoothingAlgorithm(test.algorithm);
        terrain->ambiantOcclusionFactor = test.occlusion;
        terrain->displayShadows = test.shadows;
        const auto warmupStart = std::chrono::steady_clock::now();
        int warmup = 0;
        do {
            window.viewer->grabFramebuffer();
            ++warmup;
        } while (warmup < 3 || (warmup < 256 && std::chrono::steady_clock::now() - warmupStart < std::chrono::milliseconds(750)));
        std::vector<double> times;
        QImage frame;
        for (int sample = 0; sample < 11; ++sample) {
            const auto start = std::chrono::steady_clock::now();
            frame = window.viewer->grabFramebuffer();
            if (frame.isNull()) throw std::runtime_error("Benchmark framebuffer capture failed");
            times.push_back(std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count());
        }
        std::sort(times.begin(), times.end());
        std::cout << "BENCH " << test.name << " median_ms=" << times[times.size() / 2]
                  << " min_ms=" << times.front() << '\n' << std::flush;
        if (!imagePrefix.isEmpty() && !frame.save(imagePrefix + "-" + test.name + ".png"))
            throw std::runtime_error("Could not save benchmark image");
    }
}

int main(int argc, char** argv)
{
    QSurfaceFormat format;
    format.setRenderableType(QSurfaceFormat::OpenGL);
    format.setVersion(4, 5);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    format.setDepthBufferSize(24);
    format.setSwapInterval(0);
    QSurfaceFormat::setDefaultFormat(format);
    QApplication app(argc, argv);
    QTemporaryDir directory;
    QImage heights(16, 16, QImage::Format_Grayscale8);
    heights.fill(120);
    const auto filename = directory.filePath("height.png");
    if (!heights.save(filename)) return 1;
    ViewerInterface window(filename.toStdString(), GRID_MODE);
    window.viewer->setStateFileName(directory.filePath("viewer.xml"));
    window.resize(900, 700);
    bool started = false;
    int result = 1;
    QObject::connect(window.viewer, &Viewer::viewerInitialized, &app, [&]() {
        if (started) return;
        started = true;
        QTimer::singleShot(100, &app, [&]() {
            try {
                benchmark(window, argc > 1 ? QString::fromLocal8Bit(argv[1]) : QString());
                result = 0;
            } catch (const std::exception& error) {
                std::cerr << "FAIL: " << error.what() << '\n';
            }
            app.quit();
        });
    });
    QTimer::singleShot(55000, &app, &QApplication::quit);
    window.show();
    app.exec();
    return result;
}
