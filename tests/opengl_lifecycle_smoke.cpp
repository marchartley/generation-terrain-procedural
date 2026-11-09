#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "Interface/Interface.h"
#include "Interface/TerrainGenerationInterface.h"
#include "Graphics/ComputeShader.h"
#include <QTemporaryDir>
#include <QTimer>
#include <QFile>
#include <array>
#include <stdexcept>
#include <set>
#include <future>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

int main(int argc, char** argv)
{
    // Construction must work before even QApplication or a context exists.
    Shader early("src/Shaders/no_shader.vert", "src/Shaders/no_shader.frag");
    require(early.programID == 0, "Shader constructor performed GPU work");
    GizmoRenderer earlyGizmo; // Must also be safe to construct before a context exists.
    std::vector<std::future<std::unique_ptr<GizmoRenderer>>> producers;
    for (int i = 0; i < 8; ++i) {
        producers.push_back(std::async(std::launch::async, []() {
            // Concurrent construction/destruction must not mutate global mesh lists.
            for (int n = 0; n < 100; ++n) {
                GizmoRenderer temporary;
            }
            return std::make_unique<GizmoRenderer>();
        }));
    }
    std::vector<std::unique_ptr<GizmoRenderer>> preparedGizmos;
    for (auto& producer : producers) preparedGizmos.push_back(producer.get());
    QSurfaceFormat format;
    format.setRenderableType(QSurfaceFormat::OpenGL);
    format.setVersion(4, 5);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    format.setDepthBufferSize(24);
    QSurfaceFormat::setDefaultFormat(format);
    qputenv("TERRAIN_GL_DEBUG", "1");
    QApplication app(argc, argv);
    QTemporaryDir directory;
    QImage heights(32, 32, QImage::Format_Grayscale8);
    for (int y = 0; y < heights.height(); ++y) {
        auto* row = heights.scanLine(y);
        for (int x = 0; x < heights.width(); ++x) {
            const float dx = x - 16.f, dy = y - 16.f;
            row[x] = static_cast<uchar>(40.f + 180.f * std::exp(-(dx*dx + dy*dy) / 120.f));
        }
    }
    const auto terrainFile = directory.filePath("height.png");
    require(heights.save(terrainFile), "Could not create test heightmap");

    auto* window = new ViewerInterface(terrainFile.toStdString(), GRID_MODE);
    window->viewer->setStateFileName(directory.filePath("viewer.xml"));
    window->resize(800, 600);
    int result = 1;
    bool ran = false;
    bool scheduled = false;
    QObject::connect(window->viewer, &Viewer::viewerInitialized, &app, [&]() {
        if (scheduled) return;
        scheduled = true;
        // Exercise a GUI callback after initializeGL has returned.
        QTimer::singleShot(500, &app, [&]() {
            try {
                require(GlobalsGL::isReady(), "Viewer never became ready");
                window->viewer->doneCurrent();
                require(!QOpenGLContext::currentContext(), "Test requires no current context");
                early.setVector("color", Vector4(1.f, 0.f, 0.f, 1.f));
                require(!QOpenGLContext::currentContext(), "Shader scope did not restore context");
                {
                    GlobalsGL::ContextScope scope;
                    require(QOpenGLContext::currentContext() == window->viewer->context(),
                            "Scope acquired the wrong context");
                    GLint linked = GL_FALSE;
                    GlobalsGL::f()->glGetProgramiv(early.programID, GL_LINK_STATUS, &linked);
                    require(linked == GL_TRUE, "Deferred shader failed to link");
                }
                {
                    GlobalsGL::ContextScope scope;
                    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
                }
                early.setTexture3D("signedFieldRegression", 10, GridF(2, 2, 2, -.75f));
                {
                    GlobalsGL::ContextScope scope;
                    GLint activeUnit = 0, oldBinding = 0;
                    glGetIntegerv(GL_ACTIVE_TEXTURE, &activeUnit);
                    require(activeUnit == GL_TEXTURE3, "Texture upload changed the active unit");
                    glGetIntegerv(GL_TEXTURE_BINDING_3D, &oldBinding);
                    glBindTexture(GL_TEXTURE_3D, early.textureSlotIndices.at(GL_TEXTURE10));
                    float samples[8] = {};
                    glGetTexImage(GL_TEXTURE_3D, 0, GL_ALPHA, GL_FLOAT, samples);
                    for (float value : samples)
                        require(std::abs(value + .75f) < 1e-6f, "Texture upload clamped a negative field value");
                    glBindTexture(GL_TEXTURE_3D, oldBinding);
                    GlobalsGL::f()->glActiveTexture(GL_TEXTURE0);
                }
                const auto uploads = early.volumeUploadCount;
                const auto allocations = early.volumeAllocationCount;
                early.setTexture3D("signedFieldRegression", 10, GridF(2, 2, 2, -.75f));
                require(early.volumeUploadCount == uploads, "Unchanged volume was uploaded again");
                early.setTexture3D("signedFieldRegression", 10, GridF(2, 2, 2, -.5f));
                require(early.volumeUploadCount == uploads + 1 && early.volumeAllocationCount == allocations,
                        "Same-size volume edit reallocated its texture");
                early.setTexture3D("signedFieldRegression", 10, GridF(3, 2, 2, -.5f));
                require(early.volumeAllocationCount == allocations + 1, "Resized volume kept old texture storage");
                {
                    GlobalsGL::ContextScope scope;
                    glBindTexture(GL_TEXTURE_3D, early.textureSlotIndices.at(GL_TEXTURE10));
                    GLint width = 0;
                    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &width);
                    require(width == 3, "Resized volume has wrong GPU dimensions");
                    float samples[12] = {};
                    glGetTexImage(GL_TEXTURE_3D, 0, GL_ALPHA, GL_FLOAT, samples);
                    for (float value : samples) require(value == -.5f, "Edited volume upload was stale");
                }
                // Direct-state uniform updates must target their own program,
                // without changing whichever program the renderer has bound.
                Shader other("src/Shaders/no_shader.vert", "src/Shaders/no_shader.frag");
                other.use();
                early.setVector("color", Vector4(.2f, .4f, .6f, 1.f));
                {
                    GlobalsGL::ContextScope scope;
                    GLint current = 0;
                    glGetIntegerv(GL_CURRENT_PROGRAM, &current);
                    require(GLuint(current) == other.programID, "Uniform update switched the current program");
                    float color[4] = {};
                    GlobalsGL::f()->glGetUniformfv(early.programID, early.uniformLocation("color"), color);
                    require(std::abs(color[1] - .4f) < 1e-6f, "Uniform cache skipped a changed value");
                }
                early.initializeGL();
                early.setVector("color", Vector4(.2f, .4f, .6f, 1.f));
                {
                    GlobalsGL::ContextScope scope;
                    float color[4] = {};
                    GlobalsGL::f()->glGetUniformfv(early.programID, early.uniformLocation("color"), color);
                    require(std::abs(color[2] - .6f) < 1e-6f, "Uniform cache survived program relinking");
                }
                const std::vector<std::array<std::string, 3>> shaderFamilies = {
                    {"MarchingCubes.vert", "MarchingCubes.frag", "MarchingCubes.geom"},
                    {"MarchingCubes.vert", "no_shader.frag", "MarchingCubes.geom"},
                    {"MarchingCubes.vert", "TerrainChanges.frag", "MarchingCubes.geom"},
                    {"voxels.vert", "voxels.frag", ""},
                    {"MarchingCubes.vert", "MarchingCubes.frag", "grid.geom"},
                    {"layer_based.vert", "MarchingCubes.frag", "layer_based.geom"},
                    {"rockShader.vert", "rockShader.frag", ""},
                    {"meshInstancesShader.vert", "meshInstancesShader.frag", ""},
                    {"gouraud.vert", "gouraud.frag", ""},
                    {"gouraud_120.vert", "gouraud_120.frag", ""},
                    {"grid.vert", "grid.frag", ""},
                    {"no_shader.vert", "no_shader.frag", ""},
                    {"grabber.vert", "grabber.frag", ""},
                    {"particle.vert", "particle.frag", "particle.geom"}
                };
                for (const auto& paths : shaderFamilies) {
                    Shader shader("src/Shaders/" + paths[0], "src/Shaders/" + paths[1],
                                  paths[2].empty() ? "" : "src/Shaders/" + paths[2]);
                    require(shader.use(), "An optimized shader family failed to link");
                }
                {
                    QFile file(directory.filePath("uniform-cache.comp"));
                    require(file.open(QIODevice::WriteOnly), "Could not create compute regression shader");
                    file.write("#version 430\nlayout(local_size_x=1) in; uniform float value; layout(std430,binding=0) buffer Result { float samples[]; }; void main(){samples[gl_GlobalInvocationID.x]=value;}");
                    file.close();
                    ComputeShader compute(file.fileName().toStdString());
                    compute.setFloat("value", .75f);
                    GlobalsGL::ContextScope scope;
                    GLuint buffer = 0;
                    GlobalsGL::f()->glGenBuffers(1, &buffer);
                    GlobalsGL::f()->glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer);
                    const float empty[2] = {};
                    GlobalsGL::f()->glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(empty), empty, GL_DYNAMIC_READ);
                    GlobalsGL::f()->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, buffer);
                    compute.apply(2, 1, 1);
                    GlobalsGL::f()->glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);
                    float samples[2] = {};
                    GlobalsGL::f()->glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(samples), samples);
                    GlobalsGL::f()->glDeleteBuffers(1, &buffer);
                    require(samples[0] == .75f && samples[1] == .75f, "Compute shader uniform cache failed");
                }
                GizmoState gizmoState;
                gizmoState.visual.body_color = Vector3::red;
                gizmoState.capabilities.allowTranslateX = true;
                gizmoState.capabilities.allowTranslateY = true;
                gizmoState.capabilities.allowTranslateZ = true;
                gizmoState.capabilities.allowRotateX = true;
                gizmoState.capabilities.allowRotateY = true;
                gizmoState.capabilities.allowRotateZ = true;
                // Reproduce the first-visible-frame path, without any worker threads.
                earlyGizmo.render(gizmoState);
                for (auto& gizmo : preparedGizmos) {
                    gizmoState.dirty = true;
                    gizmo->render(gizmoState);
                    gizmoState.settings.visible = false;
                    gizmo->render(gizmoState);
                    gizmoState.settings.visible = true;
                }
                // Exercise the main mouse gizmo on its first frame while pressed.
                const auto center = window->viewer->heightmap->getDimensions() * .5f;
                window->viewer->mainGrabber->move(center);
                window->viewer->mainGrabber->show();
                window->viewer->mouseDown = true;
                require(!window->viewer->grabFramebuffer().isNull(), "First mouse gizmo frame failed");
                window->viewer->mouseDown = false;
                const auto pixel = window->viewer->camera()->projectedCoordinatesOf(center);
                const QPointF local(pixel.x, pixel.y);
                window->viewer->mousePos = local.toPoint();
                QMouseEvent move(QEvent::MouseMove, local, local, Qt::NoButton,
                                 Qt::NoButton, Qt::NoModifier);
                QApplication::sendEvent(window->viewer, &move);
                QMouseEvent press(QEvent::MouseButtonPress, local, local, Qt::LeftButton,
                                  Qt::LeftButton, Qt::NoModifier);
                QApplication::sendEvent(window->viewer, &press);
                require(!window->viewer->grabFramebuffer().isNull(), "Pressed mouse frame failed");
                QMouseEvent release(QEvent::MouseButtonRelease, local, local, Qt::LeftButton,
                                    Qt::NoButton, Qt::NoModifier);
                QApplication::sendEvent(window->viewer, &release);
                auto shader = std::make_shared<Shader>(early);
                Mesh original(std::vector<Vector3>{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, shader);
                original.update();
                Mesh copied = original;
                copied.translate(Vector3(1, 0, 0));
                copied.update();
                require(copied.vao != original.vao, "Mesh copy reused writable GPU buffers");
                std::weak_ptr<GlobalsGL::Resources> discardedProgram;
                {
                    Mesh temporary(std::make_shared<Shader>(early));
                    temporary.update();
                    discardedProgram = temporary.shader->gpuResources;
                }
                Shader::applyToAllShaders([](const auto&) {});
                require(discardedProgram.expired(), "Shader registry retained a discarded mesh program");
                window->viewer->reloadAllShaders();
                auto terrain = std::dynamic_pointer_cast<TerrainGenerationInterface>(
                    window->viewer->interfaces.at("terraingeneration"));
                terrain->heightmap->heights = GridF(12, 10, 1);
                terrain->heightmap->heights.iterate([&](int x, int y, int) {
                    terrain->heightmap->heights(x, y) = 2.25f + .13f * x + .07f * y;
                });
                terrain->voxelGrid->setVoxelValues(GridF(12, 10, 8, -1.f));
                const auto expectedHeights = terrain->heightmap->heights;
                auto requireSurface = [&]() {
                    require(terrain->heightmap->heights.getDimensions() == expectedHeights.getDimensions(), "UI conversion changed dimensions");
                    for (size_t i = 0; i < expectedHeights.size(); ++i)
                        require(std::abs(terrain->heightmap->heights[i] - expectedHeights[i]) < 1e-4f, "UI conversion moved the surface");
                };
                terrain->heightmapToAll();
                terrain->voxelsToAll();
                requireSurface();
                terrain->layersToAll();
                requireSurface();
                terrain->implicitToAll();
                requireSurface();
                terrain->heightmapToAll();
                window->viewer->setSmoothingAlgorithm(MARCHING_CUBES);
                window->viewer->setSceneCenter(Vector3(6, 5, 2));
                terrain->setWaterLevel(0.25f);
                window->viewer->setSceneRadius(window->viewer->heightmap->getDimensions().maxComp());
                window->viewer->showEntireScene();
                QImage referenceVolume;
                for (auto mode : {GRID_MODE, VOXEL_MODE, LAYER_MODE, IMPLICIT_MODE}) {
                    window->viewer->setMapMode(mode);
                    const auto rendered = window->viewer->grabFramebuffer();
                    require(!rendered.isNull(), "Terrain mode failed to render");
                    if (mode != GRID_MODE) {
                        std::set<QRgb> volumeColors;
                        for (int y = 0; y < rendered.height(); y += 4)
                            for (int x = 0; x < rendered.width(); x += 4)
                                volumeColors.insert(rendered.pixel(x, y));
                        require(volumeColors.size() > 8, "Volumetric terrain frame is blank");
                    }
                    if (argc > 1)
                        require(rendered.save(QString::fromLocal8Bit(argv[1]) + QString(".mode-%1.png").arg(int(mode))), "Could not save conversion comparison");
                    if (mode == VOXEL_MODE) referenceVolume = rendered;
                    if (mode == LAYER_MODE || mode == IMPLICIT_MODE) {
                        int different = 0;
                        for (int y = 0; y < rendered.height(); ++y) {
                            for (int x = 0; x < rendered.width(); ++x) {
                                const auto a = rendered.pixelColor(x, y), b = referenceVolume.pixelColor(x, y);
                                if (std::abs(a.red() - b.red()) + std::abs(a.green() - b.green()) + std::abs(a.blue() - b.blue()) > 30) ++different;
                            }
                        }
                        qInfo() << "Conversion image comparison, mode" << int(mode) << ":" << different << "of" << rendered.width() * rendered.height() << "pixels differ significantly";
                        require(different < rendered.width() * rendered.height() / 50, "Converted volumetric render changed visibly");
                    }
                }
                window->viewer->setMapMode(GRID_MODE);
                Shader::applyToAllShaders([](const auto& shader) { shader->setBool("adaptiveTessellation", false); });
                const auto detailedFrame = window->viewer->grabFramebuffer();
                Shader::applyToAllShaders([](const auto& shader) { shader->setBool("adaptiveTessellation", true); });
                const auto adaptiveFrame = window->viewer->grabFramebuffer();
                double adaptiveDifference = 0.0;
                for (int y = 0; y < adaptiveFrame.height(); ++y) {
                    for (int x = 0; x < adaptiveFrame.width(); ++x) {
                        const auto a = adaptiveFrame.pixelColor(x, y), b = detailedFrame.pixelColor(x, y);
                        adaptiveDifference += std::abs(a.red() - b.red()) + std::abs(a.green() - b.green()) + std::abs(a.blue() - b.blue());
                    }
                }
                adaptiveDifference /= (3.0 * adaptiveFrame.width() * adaptiveFrame.height());
                qInfo() << "Adaptive heightmap mean channel difference:" << adaptiveDifference;
                require(adaptiveDifference < 2.0, "Adaptive tessellation changed the heightmap visibly");
                const auto frame = window->viewer->grabFramebuffer();
                require(!frame.isNull(), "Framebuffer capture failed");
                std::set<QRgb> colors;
                for (int y = 0; y < frame.height(); y += 4)
                    for (int x = 0; x < frame.width(); x += 4)
                        colors.insert(frame.pixel(x, y));
                require(colors.size() > 8, "Terrain lost its shader/texture bindings");
                const auto screenshot = argc > 1 ? QString::fromLocal8Bit(argv[1])
                                                : directory.filePath("frame.png");
                require(frame.save(screenshot), "Screenshot save failed");
                auto cpuTerrain = window->viewer->heightmap;
                auto previousResources = early.gpuResources;
                {
                    GlobalsGL::ContextScope scope;
                    // Simulate losing all context-owned objects, then rerun GL setup.
                    GlobalsGL::Resources::releaseContext(window->viewer->context());
                    window->viewer->init();
                }
                require(window->viewer->heightmap == cpuTerrain, "GPU rebuild replaced CPU terrain");
                require(previousResources->released, "GPU reset left old resources active");
                early.setVector("color", Vector4(0.f, 1.f, 0.f, 1.f));
                require(early.gpuResources != previousResources, "Shader reused released resources");
                require(!window->viewer->grabFramebuffer().isNull(), "GPU rebuild failed to render");
                auto lease = early.gpuResources;
                delete window;
                window = nullptr;
                require(lease->released, "GPU resources survived viewer teardown");
                require(!GlobalsGL::isReady(), "Destroyed viewer remained registered");
                qInfo() << "PASS: concurrent gizmo preparation, first render, mouse press/release, deferred shaders, GUI callbacks, copied meshes, shader reload,"
                           " water level, terrain conversion actions, equivalent volumetric renders, four terrain modes, screenshot, GPU rebuild and context cleanup";
                result = 0;
            } catch (const std::exception& e) {
                qCritical() << "FAIL:" << e.what();
            }
            ran = true;
            app.quit();
        });
    });
    QTimer::singleShot(45000, &app, [&]() {
        if (!ran) qCritical() << "FAIL: viewer initialization timed out";
        app.quit();
    });
    window->show();
    app.exec();
    delete window;
    return result;
}
