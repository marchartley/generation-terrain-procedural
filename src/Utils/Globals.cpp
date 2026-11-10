#include "Utils/Globals.h"
#include "Utils/Random.h"

#include <string>
#include <fstream>
#include <algorithm>

#define UNUSED(expr) (void)(expr)

std::default_random_engine random_gen::random_generator;
FastNoiseLite random_gen::perlinNoise;


QOpenGLContext* GlobalsGL::_context;
QOpenGLFunctions* GlobalsGL::_f;
QOpenGLExtraFunctions* GlobalsGL::_ef;
QOpenGLFunctions_4_5_Compatibility* GlobalsGL::_ef45;
//QOpenGLFunctions_4_6_Compatibility* GlobalsGL::_ef46;
GLuint GlobalsGL::_renderingProgram;
GLuint GlobalsGL::vao[numVAOs];
GLuint GlobalsGL::vbo[numVBOs];
GLuint GlobalsGL::currentBufferId = 0;
bool GlobalsGL::buffersGenerated = false;

QPointer<QOpenGLWidget> GlobalsGL::viewer;

namespace { std::vector<std::weak_ptr<GlobalsGL::Resources>> resourceLeases; }

void GlobalsGL::Resources::releaseContext(QOpenGLContext* context) {
    for (auto it = resourceLeases.begin(); it != resourceLeases.end();) {
        if (auto lease = it->lock()) {
            if (lease->owner == context) lease->release();
            ++it;
        } else it = resourceLeases.erase(it);
    }
}

std::shared_ptr<GlobalsGL::Resources> GlobalsGL::Resources::create() {
    auto resources = std::make_shared<Resources>();
    resources->owner = GlobalsGL::context();
    resourceLeases.erase(std::remove_if(resourceLeases.begin(), resourceLeases.end(),
        [](const auto& weak) { return weak.expired(); }), resourceLeases.end());
    resourceLeases.push_back(resources);
    std::weak_ptr<Resources> weak = resources;
    resources->destructionConnection = QObject::connect(resources->owner, &QOpenGLContext::aboutToBeDestroyed,
                     resources->owner, [weak]() {
        if (auto live = weak.lock()) live->release();
    }, Qt::DirectConnection);
    return resources;
}

GlobalsGL::Resources::~Resources() {
    QObject::disconnect(destructionConnection);
    release();
}

void GlobalsGL::Resources::release() {
    if (released) return;
    // A destroyed context has already reclaimed its objects. Never delete its IDs
    // in a replacement context, where the same numbers may identify other objects.
    if (owner && (QOpenGLContext::currentContext() == owner ||
                  (GlobalsGL::isReady() && viewer->context() == owner &&
                   QThread::currentThread() == viewer->thread()))) {
        ContextScope scope;
        if (QOpenGLContext::currentContext() == owner) {
            auto* gl = GlobalsGL::f();
            for (auto program : programs) gl->glDeleteProgram(program);
            gl->glDeleteFramebuffers(GLsizei(framebuffers.size()), framebuffers.data());
            gl->glDeleteTextures(GLsizei(textures.size()), textures.data());
            gl->glDeleteBuffers(GLsizei(buffers.size()), buffers.data());
            gl->glDeleteVertexArrays(GLsizei(arrays.size()), arrays.data());
        }
    }
    programs.clear(); textures.clear(); buffers.clear(); arrays.clear(); framebuffers.clear();
    released = true;
}


void GlobalsGL::setViewer(QOpenGLWidget* widget) { viewer = widget; }
bool GlobalsGL::isReady() { return viewer && viewer->isValid(); }

GlobalsGL::ContextScope::ContextScope()
    : previous(QOpenGLContext::currentContext())
{
    if (!GlobalsGL::isReady()) {
        if (!previous) qFatal("OpenGL resources requested before viewer initialization");
        return;
    }
    if (QThread::currentThread() != viewer->thread())
        qFatal("Viewer OpenGL operations must run on the GUI thread");
    if (previous == viewer->context()) return;
    if (previous) previousSurface = previous->surface();
    viewer->makeCurrent();
    if (QOpenGLContext::currentContext() != viewer->context())
        qFatal("Could not make the viewer OpenGL context current");
    switched = true;
}

GlobalsGL::ContextScope::~ContextScope()
{
    if (!switched) return;
    if (previous && previousSurface) previous->makeCurrent(previousSurface);
    else if (auto* current = QOpenGLContext::currentContext()) current->doneCurrent();
}

QOpenGLContext* GlobalsGL::context() {
    auto* ctx = QOpenGLContext::currentContext();
    if (!ctx) qFatal("OpenGL operation requires a current context");
    return ctx;
}

QOpenGLFunctions_4_5_Compatibility* GlobalsGL::f() { return f45(); }

QOpenGLFunctions_4_5_Compatibility* GlobalsGL::f45() {
    // Qt owns the wrapper and caches it per context. Never cache a raw pointer here.
    auto* functions = QOpenGLVersionFunctionsFactory::
        get<QOpenGLFunctions_4_5_Compatibility>(context());
    if (!functions) qFatal("OpenGL 4.5 compatibility profile unavailable");
    return functions;
}

QOpenGLExtraFunctions* GlobalsGL::ef() { return context()->extraFunctions(); }

void GlobalsGL::generateBuffers()
{
#if useModernOpenGL
    if(GlobalsGL::buffersGenerated)
        return;
    GlobalsGL::f()->glGenVertexArrays(numVAOs, GlobalsGL::vao);
    GlobalsGL::f()->glBindVertexArray(GlobalsGL::vao[0]);
    GlobalsGL::f()->glGenBuffers(numVBOs, GlobalsGL::vbo);
#endif
    GlobalsGL::buffersGenerated = true;
}
GLuint GlobalsGL::newBufferId()
{
    return GlobalsGL::currentBufferId += 1; // Gives space for vertex, texture and normals
}
bool GlobalsGL::checkOpenGLError()
{
    bool error = false;
    int glErr = glGetError();
    while(glErr != GL_NO_ERROR)
    {
        std::cout << "[OpenGL] Error: " << glErr << std::endl;
        error = true;
        glErr = glGetError();
    }
    return !error;
}

bool GlobalsGL::printShaderErrors(GLuint shader)
{
#if useModernOpenGL
    int state = 0;
    GlobalsGL::f()->glGetShaderiv(shader, GL_COMPILE_STATUS, &state);
    if (state == 1)
        return true;
    int len = 0;
    int chWritten = 0;
    char* log;
    GlobalsGL::f()->glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
    if (len > 0)
    {
        log = (char*)malloc(len);
        GlobalsGL::f()->glGetShaderInfoLog(shader, len, &chWritten, log);
        std::cout << "[OpenGL] Shader error: " << log << std::endl;
        free(log);
    }
#endif
    return false;
}
bool GlobalsGL::printProgramErrors(int program)
{

#if useModernOpenGL
    int state = 0;
    GlobalsGL::f()->glGetProgramiv(program, GL_LINK_STATUS, &state);
    if (state == 1)
        return true;
    int len = 0;
    int chWritten = 0;
    char* log;
    GlobalsGL::f()->glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
    if (len > 0)
    {
        log = (char*)malloc(len);
        GlobalsGL::f()->glGetProgramInfoLog(program, len, &chWritten, log);
        std::cout << "[OpenGL] Program error: " << log << std::endl;
        free(log);
    }
#endif
    return false;
}
/*
void GLAPIENTRY GlobalsGL::MessageCallback( GLenum source, GLenum type,
                                            GLuint id, GLenum severity,
                                            GLsizei length, const GLchar* message,
                                            const void* userParam )
{
    UNUSED(source);
    UNUSED(type);
    UNUSED(length);
    UNUSED(userParam);
    return;
    if (id == 131154) return; // Ignore "Pixel-path performance warning: Pixel transfer is synchronized with 3D rendering." due to screenshots
    if (severity == GL_DEBUG_SEVERITY_HIGH || severity == GL_DEBUG_SEVERITY_MEDIUM || severity == GL_DEBUG_SEVERITY_LOW) {
        std::string s_severity = (severity == GL_DEBUG_SEVERITY_HIGH ? "High" : severity == GL_DEBUG_SEVERITY_MEDIUM ? "Medium" : "Low");
        std::cout << "Error " << id << " [severity=" << s_severity << "]: " << message << std::endl;
    }
}*/
