#ifndef GLOBALS_H
#define GLOBALS_H

#define useModernOpenGL 1

#include <QOpenGLFunctions_4_5_Compatibility>
#include <QOpenGLExtraFunctions>
#include <QOpenGLWidget>
#include <QPointer>
#include <memory>
#include <vector>
#include <string>

#define numVAOs 1000
#define numVBOs 4000

class GlobalsGL {
public:
    // GPU entry points use this scope; construction and CPU work need no context.
    class ContextScope {
    public:
        ContextScope();
        ~ContextScope();
        ContextScope(const ContextScope&) = delete;
        ContextScope& operator=(const ContextScope&) = delete;
    private:
        QPointer<QOpenGLContext> previous;
        QSurface* previousSurface = nullptr;
        bool switched = false;
    };
    // Shared leases prevent copied meshes from double-deleting GPU objects.
    // Each lease is also released before its context is destroyed.
    struct Resources {
        QPointer<QOpenGLContext> owner;
        bool released = false;
        QMetaObject::Connection destructionConnection;
        std::vector<GLuint> programs, textures, buffers, arrays, framebuffers;
        ~Resources();
        void release();
        static std::shared_ptr<Resources> create();
        static void releaseContext(QOpenGLContext* context);
    };
    static void setViewer(QOpenGLWidget* viewer);
    static bool isReady();
    static QPointer<QOpenGLWidget> viewer;
    static QOpenGLContext *_context;
    static QOpenGLContext *context();
    static QOpenGLFunctions* _f;
//    static QOpenGLFunctions_4_6_Compatibility *f(); // Alias for ef()
    // static QOpenGLFunctions_4_5_Compatibility *f(); // Alias for ef()
    static QOpenGLFunctions_4_5_Compatibility* f();
    static QOpenGLExtraFunctions* _ef;
    static QOpenGLFunctions_4_5_Compatibility* f45();
    static QOpenGLFunctions_4_5_Compatibility* _ef45;
//    static QOpenGLFunctions_4_6_Compatibility* f46();
//    static QOpenGLFunctions_4_6_Compatibility* _ef46;
    static QOpenGLExtraFunctions *ef();
    static GLuint _renderingProgram;
    static GLuint renderingProgram();

    static GLuint vao[numVAOs];
    static GLuint vbo[numVBOs];
    static bool buffersGenerated;

    static GLuint currentBufferId;

    static std::string readShaderSource(const std::string& filename);


    static GLuint createShaderProgram(const std::string& vertexShaderFile = "", std::string fragmentShaderFile = "");

    static GLuint newBufferId();

    static void generateBuffers();

    static bool checkOpenGLError();
    static bool printShaderErrors(GLuint shader);
    static bool printProgramErrors(int program);

//    static void GLAPIENTRY MessageCallback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam );
};
#endif // GLOBALS_H
