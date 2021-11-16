QT += quick opengl xml

unix {
    INCLUDEPATH *= /home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2 /home/simulateurrsm/Documents/eigen
    LIBS *= -L/home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2/QGLViewer -lQGLViewer-qt5
}
win32 {
    INCLUDEPATH *= "C:\codes\CPP\glew-2.1.0\include" "C:\Qt\libQGLViewer-2.7.2" C:\codes\CPP\eigen
    LIBS *= -L"C:\codes\CPP\glew-2.1.0\lib\Release\x64\glew32.lib" -L"C:\Qt\libQGLViewer-2.7.2\QGLViewer" -lQGLViewer2 -lOpengl32
}
CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        Globals.cpp \
        Grid.cpp \
        MarchingCubes.cpp \
        Matrix.cpp \
        Mesh.cpp \
        RockErosion.cpp \
        Shader.cpp \
        ShaderElement.cpp \
        UnderwaterErosion.cpp \
        Vector3.cpp \
        Vertex.cpp \
        Viewer.cpp \
        Voxel.cpp \
        VoxelChunk.cpp \
        VoxelGrid.cpp \
        main.cpp

RESOURCES += qml.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Additional import path used to resolve QML modules just for Qt Quick Designer
QML_DESIGNER_IMPORT_PATH =

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    FastNoiseLit.h \
    Globals.h \
    Grid.h \
    MarchingCubes.h \
    Matrix.h \
    Mesh.h \
    RockErosion.h \
    Shader.h \
    ShaderElement.h \
    Shader_templates.h \
    UnderwaterErosion.h \
    Vector3.h \
    Vertex.h \
    Viewer.h \
    Voxel.h \
    VoxelChunk.h \
    VoxelGrid.h

RESOURCES += \
fragment_shader_gouraud.glsl \
vertex_shader_gouraud.glsl \
fragment_shader_blinn_phong.glsl \
vertex_shader_blinn_phong.glsl

DISTFILES += \
    noise.glsl
