QT *= quick opengl xml widgets gui
CONFIG += qt opengl warn_on thread rtti console embed_manifest_exe no_keywords

unix {
    INCLUDEPATH *= /home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2 /home/simulateurrsm/Documents/eigen #"/home/simulateurrsm/Documents/App downloads/tbb/include"
    LIBS *= -L/home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2/QGLViewer -lQGLViewer-qt5 -ltbb -ltbbmalloc
}
win32 {
    # I installed the sources of:
    # glew : https://github.com/nigels-com/glew
    # libqglviewer : http://www.libqglviewer.com/src/libQGLViewer-2.7.2.zip
    # eigen : https://gitlab.com/libeigen/eigen
    # OpenVDB : https://github.com/AcademySoftwareFoundation/openvdb
    # Boost and TBB are installed by VCPKG (https://github.com/microsoft/vcpkg) at the location C:\Programs_installations\vcpkg
    # For TBB on Windows, do everything on Release...
    INCLUDEPATH *= "C:\codes\CPP\glew-2.1.0\include" "C:\Qt\libQGLViewer-2.7.2" C:\codes\CPP\eigen "C:/Program Files/OpenVDB/include" C:\Programs_installations\vcpkg\installed\x64-windows\include
    LIBS *= -L"C:\codes\CPP\glew-2.1.0\lib\Release\x64\glew32.lib" -L"C:\Qt\libQGLViewer-2.7.2\QGLViewer" -lQGLViewer2 -lOpengl32
    LIBS *= -LC:\Programs_installations\vcpkg\installed\x64-windows\lib -ltbb -ltbbmalloc -LC:\Programs_installations\vcpkg\installed\x64-windows\debug\bin -ltbb_debug
    LIBS *= -L"C:\Program Files\OpenVDB\bin" -lopenvdb
    DEFINES += "OPENVDB_DLL"
}
CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
TEMPLATE = app
TARGET = interface

SOURCES += \
        BSpline.cpp \
        FluidSimulation.cpp \
        Globals.cpp \
        Grid.cpp \
        Interface.cpp \
        LayerBasedGrid.cpp \
        MarchingCubes.cpp \
        Matrix.cpp \
        Matrix3.cpp \
        Mesh.cpp \
        RockErosion.cpp \
        Shader.cpp \
        ShaderElement.cpp \
        Sphere.cpp \
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
    BSpline.h \
    FastNoiseLit.h \
    FluidSimulation.h \
    Globals.h \
    Grid.h \
    Interface.h \
    LayerBasedGrid.h \
    MarchingCubes.h \
    Matrix.h \
    Matrix3.h \
    Mesh.h \
    RockErosion.h \
    Shader.h \
    ShaderElement.h \
    Sphere.h \
    UnderwaterErosion.h \
    Vector3.h \
    Vertex.h \
    Viewer.h \
    Voxel.h \
    VoxelChunk.h \
    VoxelGrid.h

RESOURCES += \
fragment_shader_gouraud.glsl \
    grid_fragment_shader_blinn_phong.glsl \
    grid_vertex_shader_blinn_phong.glsl \
vertex_shader_gouraud.glsl


DISTFILES += \
    grabber_fragment_shader.glsl \
    grabber_vertex_shader.glsl \
    layer_based_fragment_shader.glsl \
    layer_based_vertex_shader.glsl \
    no_fragment_shader.glsl \
    no_vertex_shader.glsl \
    noise.glsl \
    voxels_fragment_shader_blinn_phong.glsl \
    voxels_vertex_shader_blinn_phong.glsl

FORMS +=
