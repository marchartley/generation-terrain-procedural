QT += quick opengl xml

INCLUDEPATH *= /home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2 /home/simulateurrsm/Documents/eigen
LIBS *= -L/home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2/QGLViewer -lQGLViewer-qt5

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        Globals.cpp \
        Grid.cpp \
        MarchingCubes.cpp \
        RockErosion.cpp \
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
    RockErosion.h \
    UnderwaterErosion.h \
    Vector3.h \
    Vertex.h \
    Viewer.h \
    Voxel.h \
    VoxelChunk.h \
    VoxelGrid.h
