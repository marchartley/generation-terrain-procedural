QT *= quick opengl xml widgets gui charts
CONFIG += qt opengl warn_on thread rtti console embed_manifest_exe no_keywords

INCLUDEPATH *= src/

unix {
    INCLUDEPATH *= /home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2 /home/simulateurrsm/Documents/eigen #"/home/simulateurrsm/Documents/App downloads/tbb/include"
    INCLUDEPATH *= src/third-party/boost_1_79_0/boost
    LIBS *= -L/home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2/QGLViewer -lQGLViewer-qt5 #-ltbb -ltbbmalloc
}
win32 {
    # I installed the sources of:
    # glew : https://github.com/nigels-com/glew
    # libqglviewer : http://www.libqglviewer.com/src/libQGLViewer-2.7.2.zip
    # eigen : https://gitlab.com/libeigen/eigen
    # OpenVDB : https://github.com/AcademySoftwareFoundation/openvdb
    # Boost and TBB are installed by VCPKG (https://github.com/microsoft/vcpkg) at the location C:\Programs_installations\vcpkg
    # For TBB on Windows, do everything on Release...
    INCLUDEPATH *= "C:\codes\CPP\boost_1_66_0" "C:\Program Files\Python39\include" "C:\codes\CPP\glew-2.1.0\include" "C:\Qt\libQGLViewer-2.7.2" C:\codes\CPP\eigen "C:/Program Files/OpenVDB/include" # C:\Programs_installations\vcpkg\installed\x64-windows\include
    LIBS *= -L"C:\Program Files\Python39\libs\python3.lib" -L"C:\codes\CPP\glew-2.1.0\lib\Release\x64\glew32.lib" -L"C:\Qt\libQGLViewer-2.7.2\QGLViewer" -lQGLViewer2 -lOpengl32
#    LIBS *= -LC:\Programs_installations\vcpkg\installed\x64-windows\lib -ltbb -ltbbmalloc -LC:\Programs_installations\vcpkg\installed\x64-windows\debug\bin -ltbb_debug
#    LIBS *= -L"C:\Program Files\OpenVDB\bin" -lopenvdb
#    DEFINES += "OPENVDB_DLL"
}
CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
TEMPLATE = app
TARGET = interface

SOURCES += \
    src/Biomes/BiomeInstance.cpp \
    src/Biomes/BiomeModel.cpp \
    src/Biomes/BiomeUtils.cpp \
    src/Biomes/InstancesTree.cpp \
    src/Biomes/ModelsTree.cpp \
        src/DataStructure/Matrix.cpp \
        src/DataStructure/Matrix3.cpp \
    src/DataStructure/Tree.cpp \
        src/DataStructure/Vector3.cpp \
        src/DataStructure/Vertex.cpp \
        src/DataStructure/Voxel.cpp \
        src/FastWFC/propagator.cpp \
        src/FastWFC/wave.cpp \
        src/FastWFC/wfc.cpp \
        src/FluidSimulation/FluidSimulation.cpp \
        src/Graph/Graph.cpp \
        src/Graph/FastPoissonGraph.cpp \
        src/Graph/GraphNode.cpp \
        src/Graph/Matrix3Graph.cpp \
        src/Graph/Pathfinding.cpp \
    src/Graph/RegularSimplicialComplex.cpp \
    src/Graph/TopoMap.cpp \
        src/Graph/WaveFunctionCollapse.cpp \
        src/Graphics/CubeMesh.cpp \
        src/Graphics/DebugShader.cpp \
    src/Graphics/DisplayGraphics.cpp \
        src/Graphics/MarchingCubes.cpp \
        src/Graphics/Shader.cpp \
        src/Graphics/Mesh.cpp \
        src/Graphics/ShaderElement.cpp \
        src/Graphics/Sphere.cpp \
        src/Interface/ActionInterface.cpp \
    src/Interface/BiomeInterface.cpp \
        src/Interface/ControlPoint.cpp \
        src/Interface/CustomInteractiveObject.cpp \
        src/Interface/ErosionInterface.cpp \
        src/Interface/FancySlider.cpp \
        src/Interface/FaultSlipInterface.cpp \
        src/Interface/FlowFieldInterface.cpp \
        src/Interface/GravityInterface.cpp \
        src/Interface/HeightmapErosionInterface.cpp \
    src/Interface/HierarchicalListWidget.cpp \
        src/Interface/InteractiveVector.cpp \
        src/Interface/Interface.cpp \
        src/Interface/InterfaceUtils.cpp \
        src/Interface/KarstPathGenerationInterface.cpp \
        src/Interface/ManualEditionInterface.cpp \
        src/Interface/PathCameraConstraint.cpp \
        src/Interface/RangeSlider.cpp \
        src/Interface/Slider3D.cpp \
        src/Interface/SpaceColonizationInterface.cpp \
        src/Interface/Spoiler.cpp \
        src/Interface/StickyFrame.cpp \
        src/Interface/TerrainGenerationInterface.cpp \
        src/Interface/TunnelInterface.cpp \
        src/Interface/UndoRedoInterface.cpp \
        src/Interface/Viewer.cpp \
        src/Interface/VisitingCamera.cpp \
        src/Karst/KarstHole.cpp \
        src/Karst/KarstHoleProfile.cpp \
        src/Karst/KarstPathsGeneration.cpp \
        src/TerrainGen/Grid.cpp \
        src/TerrainGen/LayerBasedGrid.cpp \
        src/TerrainGen/VoxelChunk.cpp \
        src/TerrainGen/VoxelGrid.cpp \
        src/TerrainModification/FaultSlip.cpp \
        src/TerrainModification/RockErosion.cpp \
        src/TerrainModification/TerrainAction.cpp \
        src/TerrainModification/UnderwaterErosion.cpp \
        src/TreeColonisation/TreeColonisation.cpp \
    src/Utils/AdjencySolver.cpp \
        src/Utils/BSpline.cpp \
        src/Utils/Collisions.cpp \
    src/Utils/ConstraintsSolver.cpp \
        src/Utils/Globals.cpp \
        src/Utils/RewritableFile.cpp \
    src/Utils/ShapeCurve.cpp \
        src/Utils/Utils.cpp \
    src/Utils/Voronoi.cpp \
        src/main.cpp
        #src/sim-fluid-ethanjli/fluidsystem.cpp \
        #src/sim-fluid-ethanjli/math.cpp

RESOURCES += qml.qrc \
    icons.qrc \
    models_3d.qrc \
    terrain_textures.qrc \
    tunnels_icons.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Additional import path used to resolve QML modules just for Qt Quick Designer
QML_DESIGNER_IMPORT_PATH =

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    src/Biomes/BiomeInstance.h \
    src/Biomes/BiomeModel.h \
    src/Biomes/BiomeUtils.h \
    src/Biomes/InstancesTree.h \
    src/Biomes/ModelsTree.h \
    src/DataStructure/Matrix.h \
    src/DataStructure/Matrix3.h \
    src/DataStructure/Tree.h \
    src/DataStructure/Vector3.h \
    src/DataStructure/Vertex.h \
    src/DataStructure/Voxel.h \
    src/FastWFC/color.hpp \
    src/FastWFC/direction.hpp \
    src/FastWFC/external/rapidxml.hpp \
    src/FastWFC/image.hpp \
    src/FastWFC/overlapping_wfc.hpp \
    src/FastWFC/propagator.hpp \
    src/FastWFC/rapidxml_utils.hpp \
    src/FastWFC/tiling_wfc.hpp \
    src/FastWFC/utils.hpp \
    src/FastWFC/utils/array2D.hpp \
    src/FastWFC/utils/array3D.hpp \
    src/FastWFC/wave.hpp \
    src/FastWFC/wfc.hpp \
    src/FluidSimulation/FluidSimulation.h \
    src/Graph/Graph.h \
    src/Graph/FastPoissonGraph.h \
    src/Graph/GraphNode.h \
    src/Graph/Matrix3Graph.h \
    src/Graph/Pathfinding.h \
    src/Graph/RegularSimplicialComplex.h \
    src/Graph/TopoMap.h \
    src/Graph/WaveFunctionCollapse.h \
    src/Graphics/CubeMesh.h \
    src/Graphics/DebugShader.h \
    src/Graphics/DisplayGraphics.h \
    src/Graphics/MarchingCubes.h \
    src/Graphics/Mesh.h \
    src/Graphics/Shader.h \
    src/Graphics/ShaderElement.h \
    src/Graphics/Sphere.h \
    src/Interface/ActionInterface.h \
    src/Interface/BiomeInterface.h \
    src/Interface/ControlPoint.h \
    src/Interface/CustomInteractiveObject.h \
    src/Interface/ErosionInterface.h \
    src/Interface/FancySlider.h \
    src/Interface/FaultSlipInterface.h \
    src/Interface/FlowFieldInterface.h \
    src/Interface/GravityInterface.h \
    src/Interface/HeightmapErosionInterface.h \
    src/Interface/HierarchicalListWidget.h \
    src/Interface/InteractiveVector.h \
    src/Interface/Interface.h \
    src/Interface/InterfaceUtils.h \
    src/Interface/KarstPathGenerationInterface.h \
    src/Interface/ManualEditionInterface.h \
    src/Interface/PathCameraConstraint.h \
    src/Interface/RangeSlider.h \
    src/Interface/Slider3D.h \
    src/Interface/SpaceColonizationInterface.h \
    src/Interface/Spoiler.h \
    src/Interface/StickyFrame.h \
    src/Interface/TerrainGenerationInterface.h \
    src/Interface/TunnelInterface.h \
    src/Interface/UndoRedoInterface.h \
    src/Interface/Viewer.h \
    src/Interface/VisitingCamera.h \
    src/Karst/KarstHole.h \
    src/Karst/KarstHoleProfile.h \
    src/TerrainGen/Grid.h \
    src/TerrainGen/LayerBasedGrid.h \
    src/TerrainGen/VoxelChunk.h \
    src/TerrainGen/VoxelGrid.h \
    src/TerrainModification/FaultSlip.h \
    src/TerrainModification/RockErosion.h \
    src/TerrainModification/TerrainAction.h \
    src/TerrainModification/UnderwaterErosion.h \
    src/TreeColonisation/TreeColonisation.h \
    src/Utils/AdjencySolver.h \
    src/Utils/BSpline.h \
    src/Utils/Collisions.h \
    src/Utils/ConstraintsSolver.h \
    src/Utils/FastNoiseLit.h \
    src/Utils/Globals.h \
    src/Karst/KarstPathsGeneration.h \
    src/Utils/RewritableFile.h \
    src/Utils/ShapeCurve.h \
    src/Utils/Utils.h \
    src/Utils/Voronoi.h \
    src/Utils/gnuplot-iostream.h \
    src/Utils/jc_voronoi.h \
    src/Utils/jc_voronoi_clip.h \
    src/Utils/json.h \
    src/Utils/matplotlibcpp.h \
    src/Utils/stb_image.h \
    src/Utils/stb_image_write.h \
    src/Utils/stl_reader.h
#    src/sim-fluid-ethanjli/fluidsystem.h \
#    src/sim-fluid-ethanjli/fluidsystem.tpp \
#    src/sim-fluid-ethanjli/math.h \
#    src/sim-fluid-ethanjli/math.tpp \
#    src/sim-fluid-ethanjli/vectorfield.h \
#    src/sim-fluid-ethanjli/vectorfield.tpp

RESOURCES +=\
    src/Shaders/MarchingCubes.frag \
    src/Shaders/MarchingCubes.geom \
    src/Shaders/MarchingCubes.vert \
    src/Shaders/gouraud.frag \
    src/Shaders/gouraud.vert \
    src/Shaders/grabber.frag \
    src/Shaders/grabber.vert \
    src/Shaders/grid.frag \
    src/Shaders/grid.geom \
    src/Shaders/grid.vert \
    src/Shaders/layer_based.frag \
    src/Shaders/layer_based.vert \
    src/Shaders/no_shader.frag \
    src/Shaders/no_shader.vert \
    src/Shaders/particle.frag \
    src/Shaders/noise.glsl \
    src/Shaders/particle.vert \
    src/Shaders/particle.geom \
    src/Shaders/voxels.frag \
    src/Shaders/voxels.vert \
    src/assets/manual-edit_button.png\
    src/Shaders/rockShader.frag \
    src/Shaders/rockShader.vert \
    src/assets/handle.png \
    src/assets/fault-slip_button.png \
    src/assets/flowfield_button.png \
    src/assets/gravity_button.png \
    src/assets/karst_button.png \
    src/assets/open_button.png \
    src/assets/recording_button.png \
    src/assets/save_button.png \
    src/assets/tunnel_button.png


#DISTFILES +=

FORMS +=

DISTFILES +=
