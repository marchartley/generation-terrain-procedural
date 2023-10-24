QT *= quick opengl xml widgets gui charts
CONFIG += qt opengl warn_on thread rtti console embed_manifest_exe no_keywords

INCLUDEPATH *= src/

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

#LIBS *= -L"C:\Program Files\ESI-OpenCFD\OpenFOAM\v2112\msys64\home\ofuser\OpenFOAM\OpenFOAM-v2112\platforms\win64MingwDPInt32Opt\bin"
#INCLUDEPATH *= src/third-party/openFoam/inInclude

#DEFINES += "WM_LINK_LANGUAGE=c++"
#DEFINES += "WM_ARCH=linux64"
#DEFINES += "WM_COMPILER_TYPE=system"
#DEFINES += "WM_OSTYPE=POSIX"
#DEFINES += "WM_THIRD_PARTY_DIR=/home/marzieh/ThirdParty-7"
#DEFINES += "WM_CXXFLAGS=-m64 -fPIC -std=c++0x"
#DEFINES += "WM_CFLAGS=-m64 -fPIC"
#DEFINES += "WM_PROJECT_VERSION=7"
#DEFINES += "WM_COMPILER_LIB_ARCH=64"
#DEFINES += "WM_PROJECT_INST_DIR=/home/marzieh"
#DEFINES += "WM_CXX=g++"
#DEFINES += "WM_PROJECT_DIR=/home/marzieh/OpenFOAM-7"
#DEFINES += "WM_LABEL_OPTION=Int32"
#DEFINES += "WM_PROJECT=OpenFOAM"
#DEFINES += "WM_LDFLAGS=-m64"
#DEFINES += "WM_COMPILER=Gcc"
#DEFINES += "WM_MPLIB=SYSTEMOPENMPI"
#DEFINES += "WM_CC=gcc"
#DEFINES += "WM_COMPILE_OPTION=Opt"
#DEFINES += "WM_DIR=/home/marzieh/OpenFOAM-7/wmake"
DEFINES += "WM_LABEL_SIZE=32"
#DEFINES += "WM_PROJECT_USER_DIR=/home/marzieh/OpenFOAM/marzieh-7"
#DEFINES += "WM_OPTIONS=linux64GccDPInt32Opt"
#DEFINES += "WM_PRECISION_OPTION=DP"
#DEFINES += "WM_ARCH_OPTION=64"

unix {
#    INCLUDEPATH *= /export/home/scharf/mhartley/codes/libQGLViewer-2.9.1 #/home/simulateurrsm/Documents/eigen #"/home/simulateurrsm/Documents/App downloads/tbb/include"
    #INCLUDEPATH *= src/third-party/boost_1_79_0/boost
#    LIBS *= -L/export/home/scharf/mhartley/codes/libQGLViewer-2.9.1/QGLViewer
    LIBS *= -lQGLViewer-qt5 #-ltbb -ltbbmalloc
    LIBS *= -lpng
#    INCLUDEPATH *= /home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2 /home/simulateurrsm/Documents/eigen #"/home/simulateurrsm/Documents/App downloads/tbb/include"
#    INCLUDEPATH *= src/third-party/boost_1_79_0/boost
#    LIBS *= -L/home/simulateurrsm/Documents/libqglviewer/libQGLViewer-2.7.2/QGLViewer -lQGLViewer-qt5 #-ltbb -ltbbmalloc
    INCLUDEPATH *= src/third-party/glad/include
    INCLUDEPATH *= src/third-party/glfw/include
    INCLUDEPATH *= src/third-party/glm
}
win32 {
    # I installed the sources of:
    # glew : https://github.com/nigels-com/glew
    # libqglviewer : http://www.libqglviewer.com/src/libQGLViewer-2.7.2.zip
    # eigen : https://gitlab.com/libeigen/eigen
    # OpenVDB : https://github.com/AcademySoftwareFoundation/openvdb
    # Boost and TBB are installed by VCPKG (https://github.com/microsoft/vcpkg) at the location C:\Programs_installations\vcpkg
    # For TBB on Windows, do everything on Release...
    INCLUDEPATH *= "C:\codes\CPP\boost_1_66_0"
#    INCLUDEPATH *= "C:\Program Files\Python39\include"
    INCLUDEPATH *= "C:\codes\CPP\glew-2.1.0\include"
    INCLUDEPATH *= "D:\code\Qt\libQGLViewer-2.7.2"
    INCLUDEPATH *= "src\third-party\glad\include"
    INCLUDEPATH *= "src\third-party\glfw\include"
    INCLUDEPATH *= "src\third-party\glm"
    INCLUDEPATH *= "src\third-party\boost_1_66_0"
    INCLUDEPATH *= "src\third-party"

    INCLUDEPATH *= "C:\Users\Marc\Downloads\libpng-1.2.37-lib\include"
#    INCLUDEPATH *= C:\codes\CPP\eigen
#    INCLUDEPATH *= "C:/Program Files/OpenVDB/include"
    # INCLUDEPATH *= C:\Programs_installations\vcpkg\installed\x64-windows\include
#    LIBS *= -L"C:\Program Files\Python39\libs\python3.lib"
    LIBS *= -L"C:\codes\CPP\glew-2.1.0\lib\Release\x64\glew32.lib"
    LIBS *= -L"C:\Qt\libQGLViewer-2.7.2\QGLViewer"

    LIBS *= -L"C:\Users\Marc\Downloads\libpng-1.2.37-lib\lib\libpng.lib"

    LIBS *= -lQGLViewer2 -lOpengl32
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
#    src/Biomes/InstancesTree.cpp \
#    src/Biomes/ModelsTree.cpp \
    src/DataStructure/BVH.cpp \
    src/DataStructure/Image.cpp \
    src/DataStructure/KDTree.cpp \
        src/DataStructure/Matrix.cpp \
        src/DataStructure/Matrix3.cpp \
    src/DataStructure/MemoryPool.cpp \
    src/DataStructure/Octree.cpp \
    src/DataStructure/Particle.cpp \
    src/DataStructure/Quaternion.cpp \
    src/DataStructure/SpacePartitioning.cpp \
    src/DataStructure/Tree.cpp \
    src/DataStructure/Triangle.cpp \
        src/DataStructure/Vector3.cpp \
    src/DataStructure/Vector4.cpp \
        src/DataStructure/Vertex.cpp \
        src/DataStructure/Voxel.cpp \
    src/EnvObject/EnvObject.cpp \
#        src/FastWFC/propagator.cpp \
#        src/FastWFC/wave.cpp \
#        src/FastWFC/wfc.cpp \
    src/EnvObject/ExpressionParser.cpp \
    src/FluidSimulation/FluidSimulation.cpp \
    src/FluidSimulation/FLIPSimulation.cpp \
    src/FluidSimulation/LBMFluidSimulation.cpp \
    src/FluidSimulation/OpenFoamParser.cpp \
    src/FluidSimulation/SPHSimulation.cpp \
    src/FluidSimulation/ShallowWaterSimulation.cpp \
    src/FluidSimulation/StableFluidsFluidSimulation.cpp \
    src/FluidSimulation/WarpedFluidSimulation.cpp \
        src/Graph/Graph.cpp \
        src/Graph/FastPoissonGraph.cpp \
        src/Graph/GraphNode.cpp \
        src/Graph/Matrix3Graph.cpp \
        src/Graph/Pathfinding.cpp \
    src/Graph/RegularSimplicialComplex.cpp \
    src/Graph/TopoMap.cpp \
        src/Graph/WaveFunctionCollapse.cpp \
    src/Graphics/ComputeShader.cpp \
        src/Graphics/CubeMesh.cpp \
        src/Graphics/DebugShader.cpp \
    src/Graphics/DisplayGraphics.cpp \
        src/Graphics/MarchingCubes.cpp \
    src/Graphics/RayMarching.cpp \
        src/Graphics/Shader.cpp \
        src/Graphics/Mesh.cpp \
        src/Graphics/ShaderElement.cpp \
        src/Graphics/Sphere.cpp \
    src/Graphics/miniz.c \
    src/Graphics/ofbx.cpp \
    src/Interface/AbstractFluidSimulationInterface.cpp \
        src/Interface/ActionInterface.cpp \
    src/Interface/BiomeInterface.cpp \
    src/Interface/CommonInterface.cpp \
        src/Interface/ControlPoint.cpp \
    src/Interface/CoralIslandGeneratorInterface.cpp \
        src/Interface/CustomInteractiveObject.cpp \
    src/Interface/EnvObjsInterface.cpp \
        src/Interface/ErosionInterface.cpp \
    src/Interface/FLIPSimulationInterface.cpp \
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
    src/Interface/LBMFluidSimulationInterface.cpp \
        src/Interface/ManualEditionInterface.cpp \
    src/Interface/MeshInstanceAmplificationInterface.cpp \
        src/Interface/PathCameraConstraint.cpp \
    src/Interface/PrimitivePatchesInterface.cpp \
        src/Interface/RangeSlider.cpp \
    src/Interface/SPHSimulationInterface.cpp \
        src/Interface/Slider3D.cpp \
    src/Interface/SmoothInterface.cpp \
        src/Interface/SpaceColonizationInterface.cpp \
    src/Interface/SpheroidalErosionInterface.cpp \
        src/Interface/Spoiler.cpp \
        src/Interface/StickyFrame.cpp \
        src/Interface/TerrainGenerationInterface.cpp \
        src/Interface/TerrainSavingInterface.cpp \
        src/Interface/TunnelInterface.cpp \
        src/Interface/UndoRedoInterface.cpp \
        src/Interface/Viewer.cpp \
        src/Interface/VisitingCamera.cpp \
    src/Interface/WarpFluidSimulationInterface.cpp \
        src/Karst/KarstHole.cpp \
        src/Karst/KarstHoleProfile.cpp \
        src/Karst/KarstPathsGeneration.cpp \
    src/TerrainGen/GlobalTerrainProperties.cpp \
    src/TerrainGen/Heightmap.cpp \
    src/TerrainGen/ImplicitPatch.cpp \
        src/TerrainGen/LayerBasedGrid.cpp \
    src/TerrainGen/TerrainModel.cpp \
        src/TerrainGen/VoxelChunk.cpp \
        src/TerrainGen/VoxelGrid.cpp \
    src/TerrainModification/CoralGrowth.cpp \
    src/TerrainModification/CoralIslandGenerator.cpp \
        src/TerrainModification/FaultSlip.cpp \
        src/TerrainModification/RockErosion.cpp \
    src/TerrainModification/SpheroidalWeathering.cpp \
        src/TerrainModification/TerrainAction.cpp \
        src/TerrainModification/UnderwaterErosion.cpp \
        src/TreeColonisation/TreeColonisation.cpp \
    src/Utils/AdjencySolver.cpp \
        src/Utils/BSpline.cpp \
        src/Utils/Collisions.cpp \
    src/Utils/ConstraintsSolver.cpp \
    src/Utils/Curve1D.cpp \
        src/Utils/Globals.cpp \
    src/Utils/PbmReader.cpp \
        src/Utils/RewritableFile.cpp \
    src/Utils/ShapeCurve.cpp \
    src/Utils/Skeletonize.cpp \
    src/Utils/Table.cpp \
        src/Utils/Utils.cpp \
    src/Utils/Voronoi.cpp \
#    src/libpng/example.c \
#    src/libpng/png.c \
#    src/libpng/pngerror.c \
#    src/libpng/pngget.c \
#    src/libpng/pngmem.c \
#    src/libpng/pngpread.c \
#    src/libpng/pngread.c \
#    src/libpng/pngrio.c \
#    src/libpng/pngrtran.c \
#    src/libpng/pngrutil.c \
#    src/libpng/pngset.c \
#    src/libpng/pngtest.c \
#    src/libpng/pngtrans.c \
#    src/libpng/pngwio.c \
#    src/libpng/pngwrite.c \
#    src/libpng/pngwtran.c \
#    src/libpng/pngwutil.c \
        src/main.cpp \
#    src/pngpp/pngpptest.cpp

RESOURCES += qml.qrc \
    icons.qrc \
#    models_3d.qrc \
#    terrain_textures.qrc \
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
    src/DataStructure/BVH.h \
    src/DataStructure/Image.h \
    src/DataStructure/KDTree.h \
    src/DataStructure/Matrix.h \
    src/DataStructure/Matrix3.h \
    src/DataStructure/MemoryPool.h \
    src/DataStructure/Octree.h \
    src/DataStructure/Particle.h \
    src/DataStructure/Quaternion.h \
    src/DataStructure/SpacePartitioning.h \
    src/DataStructure/Tree.h \
    src/DataStructure/Triangle.h \
    src/DataStructure/Vector3.h \
    src/DataStructure/Vector4.h \
    src/DataStructure/Vertex.h \
    src/DataStructure/Voxel.h \
    src/EnvObject/EnvObject.h \
#    src/FastWFC/color.hpp \
#    src/FastWFC/direction.hpp \
#    src/FastWFC/external/rapidxml.hpp \
#    src/FastWFC/image.hpp \
#    src/FastWFC/overlapping_wfc.hpp \
#    src/FastWFC/propagator.hpp \
#    src/FastWFC/rapidxml_utils.hpp \
#    src/FastWFC/tiling_wfc.hpp \
#    src/FastWFC/utils.hpp \
#    src/FastWFC/utils/array2D.hpp \
#    src/FastWFC/utils/array3D.hpp \
#    src/FastWFC/wave.hpp \
#    src/FastWFC/wfc.hpp \
    src/EnvObject/ExpressionParser.h \
    src/FluidSimulation/FluidSimulation.h \
    src/FluidSimulation/FLIPSimulation.h \
    src/FluidSimulation/LBMFluidSimulation.h \
    src/FluidSimulation/OpenFoamParser.h \
    src/FluidSimulation/SPHSimulation.h \
    src/FluidSimulation/ShallowWaterSimulation.h \
    src/FluidSimulation/StableFluidsFluidSimulation.h \
    src/FluidSimulation/WarpedFluidSimulation.h \
    src/Graph/Graph.h \
    src/Graph/FastPoissonGraph.h \
    src/Graph/GraphNode.h \
    src/Graph/Matrix3Graph.h \
    src/Graph/Pathfinding.h \
    src/Graph/RegularSimplicialComplex.h \
    src/Graph/TopoMap.h \
    src/Graph/WaveFunctionCollapse.h \
    src/Graphics/ComputeShader.h \
    src/Graphics/CubeMesh.h \
    src/Graphics/DebugShader.h \
    src/Graphics/DisplayGraphics.h \
    src/Graphics/MarchingCubes.h \
    src/Graphics/Mesh.h \
    src/Graphics/RayMarching.h \
    src/Graphics/Shader.h \
    src/Graphics/ShaderElement.h \
    src/Graphics/Sphere.h \
    src/Graphics/miniz.h \
    src/Graphics/ofbx.h \
    src/Interface/AbstractFluidSimulationInterface.h \
    src/Interface/ActionInterface.h \
    src/Interface/BiomeInterface.h \
    src/Interface/CommonInterface.h \
    src/Interface/ControlPoint.h \
    src/Interface/CoralIslandGeneratorInterface.h \
    src/Interface/CustomInteractiveObject.h \
    src/Interface/EnvObjsInterface.h \
    src/Interface/ErosionInterface.h \
    src/Interface/FLIPSimulationInterface.h \
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
    src/Interface/LBMFluidSimulationInterface.h \
    src/Interface/ManualEditionInterface.h \
    src/Interface/MeshInstanceAmplificationInterface.h \
    src/Interface/PathCameraConstraint.h \
    src/Interface/PrimitivePatchesInterface.h \
    src/Interface/RangeSlider.h \
    src/Interface/SPHSimulationInterface.h \
    src/Interface/Slider3D.h \
    src/Interface/SmoothInterface.h \
    src/Interface/SpaceColonizationInterface.h \
    src/Interface/SpheroidalErosionInterface.h \
    src/Interface/Spoiler.h \
    src/Interface/StickyFrame.h \
    src/Interface/TerrainGenerationInterface.h \
    src/Interface/TerrainSavingInterface.h \
    src/Interface/TunnelInterface.h \
    src/Interface/UndoRedoInterface.h \
    src/Interface/Viewer.h \
    src/Interface/VisitingCamera.h \
    src/Interface/WarpFluidSimulationInterface.h \
    src/Karst/KarstHole.h \
    src/Karst/KarstHoleProfile.h \
    src/TerrainGen/GlobalTerrainProperties.h \
    src/TerrainGen/Heightmap.h \
    src/TerrainGen/ImplicitPatch.h \
    src/TerrainGen/LayerBasedGrid.h \
    src/TerrainGen/TerrainModel.h \
    src/TerrainGen/VoxelChunk.h \
    src/TerrainGen/VoxelGrid.h \
    src/TerrainModification/CoralGrowth.h \
    src/TerrainModification/CoralIslandGenerator.h \
    src/TerrainModification/FaultSlip.h \
    src/TerrainModification/RockErosion.h \
    src/TerrainModification/SpheroidalWeathering.h \
    src/TerrainModification/TerrainAction.h \
    src/TerrainModification/UnderwaterErosion.h \
    src/TreeColonisation/TreeColonisation.h \
    src/Utils/AdjencySolver.h \
    src/Utils/BSpline.h \
    src/Utils/Collisions.h \
    src/Utils/ConstraintsSolver.h \
    src/Utils/Curve1D.h \
    src/Utils/FastNoiseLit.h \
    src/Utils/Globals.h \
    src/Karst/KarstPathsGeneration.h \
    src/Utils/PbmReader.h \
    src/Utils/RewritableFile.h \
    src/Utils/ShapeCurve.h \
    src/Utils/Skeletonize.h \
    src/Utils/Table.h \
    src/Utils/Utils.h \
    src/Utils/Voronoi.h \
    src/Utils/gnuplot-iostream.h \
    src/Utils/jc_voronoi.h \
    src/Utils/jc_voronoi_clip.h \
    src/Utils/json.h \
    src/Utils/matplotlibcpp.h \
    src/Utils/stb_image.h \
    src/Utils/stb_image_write.h \
    src/Utils/stb_image_resize.h \
    src/Utils/stl_reader.h \
#    src/libpng/png.h \
#    src/libpng/pngconf.h \
#    src/libpng/pngdebug.h \
#    src/libpng/pnginfo.h \
#    src/libpng/pngpriv.h \
#    src/libpng/pngstruct.h \
#    src/pngpp/convert_color_space.hpp \
#    src/pngpp/end_info.hpp \
#    src/pngpp/error.hpp \
#    src/pngpp/image.hpp \
#    src/pngpp/info.hpp \
#    src/pngpp/info_base.hpp \
#    src/pngpp/io_base.hpp \
#    src/pngpp/io_transform.hpp \
#    src/pngpp/pixel_buffer.hpp \
#    src/pngpp/pixel_traits.hpp \
#    src/pngpp/png.hpp \
#    src/pngpp/reader.hpp \
#    src/pngpp/require_color_space.hpp \
#    src/pngpp/rgb_pixel.hpp \
#    src/pngpp/rgba_pixel.hpp \
#    src/pngpp/types.hpp \
#    src/pngpp/writer.hpp \
    src/sim-fluid-loganzartman/Box.hpp \
    src/sim-fluid-loganzartman/DebugLine.hpp \
    src/sim-fluid-loganzartman/Fluid.hpp \
    src/sim-fluid-loganzartman/Game.hpp \
    src/sim-fluid-loganzartman/GridCell.hpp \
    src/sim-fluid-loganzartman/P2GTransfer.hpp \
    src/sim-fluid-loganzartman/Particle.hpp \
    src/sim-fluid-loganzartman/Quad.hpp \
    src/sim-fluid-loganzartman/SSFBufferElement.hpp \
    src/sim-fluid-loganzartman/SSFRenderTexture.hpp \
    src/sim-fluid-loganzartman/gfx/object.hpp \
    src/sim-fluid-loganzartman/gfx/program.hpp \
    src/sim-fluid-loganzartman/gfx/rendertexture.hpp \
    src/sim-fluid-loganzartman/util.hpp



DISTFILES += \
    src/Shaders/computeMC.comp \
    src/Shaders/gouraud_120.frag \
    src/Shaders/gouraud_120.vert \
    src/Shaders/layer_based.geom \
    src/Shaders/test_raymarching_voxels.frag \
    src/Shaders/test_raymarching_voxels.vert \
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
    src/Shaders/rockShader.frag \
    src/Shaders/rockShader.vert
