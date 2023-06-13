QT *= quick opengl xml widgets gui charts
CONFIG += qt opengl warn_on thread rtti console embed_manifest_exe no_keywords

INCLUDEPATH *= src/

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

unix {
#    INCLUDEPATH *= /export/home/scharf/mhartley/codes/libQGLViewer-2.9.1 #/home/simulateurrsm/Documents/eigen #"/home/simulateurrsm/Documents/App downloads/tbb/include"
    #INCLUDEPATH *= src/third-party/boost_1_79_0/boost
#    LIBS *= -L/export/home/scharf/mhartley/codes/libQGLViewer-2.9.1/QGLViewer
    LIBS *= -lQGLViewer-qt5 #-ltbb -ltbbmalloc
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
    INCLUDEPATH *= "C:\Qt\libQGLViewer-2.7.2"
    INCLUDEPATH *= "src\third-party\glad\include"
    INCLUDEPATH *= "src\third-party\glfw\include"
    INCLUDEPATH *= "src\third-party\glm"
    INCLUDEPATH *= "src\third-party"
#    INCLUDEPATH *= C:\codes\CPP\eigen
#    INCLUDEPATH *= "C:/Program Files/OpenVDB/include"
    # INCLUDEPATH *= C:\Programs_installations\vcpkg\installed\x64-windows\include
#    LIBS *= -L"C:\Program Files\Python39\libs\python3.lib"
    LIBS *= -L"C:\codes\CPP\glew-2.1.0\lib\Release\x64\glew32.lib"
    LIBS *= -L"C:\Qt\libQGLViewer-2.7.2\QGLViewer"

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
    src/FluidSimulation_FLIP/FLIPSimulation.cpp \
    src/FluidSimulation_SPH/SPHSimulation.cpp \
    src/FluidSimulation_ShallowWater/ShallowWaterSimulation.cpp \
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
    src/Interface/MeshInstanceAmplificationInterface.cpp \
        src/Interface/PathCameraConstraint.cpp \
    src/Interface/PrimitivePatchesInterface.cpp \
        src/Interface/RangeSlider.cpp \
        src/Interface/Slider3D.cpp \
    src/Interface/SmoothInterface.cpp \
        src/Interface/SpaceColonizationInterface.cpp \
        src/Interface/Spoiler.cpp \
        src/Interface/StickyFrame.cpp \
        src/Interface/TerrainGenerationInterface.cpp \
        src/Interface/TerrainSavingInterface.cpp \
        src/Interface/TunnelInterface.cpp \
        src/Interface/UndoRedoInterface.cpp \
        src/Interface/Viewer.cpp \
        src/Interface/VisitingCamera.cpp \
        src/Karst/KarstHole.cpp \
        src/Karst/KarstHoleProfile.cpp \
        src/Karst/KarstPathsGeneration.cpp \
    src/TerrainGen/Heightmap.cpp \
    src/TerrainGen/ImplicitPatch.cpp \
        src/TerrainGen/LayerBasedGrid.cpp \
    src/TerrainGen/TerrainModel.cpp \
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
    src/Utils/PbmReader.cpp \
        src/Utils/RewritableFile.cpp \
    src/Utils/ShapeCurve.cpp \
        src/Utils/Utils.cpp \
    src/Utils/Voronoi.cpp \
        src/main.cpp \
#    src/third-party/GridFluidSim3D/src/aabb.cpp \
#    src/third-party/GridFluidSim3D/src/anisotropicparticlemesher.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/cbindings.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/config_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/cuboidfluidsource_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/fluidsimulation_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/fluidsimulationsavestate_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/fluidsource_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/sphericalfluidsource_c.cpp \
#    src/third-party/GridFluidSim3D/src/c_bindings/utils_c.cpp \
#    src/third-party/GridFluidSim3D/src/clscalarfield.cpp \
#    src/third-party/GridFluidSim3D/src/collision.cpp \
#    src/third-party/GridFluidSim3D/src/config.cpp \
#    src/third-party/GridFluidSim3D/src/cuboidfluidsource.cpp \
#    src/third-party/GridFluidSim3D/src/diffuseparticlesimulation.cpp \
#    src/third-party/GridFluidSim3D/src/fluidbrickgrid.cpp \
#    src/third-party/GridFluidSim3D/src/fluidbrickgridsavestate.cpp \
#    src/third-party/GridFluidSim3D/src/fluidmaterialgrid.cpp \
#    src/third-party/GridFluidSim3D/src/fluidsimulation.cpp \
#    src/third-party/GridFluidSim3D/src/fluidsimulationsavestate.cpp \
#    src/third-party/GridFluidSim3D/src/fluidsource.cpp \
#    src/third-party/GridFluidSim3D/src/gridindexkeymap.cpp \
#    src/third-party/GridFluidSim3D/src/gridindexvector.cpp \
#    src/third-party/GridFluidSim3D/src/implicitpointprimitive.cpp \
#    src/third-party/GridFluidSim3D/src/interpolation.cpp \
#    src/third-party/GridFluidSim3D/src/isotropicparticlemesher.cpp \
#    src/third-party/GridFluidSim3D/src/kernels/kernels.cpp \
#    src/third-party/GridFluidSim3D/src/levelset.cpp \
#    src/third-party/GridFluidSim3D/src/logfile.cpp \
#    src/third-party/GridFluidSim3D/src/macvelocityfield.cpp \
#    src/third-party/GridFluidSim3D/src/main.cpp \
#    src/third-party/GridFluidSim3D/src/particleadvector.cpp \
#    src/third-party/GridFluidSim3D/src/polygonizer3d.cpp \
#    src/third-party/GridFluidSim3D/src/pressuresolver.cpp \
#    src/third-party/GridFluidSim3D/src/render/brick_texture_packer/brick_texture_packer.cpp \
#    src/third-party/GridFluidSim3D/src/render/brick_texture_packer/lodepng/lodepng.cpp \
#    src/third-party/GridFluidSim3D/src/scalarfield.cpp \
#    src/third-party/GridFluidSim3D/src/spatialpointgrid.cpp \
#    src/third-party/GridFluidSim3D/src/sphericalfluidsource.cpp \
#    src/third-party/GridFluidSim3D/src/stopwatch.cpp \
#    src/third-party/GridFluidSim3D/src/trianglemesh.cpp \
#    src/third-party/GridFluidSim3D/src/turbulencefield.cpp \
#    src/third-party/GridFluidSim3D/src/utils.cpp \
#    src/third-party/GridFluidSim3D/src/vmath.cpp

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
    src/FluidSimulation_FLIP/FLIPSimulation.h \
    src/FluidSimulation_SPH/SPHSimulation.h \
    src/FluidSimulation_ShallowWater/ShallowWaterSimulation.h \
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
    src/Interface/MeshInstanceAmplificationInterface.h \
    src/Interface/PathCameraConstraint.h \
    src/Interface/PrimitivePatchesInterface.h \
    src/Interface/RangeSlider.h \
    src/Interface/Slider3D.h \
    src/Interface/SmoothInterface.h \
    src/Interface/SpaceColonizationInterface.h \
    src/Interface/Spoiler.h \
    src/Interface/StickyFrame.h \
    src/Interface/TerrainGenerationInterface.h \
    src/Interface/TerrainSavingInterface.h \
    src/Interface/TunnelInterface.h \
    src/Interface/UndoRedoInterface.h \
    src/Interface/Viewer.h \
    src/Interface/VisitingCamera.h \
    src/Karst/KarstHole.h \
    src/Karst/KarstHoleProfile.h \
    src/TerrainGen/Heightmap.h \
    src/TerrainGen/ImplicitPatch.h \
    src/TerrainGen/LayerBasedGrid.h \
    src/TerrainGen/TerrainModel.h \
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
    src/Utils/PbmReader.h \
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
    src/Utils/stl_reader.h \
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
    src/sim-fluid-loganzartman/util.hpp \
#    src/third-party/GridFluidSim3D/src/aabb.h \
#    src/third-party/GridFluidSim3D/src/anisotropicparticlemesher.h \
#    src/third-party/GridFluidSim3D/src/array3d.h \
#    src/third-party/GridFluidSim3D/src/arrayview3d.h \
#    src/third-party/GridFluidSim3D/src/brick.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/aabb_c.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/cbindings.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/diffuseparticle_c.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/gridindex_c.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/markerparticle_c.h \
#    src/third-party/GridFluidSim3D/src/c_bindings/vector3_c.h \
#    src/third-party/GridFluidSim3D/src/clscalarfield.h \
#    src/third-party/GridFluidSim3D/src/collision.h \
#    src/third-party/GridFluidSim3D/src/config.h.in \
#    src/third-party/GridFluidSim3D/src/cuboidfluidsource.h \
#    src/third-party/GridFluidSim3D/src/diffuseparticle.h \
#    src/third-party/GridFluidSim3D/src/diffuseparticlesimulation.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_dambreak.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_diffuse_inflow.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_export_particle_positions.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_gravity_field.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_hello_world.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_import_mesh_obstacle.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_inflow_outflow.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_lego_sphere_drop.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_save_load_state.h \
#    src/third-party/GridFluidSim3D/src/examples/cpp/example_sphere_drop.h \
#    src/third-party/GridFluidSim3D/src/fluidbrickgrid.h \
#    src/third-party/GridFluidSim3D/src/fluidbrickgridsavestate.h \
#    src/third-party/GridFluidSim3D/src/fluidmaterialgrid.h \
#    src/third-party/GridFluidSim3D/src/fluidsimassert.h \
#    src/third-party/GridFluidSim3D/src/fluidsimulation.h \
#    src/third-party/GridFluidSim3D/src/fluidsimulationsavestate.h \
#    src/third-party/GridFluidSim3D/src/fluidsource.h \
#    src/third-party/GridFluidSim3D/src/fragmentedvector.h \
#    src/third-party/GridFluidSim3D/src/grid3d.h \
#    src/third-party/GridFluidSim3D/src/gridindexkeymap.h \
#    src/third-party/GridFluidSim3D/src/gridindexvector.h \
#    src/third-party/GridFluidSim3D/src/implicitpointprimitive.h \
#    src/third-party/GridFluidSim3D/src/interpolation.h \
#    src/third-party/GridFluidSim3D/src/isotropicparticlemesher.h \
#    src/third-party/GridFluidSim3D/src/kernels/kernels.cpp.in \
#    src/third-party/GridFluidSim3D/src/kernels/kernels.h \
#    src/third-party/GridFluidSim3D/src/levelset.h \
#    src/third-party/GridFluidSim3D/src/logfile.h \
#    src/third-party/GridFluidSim3D/src/macvelocityfield.h \
#    src/third-party/GridFluidSim3D/src/main.h \
#    src/third-party/GridFluidSim3D/src/markerparticle.h \
#    src/third-party/GridFluidSim3D/src/mortonarray3d.h \
#    src/third-party/GridFluidSim3D/src/particleadvector.h \
#    src/third-party/GridFluidSim3D/src/polygonizer3d.h \
#    src/third-party/GridFluidSim3D/src/pressuresolver.h \
#    src/third-party/GridFluidSim3D/src/render/brick_texture_packer/lodepng/lodepng.h \
#    src/third-party/GridFluidSim3D/src/scalarfield.h \
#    src/third-party/GridFluidSim3D/src/spatialpointgrid.h \
#    src/third-party/GridFluidSim3D/src/sphericalfluidsource.h \
#    src/third-party/GridFluidSim3D/src/stopwatch.h \
#    src/third-party/GridFluidSim3D/src/subdividedarray3d.h \
#    src/third-party/GridFluidSim3D/src/triangle.h \
#    src/third-party/GridFluidSim3D/src/trianglemesh.h \
#    src/third-party/GridFluidSim3D/src/turbulencefield.h \
#    src/third-party/GridFluidSim3D/src/utils.h \
#    src/third-party/GridFluidSim3D/src/vmath.h
#    src/sim-fluid-ethanjli/fluidsystem.h \
#    src/sim-fluid-ethanjli/fluidsystem.tpp \
#    src/sim-fluid-ethanjli/math.h \
#    src/sim-fluid-ethanjli/math.tpp \
#    src/sim-fluid-ethanjli/vectorfield.h \
#    src/sim-fluid-ethanjli/vectorfield.tpp


#RESOURCES +=\
#    src/assets/handle.png \
#    src/assets/fault-slip_button.png \
#    src/assets/flowfield_button.png \
#    src/assets/gravity_button.png \
#    src/assets/karst_button.png \
#    src/assets/open_button.png \
#    src/assets/recording_button.png \
#    src/assets/save_button.png \
#    src/assets/tunnel_button.png


##DISTFILES +=

#FORMS +=

DISTFILES += \
    src/Shaders/computeMC.comp \
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
    src/Shaders/rockShader.vert \
 \#    src/sim-fluid-loganzartman/CMakeLists.txt
    src/third-party/GridFluidSim3D/.git/HEAD \
    src/third-party/GridFluidSim3D/.git/config \
    src/third-party/GridFluidSim3D/.git/description \
    src/third-party/GridFluidSim3D/.git/hooks/applypatch-msg.sample \
    src/third-party/GridFluidSim3D/.git/hooks/commit-msg.sample \
    src/third-party/GridFluidSim3D/.git/hooks/fsmonitor-watchman.sample \
    src/third-party/GridFluidSim3D/.git/hooks/post-update.sample \
    src/third-party/GridFluidSim3D/.git/hooks/pre-applypatch.sample \
    src/third-party/GridFluidSim3D/.git/hooks/pre-commit.sample \
    src/third-party/GridFluidSim3D/.git/hooks/pre-push.sample \
    src/third-party/GridFluidSim3D/.git/hooks/pre-rebase.sample \
    src/third-party/GridFluidSim3D/.git/hooks/pre-receive.sample \
    src/third-party/GridFluidSim3D/.git/hooks/prepare-commit-msg.sample \
    src/third-party/GridFluidSim3D/.git/hooks/update.sample \
    src/third-party/GridFluidSim3D/.git/index \
    src/third-party/GridFluidSim3D/.git/info/exclude \
    src/third-party/GridFluidSim3D/.git/logs/HEAD \
    src/third-party/GridFluidSim3D/.git/logs/refs/heads/master \
    src/third-party/GridFluidSim3D/.git/logs/refs/remotes/origin/HEAD \
    src/third-party/GridFluidSim3D/.git/objects/pack/pack-92fd604aa35276ee3371ff5987591a000327f956.idx \
    src/third-party/GridFluidSim3D/.git/objects/pack/pack-92fd604aa35276ee3371ff5987591a000327f956.pack \
    src/third-party/GridFluidSim3D/.git/packed-refs \
    src/third-party/GridFluidSim3D/.git/refs/heads/master \
    src/third-party/GridFluidSim3D/.git/refs/remotes/origin/HEAD \
    src/third-party/GridFluidSim3D/.gitignore \
    src/third-party/GridFluidSim3D/CMakeLists.txt \
    src/third-party/GridFluidSim3D/LICENSE.md \
    src/third-party/GridFluidSim3D/README.md \
    src/third-party/GridFluidSim3D/src/examples/python/example_dambreak.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_diffuse_inflow.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_hello_world.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_inflow_outflow.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_lego_sphere_drop.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_save_load_state.py \
    src/third-party/GridFluidSim3D/src/examples/python/example_sphere_drop.py \
    src/third-party/GridFluidSim3D/src/kernels/scalarfield.cl \
    src/third-party/GridFluidSim3D/src/kernels/tricubicinterpolate.cl \
    src/third-party/GridFluidSim3D/src/pyfluid/__init__.py \
    src/third-party/GridFluidSim3D/src/pyfluid/aabb.py \
    src/third-party/GridFluidSim3D/src/pyfluid/array3d.py \
    src/third-party/GridFluidSim3D/src/pyfluid/config.py \
    src/third-party/GridFluidSim3D/src/pyfluid/fluidsimulation.py \
    src/third-party/GridFluidSim3D/src/pyfluid/fluidsimulationsavestate.py \
    src/third-party/GridFluidSim3D/src/pyfluid/fluidsource.py \
    src/third-party/GridFluidSim3D/src/pyfluid/gridindex.py \
    src/third-party/GridFluidSim3D/src/pyfluid/method_decorators.py \
    src/third-party/GridFluidSim3D/src/pyfluid/pybindings.py \
    src/third-party/GridFluidSim3D/src/pyfluid/pyfluid.py \
    src/third-party/GridFluidSim3D/src/pyfluid/utils.py \
    src/third-party/GridFluidSim3D/src/pyfluid/vector3.py \
    src/third-party/GridFluidSim3D/src/render/blender_scripts/brick_texture_lookup_shader.osl \
    src/third-party/GridFluidSim3D/src/render/blender_scripts/import_animation-basic.py \
    src/third-party/GridFluidSim3D/src/render/blender_scripts/import_animation-brick.py \
    src/third-party/GridFluidSim3D/src/render/blender_scripts/import_animation-diffuse.py \
    src/third-party/GridFluidSim3D/src/render/blender_scripts/import_animation-diffuse_particles_only.py \
    src/third-party/GridFluidSim3D/src/render/brick_texture_packer/README.md \
    src/third-party/GridFluidSim3D/src/render/brick_texture_packer/lodepng/README.md
