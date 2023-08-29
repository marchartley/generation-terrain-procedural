#include "FlowFieldInterface.h"

FlowFieldInterface::FlowFieldInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("stablefluidssimulation", "Stable Fluids Simulation", "physics", "Stable fluids simulation", "flowfield_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[STABLE]; // = dynamic_cast<StableFluidsSimulation*>(_simulation);
//    *_simulation = StableFluidsSimulation();
    this->nbComputationsPerFrame = 5;
}


/*
FlowFieldInterface::FlowFieldInterface(QWidget *parent) : ActionInterface("flowfield", parent)
{

}

void FlowFieldInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    this->displayedFlowfields = std::vector<bool>(voxelGrid->multipleFlowFields.size(), false);
    if (!this->displayedFlowfields.empty()) this->displayedFlowfields[0] = true;
    this->flowMeshes = std::vector<Mesh>(voxelGrid->multipleFlowFields.size());

    // Waiting for OpenGL to be available to create the shaders...
    std::string shaderPath = "src/Shaders/";
    std::string vertexFilename = shaderPath + "MarchingCubes.vert";
    std::string geometryFilename = shaderPath + "MarchingCubes.geom";
    std::string fragmentFilename = shaderPath + "no_shader.frag";

    std::string noShaderVertexFilename = shaderPath + "no_shader.vert";
    std::string noShaderFragmentFilename = shaderPath + "no_shader.frag";

    pressureDensityMesh.shader = std::make_shared<Shader>(vertexFilename, fragmentFilename, geometryFilename);
    pressureDensityMesh.useIndices = false;
    pressureDensityMesh.cullFace = true;

    sumOfFlowsMesh.shader = std::make_shared<Shader>(noShaderVertexFilename, noShaderFragmentFilename);
    sumOfFlowsMesh.useIndices = false;
    sumOfFlowsMesh.cullFace = true;

    for (auto& flowMesh : this->flowMeshes) {
        flowMesh.shader = std::make_shared<Shader>(noShaderVertexFilename, noShaderFragmentFilename);
        flowMesh.useIndices = false;
        flowMesh.cullFace = true;
    }
}

void FlowFieldInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    this->displayPressureDensities(camPos);
    this->displayFlows(camPos);
    this->displaySumOfFlows(camPos);
    if (this->computeAtEachFrame)
        this->recomputeFlowfield(10);
}

void FlowFieldInterface::displayPressureDensities(const Vector3& camPos)
{
    if (!this->displayingPressure)
        return;

    std::vector<float> isoValues = {0.5f}; // {0.9f, 0.5f, 0.2f};
    pressureDensityMesh.displayAsScalarField(pressureDensityVoxels, camPos, isoValues);


//    std::vector<Vector3> positions(pressureDensityVoxels.size());
//    for (size_t i = 0; i < positions.size(); i++) {
//        positions[i] = pressureDensityVoxels.getCoordAsVector3(i);
//    }
//    pressureDensityMesh.fromArray(positions);
//    pressureDensityMesh.update();
//    GlobalsGL::f()->glBindVertexArray(pressureDensityMesh.vao);
//    pressureDensityMesh.shader->setTexture3D("dataFieldTex", 0, pressureDensityVoxels + .5f);
//    pressureDensityMesh.shader->setInt("dataFieldTex", 0);
//    pressureDensityMesh.shader->setInt("edgeTableTex", 1);
//    pressureDensityMesh.shader->setInt("triTableTex", 2);
//    pressureDensityMesh.shader->setFloat("isolevel", 0.f);
//    pressureDensityMesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
//    pressureDensityMesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
//    pressureDensityMesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
//    pressureDensityMesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
//    pressureDensityMesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
//    pressureDensityMesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
//    pressureDensityMesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
//    pressureDensityMesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
//    pressureDensityMesh.shader->setBool("useMarchingCubes", true);
//    //Edge Table texture//
//    //This texture store the 256 different configurations of a marching cube.
//    //This is a table accessed with a bitfield of the 8 cube edges states
//    //(edge cut by isosurface or totally in or out).
//    //(cf. MarchingCubes.cpp)

//    GLuint edgeTableTex, triTableTex;
//    GlobalsGL::f()->glGenTextures(1, &edgeTableTex);
//    GlobalsGL::f()->glActiveTexture(GL_TEXTURE1);
//    glEnable(GL_TEXTURE_2D);
//    glBindTexture(GL_TEXTURE_2D, edgeTableTex);
//    //Integer textures must use nearest filtering mode

//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
//    //We create an integer texture with new GL_EXT_texture_integer formats
//    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 256, 1, 0,
//    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::cubeEdges));

//    //Triangle Table texture//
//    //This texture store the vertex index list for
//    //generating the triangles of each configurations.
//    //(cf. MarchingCubes.cpp)

//    glGenTextures(1, &triTableTex);
//    GlobalsGL::f()->glActiveTexture(GL_TEXTURE2);
//    glEnable(GL_TEXTURE_2D);
//    glBindTexture(GL_TEXTURE_2D, triTableTex);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
//    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 16, 256, 0,
//    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::triangleTable));


//    GlobalsGL::f()->glActiveTexture(GL_TEXTURE0);

//    // Ignore parameters to hide some voxels
////    pressureDensityMesh.shader->setVector("min_vertice_positions", Vector3::min());
////    pressureDensityMesh.shader->setVector("max_vertice_positions", Vector3::max());
////    pressureDensityMesh.shader->setFloat("min_isolevel", -1000.f); // 3.5f);
////    pressureDensityMesh.shader->setFloat("max_isolevel", 1000.f);

//    for (size_t i = 0; i < isoValues.size(); i++) {
//        float iso = isoValues[i];
//        Vector3 color = HSVtoRGB(i / float(isoValues.size()), 1.f, 1.f);
//        pressureDensityMesh.shader->setVector("color", std::vector<float> {color.x, color.y, color.z, .3f});
//        pressureDensityMesh.shader->setFloat("isolevel", iso);

//        // display the mesh
//        pressureDensityMesh.reorderVertices(camPos);
//        pressureDensityMesh.display(GL_POINTS);
//    }
}

void FlowFieldInterface::displayFlows(const Vector3& camPos)
{
    for (size_t i = 0; i < this->displayedFlowfields.size(); i++) {
        if (!this->displayedFlowfields[i])
            continue;

        Vector3 color = HSVtoRGB(float(i) / float(displayedFlowfields.size() - 1), 1.f, 1.f);
        Mesh& flowMesh = this->flowMeshes[i];
        if (flowMesh.shader != nullptr)
            flowMesh.shader->setVector("color", std::vector<float>{color.x, color.y, color.z, .5f}); //std::vector<float>({0.f, 0.f, 1.f, .4f}));
        flowMesh.reorderLines(camPos);
        flowMesh.display(GL_LINES, 5.f);
    }
}

void FlowFieldInterface::displaySumOfFlows(const Vector3& camPos)
{
    if (!this->displayingSumOfFlows)
        return;
    Vector3 color = HSVtoRGB(1.f, 1.f, 1.f);
    if (this->sumOfFlowsMesh.shader != nullptr)
        this->sumOfFlowsMesh.shader->setVector("color", std::vector<float>{color.x, color.y, color.z, .5f}); //std::vector<float>({0.f, 0.f, 1.f, .4f}));
    this->sumOfFlowsMesh.reorderLines(camPos);
    this->sumOfFlowsMesh.display(GL_LINES, 5.f);
}

void FlowFieldInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        this->voxelGrid->computeFlowfield(LBM);
    }
}

void FlowFieldInterface::recomputeFlowfield(int steps)
{
    this->voxelGrid->computeMultipleFlowfields(LBM, steps, this->implicitTerrain.get());

    totalFlow = voxelGrid->getFlowfield(LBM);

//    this->totalFlow = this->voxelGrid->getFlowfield(0);
//    for (size_t i = 1; i < this->voxelGrid->multipleFlowFields.size(); i++)
//        this->totalFlow += this->voxelGrid->getFlowfield(i);

    pressureDensityVoxels = -this->totalFlow.divergence();
//    }
    pressureDensityVoxels.normalize();

    this->updateFlowfieldDebugMesh();

    this->addTerrainAction(nlohmann::json({
                                              {}
                                          }));
    if (steps == 30)
        Q_EMIT this->updated();
}

void FlowFieldInterface::updateFlowfieldDebugMesh()
{
//    float maxPressure = this->pressureDensityVoxels.max();
    float maxFlow = 0.f;
    for (size_t i = 0; i < this->flowMeshes.size(); i++) {
//        Mesh::createVectorField(this->voxelGrid->getFlowfield(i).resize(this->voxelGrid->multipleFluidSimulations[0].dimensions) / 100.f, this->voxelGrid->getDimensions(), &flowMeshes[i]);
//        maxFlow = std::max(maxFlow, this->voxelGrid->getFlowfield(i).resize(this->voxelGrid->multipleFluidSimulations[0].dimensions).max().norm());
        Mesh::createVectorField(this->voxelGrid->getFlowfield(LBM).resize(this->voxelGrid->multipleFluidSimulations[0].dimensions) / 100.f, this->voxelGrid->getDimensions(), &flowMeshes[i]);
        maxFlow = std::max(maxFlow, this->voxelGrid->getFlowfield(LBM).resize(this->voxelGrid->multipleFluidSimulations[0].dimensions).max().norm());
    }
    std::cout << "Max flow: " << maxFlow << std::endl;

    Mesh::createVectorField(this->totalFlow.resize(1.f / this->voxelGrid->fluidSimRescale) / 100.f, this->voxelGrid->getDimensions(), &sumOfFlowsMesh);
}

void FlowFieldInterface::hide()
{
//    this->flowMesh.hide();
    CustomInteractiveObject::hide();
}

void FlowFieldInterface::show()
{
//    this->flowMesh.show();
    CustomInteractiveObject::show();
}

QLayout* FlowFieldInterface::createGUI()
{
    QVBoxLayout* flowFieldLayout = new QVBoxLayout;

    QPushButton* flowFieldComputeButton = new QPushButton("Calculer");
    QCheckBox* flowFieldDisplayPressureButton = new QCheckBox("Display pressure");
    QCheckBox* flowFieldDisplayAllFlowsButton = new QCheckBox("Display sum");
    QCheckBox* autoComputeFlowsButton = new QCheckBox("Compute each frame");

    for (size_t i = 0; i < this->displayedFlowfields.size(); i++) {
        QCheckBox* displayFlow = new QCheckBox(QString::fromStdString("Flow #" + std::to_string(i + 1)));
        displayFlow->setChecked(this->displayedFlowfields[i]);
        flowFieldLayout->addWidget(displayFlow);

        QObject::connect(displayFlow, &QCheckBox::toggled, this, [=](bool checked) { this->displayedFlowfields[i] = checked; } );
    }


    flowFieldLayout->addWidget(flowFieldDisplayPressureButton);
    flowFieldLayout->addWidget(flowFieldDisplayAllFlowsButton);
    flowFieldLayout->addWidget(flowFieldComputeButton);
    flowFieldLayout->addWidget(autoComputeFlowsButton);

    flowFieldDisplayPressureButton->setChecked(this->displayingPressure);
    flowFieldDisplayAllFlowsButton->setChecked(this->displayingSumOfFlows);
    autoComputeFlowsButton->setChecked(this->computeAtEachFrame);
    QObject::connect(flowFieldDisplayPressureButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayingPressure = checked; } );
    QObject::connect(flowFieldDisplayAllFlowsButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayingSumOfFlows = checked; } );
    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=]() { this->recomputeFlowfield(30); } );
    QObject::connect(autoComputeFlowsButton, &QCheckBox::toggled, this, [=](bool checked) { this->computeAtEachFrame = checked; } );

    return flowFieldLayout;
}
*/
