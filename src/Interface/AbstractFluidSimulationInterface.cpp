#include "AbstractFluidSimulationInterface.h"
#include "Interface/InterfaceUtils.h"

AbstractFluidSimulationInterface::AbstractFluidSimulationInterface(std::string actionTypeName, std::string interfaceName, std::string interfaceType, std::string mainActionDescription, std::string mainActionButtonLogo, QWidget *parent)
    : ActionInterface(actionTypeName, interfaceName, interfaceType, mainActionDescription, mainActionButtonLogo, parent)
{
}

void AbstractFluidSimulationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    const char* vParticleShader = "src/Shaders/particle.vert";
    const char* gParticleShader = "src/Shaders/particle.geom";
    const char* fParticleShader = "src/Shaders/particle.frag";
    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";
    const char* vMCShader = "src/Shaders/MarchingCubes.vert";
    const char* fMCShader = "src/Shaders/no_shader.frag";
    const char* gMCShader = "src/Shaders/MarchingCubes.geom";

    this->particlesMesh = Mesh(std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader));
    this->particlesMesh.useIndices = false;
    this->vectorsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->vectorsMesh.useIndices = false;
//    this->boundariesMesh = Mesh(std::make_shared<Shader>(vMCShader, fMCShader, gMCShader));
    this->boundariesMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->boundariesMesh.useIndices = false;
    this->gridBoundaryMesh = Mesh(std::make_shared<Shader>(vMCShader, fMCShader, gMCShader));
    this->gridBoundaryMesh.useIndices = false;
    this->otherMeshToDisplay = Mesh(std::make_shared<Shader>(vNoShader, fNoShader)); //(vMCShader, fMCShader, gMCShader));
    this->otherMeshToDisplay.useIndices = false;

    this->updateBoundariesMesh();
    this->updateSimulationMeshes();
}

void AbstractFluidSimulationInterface::display(const Vector3& camPos)
{
    if (computeAtEachFrame) {
        this->computeSimulation(nbComputationsPerFrame);
        updateSimulationMeshes();
    }

    if (!this->isVisible())
        return;
//    otherMeshToDisplay.displayWithOutlines(std::vector<float>{0.f, 1.f, 0.f, .4f}, GL_TRIANGLES);
    if (displayBoundaries) {
        boundariesMesh.displayWithOutlines(std::vector<float>{0.f, 1.f, 0.f, .4f}, GL_TRIANGLES);
    }
    if (displayGridBoundaries) {
        Mesh::displayScalarField(this->_simulation->obstacleGrid, gridBoundaryMesh, camPos, {.01f, .2f, .5f, .7f, .99f});
    }
    if (displayParticles) {
        particlesMesh.shader->setVector("color", std::vector<float>{0.f, .5f, 1.f, .4f});
        particlesMesh.reorderVertices(camPos);
        particlesMesh.display(GL_POINTS);
    }
    if (displayVectors) {
        vectorsMesh.shader->setVector("color", std::vector<float>{0.f, .5f, 1.f, .8f});
//        vectorsMesh.reorderLines(camPos);
        vectorsMesh.display(GL_LINES, 3);
    }

//    otherMeshToDisplay.display();
}

void AbstractFluidSimulationInterface::replay(nlohmann::json action)
{
    // ActionInterface::replay(action);
}

QLayout *AbstractFluidSimulationInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    CheckboxElement* onlyAtSurfaceButton = new CheckboxElement("At surface", displayOnlyAtSurface);
    CheckboxElement* displayBoundariesButton = new CheckboxElement("Display boundaries", displayBoundaries);
    CheckboxElement* displayGridBoundariesButton = new CheckboxElement("Display boundaries (grid)", displayGridBoundaries);
    CheckboxElement* displayParticlesButton = new CheckboxElement("Display particles", displayParticles);
    CheckboxElement* displayVectorsButton = new CheckboxElement("Display vectors", displayVectors);
    CheckboxElement* autoComputeButton = new CheckboxElement("Compute at each frame", computeAtEachFrame);
    QPushButton* computeButton = new QPushButton("Compute");
    QPushButton* updateMeshButton = new QPushButton("Update terrain");

    layout->addWidget(onlyAtSurfaceButton->get());
    layout->addWidget(displayBoundariesButton->get());
    layout->addWidget(displayGridBoundariesButton->get());
    layout->addWidget(displayParticlesButton->get());
    layout->addWidget(displayVectorsButton->get());
    layout->addWidget(autoComputeButton->get());
    layout->addWidget(createHorizontalGroup({
                                                computeButton,
                                                updateMeshButton
                                            }));

    /*displayBoundariesButton->setChecked(this->displayBoundaries);
    displayParticlesButton->setChecked(this->displayParticles);
    displayVectorsButton->setChecked(this->displayVectors);
    autoComputeButton->setChecked(this->computeAtEachFrame);

    QObject::connect(displayBoundariesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayBoundaries = checked; });
    QObject::connect(displayGridBoundariesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayGridBoundaries = checked; });
    QObject::connect(displayParticlesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayParticles = checked; });
    QObject::connect(displayVectorsButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayVectors = checked; });
    QObject::connect(autoComputeButton, &QCheckBox::toggled, this, [=](bool checked) { this->computeAtEachFrame = checked; });*/
    QObject::connect(computeButton, &QPushButton::pressed, this, [=]() { computeSimulation(this->nbComputationsPerFrame); });
    QObject::connect(updateMeshButton, &QPushButton::pressed, this, &AbstractFluidSimulationInterface::updateBoundariesMesh);

    return layout;
}

void AbstractFluidSimulationInterface::updateVectorsMesh()
{
    Vector3 resolution(25, 25, 10);
    GridV3 velocities = _simulation->getVelocities(resolution);
    if (displayOnlyAtSurface) {
        GridF surface = voxelGrid->getVoxelValues().resize(resolution).binarize().dilate();
        for (int x = 0; x < resolution.x; x++)
            for (int y = 0; y < resolution.y; y++)
                for (int z = 0; z < 3; z++)
                    surface.at(x, y, z) = 1;
        velocities *= surface;
    }
    velocities.iterateParallel([&](size_t i) { velocities[i].normalize(); });
    Mesh::createVectorField(velocities, this->voxelGrid->getDimensions(), &vectorsMesh, 1.f, false, true);
}

void AbstractFluidSimulationInterface::updateSimulationMeshes()
{
    float renderTime = timeIt([=]() {
        this->updateVectorsMesh();
        this->updateParticlesMesh();
//        this->updateBoundariesMesh();
    });
//    std::cout << showTime(renderTime) << std::endl;
}

void AbstractFluidSimulationInterface::show()
{
    ActionInterface::show();
}

void AbstractFluidSimulationInterface::hide()
{
    ActionInterface::hide();
}

void AbstractFluidSimulationInterface::afterTerrainUpdated()
{
    this->updateBoundariesMesh();
}

void AbstractFluidSimulationInterface::updateParticlesMesh()
{

}

void AbstractFluidSimulationInterface::updateBoundariesMesh()
{
    if (voxelGrid == nullptr) return;
    Vector3 finalDimensions = voxelGrid->getDimensions();

    GridF bigValues = voxelGrid->getVoxelValues();
    GridF values = bigValues; //.resize(20, 20, 10); //.meanSmooth(5, 5, 5); //.resize(100, 100, 10).meanSmooth();

    values.iterateParallel([&](const Vector3& p) {
        values(p) = (values(p) > 0 ? values(p) : Vector3::isInBox(p, Vector3(-1, 10, 10), voxelGrid->getDimensions() - Vector3(-1, 10, -1)) ? -1.f : 1.f);
    });

    AABBox box;
    auto triangles = Mesh::applyMarchingCubes(values).getTriangles();
    for (auto& tri : triangles) {
        for (auto& p : tri) {
            p = p *  _simulation->dimensions / values.getDimensions();
            box.expand(p);
        }
    }

//    std::cout << this->interfaceName << ": From ----- " << box.min() << " ----- to ----- " << box.max() << std::endl;

//    otherMeshToDisplay.fromArray(flattenArray(Mesh::applyMarchingCubes(values.binarize()).getTriangles()));

    _simulation->setObstacles(triangles);
    _simulation->setObstacles(values.binarize());

    auto usedTrianglesIndices = _simulation->obstacleTriangleTree.getAllStoredTrianglesIndices();
//    std::vector<size_t> usedTrianglesIndices(triangles.size());
//    for (size_t i = 0; i < usedTrianglesIndices.size(); i++)
//        usedTrianglesIndices[i] = i;
    std::vector<Vector3> allVertices;
    allVertices.reserve(usedTrianglesIndices.size() * 3);
    for (auto iTriangle : usedTrianglesIndices) {
        const auto& triangle = triangles[iTriangle];
        for (const auto& vertex : triangle) {
            allVertices.push_back(vertex * (finalDimensions / _simulation->dimensions));
        }
    }
    boundariesMesh.fromArray(allVertices);
}

void AbstractFluidSimulationInterface::computeSimulation(int nbSteps)
{
//    this->voxelGrid->computeFlowfield(LBM, nbSteps, this->implicitTerrain.get());
    for (int x = 0; x < _simulation->dimensions.x; x++) {
        for (int y = 0; y < _simulation->dimensions.y; y++) {
            for (int z = 0; z < _simulation->dimensions.z; z++) {
                if (x == 0) {
                    _simulation->addVelocity(x, y, z, Vector3(1, 0, 0));
                }
            }
        }
    }
    for (int i = 0; i < nbSteps; i++) {
        _simulation->step();
    }
}
