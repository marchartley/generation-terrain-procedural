#include "EnvObjectFluidSimulation.h"
#include "Interface/InterfaceUtils.h"

EnvObjectFluidSimulation::EnvObjectFluidSimulation(QWidget *parent)
    : ActionInterface("envbjfluids", "Environmental Obj. fluids", "physics", "Display fluid simulation based on Environmental Obj.", "envobjsfluids_button.png", parent)
{
}

void EnvObjectFluidSimulation::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
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

void EnvObjectFluidSimulation::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;
    if (computeAtEachFrame) {
        this->computeSimulation(nbComputationsPerFrame);
        updateSimulationMeshes();
    }

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

void EnvObjectFluidSimulation::replay(nlohmann::json action)
{
    // ActionInterface::replay(action);
}

QLayout *EnvObjectFluidSimulation::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    QCheckBox* displayBoundariesButton = new QCheckBox("Display boundaries");
    QCheckBox* displayGridBoundariesButton = new QCheckBox("Display boundaries (grid)");
    QCheckBox* displayParticlesButton = new QCheckBox("Display particles");
    QCheckBox* displayVectorsButton = new QCheckBox("Display vectors");
    QCheckBox* autoComputeButton = new QCheckBox("Compute at each frame");
    QPushButton* computeButton = new QPushButton("Compute");
    QPushButton* updateMeshButton = new QPushButton("Update terrain");

    layout->addWidget(displayBoundariesButton);
    layout->addWidget(displayGridBoundariesButton);
    layout->addWidget(displayParticlesButton);
    layout->addWidget(displayVectorsButton);
    layout->addWidget(autoComputeButton);
    layout->addWidget(createHorizontalGroup({
                                                computeButton,
                                                updateMeshButton
                                            }));

    displayBoundariesButton->setChecked(this->displayBoundaries);
    displayParticlesButton->setChecked(this->displayParticles);
    displayVectorsButton->setChecked(this->displayVectors);
    autoComputeButton->setChecked(this->computeAtEachFrame);

    QObject::connect(displayBoundariesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayBoundaries = checked; });
    QObject::connect(displayGridBoundariesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayGridBoundaries = checked; });
    QObject::connect(displayParticlesButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayParticles = checked; });
    QObject::connect(displayVectorsButton, &QCheckBox::toggled, this, [=](bool checked) { this->displayVectors = checked; });
    QObject::connect(autoComputeButton, &QCheckBox::toggled, this, [=](bool checked) { this->computeAtEachFrame = checked; });
    QObject::connect(computeButton, &QPushButton::pressed, this, [=]() { computeSimulation(this->nbComputationsPerFrame); });
    QObject::connect(updateMeshButton, &QPushButton::pressed, this, &EnvObjectFluidSimulation::updateBoundariesMesh);

    return layout;
}

void EnvObjectFluidSimulation::updateVectorsMesh()
{
    GridV3 velocities = EnvObject::flowfield.resize(Vector3(20, 20, 1));
    Mesh::createVectorField(velocities, this->voxelGrid->getDimensions(), &vectorsMesh, 1.f, false, true);
}

void EnvObjectFluidSimulation::updateSimulationMeshes()
{
    /*std::cout << */timeIt([=]() {
        this->updateVectorsMesh();
        this->updateParticlesMesh();
//        this->updateBoundariesMesh();
    }); /* << "ms render" << std::endl;*/
}

void EnvObjectFluidSimulation::show()
{
    ActionInterface::show();
}

void EnvObjectFluidSimulation::hide()
{
    ActionInterface::hide();
}

void EnvObjectFluidSimulation::afterTerrainUpdated()
{
    this->updateBoundariesMesh();
}

void EnvObjectFluidSimulation::updateParticlesMesh()
{

}

void EnvObjectFluidSimulation::updateBoundariesMesh()
{
}

void EnvObjectFluidSimulation::computeSimulation(int nbSteps)
{
}
