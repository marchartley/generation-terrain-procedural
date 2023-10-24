#include "SPHSimulationInterface.h"

SPHSimulationInterface::SPHSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("sphsimulation", "SPH Fluid simulation", "physics", "SPH Simulation", "sph_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[SPH];// = dynamic_cast<SPHSimulation*>(_simulation);
//    _simulation = new SPHSimulation();
}
/*
void SPHSimulationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    this->updateParticles();

    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";

    const char* vParticleShader = "src/Shaders/particle.vert";
    const char* gParticleShader = "src/Shaders/particle.geom";
    const char* fParticleShader = "src/Shaders/particle.frag";

//    this->particlesMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
//    this->particlesMesh.useIndices = false;
//    this->ghostsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
//    this->ghostsMesh.useIndices = false;
//    this->connectionsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
//    this->connectionsMesh.useIndices = false;

    this->particlesMesh = Mesh(std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader));
    this->particlesMesh.useIndices = false;
    this->ghostsMesh = Mesh(std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader));
    this->ghostsMesh.useIndices = false;
//    this->connectionsMesh = Mesh(std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader));
//    this->connectionsMesh.useIndices = false;
}

void SPHSimulationInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    for (int _ = 0; _ < 1; _++)
        simulation.step();
    updateSimulationMesh();

//        ghostsMesh.shader->setVector("color", std::vector<float>{1.f, .5f, 0.f, .4f});
//        ghostsMesh.display(GL_POINTS); //GL_TRIANGLES);
    particlesMesh.shader->setVector("color", std::vector<float>{0.f, .5f, 1.f, .4f});
    particlesMesh.reorderVertices(camPos);
    particlesMesh.display(GL_POINTS); //GL_TRIANGLES);
//        connectionsMesh.shader->setVector("color", std::vector<float>{0.f, 1.f, 0.f, .4f});
//        connectionsMesh.display(GL_LINES);
}

void SPHSimulationInterface::replay(nlohmann::json action)
{
    // ActionInterface::replay(action);
}

QLayout *SPHSimulationInterface::createGUI()
{
    QHBoxLayout* layout = new QHBoxLayout();

    return layout;
}

void SPHSimulationInterface::updateSimulationMesh()
{
    std::cout << timeIt([=]() {
        Vector3 scale = voxelGrid->getDimensions() / simulation.dimensions;

        std::vector<Vector3> trianglesParticles;
        std::vector<Vector3> trianglesGhosts;
//        trianglesParticles.reserve(simulation.nbParticles * CubeMesh::cubesVertices.size());
//        trianglesGhosts.reserve((simulation.particles.size() - simulation.nbParticles) * CubeMesh::cubesVertices.size());
        for (size_t i = 0; i < simulation.particles.size(); i++) {
//            auto cubeTriangles = CubeMesh::cubesVertices;
//            for (auto& vert : cubeTriangles)
//                vert = vert * .5f + simulation.particles[i].position * scale;
            if (simulation.particles[i].isGhost) {
//                trianglesGhosts.insert(trianglesGhosts.end(), cubeTriangles.begin(), cubeTriangles.end());
                trianglesGhosts.push_back(simulation.particles[i].position * scale);
            } else {
//                trianglesParticles.insert(trianglesParticles.end(), cubeTriangles.begin(), cubeTriangles.end());
                trianglesParticles.push_back(simulation.particles[i].position * scale);
            }
        }
        particlesMesh.fromArray(trianglesParticles);
        ghostsMesh.fromArray(trianglesGhosts);
    }) << "ms render" << std::endl;
}

void SPHSimulationInterface::show()
{
    ActionInterface::show();
}

void SPHSimulationInterface::hide()
{
    ActionInterface::hide();
}

void SPHSimulationInterface::afterTerrainUpdated()
{

}

void SPHSimulationInterface::updateParticles()
{
    if (voxelGrid == nullptr) return;

    GridF values = voxelGrid->getVoxelValues();
    values = values.resize(20, 20, 10);
    Mesh m;
    auto triangles = m.applyMarchingCubes(values).getTriangles();
    for (auto& tri : triangles) {
        for (auto& p : tri) {
            p *= simulation.dimensions / values.getDimensions();
        }
    }
    simulation.initialize(triangles);
}
*/
