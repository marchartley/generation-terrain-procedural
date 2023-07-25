#include "FLIPSimulationInterface.h"

using namespace FLIP;

FLIPSimulationInterface::FLIPSimulationInterface(QWidget *parent) : ActionInterface("FLIPSimulation", parent)
{
    simulation = FLIPSimulation(1000, 20, 20, 1, 1, .2f, 2000, 0.1f);
}

void FLIPSimulationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);


    const char* vParticleShader = "src/Shaders/particle.vert";
    const char* gParticleShader = "src/Shaders/particle.geom";
    const char* fParticleShader = "src/Shaders/particle.frag";
    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";

    this->particlesMesh = Mesh(std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader));
    this->particlesMesh.useIndices = false;
    this->debugObstacleMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->debugObstacleMesh.useIndices = true;
    this->debugObstacleMesh.cullFace = true;
    this->updateParticles();
}

void FLIPSimulationInterface::display(Vector3 camPos)
{
    if (!this->isVisible())
        return;
    for (int _ = 0; _ < 1; _++)
        simulation.step();
    updateSimulationMesh();

    debugObstacleMesh.shader->setVector("color", std::vector<float>{0.f, 1.f, 0.f, .4f});
    debugObstacleMesh.display(GL_TRIANGLES);
    debugObstacleMesh.shader->setVector("color", std::vector<float>{0.f, 0.f, 0.f, .4f});
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    debugObstacleMesh.display(GL_TRIANGLES);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    particlesMesh.shader->setVector("color", std::vector<float>{0.f, .5f, 1.f, .4f});
    particlesMesh.reorderVertices(camPos);
    particlesMesh.display(GL_POINTS); //GL_TRIANGLES);
}

void FLIPSimulationInterface::replay(nlohmann::json action)
{
    // ActionInterface::replay(action);
}

QLayout *FLIPSimulationInterface::createGUI()
{
    QHBoxLayout* layout = new QHBoxLayout();

    return layout;
}

void FLIPSimulationInterface::updateSimulationMesh()
{
    std::cout << timeIt([=]() {
        Vector3 scale = voxelGrid->getDimensions() / simulation.dimensions;
        auto& particles = simulation.particles;
        std::vector<Vector3> trianglesParticles;
        std::vector<Vector3> colors;

        Matrix3<float> pressure = Matrix3<float>(simulation.dimensions);
        pressure.data = simulation.p;
        pressure.normalize();

        std::vector<Vector3> positions(simulation.particles.size());
        std::vector<Vector3> velocities(simulation.particles.size());
        Vector3 minVel = Vector3::max(), maxVel = Vector3::min();
        for (size_t i = 0; i < simulation.particles.size(); i++) {
            velocities[i] = simulation.particles[i].velocity;
            minVel = std::min(minVel, velocities[i]);
            maxVel = std::max(maxVel, velocities[i]);
            positions[i] = simulation.particles[i].position / simulation.dimensions;
        }
        for (size_t i = 0; i < simulation.particles.size(); i++)
            velocities[i] = (velocities[i] - minVel) / (maxVel - minVel);


        for (size_t i = 0; i < simulation.particles.size(); i++) {
                trianglesParticles.push_back(particles[i].position * scale);
                colors.push_back(Vector3(particles[i].velocity.norm(), 0, 1));
        }
        particlesMesh.colorsArray = colors;
        particlesMesh.fromArray(trianglesParticles);
    }) << "ms render" << std::endl;
}

void FLIPSimulationInterface::show()
{
    ActionInterface::show();
}

void FLIPSimulationInterface::hide()
{
    ActionInterface::hide();
}

void FLIPSimulationInterface::afterTerrainUpdated()
{

}

void FLIPSimulationInterface::updateParticles()
{
    if (voxelGrid == nullptr) return;

    Matrix3<float> bigValues = voxelGrid->getVoxelValues();
    Matrix3<float> values = bigValues.resize(20, 20, 10);
    Vector3 simDims = simulation.dimensions.xy() + Vector3(0, 0, 1);
    Mesh m;
    auto triangles = m.applyMarchingCubes(values).getTriangles();
    for (auto& tri : triangles) {
        for (auto& p : tri) {
            p *= simDims / values.getDimensions();
        }
    }

    simulation.setObstacles(triangles);

    std::vector<Vector3> allVertices;
    allVertices.reserve(triangles.size() * 3);
    for (size_t iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        const auto& triangle = triangles[iTriangle];
        for (const auto& vertex : triangle) {
            allVertices.push_back(vertex * (bigValues.getDimensions() / simDims));
        }
    }
    debugObstacleMesh.fromArray(allVertices);
//    debugObstacleMesh.colorsArray = std::vector<Vector3>(allVertices.size(), Vector3(0, 1.f, 0));
//    debugObstacleMesh.computeNormals();
}
