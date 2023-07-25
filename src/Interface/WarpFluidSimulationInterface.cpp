#include "WarpFluidSimulationInterface.h"

WarpFluidSimulationInterface::WarpFluidSimulationInterface(QWidget *parent) : AbstractFluidSimulationInterface("WarpFluidSimulation", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[WARP]; // = dynamic_cast<WarpedFluidSimulation*>(_simulation);
//    _simulation = new WarpedFluidSimulation(20, 20, 1);
}

void WarpFluidSimulationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    AbstractFluidSimulationInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    this->computeFromTerrain(voxelGrid.get());
}
    /*

    const char* vParticleShader = "src/Shaders/particle.vert";
    const char* gParticleShader = "src/Shaders/particle.geom";
    const char* fParticleShader = "src/Shaders/particle.frag";
    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";

    this->particlesMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->particlesMesh.useIndices = false;
//    this->debugObstacleMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
//    this->debugObstacleMesh.useIndices = true;
//    this->debugObstacleMesh.cullFace = true;
//    this->updateParticles();
    this->computeFromTerrain(voxelGrid.get());
    this->updateSimulationMesh();
}

void WarpFluidSimulationInterface::display(Vector3 camPos)
{
    if (!this->isVisible())
        return;
//    if (computeAtEachFrame) {
//        for (int _ = 0; _ < 1; _++)
//            simulation.step();
//        updateSimulationMesh();
//    }

//    if (displayBoundaries) {
//        debugObstacleMesh.shader->setVector("color", std::vector<float>{0.f, 1.f, 0.f, .4f});
//        debugObstacleMesh.display(GL_TRIANGLES);
//        debugObstacleMesh.shader->setVector("color", std::vector<float>{0.f, 0.f, 0.f, .4f});
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//        debugObstacleMesh.display(GL_TRIANGLES);
//        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//    }
    particlesMesh.shader->setVector("color", std::vector<float>{0.f, .5f, 1.f, .8f});
    particlesMesh.reorderLines(camPos);
    particlesMesh.display(GL_LINES, 3);
}

void WarpFluidSimulationInterface::replay(nlohmann::json action)
{
    // ActionInterface::replay(action);
}

QLayout *WarpFluidSimulationInterface::createGUI()
{
    QHBoxLayout* layout = new QHBoxLayout();

    QCheckBox* displayBoundariesButton = new QCheckBox("Display boundaries");
//    QCheckBox* autoComputeButton = new QCheckBox("Compute at each frame");

    layout->addWidget(displayBoundariesButton);
//    layout->addWidget(autoComputeButton);

    displayBoundariesButton->setChecked(this->displayBoundaries);
//    autoComputeButton->setChecked(this->computeAtEachFrame);

    QObject::connect(displayBoundariesButton, &QCheckBox::toggled, this, [=](bool cheecked) { this->displayBoundaries = cheecked; });
//    QObject::connect(autoComputeButton, &QCheckBox::toggled, this, [=](bool cheecked) { this->computeAtEachFrame = cheecked; });

    return layout;
}

void WarpFluidSimulationInterface::updateSimulationMesh()
{
    particlesMesh = Mesh::createVectorField(this->velocity.resize(20, 20, 20), this->voxelGrid->getDimensions(), &particlesMesh, 2.f);
    /*
    std::cout << timeIt([=]() {
//        Vector3 scale = voxelGrid->getDimensions() / simulation.dimensions;
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
                trianglesParticles.push_back(positions[i] * voxelGrid->getDimensions());
                colors.push_back(Vector3(particles[i].velocity.norm(), 0, 1));
        }
        particlesMesh.colorsArray = colors;
        particlesMesh.fromArray(trianglesParticles);
        std::cout << Vector3::max(positions) << std::endl;
    }) << "ms render" << std::endl;
    *//*
}*/

void WarpFluidSimulationInterface::computeFromTerrain(TerrainModel *terrain)
{
    auto simulation = dynamic_cast<WarpedFluidSimulation*>(_simulation);
    auto asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    auto asHeightmap = dynamic_cast<Heightmap*>(terrain);
    auto asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    auto asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    Matrix3<float> heights;
    if (asVoxels) {
        heights = Heightmap().fromVoxelGrid(*asVoxels).heights;
    } else if (asHeightmap) {
        heights = asHeightmap->heights;
    } else if (asImplicit) {
        heights = Heightmap().fromImplicit(asImplicit).heights;
    } else if (asLayers) {
        heights = Heightmap().fromLayerGrid(*asLayers).heights;
    } else {
        throw std::domain_error("Was not able to compute WarpFluid because we couldn't cast the terrain model");
    }

    Matrix3<Vector3> gradients = heights.gradient();
    Vector3 flowDir = simulation->mainDirection.normalized();
    auto velocity = Matrix3<Vector3>(gradients.getDimensions());
    for (size_t i = 0; i < gradients.size(); i++) {
        auto& g = gradients[i];
        if (g != Vector3()) {
            float similarity = clamp(g.normalized().dot(flowDir), 0.f, 1.f);
            float slope = clamp(g.norm(), 0.f, 1.f);
            float t = (similarity + slope) / 2.f;
            velocity[i] = lerp(t, flowDir / 2.f, flowDir) - (g.norm2() > 1.f ? g.normalized() : g);
        }
    }
    simulation->velocities = velocity;
//    this->velocity = gradients + simulation.mainDirection;
    std::cout << "Recomputed!" << std::endl;
}
/*
void WarpFluidSimulationInterface::show()
{
    ActionInterface::show();
}

void WarpFluidSimulationInterface::hide()
{
    ActionInterface::hide();
}
*/
void WarpFluidSimulationInterface::afterTerrainUpdated()
{
    this->computeFromTerrain(voxelGrid.get());
}
/*
void WarpFluidSimulationInterface::updateParticles()
{
}
*/
