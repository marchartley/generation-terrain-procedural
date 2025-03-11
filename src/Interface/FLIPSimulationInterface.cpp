#include "FLIPSimulationInterface.h"

FLIPSimulationInterface::FLIPSimulationInterface(QWidget *parent)
    : AbstractFluidSimulationInterface("flipsimulation", "FLIP simulation", "physics", "FLIP Simulation", "flip_button.png", parent)
{
    _simulation = GlobalTerrainProperties::get()->simulations[FLIP];
    simulation = dynamic_cast<FLIPSimulation*>(_simulation);
//    *_simulation = FLIPSimulation(1000, 20, 20, 1, 1, .2f, 20000, 0.1f);
    nbComputationsPerFrame = 1;
    computeAtEachFrame = false;
    displayParticles = true;
    displayBoundaries = false;
    displayVectors = false;
    displayOnlyAtSurface = false;
}

QLayout *FLIPSimulationInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout();

    InterfaceUI* ui = new InterfaceUI(new QVBoxLayout, "FLIP");


    auto resetButton = new ButtonElement("Reset", [&](){ this->resetParticles(); });
    auto particleCountSlider = new SliderElement("# particles", 5, 10000, 1);
    particleCountSlider->slider()->setfValue(this->simulation->maxParticles);
    particleCountSlider->setOnValueChanged([&](float newValue) { this->simulation->maxParticles = newValue; });
    auto particleRadiusSlider = new SliderElement("Particle radius", .1f, 3.f, .01f);
    particleRadiusSlider->slider()->setfValue(this->simulation->particleRadius);
    particleRadiusSlider->setOnValueChanged([&](float newValue) { this->simulation->particleRadius = newValue; });
    auto particleDensitySlider = new SliderElement("Particle density", 1.f, 2000.f, 1.f);
    particleDensitySlider->slider()->setfValue(this->simulation->density);
    particleDensitySlider->setOnValueChanged([&](float newValue) { this->simulation->density = newValue; });
    auto particleFLIPRatioSlider = new SliderElement("FLIP ratio", 0.f, 1.f, .05f);
    particleFLIPRatioSlider->slider()->setfValue(this->simulation->flipRatio);
    particleFLIPRatioSlider->setOnValueChanged([&](float newValue) { this->simulation->flipRatio = newValue; });
    auto averagingSlider = new SliderElement("Averaging", 0.f, 100.f, 1.f);
    averagingSlider->slider()->setfValue(this->simulation->averaging);
    averagingSlider->setOnValueChanged([&](float newValue) { this->simulation->averaging = newValue; });
    /*auto spacingSlider = new SliderElement("Spacing", 0.1f, 1.f, .01f);
    spacingSlider->slider()->setfValue(1.f / this->simulation->fInvSpacing);
    spacingSlider->setOnValueChanged([&](float newValue) { this->simulation->fInvSpacing = 1.f / newValue; });*/
    auto driftCheckbox = new CheckboxElement("Compensate drift?", simulation->compensateDrift);
    auto overrelaxationSlider = new SliderElement("Overrelaxation", 0.01f, 2.f, 0.01f, simulation->overRelaxation);
    auto numIterationsSlider = new SliderElement("Iterations", 0, 30, 1);
    numIterationsSlider->slider()->setValue(simulation->numIterations);
    numIterationsSlider->setOnValueChanged([&](float newVal) { simulation->numIterations = newVal; });
    ui->add({
                resetButton,
                particleCountSlider,
                particleRadiusSlider,
                particleDensitySlider,
                particleFLIPRatioSlider,
                averagingSlider,
//                spacingSlider,
                driftCheckbox,
                overrelaxationSlider,
                numIterationsSlider
            });
    ui->add(AbstractFluidSimulationInterface::createGUI()); // Still use the default system

    layout->addWidget(ui->get());

    return layout;
}

void FLIPSimulationInterface::updateParticlesMesh()
{
    auto simulation = dynamic_cast<FLIPSimulation*>(_simulation);
    auto& particles = simulation->particles;
    std::vector<Vector3> particlePositions;
    std::vector<Vector3> colors;

//    GridF pressure = GridF(simulation->dimensions);
    GridF pressure = simulation->p;
    pressure.normalize();

    std::vector<Vector3> positions(simulation->particles.size());
    std::vector<Vector3> velocities(simulation->particles.size());
    Vector3 minVel = Vector3::max(), maxVel = Vector3::min();
    for (size_t i = 0; i < simulation->particles.size(); i++) {
        velocities[i] = simulation->particles[i].velocity;
        minVel = std::min(minVel, velocities[i]);
        maxVel = std::max(maxVel, velocities[i]);
        positions[i] = simulation->particles[i].position / simulation->dimensions;
    }
    for (size_t i = 0; i < simulation->particles.size(); i++)
        velocities[i] = (velocities[i] - minVel) / (maxVel - minVel);


    for (size_t i = 0; i < simulation->particles.size(); i++) {
            particlePositions.push_back(positions[i] * voxelGrid->getDimensions());
            if (simulation->savedState.size() > i)
                colors.push_back(Vector3((particles[i].position - simulation->savedState[i].position).norm() / simulation->dt, 0, 1));
    }
    particlesMesh.colorsArray = colors;
    particlesMesh.fromArray(particlePositions);
}

void FLIPSimulationInterface::resetParticles()
{
    this->simulation->reset();
}

void FLIPSimulationInterface::display(const Vector3 &camPos)
{
    AbstractFluidSimulationInterface::display(camPos);
    if (!this->isVisible()) return;

    gridBoundaryMesh.shader->setVector("scale", voxelGrid->getDimensions() / simulation->particleDensity.getDimensions());
//    Mesh::displayScalarField((simulation->p / simulation->p.max()).resize(voxelGrid->getDimensions() / 2.f), gridBoundaryMesh, camPos, {0.f, .1f, .2f, .3f, .4f, .5f, .6f, .7f, .8f, .9f, 1.f}, voxelGrid->getDimensions());
    GridV3 velocitiesFromUVW(simulation->dimensions);
    velocitiesFromUVW.iterateParallel([&](const Vector3& pos) {
        velocitiesFromUVW(pos) = Vector3(simulation->u.interpolate(pos + Vector3(.5f, .0f, .0f)), simulation->v.interpolate(pos + Vector3(.0f, .5f, .0f)), simulation->w.interpolate(pos + Vector3(.0f, .0f, .5f)));
    });
    Mesh::createVectorField(velocitiesFromUVW, velocitiesFromUVW.getDimensions(), &vectorsMesh, -1, false, true);

//    vectorsMesh.clear();
    std::vector<Vector3> velocitiesFromParticles = vectorsMesh.vertexArray;
    for (size_t i = 0; i < simulation->particles.size(); i++) {
        Vector3 before = simulation->savedState[i].position;
        Vector3 after = simulation->particles[i].position;
        auto points = Mesh::getPointsForArrow(before, before + (after - before).normalized());
        velocitiesFromParticles.insert(velocitiesFromParticles.end(), points.begin(), points.end());
    }
    vectorsMesh.fromArray(velocitiesFromParticles);
    vectorsMesh.scale(voxelGrid->getDimensions() / simulation->dimensions);
    vectorsMesh.display(GL_LINES);
//    std::cout << simulation->particleDensity.min() << " -- " << simulation->particleDensity.max() << " -- " << simulation->particleDensity.sum() / float(simulation->particleDensity.size()) << std::endl;
}


void FLIPSimulationInterface::updateParticles()
{
}

void FLIPSimulationInterface::computeSimulation(int nbSteps)
{
    Vector3 errorsDimensions(30, 30, 30);
    this->_currentVelocities = _simulation->getVelocities(errorsDimensions);
    GridV3 oldVelocities;
    float MSE;
    for (int i = 0; i < nbSteps; i++) {
        oldVelocities = _currentVelocities;
        simulation->step();
        this->_currentVelocities = _simulation->getVelocities(errorsDimensions);
        MSE = 0.f;
        oldVelocities.iterate([&](size_t i) {
            MSE += (oldVelocities[i] - _currentVelocities[i]).norm2();
        });
        if (_simulation->currentStep > 5) {
            changesHistory.push_back(MSE / float(oldVelocities.size()));
        }
    }
}

void FLIPSimulationInterface::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_V) {
        GridF data(changesHistory.size(), 1, 1);
        data.data = changesHistory;
        Plotter::get()->reset();
        Plotter::get()->addPlot(data.medianBlur(20, 3, 3, true).data);
        Plotter::get()->show();
        e->accept();
    }
    return AbstractFluidSimulationInterface::keyPressEvent(e);
}
