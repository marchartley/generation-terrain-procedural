#include "EnvObjsInterface.h"

#include "Interface/HierarchicalListWidget.h"
#include "Interface/InterfaceUtils.h"
#include "TerrainModification/UnderwaterErosion.h"

EnvObjsInterface::EnvObjsInterface(QWidget *parent)
    : ActionInterface("envobjects", "Environmental Objects", "model", "Management of environmental objects generation", "envobjs_button.png", parent)
{

}

void EnvObjsInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";
    const char* vMCShader = "src/Shaders/MarchingCubes.vert";
    const char* fMCShader = "src/Shaders/no_shader.frag";
    const char* gMCShader = "src/Shaders/MarchingCubes.geom";

    this->velocitiesMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->velocitiesMesh.useIndices = false;
    this->highErosionsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->highErosionsMesh.useIndices = false;
    this->highDepositionMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->highDepositionMesh.useIndices = false;

    createEnvObjectsFromImplicitTerrain();
//    recomputeErosionValues();
}

void EnvObjsInterface::display(const Vector3 &camPos)
{
    if (!this->visible)
        return;

    if (displayVelocities) {
        velocitiesMesh.shader->setVector("color", std::vector<float>{.2f, .2f, .8f, .5f});
        velocitiesMesh.reorderLines(camPos);
        velocitiesMesh.display(GL_LINES, 5.f);
    }
    if (displayHighErosions) {
        highErosionsMesh.shader->setVector("color", std::vector<float>{.8f, .2f, .2f, .5f});
        highErosionsMesh.reorderTriangles(camPos);
        highErosionsMesh.display();

        highDepositionMesh.shader->setVector("color", std::vector<float>{.2f, .8f, .2f, .5f});
        highDepositionMesh.reorderTriangles(camPos);
        highDepositionMesh.display();
    }
}

void EnvObjsInterface::replay(nlohmann::json action)
{

}

QLayout *EnvObjsInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout;

    QCheckBox* displayCurrentsButton = new QCheckBox("Display currents");
    // TODO: Simulation type
    QCheckBox* displaySedimentsButton = new QCheckBox("Display sediments");
    QCheckBox* displayHighCurrentsButton = new QCheckBox("Display high currents");
    QCheckBox* displayErosionsButton = new QCheckBox("Display erosion");

    QPushButton* instantiateButton = new QPushButton("Instantiate");
    QPushButton* instantiatePassButton = new QPushButton("Passe");
    QPushButton* instantiateReefButton = new QPushButton("Reef");
    QPushButton* instantiateIslandButton = new QPushButton("Island");
    QPushButton* instantiateBoulderButton = new QPushButton("Boulder");
    QPushButton* instantiateDeltaButton = new QPushButton("Delta");

    QPushButton* recomputeErosionButton = new QPushButton("Erosion values");

    HierarchicalListWidget* objectsListWidget = new HierarchicalListWidget;

    layout->addWidget(displayCurrentsButton);
    layout->addWidget(displaySedimentsButton);
    layout->addWidget(displayHighCurrentsButton);
    layout->addWidget(displayErosionsButton);

    layout->addWidget(instantiateButton);
    layout->addWidget(createVerticalGroup({
                                              instantiatePassButton,
                                              instantiateReefButton,
                                              instantiateIslandButton,
                                              instantiateBoulderButton,
                                              instantiateDeltaButton
                                          }));
    layout->addWidget(objectsListWidget);
    layout->addWidget(recomputeErosionButton);

    QObject::connect(displayCurrentsButton, &QCheckBox::toggled, this, [=](bool checked) { displayVelocities = checked; });
    QObject::connect(displayHighCurrentsButton, &QCheckBox::toggled, this, [=](bool checked) { displayHighCurrents = checked; });
    QObject::connect(displaySedimentsButton, &QCheckBox::toggled, this, [=](bool checked) { displaySediments = checked; });
    QObject::connect(displayErosionsButton, &QCheckBox::toggled, this, [=](bool checked) { displayHighErosions = checked; });
    QObject::connect(recomputeErosionButton, &QPushButton::pressed, this, &EnvObjsInterface::recomputeErosionValues);

    displayCurrentsButton->setChecked(displayVelocities);
    displayHighCurrentsButton->setChecked(displayHighCurrents);
    displaySedimentsButton->setChecked(displaySediments);
    displayErosionsButton->setChecked(displayHighErosions);

    return layout;
}

std::tuple<GridF, GridV3> EnvObjsInterface::extractErosionDataOnTerrain()
{

    TerrainModel *terrain = voxelGrid.get();
    BVHTree boundariesTree;
    GridF densityField;
    GridV3 waterFlowfield = GridV3();
    GridV3 airFlowfield = GridV3();

    int nbPos, nbErosions;
    float particleSimulationTime, terrainModifTime;
    Vector3 terrainDims = terrain->getDimensions();
    std::vector<std::vector<Vector3>> triangles;
    Vector3 geomSize = Vector3::min(terrainDims, Vector3(100, 100, 50));
    triangles = terrain->getGeometry(geomSize).getTriangles();

    boundariesTree = BVHTree();
    boundariesTree.build(Triangle::vectorsToTriangles(triangles));


    std::vector<std::vector<std::pair<float, Vector3>>> allErosions;
    std::vector<BSpline> lastRocksLaunched;

    float erosionSize = 8.f;
    float erosionStrength = .5; // .35f;
    int erosionQtt = 1000;
    float gravity = .981f;
    float bouncingCoefficient = 0.15f; // 1.f;
    float bounciness = 1.f;
    float minSpeed = .1f;
    float maxSpeed = 5.f;
    float maxCapacityFactor = 1.f;
    float erosionFactor = 1.f;
    float depositFactor = 1.f;
    float matterDensity = 500.f;
    float materialImpact = 1.f;

    float airFlowfieldRotation = 270.f;
    float waterFlowfieldRotation = 90.f;

    float airForce = .5f;
    float waterForce = 0.f;

    float dt = 1.f;

    float shearingStressConstantK = 1.f;
    float shearingRatePower = .5f;
    float erosionPowerValue = 1.f;
    float criticalShearStress = .8f;

    bool wrapParticles = true;

    float initialCapacity = .0f;

    FluidSimType selectedSimulationType = WARP;
    dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[selectedSimulationType])->recomputeVelocities();

    std::vector<std::pair<Vector3, Vector3>> initialPositionsAndDirections(erosionQtt);
    for (auto& [pos, dir] : initialPositionsAndDirections) {
        pos = Vector3::random(Vector3(), terrainDims);
    }

    UnderwaterErosion::EROSION_APPLIED applyOn = UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS;
    UnderwaterErosion::FLOWFIELD_TYPE flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS;
    UnderwaterErosion::DENSITY_TYPE densityUsed = UnderwaterErosion::DENSITY_TYPE::NATIVE;

    UnderwaterErosion erod(voxelGrid.get(), erosionSize, erosionStrength, erosionQtt);
    std::tie(lastRocksLaunched, nbPos, nbErosions, allErosions) = erod.Apply(applyOn,
                                                                terrain,
                                                                boundariesTree,
                                                                particleSimulationTime, terrainModifTime,
                                                                Vector3(false),
                                                                Vector3(false),
                                                                0.f,
                                                                true,
                                                                gravity,
                                                                bouncingCoefficient,
                                                                bounciness,
                                                                minSpeed,
                                                                maxSpeed,
                                                                maxCapacityFactor,
                                                                erosionFactor,
                                                                depositFactor,
                                                                matterDensity, // + .1f,
                                                                materialImpact,
                                                                airFlowfieldRotation,
                                                                waterFlowfieldRotation,
                                                                airForce,
                                                                waterForce,
                                                                dt,
                                                                shearingStressConstantK,
                                                                shearingRatePower,
                                                                erosionPowerValue,
                                                                criticalShearStress,
                                                                initialPositionsAndDirections,
                                                                flowfieldUsed,
                                                                waterFlowfield,
                                                                airFlowfield,
                                                                densityUsed,
                                                                densityField,
                                                                initialCapacity,
                                                                selectedSimulationType,
                                                                wrapParticles,
                                                                false
                                                                );

    GridF erosionsAmount(terrainDims);
    std::cout << "Just for the erosion part: " << showTime(timeIt([&]() {
        float size = erosionSize;
        std::vector<GridF> suberosions(allErosions.size(), erosionsAmount);
        #pragma omp parallel for
        for (size_t i = 0; i < suberosions.size(); i++) {
            for (auto& [val, pos] : allErosions[i]) {
                RockErosion(size, val).computeErosionMatrix(suberosions[i], pos - Vector3(.5f, .5f, .5f));
            }
        }
        for (auto& sub : suberosions)
            erosionsAmount += sub;
    })) << std::endl;

    GridV3 velocities(terrainDims);
    GridF evaluationAmounts(terrainDims);
    for (const auto& path : lastRocksLaunched) {
        for (int i = 1; i < int(path.points.size()) - 1; i++) {
            auto& pos = path.points[i];
            auto& pPrev = path.points[i - 1];
            auto& pNext = path.points[i + 1];
            auto velocity = (pNext - pPrev).normalize();
            if (velocity.norm2() > 50 * 50) continue; // When a particle is wraped from one side to the other of the terrain
            velocities(pos) += velocity * .5f;
            evaluationAmounts(pos) ++;
            velocities(pPrev) += velocity * .5f;
            evaluationAmounts(pPrev) ++;
            velocities(pNext) += velocity * .5f;
            evaluationAmounts(pNext) ++;
        }
    }
    for (size_t i = 0; i < velocities.size(); i++) {
        if (evaluationAmounts[i] == 0) continue;
        velocities[i] /= evaluationAmounts[i];
    }
    return {erosionsAmount, velocities};
}

void EnvObjsInterface::createEnvObjectsFromImplicitTerrain()
{
    EnvObject::removeAllObjects();
    auto tunnelsPatches = implicitTerrain->findAll(ImplicitPatch::ParametricTunnel);
    for (auto& tunnelPatch : tunnelsPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(tunnelPatch);
        if (asPrimitive && asPrimitive->material == WATER) {
            BSpline curve = asPrimitive->optionalCurve;
            for (auto& p : curve) {
                p = asPrimitive->getGlobalPositionOf(p);
            }
            EnvCurve* passe = dynamic_cast<EnvCurve*>(EnvObject::instantiate("passe"));
            passe->curve = curve;
        }
    }

    auto reefPatches = implicitTerrain->findAll(ImplicitPatch::MountainChain);
    for (auto& reefPatch : reefPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(reefPatch);
        if (asPrimitive) {
            BSpline curve = asPrimitive->optionalCurve;
            for (auto& p : curve) {
                p = asPrimitive->getGlobalPositionOf(p);
            }
//            EnvArea* reef = dynamic_cast<EnvArea*>(EnvObject::instantiate("reef"));
//            reef->area = curve;
        }
    }

    auto lagoonPatches = implicitTerrain->findAll(ImplicitPatch::Polygon);
    for (auto& lagoonPatch : lagoonPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(lagoonPatch);
        if (asPrimitive /* && asPrimitive->material == WATER*/) {
            ShapeCurve curve = asPrimitive->optionalCurve;
            for (auto& p : curve) {
                p = asPrimitive->getGlobalPositionOf(p);
            }
            EnvArea* lagoon = dynamic_cast<EnvArea*>(EnvObject::instantiate("lagoon"));
            lagoon->area = curve;
        }
    }

}

void EnvObjsInterface::show()
{
    ActionInterface::show();
}

void EnvObjsInterface::hide()
{
    ActionInterface::hide();
}

void EnvObjsInterface::afterTerrainUpdated()
{

}

void EnvObjsInterface::instantiateObject()
{
    /*
    voxelGrid->computeFlowfield(LBM, 1, implicitTerrain.get());
    GridF terrainSurface = voxelGrid->getVoxelValues().binarize();
    GridV3 waterFlow = voxelGrid->getFlowfield().resize(Vector3(10, 10, 10), RESIZE_MODE::MAX_VAL).resize(terrainSurface.getDimensions());
    GridV3 surfaceNormals = terrainSurface.toDistanceMap().gradient();
    for (size_t i = 0; i < surfaceNormals.size(); i++)
        surfaceNormals[i].normalize();

    terrainSurface = terrainSurface - terrainSurface.erode();
    std::vector<Vector3> availablePositions;
    availablePositions.reserve(terrainSurface.sizeX * terrainSurface.sizeY * 2);
    for (size_t i = 0; i < terrainSurface.size(); i++)
        if (terrainSurface[i]) availablePositions.push_back(terrainSurface.getCoordAsVector3(i));

    // Get 1000 random points
    std::shuffle(availablePositions.begin(), availablePositions.end(), random_gen::random_generator);
    availablePositions.resize(std::min(1000, int(availablePositions.size())));


    if (this->implicitTerrain->findAll(ImplicitPatch::ParametricTunnel).empty()) {
        this->autoGeneratePasse(terrainSurface, waterFlow, surfaceNormals, availablePositions);
    } else if (this->implicitTerrain->findAll(ImplicitPatch::ParametricTunnel).size() == 1) {
        this->autoGenerateDelta(terrainSurface, waterFlow, surfaceNormals, availablePositions);
    } else {
        this->autoGenerateMotu(terrainSurface, waterFlow, surfaceNormals, availablePositions);
    }*/
}

void EnvObjsInterface::recomputeErosionValues()
{
    auto [erosion, velocities] = this->extractErosionDataOnTerrain();
    if  (erosionGrid.getDimensions() == erosion.getDimensions()) {
        erosionGrid = (erosionGrid + erosion) * .5f;
        velocitiesGrid = (velocitiesGrid + velocities) * .5f;
    } else {
        erosionGrid  = erosion;
        velocitiesGrid = velocities;
    }

    float iso = .01f;
//    Mesh::createVectorField(velocitiesGrid.resize(.2f), voxelGrid->getDimensions(), &velocitiesMesh);
    highErosionsMesh.fromArray(flattenArray(Mesh::applyMarchingCubes((-erosionGrid) - iso).getTriangles()));
    highDepositionMesh.fromArray(flattenArray(Mesh::applyMarchingCubes(erosionGrid - iso).getTriangles()));

    this->updateEnvironmentFromEnvObjects();
}

void EnvObjsInterface::updateEnvironmentFromEnvObjects()
{
    // Get original flowfield, do not accumulate effects (for now).
    EnvObject::flowfield = dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->getVelocities(EnvObject::flowfield.sizeX, EnvObject::flowfield.sizeY, EnvObject::flowfield.sizeZ);
    EnvObject::applyEffects();
    Mesh::createVectorField(EnvObject::flowfield, voxelGrid->getDimensions(), &velocitiesMesh, -1, false, true);
}
