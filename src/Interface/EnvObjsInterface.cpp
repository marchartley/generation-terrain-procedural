#include "EnvObjsInterface.h"

#include "Interface/InterfaceUtils.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/CommonInterface.h"

EnvObjsInterface::EnvObjsInterface(QWidget *parent)
    : ActionInterface("envobjects", "Environmental Objects", "model", "Management of environmental objects generation", "envobjs_button.png", parent)
{

    QTimer* hotreloadTimer = new QTimer(this);
    hotreloadTimer->setInterval(500);
    QObject::connect(hotreloadTimer, &QTimer::timeout, this, &EnvObjsInterface::hotReloadFile);
    hotreloadTimer->start();
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
    this->objectsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->objectsMesh.useIndices = false;

    createEnvObjectsFromImplicitTerrain();
//    recomputeErosionValues();
    EnvObject::precomputeTerrainProperties(*heightmap);
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
    objectsMesh.shader->setVector("color", std::vector<float>{.8f, .2f, .8f, .5f});
    objectsMesh.display(GL_LINES, 5.f);
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
//    QPushButton* instantiatePassButton = new QPushButton("Passe");
//    QPushButton* instantiateReefButton = new QPushButton("Reef");
//    QPushButton* instantiateIslandButton = new QPushButton("Island");
//    QPushButton* instantiateBoulderButton = new QPushButton("Boulder");
//    QPushButton* instantiateDeltaButton = new QPushButton("Delta");

    QPushButton* recomputeErosionButton = new QPushButton("Erosion values");

    ButtonElement* spendTimeButton = new ButtonElement("Wait");
    spendTimeButton->setOnClick([&]() { this->updateEnvironmentFromEnvObjects(); });

    ButtonElement* showDepositionButton = new ButtonElement("Show deposition");
    showDepositionButton->setOnClick([&]() { this->displaSedimentsDistrib(); });

//    if (objectsListWidget != nullptr)
//        delete objectsListWidget;
    objectsListWidget = new HierarchicalListWidget;
    updateObjectsList();
    QObject::connect(objectsListWidget, &HierarchicalListWidget::itemClicked, this, &EnvObjsInterface::updateObjectsListSelection);

    layout->addWidget(displayCurrentsButton);
//    layout->addWidget(displaySedimentsButton);
//    layout->addWidget(displayHighCurrentsButton);
    layout->addWidget(displayErosionsButton);

    layout->addWidget(spendTimeButton->get());
    layout->addWidget(showDepositionButton->get());

    layout->addWidget(instantiateButton);

    std::vector<ButtonElement*> probaButtons;
    for (auto& [name, obj] : EnvObject::availableObjects) {
        ButtonElement* showButton = new ButtonElement("Show " + obj->name);
        ButtonElement* forceButton = new ButtonElement("Force");
        showButton->setOnClick([&](){ this->displayProbas(name); });
        forceButton->setOnClick([&](){ this->instantiateSpecific(name); });
//        probaButtons.push_back(showButton);
        layout->addWidget(createHorizontalGroupUI({showButton, forceButton})->get());
    }
//    layout->addWidget(createVerticalGroupUI(probaButtons));
    layout->addWidget(objectsListWidget);
    layout->addWidget(recomputeErosionButton);

    QObject::connect(displayCurrentsButton, &QCheckBox::toggled, this, [=](bool checked) { displayVelocities = checked; });
    QObject::connect(displayHighCurrentsButton, &QCheckBox::toggled, this, [=](bool checked) { displayHighCurrents = checked; });
    QObject::connect(displaySedimentsButton, &QCheckBox::toggled, this, [=](bool checked) { displaySediments = checked; });
    QObject::connect(displayErosionsButton, &QCheckBox::toggled, this, [=](bool checked) { displayHighErosions = checked; });
    QObject::connect(recomputeErosionButton, &QPushButton::pressed, this, &EnvObjsInterface::recomputeErosionValues);
    QObject::connect(instantiateButton, &QPushButton::pressed, this, &EnvObjsInterface::instantiateObject);

//    QObject::connect(instantiateReefButton, &QPushButton::pressed, this, [&](){ displayProbas("reef"); });

    displayCurrentsButton->setChecked(displayVelocities);
    displayHighCurrentsButton->setChecked(displayHighCurrents);
    displaySedimentsButton->setChecked(displaySediments);
    displayErosionsButton->setChecked(displayHighErosions);

    return layout;
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

    if (objectsListWidget != nullptr)
        updateObjectsList();
}

void EnvObjsInterface::setDefinitionFile(std::string filename)
{
    this->primitiveDefinitionFile = filename;
    EnvObject::readFile(filename);
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

GridF computeScoreMap(std::string objectName, const Vector3& dimensions, bool& possible) {
    auto obj = EnvObject::availableObjects[objectName];
    GridF score = GridF(dimensions);
    score.iterateParallel([&](const Vector3& pos) {
        score(pos) = obj->evaluate(pos);
    });
    if (score.min() < 1e5 || abs(score.min() - score.max()) > 1e-5) {
        possible = true;
        float secondMax = 0.f;
        score.iterate([&](size_t i) {
            if (score[i] < 1e5 && score[i] > secondMax)
                secondMax = score[i];
        });
        score.iterateParallel([&](size_t i) {
            score[i] = std::clamp(score[i], 0.f, secondMax);
        });

        score.normalizeUsing(NORMALIZE_METHOD::NORMALIZE_MINMAX);
        score = 1.f - score; // 1 = good, 0 = bad
    } else {
        possible = false;
    }
    return score;
}

Vector3 bestPositionForInstantiation(std::string objectName, const GridF& score) {
    if (std::abs(score.max() - score.min()) < 1e-5)
        return Vector3::random(Vector3(score.sizeX, score.sizeY, 0));
    float scoreToObtain = random_gen::generate(score.sum());
    float bestScore = 10000;
    Vector3 bestPos = Vector3::invalid();
    score.iterateRandomly([&](size_t i) {
        if (scoreToObtain - score[i] > 0) {
            scoreToObtain -= score[i];
        } else {
            bestPos = score.getCoordAsVector3(i);
        }
        /*
        if (score[i] < bestScore) {
            bestScore = score[i];
            bestPos = score.getCoordAsVector3(i);
        }*/
    });
    return bestPos;
}

EnvObject* instantiateObjectAtBestPosition(std::string objectName, const Vector3& position, const GridF& score) {
    EnvObject* newObject = EnvObject::instantiate(objectName);

    auto objAsPoint = dynamic_cast<EnvPoint*>(newObject);
    auto objAsCurve = dynamic_cast<EnvCurve*>(newObject);
    auto objAsArea = dynamic_cast<EnvArea*>(newObject);

    if (objAsPoint) {
        // Nothing to do...
//        objAsPoint->position = position;
    } else if (objAsCurve) {
        Vector3 direction = -score.gradient(position).normalize();
        if (direction.norm2() < 1e-5) direction = Vector3(1, 0);
        Vector3 bidirection = (isIn(newObject->material, LayerBasedGrid::invisibleLayers) ? direction : Vector3(direction.y, -direction.x));
        objAsCurve->curve = BSpline({
            bidirection * -objAsCurve->length* .5f,
            bidirection * objAsCurve->length * .5f,
        });
    } else if (objAsArea) {
        objAsArea->area = ShapeCurve::circle(objAsArea->width, Vector3(), 8);
    }

    newObject->translate(position.xy());
    return newObject;
}

void EnvObjsInterface::instantiateObject()
{
    std::cout << "Instantiate new object : " << showTime(timeIt([&]() {
        GridF heights = heightmap->getHeights();
        /*voxelGrid->computeFlowfield(LBM, 1, implicitTerrain.get());
        GridV3 surfaceNormals = heightmap->getNormals();
        GridV3 waterFlow = heightmap->properties->simulations[LBM]->getVelocities(heights.getDimensions());
        for (size_t i = 0; i < surfaceNormals.size(); i++)
            surfaceNormals[i].normalize();
        */
        std::map<std::string, GridF> scores;
        std::map<std::string, bool> possible;
        for (auto& [name, obj] : EnvObject::availableObjects) {
            auto& score = scores[name];
            score = computeScoreMap(name, heights.getDimensions(), possible[name]);
        }

        std::vector<std::string> possibleObjects;
        for (auto& [name, isPossible] : possible)
            if (isPossible)
                possibleObjects.push_back(name);

        if (possibleObjects.size() > 0) {
            std::shuffle(possibleObjects.begin(), possibleObjects.end(), random_gen::random_generator);
            std::string name = possibleObjects[0];
            auto& score = scores[name];

            Vector3 bestPos = bestPositionForInstantiation(name, score);
            EnvObject* newObject = instantiateObjectAtBestPosition(name, bestPos, score);
            auto implicit = newObject->createImplicitPatch();
            implicitTerrain->addChild(implicit);
            implicitTerrain->updateCache();
            implicitTerrain->update();
            voxelGrid->fromImplicit(implicitTerrain.get(), 40);
            heightmap->fromVoxelGrid(*voxelGrid.get());
            std::cout << "Instantiating " << name << " at position " << bestPos << std::endl;
            EnvObject::precomputeTerrainProperties(*heightmap);
        } else {
            std::cout << "No object to instantiate..." << std::endl;
        }
    })) << std::endl;
    updateObjectsList();
}

void EnvObjsInterface::instantiateSpecific(std::string objectName)
{
    std::cout << "Instantiate new " << objectName << " (forced) : " << showTime(timeIt([&]() {
        GridF heights = heightmap->getHeights();
        /*voxelGrid->computeFlowfield(LBM, 1, implicitTerrain.get());
        GridV3 surfaceNormals = heightmap->getNormals();
        GridV3 waterFlow = heightmap->properties->simulations[LBM]->getVelocities(heights.getDimensions());
        for (size_t i = 0; i < surfaceNormals.size(); i++)
            surfaceNormals[i].normalize();
        */

        bool possible;
        GridF score = computeScoreMap(objectName, heights.getDimensions(), possible);

        if (possible) {
            Vector3 bestPos = bestPositionForInstantiation(objectName, score);
            EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPos, score);
            auto implicit = newObject->createImplicitPatch();
            implicitTerrain->addChild(implicit);
            implicitTerrain->updateCache();
            implicitTerrain->update();
            voxelGrid->fromImplicit(implicitTerrain.get(), 40);
            heightmap->fromVoxelGrid(*voxelGrid.get());
            std::cout << "Instantiating " << objectName << " at position " << bestPos << std::endl;
            EnvObject::precomputeTerrainProperties(*heightmap);
        } else {
            std::cout << "Nope, impossible to instantiate..." << std::endl;
        }
    })) << std::endl;
    updateObjectsList();
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
//    Mesh::createVectorField(EnvObject::flowfield, voxelGrid->getDimensions(), &velocitiesMesh, -1, false, true);
    //    std::cout << EnvObject::sandDeposit << " -> " << EnvObject::sandDeposit.min() << " -- " << EnvObject::sandDeposit.max() << std::endl;
}

void EnvObjsInterface::displayProbas(std::string objectName)
{
    Vector3 dimensions = Vector3(heightmap->getSizeX(), heightmap->getSizeY(), 1);
    bool possible;
    GridF score = computeScoreMap(objectName, dimensions, possible);
    if (!possible) {
        Plotter::get()->addImage(score * 0.f, false);
    } else {
        Plotter::get()->addImage(score, false);
    }
    Plotter::get()->exec();
}

void EnvObjsInterface::displaSedimentsDistrib()
{
    GridF sediments = EnvObject::sandDeposit;
    Plotter::get()->addImage(sediments);
    Plotter::get()->exec();
}

void EnvObjsInterface::updateObjectsList()
{
    objectsListWidget->clear();
    auto list = EnvObject::instantiatedObjects;

    for (auto& obj : list) {
        objectsListWidget->addItem(new HierarchicalListWidgetItem(obj->name, obj->ID, 0));
    }
}

void EnvObjsInterface::updateObjectsListSelection(QListWidgetItem *newSelectionItem)
{
    auto newSelection = dynamic_cast<HierarchicalListWidgetItem*>(newSelectionItem);
    int objID = newSelection->ID;
    EnvObject* selection = nullptr;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (obj->ID == objID) {
            selection = obj;
            break;
        }
    }
    if (!selection)
        return;

    Vector3 selectionPos;
    if (auto asPoint = dynamic_cast<EnvPoint*>(selection)) {
        selectionPos = asPoint->position;
        std::cout << "Pos " << selection->name << ": " << asPoint->position << "\n-> selection at " << selectionPos.xy() << std::endl;
    } else if (auto asCurve = dynamic_cast<EnvCurve*>(selection)) {
        selectionPos = asCurve->curve.center();
        std::cout << "Curve " << selection->name << ":" << asCurve->curve.toString() << "\n-> selection at " << selectionPos.xy() << std::endl;
    } else if (auto asArea = dynamic_cast<EnvArea*>(selection)) {
        selectionPos = asArea->area.center();
        std::cout << "Area: " << selection->name << "" << asArea->area.toString() << "\n-> selection at " << selectionPos.xy() << std::endl;
    } else {
        std::cerr << "Object #" << selection->ID << " (" << selection->name << ") could not be casted to Point, Curve or Area..." << std::endl;
        return;
    }

    selectionPos.z = voxelGrid->getHeight(selectionPos.x, selectionPos.y);
    objectsMesh.fromArray(Mesh::getPointsForArrow(selectionPos + Vector3(0, 0, 20), selectionPos));
    Q_EMIT this->update();
}

void EnvObjsInterface::hotReloadFile()
{
    bool needReload = false;
    if (this->primitiveDefinitionFile != "") {
        QFileInfo infos(QString::fromStdString(this->primitiveDefinitionFile));
        auto lastModifTime = infos.lastModified();
        if (lastModifTime > this->lastTimeFileHasBeenModified || !this->lastTimeFileHasBeenModified.isValid()) {
            needReload = true;
            this->lastTimeFileHasBeenModified = lastModifTime;
        }
    }
    if (needReload) {
        EnvObject::readFile(this->primitiveDefinitionFile); //this->loadPatchesFromFile(this->primitiveDefinitionFile);
        std::cout << "Environmental objects file has changed" << std::endl;
    }
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
    std::cout << "Flow evaluation: " << showTime(timeIt([&]() {
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
    })) << std::endl;
    return {erosionsAmount, velocities};
}
