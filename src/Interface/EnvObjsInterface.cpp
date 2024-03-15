#include "EnvObjsInterface.h"

#include "Interface/InterfaceUtils.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/CommonInterface.h"
#include "TerrainModification/CoralIslandGenerator.h"
#include "DataStructure/Image.h"
#include "EnvObject/ExpressionParser.h"

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

    this->rootPatch = new Implicit2DNary;
    this->implicitTerrain->addChild(rootPatch);
//    createEnvObjectsFromImplicitTerrain();
//    recomputeErosionValues();
    EnvObject::precomputeTerrainProperties(*heightmap);


    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, this, [&](const Vector3& clickPos, Vector3 value) {
        if (!this->isVisible()) return;
        GridV3 dataV3 = Plotter::get()->displayedImage;
        GridF data(dataV3.getDimensions());
        dataV3.iterateParallel([&](size_t i) {
            data[i] = dataV3[i].x;
        });
        GridV3 gradients = data.gradient();
        ShapeCurve isoline = this->computeNewObjectsShapeAtPosition(clickPos, gradients, 50.f, 10.f);

        int nbSamples = 500;

        GridF result(data.getDimensions());

        auto path = isoline.getPath(nbSamples);
        result.iterateParallel([&](const Vector3& pos) {
            if (isoline.containsXY(pos)) {
                result(pos) = .5f;
            }
        });
        for (size_t i = 0; i < path.size(); i++) {
            result(path[i]) = 1.f;
        }
        Plotter::get()->addImage(result);
        Plotter::get()->show();
    });
}

void EnvObjsInterface::display(const Vector3 &camPos)
{
    if (!this->visible)
        return;

    if (this->waitAtEachFrame) {
        this->updateEnvironmentFromEnvObjects(true, false);
        this->displaySedimentsDistrib();
    }

//    if (displayVelocities) {
        velocitiesMesh.shader->setVector("color", std::vector<float>{.2f, .2f, .8f, .5f});
//        velocitiesMesh.reorderLines(camPos);
        velocitiesMesh.display(GL_LINES, 5.f);
//    }
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

    ButtonElement* instantiateButton = new ButtonElement("Instantiate", [&]() { this->instantiateObject(); });
    ButtonElement* recomputeErosionButton = new ButtonElement("Erosion values", [&]() { this->recomputeErosionValues(); });
    ButtonElement* spendTimeButton = new ButtonElement("Wait", [&]() { this->updateEnvironmentFromEnvObjects(true); });
    ButtonElement* showDepositionButton = new ButtonElement("Show deposition", [&]() { this->displaySedimentsDistrib(); });
    ButtonElement* showPolypButton = new ButtonElement("Show coral seeds", [&]() { this->displayPolypDistrib(); });
    ButtonElement* showFlowfieldButton = new ButtonElement("Show flow", [&]() { this->displayFlowfieldAsImage(); });
    CheckboxElement* waitAtEachFrameButton = new CheckboxElement("Auto wait", this->waitAtEachFrame);
    ButtonElement* createFromGAN = new ButtonElement("From GAN", [&]() { this->fromGanUI(); });
    TextEditElement* testingFormula = new TextEditElement("", "Try: ");
    testingFormula->setOnTextChange([&](std::string expression) { this->evaluateAndDisplayCustomCostFormula(expression); });
    ButtonElement* testPerformancesButton = new ButtonElement("Run test", [&]() { this->runPerformanceTest(); });

    objectsListWidget = new HierarchicalListWidget;
    updateObjectsList();
//    QObject::connect(objectsListWidget, &HierarchicalListWidget::clicked, this, [&](QModelIndex item) { qDebug() << item; }); //&EnvObjsInterface::updateObjectsListSelection);
    QObject::connect(objectsListWidget, &HierarchicalListWidget::currentItemChanged, this, [&](QListWidgetItem* current, QListWidgetItem* previous) { this->updateObjectsListSelection(current); });


    std::vector<QWidget*> probaButtons;
    for (auto& [name, obj] : EnvObject::availableObjects) {
        ButtonElement* showButton = new ButtonElement("Show " + obj->name, [&](){ this->displayProbas(name); });
        ButtonElement* forceButton = new ButtonElement("Force", [&](){ this->instantiateSpecific(name); });
        probaButtons.push_back(showButton->get());
        probaButtons.push_back(forceButton->get());
    }

    layout->addWidget(spendTimeButton->get());
    layout->addWidget(waitAtEachFrameButton->get());
    layout->addWidget(showDepositionButton->get());
    layout->addWidget(showPolypButton->get());
    layout->addWidget(showFlowfieldButton->get());
//    layout->addWidget(createFromGAN->get());
//    layout->addWidget(instantiateButton->get());
    layout->addWidget(createMultiColumnGroup(probaButtons, 2));
    layout->addWidget(objectsListWidget);
//    layout->addWidget(recomputeErosionButton->get());
    layout->addWidget(testingFormula->get());
//    layout->addWidget(testPerformancesButton->get());


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
    updateObjectsListSelection(nullptr); // Hide single object's data
    ActionInterface::hide();
}

void EnvObjsInterface::afterTerrainUpdated()
{

}

void EnvObjsInterface::afterWaterLevelChanged()
{
    EnvObject::recomputeFlowAndSandProperties(*heightmap);
}

void EnvObjsInterface::mouseClickedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    if (!this->isVisible()) return;
    if (!mouseInMap) return;
    if (!event->modifiers().testFlag(Qt::KeyboardModifier::ControlModifier)) return;

    Vector3 newPos = mouseWorldPosition.xy();

    std::string name = "motu";
    EnvObject* newObject = this->instantiateObjectAtBestPosition(name, newPos, GridF());
    auto implicit = newObject->createImplicitPatch();
    this->implicitPatchesFromObjects[newObject] = implicit;
    if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
        this->implicitTerrain->addChild(this->rootPatch);
    rootPatch->addChild(implicit);
    rootPatch->reevaluateAll();
    implicitTerrain->updateCache();
    implicitTerrain->update();
    voxelGrid->fromImplicit(implicitTerrain.get(), 40);
    heightmap->fromVoxelGrid(*voxelGrid.get());
    std::cout << "Instantiating " << name << " at position " << newPos << std::endl;
    EnvObject::recomputeTerrainPropertiesForObject(*heightmap, name);
    EnvObject::recomputeFlowAndSandProperties(*heightmap);
    updateObjectsList();
}

GridF computeScoreMap(std::string objectName, const Vector3& dimensions, bool& possible, bool applyNormalization = true) {
    auto obj = EnvObject::availableObjects[objectName];
    GridF score = GridF(dimensions);
    score.iterateParallel([&](const Vector3& pos) {
        score(pos) = std::max(obj->evaluate(pos), 0.f);
//        score(pos) = obj->evaluate(pos);
    });
    if (abs(score.min() - score.max()) > 1e-5) {
        possible = true;
        /*float secondMax = 0.f;
        score.iterate([&](size_t i) {
            if (score[i] < 1e5 && score[i] > secondMax)
                secondMax = score[i];
        });
        score.iterateParallel([&](size_t i) {
            score[i] = std::max(score[i], 0.f);
//            score[i] = std::min(score[i], secondMax);
//            score[i] = std::clamp(score[i], 0.f, secondMax);
        });*/
        if (applyNormalization)
            score.normalizeUsing(NORMALIZE_METHOD::NORMALIZE_MINMAX);
//        score = 1.f - score; // 1 = good, 0 = bad
    } else {
        possible = false;
    }
    return score;
}

Vector3 bestPositionForInstantiation(std::string objectName, const GridF& score) {
    GridV3 gradients = score.gradient();
    if (std::abs(score.max() - score.min()) < 1e-5)
        return Vector3::random(Vector3(score.sizeX, score.sizeY, 0));
    float scoreToObtain = random_gen::generate(score.sum());
    float bestScore = 10000;
    Vector3 bestPos = Vector3::invalid();
//    std::vector<Vector3> possiblePositions(300);
//    for (auto& p : possiblePositions) {
    while (scoreToObtain > 0) {
        Vector3 p = Vector3::random(Vector3(), score.getDimensions().xy());
        for (int iter = 0; iter < 100; iter++) {
            p += gradients(p).normalize();
        }
        scoreToObtain -= score(p);
        bestPos = p;
    }
    /*score.iterateRandomly([&](size_t i) {
        if (scoreToObtain - score[i] > 0) {
            scoreToObtain -= score[i];
        } else {
            bestPos = score.getCoordAsVector3(i);
        }
    });*/
    return bestPos;
}

EnvObject* EnvObjsInterface::instantiateObjectAtBestPosition(std::string objectName, const Vector3& position, const GridF& score) {
    EnvObject* newObject = EnvObject::instantiate(objectName);

    auto objAsPoint = dynamic_cast<EnvPoint*>(newObject);
    auto objAsCurve = dynamic_cast<EnvCurve*>(newObject);
    auto objAsArea = dynamic_cast<EnvArea*>(newObject);

    if (objAsPoint) {
        // Nothing to do...
//        objAsPoint->position = position;
    } else if (objAsCurve) {

        GridV3 gradients = score.gradient();
        BSpline curve = this->computeNewObjectsCurveAtPosition(position, gradients, 200.f, 0.f);
        curve.translate(-position);
        objAsCurve->curve = curve.resamplePoints(5);

        /*
        Vector3 direction = -score.gradient(position).normalize();
        if (direction.norm2() < 1e-5) direction = Vector3(1, 0);
        Vector3 bidirection = Vector3(direction.y, -direction.x); // (isIn(newObject->material, LayerBasedGrid::invisibleLayers) ? direction : Vector3(direction.y, -direction.x));
        objAsCurve->curve = BSpline({
            bidirection * -objAsCurve->length* .5f,
            bidirection * objAsCurve->length * .5f,
        });*/
    } else if (objAsArea) {
        objAsArea->area = ShapeCurve::circle(objAsArea->width, Vector3(), 8);
    }

    newObject->translate(position.xy());
    return newObject;
}

void EnvObjsInterface::instantiateObject()
{
    displayProcessTime("Instantiate new object... ", [&]() {
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
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);
            rootPatch->addChild(implicit);
            rootPatch->reevaluateAll();
            implicitTerrain->updateCache();
            implicitTerrain->update();
            voxelGrid->fromImplicit(implicitTerrain.get(), 40);
            heightmap->fromVoxelGrid(*voxelGrid.get());
            std::cout << "Instantiating " << name << " at position " << bestPos << std::endl;
//            EnvObject::precomputeTerrainProperties(*heightmap);
            EnvObject::recomputeTerrainPropertiesForObject(*heightmap, name);
            EnvObject::recomputeFlowAndSandProperties(*heightmap);
        } else {
            std::cout << "No object to instantiate..." << std::endl;
        }
    });
    updateObjectsList();
}

void EnvObjsInterface::instantiateSpecific(std::string objectName)
{
    displayProcessTime("Instantiate new " + objectName + " (forced)... ", [&]() {
        GridF heights = heightmap->getHeights();

        bool possible;
        GridF score = computeScoreMap(objectName, heights.getDimensions(), possible);

        if (possible) {
            Vector3 bestPos = bestPositionForInstantiation(objectName, score);
            EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPos, score);
            auto implicit = newObject->createImplicitPatch();
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);
            rootPatch->addChild(implicit);
            EnvObject::recomputeTerrainPropertiesForObject(*heightmap, objectName);
            /*rootPatch->reevaluateAll();
            implicitTerrain->updateCache();
            implicitTerrain->update();
            voxelGrid->fromImplicit(implicitTerrain.get(), 20);
            heightmap->fromVoxelGrid(*voxelGrid.get());
            std::cout << "Instantiating " << objectName << " at position " << bestPos << std::endl;
//            EnvObject::precomputeTerrainProperties(*heightmap);
            EnvObject::recomputeFlowAndSandProperties(*heightmap);*/
            this->updateEnvironmentFromEnvObjects(true);
        } else {
            std::cout << "Nope, impossible to instantiate..." << std::endl;
        }
    });
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

//    this->updateEnvironmentFromEnvObjects();
}

void EnvObjsInterface::updateEnvironmentFromEnvObjects(bool updateImplicitTerrain, bool emitUpdateSignal)
{
    bool verbose = false;
    // Get original flowfield, do not accumulate effects (for now).
    displayProcessTime("Get velocity... ", [&]() {
        EnvObject::flowfield = dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->getVelocities(EnvObject::flowfield.sizeX, EnvObject::flowfield.sizeY, EnvObject::flowfield.sizeZ);
    }, verbose);
    displayProcessTime("Apply effects... ", [&]() {
        EnvObject::applyEffects();
    }, verbose);
    displayProcessTime("Get impacted... ", [&]() {
        EnvObject::beImpactedByEvents();
    }, verbose);
    displayProcessTime("Recompute properties... ", [&]() {
        EnvObject::recomputeFlowAndSandProperties(*heightmap);
    }, verbose);
    if (updateImplicitTerrain) {
        for (auto& [obj, implicit] : this->implicitPatchesFromObjects) {
            obj->createImplicitPatch(dynamic_cast<ImplicitPrimitive*>(implicit));
        }
        rootPatch->reevaluateAll();
        implicitTerrain->updateCache();
        implicitTerrain->update();
        displayProcessTime("Update voxels from implicit terrain... ", [&]() {
            voxelGrid->fromImplicit(implicitTerrain.get(), 20);
        }, verbose);
        displayProcessTime("Update heightmap from implicit terrain... ", [&]() {
            heightmap->fromVoxelGrid(*voxelGrid.get());
        }, verbose);
        if (emitUpdateSignal) {
            Q_EMIT this->updated();
        }
    }
}

void EnvObjsInterface::onlyUpdateFlowAndSandFromEnvObjects()
{
    EnvObject::applyEffects();
}

void EnvObjsInterface::displayProbas(std::string objectName)
{
    Vector3 dimensions = Vector3(heightmap->getSizeX(), heightmap->getSizeY(), 1);
    bool possible;
    GridF score = computeScoreMap(objectName, dimensions, possible, false);
    if (!possible) {
        Plotter::get()->addImage(score * 0.f);
    } else {
        float smallestPositive = score.max();
        score.iterate([&](size_t i) {
            if (score[i] > 0.f && score[i] < smallestPositive)
                smallestPositive = score[i];
        });
        score.iterateParallel([&](size_t i) {
            score[i] = std::max(score[i], smallestPositive);
        });
        Plotter::get()->addImage(score);
    }
    Plotter::get()->show();
}

void EnvObjsInterface::displaySedimentsDistrib()
{
    GridF sediments = EnvObject::sandDeposit;
    Plotter::get()->addImage(sediments);
    Plotter::get()->show();
}

void EnvObjsInterface::displayPolypDistrib()
{
    GridF polyp = EnvObject::polypDeposit;
    Plotter::get()->addImage(polyp);
    Plotter::get()->show();
}

void EnvObjsInterface::displayFlowfieldAsImage()
{
    GridV3 flow = EnvObject::flowfield;
    Plotter::get()->addImage(flow);
    Plotter::get()->show();
}

void EnvObjsInterface::updateObjectsList()
{
    if (!objectsListWidget) return;
    objectsListWidget->clear();
    auto list = EnvObject::instantiatedObjects;

    for (auto& obj : list) {
        objectsListWidget->addItem(new HierarchicalListWidgetItem(obj->name, obj->ID, 0));
    }
}

void EnvObjsInterface::updateObjectsListSelection(QListWidgetItem *newSelectionItem)
{
    this->velocitiesMesh.fromArray(std::vector<float>{});
    this->objectsMesh.fromArray(std::vector<float>{});

    auto newSelection = dynamic_cast<HierarchicalListWidgetItem*>(newSelectionItem);
    if (!newSelection) {
        currentSelection = nullptr;
        return;
    }
    int objID = newSelection->ID;
    EnvObject* selection = nullptr;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (obj->ID == objID) {
            selection = obj;
            break;
        }
    }
    if (selection == this->currentSelection) {
        currentSelection = nullptr;
        return;
    }
    currentSelection = selection;
    if (!selection)
        return;

    Vector3 selectionPos;
    if (auto asPoint = dynamic_cast<EnvPoint*>(selection)) {
        selectionPos = asPoint->position;
        std::cout << "Pos " << selection->name << ": " << asPoint->position << "\n-> selection at " << selectionPos.xy() << std::endl;
        selectionPos.z = voxelGrid->getHeight(selectionPos.x, selectionPos.y) + 5.f;
        objectsMesh.fromArray(Mesh::getPointsForArrow(selectionPos + Vector3(0, 0, 20), selectionPos));
    } else if (auto asCurve = dynamic_cast<EnvCurve*>(selection)) {
        selectionPos = asCurve->curve.center();
        std::cout << "Curve " << selection->name << ":" << asCurve->curve.toString() << "\n-> selection at " << selectionPos.xy() << std::endl;
        std::vector<Vector3> meshPoints;
        for (size_t i = 0; i < asCurve->curve.size() - 1; i++) {
            auto p1 = asCurve->curve[i];
            auto p2 = asCurve->curve[i + 1];
            meshPoints.push_back(p1 + Vector3(0, 0, voxelGrid->getHeight(p1.x, p1.y) + 5.f));
            meshPoints.push_back(p2 + Vector3(0, 0, voxelGrid->getHeight(p2.x, p2.y) + 5.f));
        }
        objectsMesh.fromArray(meshPoints);
    } else if (auto asArea = dynamic_cast<EnvArea*>(selection)) {
        selectionPos = asArea->area.center();
        std::cout << "Area: " << selection->name << "" << asArea->area.toString() << "\n-> selection at " << selectionPos.xy() << std::endl;
        std::vector<Vector3> meshPoints;
        for (size_t i = 0; i < asArea->area.size() - 1; i++) {
            auto p1 = asArea->area[i];
            auto p2 = asArea->area[i + 1];
            meshPoints.push_back(p1 + Vector3(0, 0, voxelGrid->getHeight(p1.x, p1.y) + 5.f));
            meshPoints.push_back(p2 + Vector3(0, 0, voxelGrid->getHeight(p2.x, p2.y) + 5.f));
        }
        objectsMesh.fromArray(meshPoints);
    } else {
        std::cerr << "Object #" << selection->ID << " (" << selection->name << ") could not be casted to Point, Curve or Area..." << std::endl;
        return;
    }


    GridV3 initialFlow = EnvObject::initialFlowfield;
    GridV3 flow;
    GridF occupancy;
    std::tie(flow, occupancy) = selection->computeFlowModification();
    float currentGrowth = selection->computeGrowingState();
    occupancy *= currentGrowth;

    flow = flow * occupancy + initialFlow * (1.f - occupancy);
    initialFlow = initialFlow * (1.f - EnvObject::flowImpactFactor) + flow * EnvObject::flowImpactFactor;
//    velocitiesMesh.fromVectorField(initialFlow.resize(30, 30, 1), voxelGrid->getDimensions());
    Mesh::createVectorField(initialFlow.resize(30, 30, 1), voxelGrid->getDimensions(), &velocitiesMesh, 1.f, false, true);
    Q_EMIT this->updated();
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
        displayProcessTime("Environmental objects file has changed. Parsing... ", [&]() {
            EnvObject::readFile(this->primitiveDefinitionFile);
        });
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomCostFormula(std::string formula) const
{
    EnvPoint fake;
    fake.s_FittingFunction = formula;
    try {
        fake.fittingFunction = EnvObject::parseFittingFunction(formula, "");

        GridF eval(EnvObject::flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3& p) {
            eval(p) = fake.evaluate(p);
        });
        Plotter::get()->addImage(eval);
        Plotter::get()->show();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}

BSpline EnvObjsInterface::computeNewObjectsCurveAtPosition(const Vector3 &seedPosition, const GridV3 &gradients, float directionLength, float widthMaxLength)
{
    Vector3 pos = seedPosition;

    /*auto followIsovalue = [gradients](const Vector3& startPoint, float maxDist) -> BSpline {
        BSpline finalPath;

        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir(1, 0, 0);
        bool didAFullCircle = false;
        float totalDistance = 0.f;
        while (maxDist > totalDistance) {
            if (path.size() > 5 && (pos - startPoint).norm2() < 3*3){
                didAFullCircle = true;
                break; // Got back close to beginning
            }
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();

            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir = newDir * (dir.dot(newDir) < 0 ? -1.f : 1.f);

            pos += dir;

            totalDistance += dir.norm();

            path.points.push_back(pos);
        }

        finalPath = path;
        if (!didAFullCircle) {
            totalDistance = 0.f;
            pos = startPoint;
            path = BSpline();
            dir = Vector3(-1, 0, 0);
            while (maxDist > totalDistance) {
                if (path.size() > 5 && (pos - startPoint).norm2() < 3*3){
                    didAFullCircle = true;
                    break; // Got back close to beginning
                }
                Vector3 gradient = gradients(pos);
                if (gradient == Vector3()) break; // Nowhere to go
                gradient.normalize();

                Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
                dir = newDir * (dir.dot(newDir) < 0 ? -1.f : 1.f);

                pos += dir;
                totalDistance += dir.norm();

                path.points.push_back(pos);
            }
            if (didAFullCircle) {
                finalPath = path; // No need to take the first half, we just did a full circle
            } else {
                std::reverse(path.points.begin(), path.points.end());
                finalPath.points.insert(finalPath.points.begin(), path.points.begin(), path.points.end());
            }
        }
        return finalPath;
    };*/

    auto followGradient = [gradients](const Vector3& startPoint, float maxDist) -> BSpline {
        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir;
        float totalDistance = 0.f;
        while (totalDistance < maxDist) {
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();
            dir = gradient;

            pos += dir;

            totalDistance += dir.norm();

            path.points.push_back(pos);

            int nbPoints = std::min(int(path.size()), 5);
            if (nbPoints > 2) {
                std::vector<Vector3> lastPositions(path.points.end() - nbPoints, path.points.end());
                Vector3 meanVel;
                for (size_t i = 0; i < nbPoints - 1; i++) {
                    meanVel += (lastPositions[i + 1] - lastPositions[i]);
                }
                if ((meanVel / float(nbPoints - 1)).norm2() < .25f) {
//                        std::cout << "Stuck in grad, stopping" << std::endl;
                    path.points.erase(path.points.end() - nbPoints, path.points.end());
                    break;
                }
            }
        }
        return path;
    };
    auto followInvGradient = [gradients](const Vector3& startPoint, float maxDistance) -> BSpline {
        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir;
        float totalDistance = 0.f;
        while (totalDistance < maxDistance) {
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();
            dir = -gradient;

            pos += dir;
            totalDistance += dir.norm();

            path.points.push_back(pos);

            int nbPoints = std::min(int(path.size()), 5);

            if (nbPoints > 2) {
                std::vector<Vector3> lastPositions(path.points.end() - nbPoints, path.points.end());
                Vector3 meanVel;
                for (size_t i = 0; i < nbPoints - 1; i++) {
                    meanVel += (lastPositions[i + 1] - lastPositions[i]);
                }
                if ((meanVel / float(nbPoints - 1)).norm2() < .25f) {
//                        std::cout << "Stuck in grad, stopping" << std::endl;
                    path.points.erase(path.points.end() - nbPoints, path.points.end());
                    break;
                }
            }
        }
        return path;
    };

    BSpline gradSpline = followGradient(pos, directionLength * .5f);
    BSpline invGradSpline = followInvGradient(pos, directionLength * .5f);
    invGradSpline.reverseVertices();

    BSpline isoline = gradSpline;
    for (auto& p : invGradSpline.points) {
        isoline.points.push_back(p);
    }

    return isoline;
}

ShapeCurve EnvObjsInterface::computeNewObjectsShapeAtPosition(const Vector3 &seedPosition, const GridV3& gradients, float directionLength, float widthMaxLength)
{
    Vector3 pos = seedPosition;
    Vector3 dir(1, 0, 0);
    BSpline path({pos});
    GridV3 result(gradients.getDimensions());

    auto followIsovalue = [gradients](const Vector3& startPoint, float maxDist) -> BSpline {
        BSpline finalPath;

        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir(1, 0, 0);
        bool didAFullCircle = false;
        float totalDistance = 0.f;
        while (maxDist > totalDistance) {
            if (path.size() > 5 && (pos - startPoint).norm2() < 3*3){
                didAFullCircle = true;
                break; // Got back close to beginning
            }
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();

            Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
            dir = newDir * (dir.dot(newDir) < 0 ? -1.f : 1.f);

            pos += dir;

            totalDistance += dir.norm();

            path.points.push_back(pos);
        }

        finalPath = path;
        if (!didAFullCircle) {
            totalDistance = 0.f;
            pos = startPoint;
            path = BSpline();
            dir = Vector3(-1, 0, 0);
            while (maxDist > totalDistance) {
                if (path.size() > 5 && (pos - startPoint).norm2() < 3*3){
                    didAFullCircle = true;
                    break; // Got back close to beginning
                }
                Vector3 gradient = gradients(pos);
                if (gradient == Vector3()) break; // Nowhere to go
                gradient.normalize();

                Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
                dir = newDir * (dir.dot(newDir) < 0 ? -1.f : 1.f);

                pos += dir;
                totalDistance += dir.norm();

                path.points.push_back(pos);
            }
            if (didAFullCircle) {
                finalPath = path; // No need to take the first half, we just did a full circle
            } else {
                std::reverse(path.points.begin(), path.points.end());
                finalPath.points.insert(finalPath.points.begin(), path.points.begin(), path.points.end());
            }
        }
        return finalPath;
    };

    /*auto followGradient = [gradients](const Vector3& startPoint, float maxDist) -> BSpline {
        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir;
        float totalDistance = 0.f;
        while (totalDistance < maxDist) {
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();
            dir = gradient;

            pos += dir;

            totalDistance += dir.norm();

            path.points.push_back(pos);

            int nbPoints = std::min(int(path.size()), 5);
            if (nbPoints > 2) {
                std::vector<Vector3> lastPositions(path.points.end() - nbPoints, path.points.end());
                Vector3 meanVel;
                for (size_t i = 0; i < nbPoints - 1; i++) {
                    meanVel += (lastPositions[i + 1] - lastPositions[i]);
                }
                if ((meanVel / float(nbPoints - 1)).norm2() < .25f) {
//                        std::cout << "Stuck in grad, stopping" << std::endl;
                    path.points.erase(path.points.end() - nbPoints, path.points.end());
                    break;
                }
            }
        }
        return path;
    };
    auto followInvGradient = [gradients](const Vector3& startPoint, float maxDistance) -> BSpline {
        Vector3 pos = startPoint;
        BSpline path({pos});
        Vector3 dir;
        float totalDistance = 0.f;
        while (totalDistance < maxDistance) {
            Vector3 gradient = gradients(pos);
            if (gradient == Vector3()) break; // Nowhere to go
            gradient.normalize();
            dir = -gradient;

            pos += dir;
            totalDistance += dir.norm();

            path.points.push_back(pos);

            int nbPoints = std::min(int(path.size()), 5);

            if (nbPoints > 2) {
                std::vector<Vector3> lastPositions(path.points.end() - nbPoints, path.points.end());
                Vector3 meanVel;
                for (size_t i = 0; i < nbPoints - 1; i++) {
                    meanVel += (lastPositions[i + 1] - lastPositions[i]);
                }
                if ((meanVel / float(nbPoints - 1)).norm2() < .25f) {
//                        std::cout << "Stuck in grad, stopping" << std::endl;
                    path.points.erase(path.points.end() - nbPoints, path.points.end());
                    break;
                }
            }
        }
        return path;
    };*/

    BSpline isoline = followIsovalue(seedPosition, directionLength);
    return isoline;
}

void EnvObjsInterface::runPerformanceTest()
{
    std::vector<float> timings;
    for (int i = 0; i < 201; i++) {
        float time = displayProcessTime("Evaluation with " + std::to_string(EnvObject::instantiatedObjects.size()) + " objects (x5) : ", [&]() {
            for (int iter = 0; iter < 5; iter++)
                this->updateEnvironmentFromEnvObjects(false, false);
        });

        Vector3 newPos = Vector3::random() * heightmap->getDimensions().xy();

        std::string name = "motu";
        EnvObject* newObject = this->instantiateObjectAtBestPosition(name, newPos, GridF());

        timings.push_back((time * 1e-6) / 5.f);
    }
    Plotter::get()->reset();
    Plotter::get()->addPlot(timings, "Timing", Qt::darkBlue);
    Plotter::get()->show();
}

void EnvObjsInterface::fromGanUI()
{
    EnvObject::reset();
    EnvObject::readFile(this->primitiveDefinitionFile);
/*
    auto gauss = dynamic_cast<EnvPoint*>(EnvObject::instantiate("motu"));
    gauss->position = Vector3(100, 50, 0);
    gauss->radius = 35;
    auto patch = gauss->createImplicitPatch();
    this->implicitPatchesFromObjects[gauss] = patch;

//    auto curve = dynamic_cast<EnvCurve*>(EnvObject::instantiate("reef"));
//    curve->curve = BSpline({
//                               Vector3(120, 0, 0),
//                               Vector3(90, 50, 0),
//                               Vector3(120, 100, 0),
//                           });
//    auto patch = curve->createImplicitPatch();
//    this->implicitPatchesFromObjects[curve] = patch;
    rootPatch->addChild(patch);
    implicitTerrain->updateCache();
    implicitTerrain->update();
    rootPatch->reevaluateAll();
    displayProcessTime("Env objects to voxels... ", [&]() {
//            voxelGrid->from2DGrid(*heightmap);
        voxelGrid->fromImplicit(rootPatch, 40);
    });
    displayProcessTime("Env objects to heightmap... ", [&]() {
        heightmap->fromVoxelGrid(*voxelGrid.get());
    });
//        implicitTerrain->addChild(obj->createImplicitPatch());
//        implicitTerrain->_cached = false;
//        voxelGrid->fromImplicit(implicitTerrain.get());

    EnvObject::precomputeTerrainProperties(*heightmap);
    this->updateEnvironmentFromEnvObjects();


    return;
    */

    std::string path = "Python_tests/test_island_heightmapfeatures/";
    QString q_filename= QString::fromStdString(path + "1.png");  //QFileDialog::getOpenFileName(this, "Open feature map", QString::fromStdString(path), "*", nullptr);
    if (!q_filename.isEmpty()) {
        std::string file = q_filename.toStdString();
        GridV3 img = Image::readFromFile(file).colorImage;

        auto envObjects = CoralIslandGenerator::envObjsFromFeatureMap(img, voxelGrid->getDimensions());
        rootPatch->deleteAllChildren();
        for (auto& newObject : envObjects) {
            auto implicit = newObject->createImplicitPatch();
            this->implicitPatchesFromObjects[newObject] = implicit;
            rootPatch->addChild(implicit);
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);
        }
        implicitTerrain->updateCache();
        implicitTerrain->update();
        rootPatch->reevaluateAll();
        std::cout << "To voxels: " << showTime(timeIt([&]() {
//            voxelGrid->from2DGrid(*heightmap);
            voxelGrid->fromImplicit(implicitTerrain.get(), 40);
        })) << std::endl;
        std::cout << "To heightmap: " << showTime(timeIt([&]() {
            heightmap->fromVoxelGrid(*voxelGrid.get());
        })) << std::endl;
//        implicitTerrain->addChild(obj->createImplicitPatch());
//        implicitTerrain->_cached = false;
//        voxelGrid->fromImplicit(implicitTerrain.get());

        EnvObject::precomputeTerrainProperties(*heightmap);
        this->updateEnvironmentFromEnvObjects();
        std::cout << "Update: " << showTime(timeIt([&]() {
            updateObjectsList();
        })) << std::endl;
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

    bool wrapParticles = false;

    float initialCapacity = .0f;

    FluidSimType selectedSimulationType = WARP;
    dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[selectedSimulationType])->recomputeVelocities();

    std::vector<std::pair<Vector3, Vector3>> initialPositionsAndDirections(erosionQtt);
    for (auto& [pos, dir] : initialPositionsAndDirections) {
        pos = Vector3::random(Vector3(), terrainDims);
    }

    EROSION_APPLIED applyOn = EROSION_APPLIED::DENSITY_VOXELS;
    FLOWFIELD_TYPE flowfieldUsed = FLOWFIELD_TYPE::BASIC; // FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS;
    DENSITY_TYPE densityUsed = DENSITY_TYPE::NATIVE;

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
    displayProcessTime("Computing erosion data... ", [&]() {
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
    });

    GridV3 velocities(terrainDims);
    GridF evaluationAmounts(terrainDims);
    displayProcessTime("Computing flow from erosion data... ", [&]() {
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
    });
    return {erosionsAmount, velocities};
}
