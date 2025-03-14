#include "EnvObjsInterface.h"

#include "Interface/InterfaceUtils.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/CommonInterface.h"
#include "TerrainModification/CoralIslandGenerator.h"
#include "DataStructure/Image.h"
#include "EnvObject/ExpressionParser.h"
#include "Utils/Voronoi.h"
#include "Interface/MeshInstanceAmplificationInterface.h"
#include "Interface/TerrainGenerationInterface.h"
#include "Utils/Delaunay.h"

EnvObjsInterface::EnvObjsInterface(QWidget *parent)
    : ActionInterface("envobjects", "Environmental Objects", "model", "Management of environmental objects generation", "envobjs_button.png", parent)
{
    primitiveDefinitionFile.onChange([&](std::string newDefinitions) { updateObjectsDefinitions(newDefinitions); });
    materialsDefinitionFile.onChange([&](std::string newDefinitions) { updateMaterialsDefinitions(newDefinitions); });
    transformationsFile.onChange([&](std::string newDefinitions) { updateMaterialsTransformationsDefinitions(newDefinitions); });
    scenarioFile.onChange([&](std::string newDefinitions) { updateScenarioDefinition(newDefinitions); });

    QTimer* hotreloadTimer = new QTimer(this);
    hotreloadTimer->setInterval(500);
    QObject::connect(hotreloadTimer, &QTimer::timeout, this, [&]() {
        materialsDefinitionFile.check();
        primitiveDefinitionFile.check();
        transformationsFile.check();
        scenarioFile.check();
    });
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
    this->selectedObjectsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->selectedObjectsMesh.useIndices = false;
    this->newObjectMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->newObjectMesh.useIndices = false;

    this->rootPatch = new ImplicitNaryOperator;
    // this->rootPatch = new Implicit2DNary;
    this->implicitTerrain->addChild(rootPatch);


    this->initialHeightmap = heightmap->heights;
    this->subsidedHeightmap = initialHeightmap;
    this->focusedArea = GridF(initialHeightmap.getDimensions(), 1.f);
    this->userFlowField = GridV3(initialHeightmap.getDimensions());
    this->simulationFlowField = GridV3(initialHeightmap.getDimensions());
    EnvObject::precomputeTerrainProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());


    // QObject::connect(Plotter::get(), &Plotter::clickedOnImage, this, [&](const Vector3& clickPos, Vector3 value) {
    QObject::connect(Plotter::get("Object Preview"), &Plotter::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* _event) {
        // if (!this->isVisible()) return;
        // if (!previewingObjectInPlotter) return;

        this->previewCurrentEnvObjectPlacement(clickPos);
    });

    QObject::connect(Plotter::get("Focus"), &Plotter::movedOnImage, this, [&](const Vector3& mousePos, const Vector3& prevPos, QMouseEvent* event) {
        // if (!this->isVisible()) return;
        // if (!this->focusAreaEditing) return;

        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;

        this->previewFocusAreaEdition(mousePos, leftPressed);
    });

    QObject::connect(Plotter::get("Flowfield"), &Plotter::movedOnImage, this, [&](const Vector3& mousePos, const Vector3& prevPos, QMouseEvent* event) {
        // if (!this->isVisible()) return;
        // if (!this->flowfieldEditing) return;

        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;

        Vector3 brushDir = (mousePos - prevPos).normalize() * .2f;
        // std::cout << prevPos << " " << brushDir << std::endl;
        this->previewFlowEdition(mousePos, brushDir);

        Q_EMIT this->updated();
    });


    /*QObject::connect(Plotter::get(), &Plotter::movedOnImage, this, [&](const Vector3& mousePos, const Vector3& prevPos, QMouseEvent* event) {

        ShapeCurve initial({
            Vector3(30, 30, 0),
            Vector3(70, 30, 0),
            Vector3(70, 70, 0),
            Vector3(50, 70),
            Vector3(30, 70, 0)
        });
        Voronoi voro(initial.points, Vector3(-100, -100, 0), Vector3(200, 200, 0));
        random_gen::random_generator.seed(0);
        std::vector<Vector3> points = initial.randomPointsInside(20);
        random_gen::random_generator.seed();
        Delaunay delaunay;
        delaunay.fromVoronoi(voro);
        auto triangles = delaunay.getTriangles();
        ShapeCurve shape(flattenArray(triangles));
        // points.insert(points.begin(), shape.points.begin(), shape.points.end());

        GridV3 img(100, 100, 1);
        for (int i = 0; i < initial.size(); i++) {
            for (auto& p : BSpline({initial[i], initial[i + 1]}).getPath(100))
                img(p) = Vector3(.5f, .5f, .5f);
        }

        std::vector<std::vector<float>> initialPos;
        for (int i = 0; i < points.size(); i++) {
            auto& p = points[i];
            auto coord = computeGreenCoordinates(p, shape);
            initialPos.push_back(coord);
        }

        for (int i = 0; i < initial.size(); i++) {
            for (auto& p : BSpline({initial[i], initial[i + 1]}).getPath(100))
                img(p) = Vector3(.5f, .5f, .5f);
        }

        delaunay.graph.nodes[3]->pos = mousePos;

        // for (auto& n : delaunay.graph.nodes)
            // if ((n->pos - mousePos).norm() < 10)
                // n->pos += mousePos - prevPos;
        triangles = delaunay.getTriangles();
        shape = ShapeCurve(flattenArray(triangles));

        for (int i = 0; i < points.size(); i++) {
            auto& coord = initialPos[i];
            auto q = computePointFromGreenCoordinates(coord, shape);
            for (auto& pp : BSpline({points[i], q}).getPath(100))
                img(pp) = Vector3(1, 1, 1);
        }

        Plotter::get()->addImage(img);
        Plotter::get()->show();
    });*/
}

void EnvObjsInterface::display(const Vector3 &camPos)
{
    if (!this->visible)
        return;

    if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
        this->implicitTerrain->addChild(this->rootPatch);

    if (this->waitAtEachFrame)
        this->updateEnvironmentFromEnvObjects(true, false);

    if (displayFlow) {
        velocitiesMesh.shader->setVector("color", std::vector<float>{.2f, .2f, .8f, .5f});
        velocitiesMesh.display(GL_LINES, 3.f);
    }

    if (displayHighErosions) {
        highErosionsMesh.shader->setVector("color", std::vector<float>{.8f, .2f, .2f, .5f});
        highErosionsMesh.reorderTriangles(camPos);
        highErosionsMesh.display();

        highDepositionMesh.shader->setVector("color", std::vector<float>{.2f, .8f, .2f, .5f});
        highDepositionMesh.reorderTriangles(camPos);
        highDepositionMesh.display();
    }
    selectedObjectsMesh.shader->setVector("color", std::vector<float>{.8f, .2f, .8f, .5f});
    selectedObjectsMesh.display(GL_LINES, 5.f);

    newObjectMesh.shader->setVector("color", std::vector<float>{.8f, .4f, .4f, .5f});
    newObjectMesh.display(GL_LINES, 5.f);
}

void EnvObjsInterface::replay(nlohmann::json action)
{

}

QLayout *EnvObjsInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout;

    InterfaceUI* ui = new InterfaceUI(layout);

    // ButtonElement* instantiateButton = new ButtonElement("Instantiate", [&]() { this->instantiateObject(); });
    // ButtonElement* recomputeErosionButton = new ButtonElement("Erosion values", [&]() { this->recomputeErosionValues(); });
    ButtonElement* spendTimeButton = new ButtonElement("Wait", [&]() {
        this->updateUntilStabilization();
        // this->updateEnvironmentFromEnvObjects(/* meh... I don't know if I should update the terrain or not */);
        this->updateSelectionMesh();

        Q_EMIT this->updated();
        //this->saveScene("testEnvObjects.json");
    });
    CheckboxElement* waitAtEachFrameButton = new CheckboxElement("Auto wait", this->waitAtEachFrame);
//    ButtonElement* createFromGAN = new ButtonElement("From GAN", [&]() { this->fromGanUI(); });
    ButtonElement* createFromFile = new ButtonElement("From file", [&]() { this->loadScene("EnvObjects/testEnvObjects.json"); });
    TextEditElement* testingFittingFormula = new TextEditElement("", "Fitting func: ");
    // testingFittingFormula->setOnTextChange([&](std::string expression) { this->evaluateAndDisplayCustomFittingFormula(expression); });
    TextEditElement* testingFitnessFormula = new TextEditElement("", "Fitness func: ");
    // testingFitnessFormula->setOnTextChange([&](std::string expression) { this->evaluateAndDisplayCustomFitnessFormula(expression); });
    testingFitnessFormula->setOnTextChange([&](std::string expression) { this->testedFitnessFunction = expression; this->evaluateAndDisplayCustomFitnessAndFittingFormula(this->testedFitnessFunction, this->testedFittingFunction); });
    testingFittingFormula->setOnTextChange([&](std::string expression) { this->testedFittingFunction = expression; this->evaluateAndDisplayCustomFitnessAndFittingFormula(this->testedFitnessFunction, this->testedFittingFunction); });
    ButtonElement* testPerformancesButton = new ButtonElement("Run test", [&]() { this->runPerformanceTest(); });
    ButtonElement* resetButton = new ButtonElement("Reset scene", [&]() { this->resetScene(); });
    CheckboxElement* addGroovesButton = new CheckboxElement("Spurs and grooves", displayGrooves);

    ButtonElement* instantiaABCbutton = new ButtonElement("Instantiate ABC", [&]() {
        this->instantiateSpecific("CoralPolypFlatA", GridF(), false);
        this->instantiateSpecific("CoralPolypFlatB", GridF(), false);
        this->instantiateSpecific("CoralPolypFlatC", GridF(), false);
    });

    // QLabel* label = new QLabel(QString::fromStdString("Objects: " + std::to_string(EnvObject::instantiatedObjects.size())));
    LabelElement* label = new LabelElement("Objects: " + std::to_string(EnvObject::instantiatedObjects.size()));

    objectsListWidget = new HierarchicalListUI;
    objectsListWidget->setSelectionMode(QAbstractItemView::SelectionMode::ExtendedSelection);
    updateObjectsList();
    objectsListWidget->setOnItemSelectionChanged([&]() { this->updateObjectsListSelection(); });
    // QObject::connect(objectsListWidget, &HierarchicalListWidget::itemSelectionChanged, this, [&]() { this->updateObjectsListSelection(); });


    std::vector<ComboboxLineElement> objectsChoices;
    for (auto& [name, obj] : EnvObject::availableObjects) {
        objectsChoices.push_back(ComboboxLineElement{name, 0});
    }
    ButtonElement* showButton = new ButtonElement("Show", [&](){
        this->displayProbas(objectCombobox->choices[objectCombobox->combobox()->currentIndex()].label);
    });
    ButtonElement* forceButton = new ButtonElement("Force", [&](){
        this->instantiateSpecific(objectCombobox->choices[objectCombobox->combobox()->currentIndex()].label, GridF(), true, true);
        // this->updateSelectionMesh();

        Q_EMIT this->updated();
    });
    objectCombobox = new ComboboxElement("Objects", objectsChoices);

    // std::vector<QWidget*> materialsButtons;
    // for (auto& [name, material] : EnvObject::materials) {
    //     ButtonElement* showButton = new ButtonElement("Show " + toCapitalize(name), [&](){ this->displayMaterialDistrib(name); });
    //     materialsButtons.push_back(showButton->get());
    // }

    std::vector<UIElement*> materialsButtons;
    for (auto& [name, material] : EnvObject::materials) {
        ButtonElement* showButton = new ButtonElement("Show " + toCapitalize(name), [&](){ this->displayMaterialDistrib(name); });
        materialsButtons.push_back(showButton);
    }

    ButtonElement* editFocusAreaButton = new ButtonElement("Edit focus", [&]() { this->manualModificationOfFocusArea(); });
    ButtonElement* editFlowfieldButton = new ButtonElement("Edit flowfield", [&]() { this->manualModificationOfFlowfield(); });
    ButtonElement* showElementsOnCanvasButton = new ButtonElement("Show all", [&]() { this->showAllElementsOnPlotter(); });

    ButtonElement* saveButton = new ButtonElement("Save", [&]() {this->saveScene("EnvObjects/testEnvObjects.json");});

    SliderElement* flowErosionSlider = new SliderElement("Erode", -10.f, 10.f, .1f);
    flowErosionSlider->setOnValueChanged([&](float newValue) {
        this->flowErosionFactor = newValue;
        this->addObjectsHeightmaps();
        this->flowErosionSimulation();

        Q_EMIT this->updated();
    });

    CheckboxElement* newObjectCreationBox = new CheckboxElement("Manual creation", [&](bool checked) {
        this->manuallyCreatingObject = checked;
        this->startNewObjectCreation();
        this->updateSelectionMesh();
        this->updateNewObjectMesh();

        Q_EMIT this->updated();
    });

    ButtonElement* nextStepButton = new ButtonElement("Step", [&]() {
        forceScenarioInterruption = true;
        this->runNextStep();
    });
    ButtonElement* runButton = new ButtonElement("Run", [&]() {
        if (!this->forceScenarioInterruption) {
            std::cout << "Stopping scenario" << std::endl;
            this->forceScenarioInterruption = true;
        } else {
            this->runScenario();
        }
    });

    CheckboxElement* displayCurrentsButton = new CheckboxElement("Flow", this->displayFlow);

//     layout->addWidget(newObjectCreationBox->get());
//     layout->addWidget(createHorizontalGroupUI({spendTimeButton, nextStepButton, runButton})->get());
//     layout->addWidget(waitAtEachFrameButton->get());
//     layout->addWidget(createMultiColumnGroup(materialsButtons, 2));
// //    layout->addWidget(showFlowfieldButton->get());
// //    layout->addWidget(createFromGAN->get());
//     layout->addWidget(flowErosionSlider->get());
//     layout->addWidget(objectCombobox->get());
//     layout->addWidget(createMultiColumnGroup({showButton->get(), forceButton->get()}, 2));
//     layout->addWidget(createHorizontalGroupUI({editFocusAreaButton, editFlowfieldButton})->get());
//     layout->addWidget(showElementsOnCanvasButton->get());
//     layout->addWidget(objectsListWidget);
//     layout->addWidget(testingFittingFormula->get());
//     layout->addWidget(createHorizontalGroupUI({instantiaABCbutton, testPerformancesButton, resetButton})->get());
//     layout->addWidget(addGroovesButton->get());
//     layout->addWidget(createHorizontalGroup({label, createFromFile->get(), saveButton->get(), displayCurrentsButton->get()}));

    ui->add({
             createHorizontalGroupUI({newObjectCreationBox, waitAtEachFrameButton}),
             createHorizontalGroupUI({spendTimeButton, nextStepButton, runButton}),
             createMultiColumnGroupUI(materialsButtons, 2),
             flowErosionSlider,
             objectCombobox,
             createMultiColumnGroupUI({showButton, forceButton}, 2),
             createHorizontalGroupUI({editFocusAreaButton, editFlowfieldButton}),
             showElementsOnCanvasButton,
             objectsListWidget,
             createVerticalGroupUI({testingFitnessFormula, testingFittingFormula}),
             createHorizontalGroupUI({/*instantiaABCbutton, testPerformancesButton, */resetButton}),
             addGroovesButton,
             createHorizontalGroupUI({label, createFromFile, saveButton, displayCurrentsButton})
    });

    return ui->get()->layout();
    // return layout;
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
            EnvCurve* passe = EnvCurve::instantiate("passe");
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
//            EnvArea* reef = EnvArea::instantiate("reef"));
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
            EnvArea* lagoon = EnvArea::instantiate("lagoon");
            lagoon->curve = curve;
        }
    }

    if (objectsListWidget != nullptr)
        updateObjectsList();
}

void EnvObjsInterface::setMaterialsDefinitionFile(std::string filename)
{
    this->materialsDefinitionFile.path = filename;
    EnvObject::readEnvMaterialsFile(filename);
}

void EnvObjsInterface::setDefinitionFile(std::string filename)
{
//    this->primitiveDefinitionFile = filename;
    this->primitiveDefinitionFile.path = filename;
    EnvObject::readEnvObjectsFile(filename);
}

void EnvObjsInterface::setTransformationsFile(std::string filename)
{
    this->transformationsFile.path = filename;
    EnvObject::readEnvMaterialsTransformationsFile(filename);
}

void EnvObjsInterface::setScenarioFile(std::string filename)
{
    this->scenarioFile.path = filename;
    EnvObject::readScenarioFile(filename);
}

void EnvObjsInterface::show()
{
    ActionInterface::show();
}

void EnvObjsInterface::hide()
{
//    updateObjectsListSelection(nullptr); // Hide single object's data
    this->manuallyCreatingObject = false;
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    ActionInterface::hide();
}

void EnvObjsInterface::afterTerrainUpdated()
{

}

void EnvObjsInterface::afterWaterLevelChanged()
{
    EnvObject::recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
}

void EnvObjsInterface::mouseClickedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    if (!this->isVisible()) return;
    if (!mouseInMap) return;

    bool moveSingleVertex = event->modifiers().testFlag(Qt::KeyboardModifier::ShiftModifier);
    bool moveWholeObject = event->modifiers().testFlag(Qt::KeyboardModifier::ControlModifier);

    if (moveSingleVertex || moveWholeObject) {
        this->startDraggingObject(mouseWorldPosition, moveSingleVertex);
    }
    else if (this->manuallyCreatingObject) {
        bool addingPoint = event->buttons().testFlag(Qt::LeftButton);
        this->addPointOnNewObjectCreation(mouseWorldPosition, addingPoint);
    }
}

void EnvObjsInterface::mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model)
{
    if (!this->isVisible()) return;
    if (!mouseWorldPosition.isValid()) return;

    this->moveDraggedObject(mouseWorldPosition);
}

void EnvObjsInterface::mouseReleasedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model)
{
    if (!this->isVisible()) return;

    bool destroyObjects = !mouseInMap;

    this->endDraggingObject(destroyObjects);
}

void EnvObjsInterface::keyPressEvent(QKeyEvent *event)
{
    if (!this->isVisible()) return;
    if (this->manuallyCreatingObject) {
        if (event->key() == Qt::Key_Enter) {
            this->endNewObjectCreation();
        }
    }
}

GridF computeScoreMap(std::string objectName, const Vector3& dimensions, bool& possible, bool applyNormalization = false) {
    auto obj = EnvObject::availableObjects[objectName];
    GridF score = GridF(dimensions);
    score.iterateParallel([&](const Vector3& pos) {
        score(pos) = std::max(obj->evaluate(pos), 0.f);
    });
    if (score.max() > 1e-5) {
        possible = true;
        if (applyNormalization)
            score.normalizeUsing(NORMALIZE_METHOD::NORMALIZE_MINMAX);
    } else {
        possible = false;
    }
    return score;
}

Vector3 bestPositionForInstantiation(std::string objectName, const GridF& score) {
    GridV3 gradients = score.gradient();
    if (std::abs(score.max() - score.min()) < 1e-5)
        return Vector3::random(Vector3(score.sizeX, score.sizeY, 0));

    Vector3 bestPos(false);
    float bestScore = std::numeric_limits<float>::lowest();
    for (int tries = 0; tries < 20; tries++) {
        Vector3 p = Vector3::random(Vector3(), score.getDimensions().xy());
        for (int iter = 0; iter < 50; iter++) {
            p += gradients(p).normalize();
        }
        float currentScore = score(p);
        if (currentScore > bestScore) {
            bestScore = currentScore;
            bestPos = p;
        }
    }
    return bestPos;
    /*
    float scoreToObtain = random_gen::generate(score.sum() / 100.f);
    float bestScore = 10000;
    Vector3 bestPos = Vector3::invalid();

    int nbIterationsToFindBestPos = 0;
    while (scoreToObtain > 0) {
        Vector3 p = Vector3::random(Vector3(), score.getDimensions().xy());
        float scoreObtained = score(p);
        if (scoreObtained < 1e-6) continue;
        scoreToObtain -= scoreObtained;
        bestPos = p;
        nbIterationsToFindBestPos++;
    }
    return bestPos;*/
}

EnvObject* EnvObjsInterface::instantiateObjectAtBestPosition(std::string objectName, Vector3 position, const GridF& score) {
    bool verbose = false;
    Vector3 initialPosition = position;
    EnvObject* newObject = EnvObject::instantiate(objectName);

    auto objAsPoint = dynamic_cast<EnvPoint*>(newObject);
    auto objAsCurve = dynamic_cast<EnvCurve*>(newObject);
    auto objAsArea = dynamic_cast<EnvArea*>(newObject);

    GridV3 gradients = score.gaussianSmooth(2.f).gradient();
    gradients.raiseErrorOnBadCoord = false;
    gradients.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;

    if (objAsPoint) {
        // position = PositionOptimizer::getHighestPosition(position, score, gradients);
        if (!position.isValid() || position == Vector3())
            return nullptr;
    } else if (objAsCurve) {
        BSpline initialCurve;
        if (objAsCurve->curveFollow == EnvCurve::SKELETON) {
            initialCurve = CurveOptimizer::getSkeletonCurve(position, score, gradients, objAsCurve->length);
        } else if (objAsCurve->curveFollow == EnvCurve::ISOVALUE) {
            initialCurve = CurveOptimizer::followIsolevel(position, score, gradients, objAsCurve->length);
        } else if (objAsCurve->curveFollow == EnvCurve::GRADIENTS) {
            initialCurve = CurveOptimizer::getExactLengthCurveFollowingGradients(position, score, gradients, objAsCurve->length);
        }
        if (initialCurve.size() == 0) {
            EnvObject::removeObject(newObject);
            delete newObject;
            return nullptr;
        }
        BSpline curve = initialCurve;
        curve.resamplePoints(10);
        position = curve[curve.size() / 2];
        curve.translate(-position);
        objAsCurve->curve = curve;
    } else if (objAsArea) {
        ShapeCurve initialCurve = AreaOptimizer::getAreaOptimizedShape(position, score, gradients, objAsArea->length * objAsArea->width);
        if (initialCurve.size() == 0) {
            EnvObject::removeObject(newObject);
            delete newObject;
            return nullptr;
        }

        position = initialCurve.centroid(); //initialCurve[0]; // The optimisation process might have moved the evaluation position greatly
        ShapeCurve curve = initialCurve.close().resamplePoints();
        curve.translate(-position);
        objAsArea->curve = curve.resamplePoints(10);
    }

    newObject->translate(position.xy());
    // newObject->recomputeEvaluationPoints();
    newObject->evaluationPositions = {initialPosition};
    int nbOutside = 0;
    for (auto& p : newObject->evaluationPositions) {
        if (!Vector3::isInBox(p, Vector3(), score.getDimensions())) {
            nbOutside++;
        }
    }
    if (nbOutside > newObject->evaluationPositions.size() / 2) {
    // if (!Vector3::isInBox(newObject->evaluationPosition, Vector3(), score.getDimensions())) {
        log("Object is outside...", verbose);
        this->destroyEnvObject(newObject, false, false);
        return nullptr;
    }
    // newObject->fitnessScoreAtCreation = newObject->evaluate(newObject->evaluationPosition);
    newObject->fitnessScoreAtCreation = newObject->evaluate();
    newObject->spawnTime = EnvObject::currentTime;
    if (newObject->fitnessScoreAtCreation <= newObject->minScore) {
        log("Object has score of 0...", verbose);
        this->destroyEnvObject(newObject, false, false);
        return nullptr;
    } else {
        this->log("Creation of obj at score = " + std::to_string(newObject->fitnessScoreAtCreation), verbose);
        return newObject;
    }
}

EnvObject *EnvObjsInterface::instantiateObjectUsingSpline(std::string objectName, const BSpline &spline)
{
    EnvObject* newObject = EnvObject::instantiate(objectName);

    auto objAsPoint = dynamic_cast<EnvPoint*>(newObject);
    auto objAsCurve = dynamic_cast<EnvCurve*>(newObject);
    auto objAsArea = dynamic_cast<EnvArea*>(newObject);

    Vector3 position;
    if (objAsPoint) {
        position = spline.points.back();
    } else if (objAsCurve) {
        BSpline curve = spline;
        curve.resamplePoints(10);
        position = curve[curve.size() / 2];
        curve.translate(-position);
        objAsCurve->curve = curve;
    } else if (objAsArea) {
        position = spline[0]; // The optimisation process might have moved the evaluation position greatly
        ShapeCurve curve = spline;
        curve.translate(-position);
        objAsArea->curve = curve.resamplePoints(10);
    }

    newObject->translate(position.xy());
    // newObject->fitnessScoreAtCreation = newObject->evaluate(newObject->evaluationPosition);
    newObject->fitnessScoreAtCreation = newObject->evaluate();
    newObject->spawnTime = EnvObject::currentTime;
    std::cout << "Manual creation of obj at score = " << newObject->fitnessScoreAtCreation << std::endl;
    return newObject;
}

EnvObject* EnvObjsInterface::instantiateSpecific(std::string objectName, GridF score, bool waitForFullyGrown, bool updateScreen)
{
    EnvObject* result = nullptr;
    objectName = toLower(objectName);
    if (EnvObject::availableObjects.count(objectName) == 0) {
        std::cerr << "No object " << objectName << " in database!" << std::endl;
        return result;
    }
    bool verbose = false;
    displayProcessTime("Instantiate new " + objectName + " (forced)... ", [&]() {
        // GridF heights = heightmap->getHeights();

        bool possible = true;
        if (score.empty())
            score = computeScoreMap(objectName, EnvObject::flowfield.getDimensions(), possible) * focusedArea;
        score.raiseErrorOnBadCoord = false;
        score.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;

        if (possible) {
            Vector3 bestPosition(false);
            float bestScore = -1;
            for (int iteration = 0; iteration < 1; iteration++) {
                Vector3 bestPos = bestPositionForInstantiation(objectName, score);
                if (!bestPos.isValid()) continue;
                EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPos, score);
                if (!newObject) continue;
                if (newObject->fitnessScoreAtCreation > bestScore) {
                    bestPosition = bestPos;
                    bestScore = newObject->fitnessScoreAtCreation;
                }
                this->destroyEnvObject(newObject, false, false);
            }
            if (!bestPosition.isValid()) {
                this->log("No valid position for object " + objectName, verbose);
                return;
            }
            EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPosition, score);
            if (!newObject) {
                this->log("Object not created", verbose);
                return;
            }
            auto implicit = newObject->createImplicitPatch(subsidedHeightmap);
            newObject->fittingFunction = EnvObject::parseFittingFunction(newObject->s_FittingFunction, newObject->name, true, newObject);
            newObject->fitnessFunction = EnvObject::parseFittingFunction(newObject->s_FitnessFunction, newObject->name, true, newObject);
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);

            if (implicit != nullptr) {
                rootPatch->addChild(implicit);
            }
            // Wait until the object is 100% grown:
            int maxIterations = 100;
            while (waitForFullyGrown && newObject && newObject->computeGrowingState() < 1.f) {
                this->updateEnvironmentFromEnvObjects(false, updateScreen);
                maxIterations--;
                if (maxIterations < 0) break;
                if (!isIn(newObject, EnvObject::instantiatedObjects)) {
                    return; // Object died in this process, stop this function now
                }
            }
            this->currentSelections = {newObject};
            EnvObject::recomputeTerrainPropertiesForObject(objectName);
            this->updateEnvironmentFromEnvObjects(implicit != nullptr, updateScreen); // If implicit is null, don't update the map
            result = newObject;
            this->materialSimulationStable = false; // We have to compute the simulation again
            return;
        } else {
            // std::cout << "Nope, impossible to instantiate..." << std::endl;
            return;
        }
    }, verbose);
    if (updateScreen)
        updateObjectsList();
    return result;
}

EnvObject *EnvObjsInterface::fakeInstantiate(std::string objectName, GridF score)
{
    EnvObject* result = nullptr;
    objectName = toLower(objectName);
    if (EnvObject::availableObjects.count(objectName) == 0) {
        this->error("No object " + objectName + " in database!");
        return result;
    }
    GridF heights = heightmap->getHeights();

    bool possible = !score.empty();
    if (score.empty())
        score = computeScoreMap(objectName, heights.getDimensions(), possible) * focusedArea;
    score.raiseErrorOnBadCoord = false;
    score.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;

    if (possible) {
        Vector3 bestPosition = bestPositionForInstantiation(objectName, score);
        EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPosition, score);
        if (!newObject) {
            return nullptr;
        }
        newObject->fittingFunction = EnvObject::parseFittingFunction(newObject->s_FittingFunction, newObject->name, true, newObject);
        newObject->fitnessFunction = EnvObject::parseFittingFunction(newObject->s_FitnessFunction, newObject->name, true, newObject);
        destroyEnvObject(newObject, false, false);
        result = newObject;
    }
    return result;
}

bool EnvObjsInterface::checkIfObjectShouldDie(EnvObject *obj, float limitFactorForDying)
{
    if (obj->createdManually) return false;
    return (obj->computeGrowingState2() <= limitFactorForDying);
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
    highErosionsMesh.fromArray(flattenArray(Mesh::applyMarchingCubes((-erosionGrid) - iso).getTriangles()));
    highDepositionMesh.fromArray(flattenArray(Mesh::applyMarchingCubes(erosionGrid - iso).getTriangles()));
}

void EnvObjsInterface::runNextStep()
{
    bool verbose = false;
    EnvObject* createdObject = nullptr;
    Scenario& scenario = EnvObject::scenario;

    for (auto& nextObject : scenario.nextObjects()) {
        bool possible;
        GridF score = computeScoreMap(nextObject.objectName, subsidedHeightmap.getDimensions(), possible) * focusedArea;
        if (possible) {
            createdObject = this->instantiateSpecific(nextObject.objectName, score, false);
        }
    }
    /*
    int maxTries = 1;
    auto nextObject = scenario.nextObject();
    bool possible;
    GridF score;
    displayProcessTime("Compute score... ", [&]() {
        score = computeScoreMap(nextObject.objectName, subsidedHeightmap.getDimensions(), possible) * focusedArea;
    }, verbose);
    if (possible) {
        while (!createdObject && maxTries > 0) {
            displayProcessTime("Instantiate " + nextObject.objectName + "... ", [&]() {
                createdObject = this->instantiateSpecific(nextObject.objectName, score, false);
            }, verbose);
            maxTries --;
        }
    }*/
    updateEnvironmentFromEnvObjects(false, false, false);
    EnvObject::currentTime += scenario.dt;
    this->log("Step: " + std::to_string(EnvObject::scenario.currentTime()));
}

void EnvObjsInterface::runScenario()
{
    updateEnvironmentFromEnvObjects(true, false, false);
    Scenario& scenario = EnvObject::scenario;
    this->forceScenarioInterruption = false;
    if (scenario.waterLevel >= 0) {
        (dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get()))->setWaterLevel(scenario.waterLevel);
    } else {
        scenario.waterLevel = heightmap->properties->waterLevel;
    }
    scenario.startTime = EnvObject::currentTime;
    while (!scenario.finished() && !forceScenarioInterruption) {
        float time = scenario.currentTime();
        float waterLevel = scenario.computeWaterLevel();
        (dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get()))->setWaterLevel(waterLevel);
        // this->log("Step: " + std::to_string(time) + (scenario.duration > 0 ? " / " + std::to_string(scenario.duration) : ""));
        runNextStep();
    }
    this->forceScenarioInterruption = true;
}

void EnvObjsInterface::updateEnvironmentFromEnvObjects(bool updateImplicitTerrain, bool emitUpdateSignal, bool killObjectsIfPossible)
{
    bool verbose = false;

    GridF subsidenceFactor = EnvObject::scenario.computeSubsidence(initialHeightmap.getDimensions());
    subsidedHeightmap = initialHeightmap * subsidenceFactor;

    std::vector<EnvObject*> immatureObjects;
    if (killObjectsIfPossible) {
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (obj->computeGrowingState() < 1.f) {
                materialSimulationStable = false;
                immatureObjects.push_back(obj);
            }
        }

        bool atLeastOneDeath = false;
        std::set<std::string> deadObjects;
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (isIn(obj, immatureObjects)) continue;
            bool shouldDie = this->checkIfObjectShouldDie(obj, .1f);
            atLeastOneDeath |= shouldDie;
            if (shouldDie) {
                float startingScore = obj->fitnessScoreAtCreation;
                // float endingScore = obj->evaluate(obj->evaluationPosition);
                float endingScore = obj->evaluate();
                this->log(obj->name + " went from " + std::to_string(startingScore) + " to " + std::to_string(endingScore) + " -> " + std::to_string(std::round(100.f * endingScore / startingScore)) + "%");
                this->destroyEnvObject(obj);
                deadObjects.insert(obj->name);
            }
        }
        if (atLeastOneDeath) {
            updateImplicitTerrain = true;
            for (auto death : deadObjects) {
                EnvObject::recomputeTerrainPropertiesForObject(death);
            }
        }
    }

    displayProcessTime("Get impacted... ", [&]() {
        EnvObject::beImpactedByEvents();
    }, verbose);

    if (!this->materialSimulationStable) { // If the simulation is stable, don't do anything
        displayProcessTime("Apply effects... ", [&]() {
            bool bigChangesInMaterials = EnvObject::applyEffects(subsidedHeightmap, userFlowField + simulationFlowField);
            //this->materialSimulationStable = !bigChangesInMaterials;
        }, verbose);
        displayProcessTime("Recompute properties... ", [&]() {
            EnvObject::recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
        }, verbose);
        // Get original flowfield, do not accumulate effects (for now).
        displayProcessTime("Get velocity... ", [&]() {
            if (!this->fluidSimulationIsStable) {
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->mainDirection = Vector3();
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->setObstacles(voxelGrid->getVoxelValues());
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->recomputeVelocities();
                this->simulationFlowField = dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->getVelocities(EnvObject::flowfield.sizeX, EnvObject::flowfield.sizeY, EnvObject::flowfield.sizeZ);
                this->simulationFlowField *= .1f;
                this->fluidSimulationIsStable = true;
            }

            updateVectorFieldVisu();
        }, verbose);


    }
    for (auto& [obj, implicit] : this->implicitPatchesFromObjects) {
        obj->createImplicitPatch(subsidedHeightmap, dynamic_cast<ImplicitPrimitive*>(implicit));
    }
    displayProcessTime("Update heightmap... ", [&]() {
        // heightmap->fromVoxelGrid(*voxelGrid.get());
        this->addObjectsHeightmaps();
        this->flowErosionSimulation();
        // this->heightmap->heights = subsidedHeightmap.gaussianSmooth(1.f, true);
        if (displayDepositionOnHeightmap) {
            for (auto& [name, material] : EnvObject::materials) {
                this->heightmap->heights += material.currentState * material.virtualHeight;
            }
        }
    }, verbose);
    EnvObject::recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());

    if (updateImplicitTerrain) {
        /*
        this->fluidSimulationIsStable = true; // Consider that the surface will change here
        for (auto& [obj, implicit] : this->implicitPatchesFromObjects) {
            obj->createImplicitPatch(subsidedHeightmap, dynamic_cast<ImplicitPrimitive*>(implicit));
        }
        // dynamic_cast<ImplicitPrimitive*>(implicitTerrain->composables[0])->cachedHeightmap *= .99f;
        // rootPatch->composables[0]->update();
        rootPatch->reevaluateAll();
        for (auto& patch : implicitTerrain->composables) {
            if (auto asHeightmap = dynamic_cast<ImplicitPrimitive*>(patch)) {
                if (asHeightmap->predefinedShape == ImplicitPatch::PredefinedShapes::ImplicitHeightmap) {
                    asHeightmap->cachedHeightmap = subsidedHeightmap * 1.f;
                }
            }
        }
        implicitTerrain->updateCache();
        implicitTerrain->update();
        displayProcessTime("Update voxels from implicit terrain... ", [&]() {
            // voxelGrid->fromImplicit(implicitTerrain.get(), 20);
        }, verbose);
        */
        displayProcessTime("Update voxels from heightmap... ", [&]() {
                voxelGrid->from2DGrid(*heightmap);
            }, verbose);
    }

    for (auto& obj : immatureObjects) {
        if (obj->computeGrowingState() >= 1.f) {
            // Got mature during this process -> now let's save the fitting score
            // obj->fitnessScoreAtCreation = obj->evaluate(obj->evaluationPosition);
            obj->fitnessScoreAtCreation = obj->evaluate();
        }
    }

    if (emitUpdateSignal) {
        Q_EMIT this->updated();
        updateObjectsList();
    }
}

void EnvObjsInterface::updateUntilStabilization()
{
    auto heights = heightmap->getHeights();
    EnvObject::stabilizeMaterials(heights);
}

void EnvObjsInterface::onlyUpdateFlowAndSandFromEnvObjects()
{
    EnvObject::applyEffects(subsidedHeightmap, userFlowField + simulationFlowField);
}

void EnvObjsInterface::destroyEnvObject(EnvObject *object, bool applyDying, bool recomputeTerrainPropertiesForObject)
{
    if (!object) return;

    if (applyDying)
        object->die();
    for (size_t i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
        if (EnvObject::instantiatedObjects[i] == object) {
            EnvObject::instantiatedObjects.erase(EnvObject::instantiatedObjects.begin() + i);
            break;
        }
    }
    if (this->implicitPatchesFromObjects.count(object) != 0) {
        for (size_t i = 0; i < rootPatch->composables.size(); i++) {
            if (rootPatch->composables[i] == this->implicitPatchesFromObjects[object])
                rootPatch->composables.erase(rootPatch->composables.begin() + i);
        }
        this->implicitPatchesFromObjects.erase(object);
    }
    if (recomputeTerrainPropertiesForObject)
        EnvObject::recomputeTerrainPropertiesForObject(object->name);
}

void EnvObjsInterface::displayProbas(std::string objectName)
{
    focusAreaEditing = false;
    flowfieldEditing = false;
    previewingObjectInPlotter = true;
    // currentlyPreviewedObject = objectName;
    Vector3 dimensions = initialHeightmap.getDimensions();
    bool possible;
    GridF score = computeScoreMap(objectName, dimensions, possible, false);
    if (!possible) {
        Plotter::get("Object Preview")->addImage(score * 0.f);
    } else {
        float smallestPositive = score.max();
        score.iterate([&](size_t i) {
            if (score[i] > 0.f && score[i] < smallestPositive)
                smallestPositive = score[i];
        });
        score.iterateParallel([&](size_t i) {
            score[i] = std::max(score[i], smallestPositive);
        });
        Plotter::get("Object Preview")->addImage(score);
    }
    Plotter::get("Object Preview")->show();
    dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(score);
    Q_EMIT updated();
}

void EnvObjsInterface::displayMaterialDistrib(std::string materialName)
{
    GridF distribution = EnvObject::materials[materialName].currentState;
    Plotter::get("Material")->addImage(distribution);
    Plotter::get("Material")->show();
    dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(distribution);
    Q_EMIT updated();
}

/*void EnvObjsInterface::displaySedimentsDistrib()
{
    GridF sediments = EnvObject::materials["sand"].currentState;
    Plotter::get()->addImage(sediments);
    Plotter::get()->show();
}

void EnvObjsInterface::displayPolypDistrib()
{
    GridF polyp = EnvObject::materials["polyp"].currentState;
    Plotter::get()->addImage(polyp);
    Plotter::get()->show();
}*/

void EnvObjsInterface::manualModificationOfFocusArea()
{
    this->focusAreaEditing = true;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = false;
    Plotter::get("Focus")->addImage(this->renderFocusArea());
    Plotter::get("Focus")->show();
}

void EnvObjsInterface::manualModificationOfFlowfield()
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = true;
    this->previewingObjectInPlotter = false;
    Plotter::get("Flowfield")->addImage(this->renderFlowfield());
    Plotter::get("Flowfield")->show();
}

void EnvObjsInterface::updateObjectsList()
{
    if (!this->isVisible()) return;
    if (!objectsListWidget) return;
    std::vector<int> currentSelectionsIDs;
    for (auto currentSelection : currentSelections) currentSelectionsIDs.push_back(currentSelection->ID);
    objectsListWidget->clear();
    auto list = EnvObject::instantiatedObjects;

    for (auto& obj : list) {

        // float startingScore = obj->fitnessScoreAtCreation;
        // float endingScore = obj->evaluate(obj->evaluationPosition);

        std::string text = obj->name;
        if (obj->createdManually)
            text += " [*]";
        else
            text += " (" + std::to_string(int(obj->computeGrowingState() * 100.f)) + "% -- " + std::to_string(int(100.f * obj->computeGrowingState2())) + "% -- " + std::to_string(obj->evaluate()) + "/" + std::to_string(obj->fitnessScoreAtCreation) + ")";
        objectsListWidget->addItem(new HierarchicalListWidgetItem(text, obj->ID, 0));
        if (isIn(obj->ID, currentSelectionsIDs))
            currentSelections.push_back(obj);
    }
    objectsListWidget->setCurrentItems(currentSelectionsIDs);
}

void EnvObjsInterface::updateObjectsListSelection(QListWidgetItem *__newSelectionItem)
{
    currentSelections.clear();

    for (auto newSelectionItem : objectsListWidget->selectedItems()) {
        auto newSelection = dynamic_cast<HierarchicalListWidgetItem*>(newSelectionItem);
        if (!newSelection) {
            continue; // Does it happen?
        }
        int objID = newSelection->ID;
        EnvObject* selection = nullptr;
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (obj->ID == objID) {
                selection = obj;
                break;
            }
        }
        currentSelections.push_back(selection);
    }

    this->updateSelectionMesh();
}

void EnvObjsInterface::updateSelectionMesh()
{
    if (currentSelections.empty()) {
        selectedObjectsMesh.clear();
        // velocitiesMesh.clear();
        return;
    }

    // this->velocitiesMesh.fromArray(std::vector<float>{});
    this->selectedObjectsMesh.fromArray(std::vector<float>{});
    Vector3 selectionPos;
    std::vector<Vector3> lines;
    std::vector<Vector3> colors;
    float offsetAbove = 5.f;
    for (auto currentSelection : currentSelections) {
        // Vector3 evalPos = currentSelection->evaluationPosition;
        for (auto evalPos : currentSelection->evaluationPositions) {
            evalPos.z = subsidedHeightmap.interpolate(evalPos.x, evalPos.y) + offsetAbove;
            std::vector<Vector3> evalLines = {evalPos - Vector3(2, 2, 0), evalPos + Vector3(2, 2, 0), evalPos - Vector3(-2, 2, 0), evalPos + Vector3(-2, 2, 0)};
            lines.insert(lines.end(), evalLines.begin(), evalLines.end());
            std::vector<Vector3> evalColors = std::vector<Vector3>(evalLines.size(), Vector3(0.5, 0.5, 1));
            colors.insert(colors.end(), evalColors.begin(), evalColors.end());
        }

        if (auto asPoint = dynamic_cast<EnvPoint*>(currentSelection)) {
            continue; // Do not display the points!!
            selectionPos = asPoint->position;
            selectionPos.z = subsidedHeightmap.interpolate(selectionPos.x, selectionPos.y) + offsetAbove;
            std::vector<Vector3> meshPoints = Mesh::getPointsForArrow(selectionPos + Vector3(0, 0, 20), selectionPos);
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(currentSelection)) {
            selectionPos = asCurve->curve.center();
            std::vector<Vector3> meshPoints;
            auto path = asCurve->curve.getPath(50);
            for (size_t i = 0; i < path.size() - 1; i++) {
                auto p1 = path[i];
                auto p2 = path[i + 1];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + offsetAbove);
                Vector3 p2leveled = p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x, p2.y) + offsetAbove);
                meshPoints.push_back(p1leveled);
                meshPoints.push_back(p2leveled);
            }
            for (int i = 0; i < asCurve->curve.size(); i++) {
                auto& p1 = asCurve->curve[i];
                auto& p2 = asCurve->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + offsetAbove);
                Vector3 perpendicular = (p2 - p1).rotate(0, 0, deg2rad(90)).normalized() * 1.f;
                meshPoints.push_back(p1leveled + perpendicular);
                meshPoints.push_back(p1leveled - perpendicular);
            }
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else if (auto asArea = dynamic_cast<EnvArea*>(currentSelection)) {
            selectionPos = asArea->curve.center();
            std::vector<Vector3> meshPoints;
            auto path = asArea->curve.getPath(20);
            for (size_t i = 0; i < path.size() - 1; i++) {
                auto p1 = path[i];
                auto p2 = path[i + 1];
                meshPoints.push_back(p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + 5.f));
                meshPoints.push_back(p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x, p2.y) + 5.f));
            }
            for (int i = 0; i < asArea->curve.size(); i++) {
                auto& p1 = asArea->curve[i];
                auto& p2 = asArea->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + offsetAbove);
                Vector3 perpendicular = (p2 - p1).rotate(0, 0, deg2rad(90)).normalized() * 1.f;
                meshPoints.push_back(p1leveled + perpendicular);
                meshPoints.push_back(p1leveled - perpendicular);
            }
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else {
            std::cerr << "Object #" << currentSelection->ID << " (" << currentSelection->name << ") could not be casted to Point, Curve or Area..." << std::endl;
//            return;
            continue;
        }

        /*
        GridV3 initialFlow = EnvObject::initialFlowfield;
        GridV3 flow;
        GridF occupancy;
        std::tie(flow, occupancy) = currentSelection->computeFlowModification();
        float currentGrowth = currentSelection->computeGrowingState();
        occupancy *= currentGrowth;

        flow = flow * occupancy + initialFlow * (1.f - occupancy);
        initialFlow = initialFlow * (1.f - EnvObject::flowImpactFactor) + flow * EnvObject::flowImpactFactor;
    //    velocitiesMesh.fromVectorField(initialFlow.resize(30, 30, 1), voxelGrid->getDimensions());
        Mesh::createVectorField(initialFlow.resize(30, 30, 1), voxelGrid->getDimensions(), &velocitiesMesh, 1.f, false, true);
        */
    }
    selectedObjectsMesh.colorsArray = colors;
    selectedObjectsMesh.fromArray(lines);
}

void EnvObjsInterface::updateNewObjectMesh()
{
    // auto objectModel = EnvObject::availableObjects[getCurrentObjectName()];

    newObjectMesh.clear();

    std::vector<Vector3> lines;
    std::vector<Vector3> colors;
    float offsetAbove = 5.f;

    std::vector<Vector3> meshPoints;
    auto path = objectSkeletonCreation.getPath(50);
    for (size_t i = 0; i < path.size() - 1; i++) {
        auto p1 = path[i];
        auto p2 = path[i + 1];
        Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + offsetAbove);
        Vector3 p2leveled = p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x, p2.y) + offsetAbove);
        meshPoints.push_back(p1leveled);
        meshPoints.push_back(p2leveled);
    }
    for (int i = 0; i < objectSkeletonCreation.size(); i++) {
        auto& p1 = objectSkeletonCreation[i];
        auto& p2 = objectSkeletonCreation[std::abs(i - 1)];
        Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x, p1.y) + offsetAbove);
        Vector3 perpendicular = (p2 - p1).rotate(0, 0, deg2rad(90)).normalized() * 1.f;
        meshPoints.push_back(p1leveled + perpendicular);
        meshPoints.push_back(p1leveled - perpendicular);
    }
    lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
    std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
    colors.insert(colors.end(), meshColors.begin(), meshColors.end());

    newObjectMesh.colorsArray = colors;
    newObjectMesh.fromArray(lines);
}

void EnvObjsInterface::updateObjectsDefinitions(const std::string &newDefinition)
{
    try {
        EnvObject::readEnvObjectsFileContent(newDefinition);
        this->previousFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << primitiveDefinitionFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousFileContent != "")
            EnvObject::readEnvObjectsFileContent(this->previousFileContent);
    }
}

void EnvObjsInterface::updateMaterialsDefinitions(const std::string &newDefinition)
{
    try {
        EnvObject::readEnvMaterialsFileContent(newDefinition);
        this->previousMaterialsFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << materialsDefinitionFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousMaterialsFileContent != "")
            EnvObject::readEnvMaterialsFileContent(this->previousMaterialsFileContent);
    }
}

void EnvObjsInterface::updateMaterialsTransformationsDefinitions(const std::string &newDefinition)
{
    try {
        EnvObject::readEnvMaterialsTransformationsFileContent(newDefinition);
        this->previousMaterialsTransformationsFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << transformationsFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousMaterialsTransformationsFileContent != "")
            EnvObject::readEnvMaterialsTransformationsFileContent(this->previousMaterialsTransformationsFileContent);
    }
}

void EnvObjsInterface::updateScenarioDefinition(const std::string &newDefinition)
{
    try {
        EnvObject::readScenarioFileContent(newDefinition);
        this->previousScenarioFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << scenarioFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousScenarioFileContent != "")
            EnvObject::readScenarioFileContent(this->previousScenarioFileContent);
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFitnessFormula(std::string formula)
{
    EnvPoint fake;
    // fake.s_FittingFunction = formula;
    fake.s_FitnessFunction = formula;
    try {
        // fake.fittingFunction = EnvObject::parseFittingFunction(formula, "");
        fake.fitnessFunction = EnvObject::parseFittingFunction(formula, "");

        GridF eval(EnvObject::flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3& p) {
            eval(p) = fake.fitnessFunction(p);
        });
        Plotter::get("Fitness Function")->addImage(eval);
        Plotter::get("Fitness Function")->show();
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFittingFormula(std::string formula)
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = true;

    EnvPoint fake;
    fake.s_FittingFunction = formula;
    // fake.s_FitnessFunction = formula;
    try {
        fake.fittingFunction = EnvObject::parseFittingFunction(formula, "");
        // fake.fitnessFunction = EnvObject::parseFittingFunction(formula, "");

        GridF eval(EnvObject::flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3& p) {
            eval(p) = fake.fittingFunction(p);
        });
        Plotter::get("Fitting Function")->addImage(eval);
        Plotter::get("Fitting Function")->show();
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFitnessAndFittingFormula(std::string fitnessFuncFormula, std::string fittingFuncFormula)
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = true;

    EnvPoint fake;
    fake.s_FittingFunction = (trim(fittingFuncFormula) == "" ? EnvObject::availableObjects[objectCombobox->getSelection().label]->s_FittingFunction : fittingFuncFormula);
    fake.s_FitnessFunction = (trim(fitnessFuncFormula) == "" ? EnvObject::availableObjects[objectCombobox->getSelection().label]->s_FitnessFunction : fitnessFuncFormula);
    try {
        fake.fittingFunction = EnvObject::parseFittingFunction(fake.s_FittingFunction, "");
        fake.fitnessFunction = EnvObject::parseFittingFunction(fake.s_FitnessFunction, "");

        GridV3 eval(EnvObject::flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3& p) {
            eval(p).x = fake.fitnessFunction(p);
            eval(p).y = fake.fittingFunction(p);
        });
        Plotter::get("Object Preview")->addImage(eval);
        Plotter::get("Object Preview")->show();
        // dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}


BSpline followGradient(const GridV3 gradients, const Vector3& startPoint, float maxDist, bool followInverse) {
    Vector3 pos = startPoint;
    BSpline path({pos});
    Vector3 dir;
    float totalDistance = 0.f;
    while (totalDistance < maxDist) {
        Vector3 gradient = gradients(pos);
        if (gradient == Vector3()) break; // Nowhere to go
        gradient.normalize();
        dir = gradient * (followInverse ? -1.f : 1.f);

        pos += dir;

        totalDistance += dir.norm();

        path.points.push_back(pos);

        int nbPoints = std::min(int(path.size()), 5);
        if (nbPoints > 4) {
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
}

std::vector<Vector3> findCandidatesPositions(const Vector3& startPosition, const Vector3& direction, float angle, float radius, int nbCandidates) {
    std::vector<Vector3> points(nbCandidates);
    float initialAngle = direction.getSignedAngleWith(Vector3(1, 0, 0));
    for (int i = 0; i < nbCandidates; i++) {
        float phi = interpolation::inv_linear(random_gen::generate(), -angle, angle);
        float r = std::sqrt(random_gen::generate()) * radius; // Use square root to svoid bias towards center of the disk

        points[i] = Vector3(r, std::sin(phi) * r).rotate(0, 0, initialAngle) + startPosition;
    }
    return points;
}

std::vector<BSpline> getCandidatesPaths(const GridV3& gradients, const std::vector<Vector3>& positions, float directionLength) {
    std::vector<BSpline> paths(positions.size());
    for (int i = 0; i < positions.size(); i++) {
        paths[i] = followGradient(gradients, positions[i], directionLength * .5f, false);
    }
    return paths;
}

BSpline getBestCandidatesPath(const GridF& score, const BSpline& initialPath, const std::vector<BSpline>& paths) {
    float longestDistance = 0.f;
    float smallestScore = std::numeric_limits<float>::max();
    int bestIndex = 0;
    for (int i = 0; i < paths.size(); i++) {
        if (paths[i].size() == 0) continue;
        float currentScore = score(paths[i].points.back()) / (initialPath.points[0] - paths[i].points.back()).norm2();
        if (currentScore < smallestScore) {
            smallestScore = currentScore;
            bestIndex = i;
        }
//            float dist = (initialPath.points[0] - paths[i].points.back()).norm2();
//            if (dist > longestDistance) {
//                longestDistance = dist;
//                bestIndex = i;
//            }
    }
    return paths[bestIndex];
}
BSpline followIsovalue(const GridF& values, const GridV3& gradients, const Vector3& startPoint, float maxDist) {
    BSpline finalPath;

    Vector3 pos = startPoint;
    float initialIsovalue = values.interpolate(pos);
    BSpline path({pos});
    Vector3 dir(0, 0, 0);
    bool didAFullCircle = false;
    float totalDistance = 0.f;
    while (maxDist > totalDistance && path.size() < 5000) {
        if (path.size() > 5 && (pos - startPoint).norm2() < 3*3){
            didAFullCircle = true;
            break; // Got back close to beginning
        }
        Vector3 gradient;
        int maxTries = 100;
        for (int iTry = 0; iTry < maxTries; iTry++) {
            Vector3 jitter = Vector3::random() * (5.f * float(iTry) / (float(maxTries)));
            if (jitter.dot(dir) < 0) continue;
            auto testPos = pos + jitter;
            gradient = gradients.interpolate(testPos);
            if (gradient.norm2() > 1e-8) {
                pos = testPos;
                break; // Nowhere to go
            }
        }
        if (gradient.norm2() < 1e-8) break; // Nowhere to go
        gradient.normalize();

        Vector3 newDir = gradient.cross(Vector3(0, 0, 1));
        dir = newDir * (dir.dot(newDir) < 0 ? -1.f : 1.f);

        if (!newDir.isValid() || !gradient.isValid()) break;

        float newVal = values.interpolate(pos + dir);

        if (std::abs(newVal - initialIsovalue) < 1e-5) {
            Vector3 newGrad = gradients.interpolate(pos + dir);
            float bestRectificationScale = 0.f;
            float closestIso = std::numeric_limits<float>::max();
            for (int i = 0; i < 10; i++) {
                float scale = float(i) / 10.f * (newVal < initialIsovalue ? 1.f : -1.f);
                float newDiff = values.interpolate(pos + dir + newGrad * scale);
                if (newDiff < closestIso) {
                    closestIso = newDiff;
                    bestRectificationScale = scale;
                }
            }
            dir += newGrad * bestRectificationScale;
        }

        pos += dir;

        totalDistance += dir.norm();

        path.points.push_back(pos);
    }

    finalPath = path;
    /*
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
            if (gradient.norm2() < 1e-8) break; // Nowhere to go
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
    }*/
    return finalPath;
}

BSpline EnvObjsInterface::computeNewObjectsCurveAtPosition(const Vector3 &seedPosition, const GridV3 &gradients, const GridF& score, float directionLength, float widthMaxLength, bool followIsolevel)
{
    Vector3 pos = seedPosition;
    BSpline isoline;
    if (followIsolevel) {
        isoline = followIsovalue(score, gradients, pos, directionLength);
    } else {
        Vector3 gradDir = gradients(pos).normalized();
        isoline = BSpline({pos - gradDir * directionLength * .5f, pos + gradDir * directionLength * .5f});
        isoline.resamplePoints(10);
        for (auto& p : isoline) {
            p += gradients(p).normalized().rotated(0, 0, deg2rad(90)) * widthMaxLength * random_gen::generate(-1, 1);
        }
    }
    return isoline;
}

ShapeCurve EnvObjsInterface::computeNewObjectsShapeAtPosition(const Vector3 &seedPosition, const GridV3& gradients, const GridF& score, float directionLength)
{
    BSpline isoline = followIsovalue(score, gradients, seedPosition, directionLength);
    return isoline;
}

ShapeCurve EnvObjsInterface::computeNewObjectsShapeAtPositionForceCircle(const Vector3 &seedPosition, const GridV3 &gradients, const GridF &score, float directionLength)
{
    ShapeCurve finalIsoline;
//    float targetArea = directionLength * _widthMaxLength;
    Vector3 pos = seedPosition;

//    int maxTries = 3;
//    float bestAreaDiff = std::numeric_limits<float>::max();
    ShapeCurve bestCurve;
    Vector3 jitterPos = pos;
//    while (maxTries > 0) {
        finalIsoline = this->computeNewObjectsShapeAtPosition(jitterPos, gradients, score, directionLength).close();
        if (finalIsoline.size() > 5 && (finalIsoline.points.front() - finalIsoline.points.back()).norm2() < 3*3) {
//            float area = finalIsoline.computeArea();
//            if (std::abs(area - targetArea) < std::abs(bestAreaDiff)) {
//                bestAreaDiff = area - targetArea;
                bestCurve = finalIsoline;
            }
//        } else {
//            jitterPos = pos + Vector3::random() * .1f;
//        }
//        maxTries--;
//    }
    return bestCurve;
}

ShapeCurve EnvObjsInterface::computeNewObjectsShapeAtPositionForceCircleOptimizedArea(const Vector3 &seedPosition, const GridV3 &gradients, const GridF &score, float directionLength, float targetArea)
{
    Vector3 currentSeedPos = seedPosition;
    float maxError = 5.f;
    int maxTries = 100;
    ShapeCurve finalCurve;
    float moveFactor = 1.f;
    bool currentlyAreaGettingSmaller = true;

    // We will move only in the direction of the gradient, since we want to optimize the isolevel.
    // And we know that higher isolevel => lower area while lower isolevel => higher area.
    // So isolevel gradient proportional to -area gradient.
    while (maxTries > 0) {
        ShapeCurve curve = computeNewObjectsShapeAtPositionForceCircle(currentSeedPos, gradients, score, directionLength);
        if (curve.size() == 0) {
            // The isocontour is too big, we didn't manage to do a full circle.
            currentSeedPos = currentSeedPos + gradients.interpolate(currentSeedPos).normalized() * 2.f;
        } else {
            float area = curve.computeArea();

            float diff = targetArea - area; // < 0 means curve too big, > 0 means curve too small
            finalCurve = curve;
            if (std::abs(diff) < maxError) break;
            currentSeedPos = currentSeedPos + gradients.interpolate(currentSeedPos).normalized() * (diff > 0 ? -1.f : 1.f) * moveFactor;

            if (currentlyAreaGettingSmaller != (diff > 0)) {
                currentlyAreaGettingSmaller = !currentlyAreaGettingSmaller;
                moveFactor *= .5f;
            }
        }
        maxTries--;
    }
    return finalCurve;
}

void EnvObjsInterface::runPerformanceTest()
{

    displayProcessTime("Benchmarking coral growth...", [&]() {
        for (int i = 0; i < 10; i++) {
            this->instantiateSpecific("coralpolyp");
        }
    });
    /*
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
    */
}

void EnvObjsInterface::resetScene()
{
    EnvObject::reset();
    this->simulationFlowField.reset();
    this->userFlowField.reset();
    this->materialSimulationStable = false; // We have to compute the simulation again
    for (auto& [obj, patch] : implicitPatchesFromObjects)
        delete patch;
    this->implicitPatchesFromObjects.clear();
//    this->rootPatch->deleteAllChildren();
    this->rootPatch->composables.clear();
    this->rootPatch->updateCache();

    this->updateEnvironmentFromEnvObjects(true);

    this->focusedArea = GridF(initialHeightmap.getDimensions(), 1.f);

    Q_EMIT this->updated();
}

void EnvObjsInterface::loadScene(std::string filename)
{
    this->resetScene();
    nlohmann::json json = nlohmann::json::parse(std::ifstream(filename));

    std::vector<nlohmann::json> allObjects = json["objects"];
    std::vector<nlohmann::json> allMaterials = json["materials"];
    if (json.contains("initialflow")) {
        std::string flowStr = json["initialflow"];
        EnvObject::initialFlowfield = loadGridV3(flowStr, false);
    }
    if (json.contains("userflow")) {
        this->userFlowField = loadGridV3(json["userflow"], false);
    }

    if (json.contains("waterlevel")) {
        float waterLevel = json["waterlevel"];
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->setWaterLevel(waterLevel);
    }
    if (json.contains("heightmap")) {
        initialHeightmap = loadGridF(json["heightmap"], false);
    }

    for (auto mat : allMaterials) {
        EnvObject::materials[mat["name"]].fromJSON(mat);
    }

    for (auto obj : allObjects) {
        std::string objectName = obj["name"];
        EnvObject* newObject = EnvObject::instantiate(objectName);
        newObject->age = obj["age"];
        newObject->fitnessScoreAtCreation = obj["fitnessScoreAtCreation"];
        // newObject->evaluationPosition = json_to_vec3(obj["evaluationPosition"]);
        /*std::vector<nlohmann::json> positions = obj["evaluationPositions"];
        for (auto position : positions) {

        }*/

        if (auto asPoint = dynamic_cast<EnvPoint*>(newObject)) {
            asPoint->position = json_to_vec3(obj["position"]);
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(newObject)) {
            asCurve->curve = json_to_bspline(obj["curve"]);
            newObject->createdManually = true;
        } else if (auto asArea = dynamic_cast<EnvArea*>(newObject)) {
            asArea->curve = json_to_bspline(obj["curve"]);
            newObject->createdManually = true;
        }
        newObject->recomputeEvaluationPoints();
        // newObject->createdManually = false;
        newObject->premature = true;
        auto implicit = newObject->createImplicitPatch(initialHeightmap);
        this->implicitPatchesFromObjects[newObject] = implicit;
        if (implicit != nullptr) {
            rootPatch->addChild(implicit);
        }
        this->currentSelections = {};
    }
    this->addObjectsHeightmaps();
    this->flowErosionSimulation();
    for (auto& [objectName, obj] : EnvObject::availableObjects) {
        EnvObject::recomputeTerrainPropertiesForObject(objectName);
    }
    EnvObject::recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
    updateEnvironmentFromEnvObjects(true, false, false);
    for (int i = 0; i < 20; i++)
        updateEnvironmentFromEnvObjects(false, false, false);
    updateEnvironmentFromEnvObjects(false, true, false);
    updateObjectsList();
    updateSelectionMesh();

    for (auto& obj : EnvObject::instantiatedObjects) {
        obj->premature = false;
        // newObject->createdManually = false;
    }
    updateEnvironmentFromEnvObjects(false, true, false);

    Q_EMIT this->updated();
}

void EnvObjsInterface::saveScene(std::string filename)
{
    nlohmann::json mainJson;
    std::vector<nlohmann::json> allObjects(EnvObject::instantiatedObjects.size());
    std::vector<nlohmann::json> allMaterials(EnvObject::materials.size());

    for (size_t i = 0; i < allObjects.size(); i++) {
        allObjects[i] = EnvObject::instantiatedObjects[i]->toJSON();
    }

    size_t i = 0;
    for (auto& [matName, material] : EnvObject::materials) {
        allMaterials[i] = material.toJSON();
        i++;
    }

    mainJson["objects"] = allObjects;
    mainJson["materials"] = allMaterials;
    mainJson["initialflow"] = stringifyGridV3(EnvObject::initialFlowfield, false);
    mainJson["userflow"] = stringifyGridV3(this->userFlowField, false);
    mainJson["waterlevel"] = heightmap->properties->waterLevel;
    mainJson["heightmap"] = stringifyGridF(initialHeightmap, false);
    std::ofstream out(filename);
    out << mainJson.dump(1, '\t');
    out.close();
}

GridV3 EnvObjsInterface::renderFocusArea() const
{
    GridV3 coloredFocus(this->focusedArea.getDimensions());
    this->focusedArea.iterateParallel([&](size_t i) {
        float value = focusedArea[i];
        coloredFocus[i] = colorPalette(value, {Vector3(1, 0, 0), Vector3(1, 1, 1), Vector3(0, 1, 0)}, {0.f, 1.f, 3.f});
    });
    return coloredFocus;
}

GridV3 EnvObjsInterface::renderFlowfield() const
{
    EnvObject::updateFlowfield(userFlowField + simulationFlowField + EnvObject::scenario.computeStorm(userFlowField.getDimensions()));
    GridV3& flow = EnvObject::flowfield;
    // return Plotter::get()->computeVectorFieldRendering(flow, 1/10.f, flow.getDimensions()  * 2.f).resize(flow.getDimensions());
    return Plotter::get("Flowfield")->computeStreamLinesRendering(flow, flow.getDimensions()  * 3.f);
}

void EnvObjsInterface::previewCurrentEnvObjectPlacement(Vector3 position)
{
    std::cout << "Preview... " << std::endl;
    GridV3 dataV3 = Plotter::get("Object Preview")->displayedImage;
    GridF fitnessScoreGrid(dataV3.getDimensions());
    GridF fittingScoreGrid(dataV3.getDimensions());
    dataV3.iterateParallel([&](size_t i) {
        fitnessScoreGrid[i] = dataV3[i].x;
        fittingScoreGrid[i] = dataV3[i].y;
    });

    auto obj = EnvObject::availableObjects[objectCombobox->getSelection().label];
    auto score = fittingScoreGrid;

    GridV3 result = GridV3(score.getDimensions(), Vector3(1, 1, 1)) * score;
    ShapeCurve isoline;
    // auto position = clickPos;
    // auto score = fitnessScoreGrid;
    // score = score.gaussianSmooth(5.f, true, true, -1);
    GridV3 gradients = score.gradient().gaussianSmooth(5.f, true, true, -1);
    fitnessScoreGrid.raiseErrorOnBadCoord = false;
    fitnessScoreGrid.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
    gradients.raiseErrorOnBadCoord = false;
    gradients.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;

    if (auto objAsPoint = dynamic_cast<EnvPoint*>(obj)) {
        //position = PositionOptimizer::getHighestPosition(position, score, gradients);
        isoline = ShapeCurve::circle(objAsPoint->radius, position, 20);
    } else if (auto objAsCurve = dynamic_cast<EnvCurve*>(obj)) {
        BSpline initialCurve;
        if (objAsCurve->curveFollow == EnvCurve::SKELETON) {
            float targetLength = objAsCurve->length;
            initialCurve = CurveOptimizer::getSkeletonCurve(position, score, gradients, targetLength);
        } else if (objAsCurve->curveFollow == EnvCurve::ISOVALUE) {
            initialCurve = CurveOptimizer::followIsolevel(position, score, gradients, objAsCurve->length);
        } else if (objAsCurve->curveFollow == EnvCurve::GRADIENTS) {
            initialCurve = CurveOptimizer::getExactLengthCurveFollowingGradients(position, score, gradients, objAsCurve->length);
        }
        isoline = initialCurve;
        isoline.closed = false;
    } else if (auto objAsArea = dynamic_cast<EnvArea*>(obj)) {
        ShapeCurve initialCurve;
        if (obj->snakeDefined) {
            float fakeRadius = std::sqrt(objAsArea->length * objAsArea->width) * .5f;
            float fakeArea = PI * fakeRadius * fakeRadius;

            ShapeCurve curve = ShapeCurve::circle(fakeRadius, position, 20);
            SnakeSegmentation s = obj->snake;
            s.contour = curve;
            s.image = score;
            s.gradientField = gradients;
            s.targetLength = 0;
            s.targetArea = fakeArea;
            s.collapseFirstAndLastPoint = true;

            int maxIterations = 10;
            for (int iteration = 0; iteration < maxIterations; iteration++) {
                initialCurve = s.runSegmentation(50);
                ShapeCurve display = initialCurve;
                display.points.push_back(display[0]);
                display.resamplePoints(100);
                if (iteration == 9) {
                    for (size_t i = 0; i < display.size(); i++) {
                        result(display[i]) = colorPalette(float(iteration) / float(maxIterations - 1));
                    }
                }
                if (iteration == 9) {
                    for (const auto& v : s.randomGreenCoords) {
                        result(computePointFromGreenCoordinates(v, ShapeCurve(s.contour))).z = 1;
                    }
                }
            }
        } else {
            initialCurve = AreaOptimizer::getAreaOptimizedShape(position, score, gradients, objAsArea->length * objAsArea->width);
        }
        isoline = initialCurve;
    }
    std::cout << isoline.toString() << std::endl;

    if (isoline.closed) {
        dataV3.iterateParallel([&](const Vector3& pos) {
            result(pos) += Vector3(.5f, .5f, .5f) * (isoline.containsXY(pos, false) ? 1.f : 0.f);
        });
    }
    int nbSamples = 500;
    auto path = isoline.getPath(nbSamples); // .resamplePoints(nbSamples).points;
    for (size_t i = 0; i < path.size(); i++) {
        result(path[i]) = Vector3(0, 0, 1.f); //colorPalette(float(i) / float(path.size() - 1));
    }
    for (size_t i = 0; i < isoline.size(); i++) {
        result(isoline[i]) = Vector3(1, 1, 1); //colorPalette(float(i) / float(path.size() - 1));
    }
    Plotter::get("Object Preview")->addImage(result);
    Plotter::get("Object Preview")->show();
    Plotter::get("Object Preview")->addImage(dataV3);
}

void EnvObjsInterface::previewFocusAreaEdition(Vector3 mousePos, bool addingFocus)
{
    //        float velocity = (prevPos - mousePos).norm(); // Typically between 0.1 to 1.0
    auto brush = GridF::normalizedGaussian(30, 30, 1, 8.f) * (addingFocus ? 1.f : -1.f) * 8.f;
    this->focusedArea.add(brush, mousePos - brush.getDimensions().xy() * .5f);

    focusedArea.iterateParallel([&](size_t i) {
        focusedArea[i] = std::clamp(focusedArea[i], 0.f, 30.f);
    });
    Plotter::get("Focus")->addImage(renderFocusArea());
    Plotter::get("Focus")->show();
}

void EnvObjsInterface::previewFlowEdition(Vector3 mousePos, Vector3 brushDir)
{

    // float velocity = (prevPos - mousePos).norm(); // Typically between 0.1 to 1.0
    GridV3 brush = GridV3(30, 30, 1, brushDir);
    brush.iterateParallel([&](const Vector3& p) {
        brush(p) *= normalizedGaussian(Vector3(30, 30, 1), p, 8.f);
    });
    // EnvObject::initialFlowfield.add(brush, (mousePos / 3.f) - brush.getDimensions().xy() * .5f);
    this->userFlowField.add(brush, (mousePos / 3.f) - brush.getDimensions().xy() * .5f);
    EnvObject::updateFlowfield(userFlowField + simulationFlowField + EnvObject::scenario.computeStorm(userFlowField.getDimensions()));
    this->updateVectorFieldVisu();

    this->addObjectsHeightmaps();
    this->flowErosionSimulation();

    Plotter::get("Flowfield")->addImage(renderFlowfield());
    Plotter::get("Flowfield")->show();
}

void EnvObjsInterface::showAllElementsOnPlotter()
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = false;

    std::map<TerrainTypes, Vector3> materialToColor = {
        {WATER, Vector3(0, 0, 1)},
        {AIR,   Vector3(0.4, 0.4, 1)},
        {SAND,  Vector3(0, 1, 1)},
        {CORAL, Vector3(0.5, 1.0, 1.0)},
        {ROCK,  Vector3(0.8, 0.8, 0.8)},
        {DIRT,  Vector3(0.7, 0.2, 0.2)}
    };
    GridV3 img(100, 100, 1);

    for (auto& obj : EnvObject::instantiatedObjects) {
        TerrainTypes material = obj->material;
        Vector3 col = materialToColor[material];
        if (auto asPoint = dynamic_cast<EnvPoint*>(obj)) {
            Vector3 pos = asPoint->position;
            for (int dx = -1; dx <= 2; dx++) {
                for (int dy = -1; dy <= 2; dy++) {
                    img(pos + Vector3(dx, dy)) = col;
                }
            }
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(obj)) {
            for (const auto& p : asCurve->curve.getPath(200)) {
                img(p) = col;
            }
        } else if (auto asArea = dynamic_cast<EnvArea*>(obj)) {
            for (const auto& p : asArea->curve.getPath(200)) {
                img(p) = col;
            }
        }
    }

    Plotter::get("Topography")->addImage(img);
    Plotter::get("Topography")->show();
}

void EnvObjsInterface::addObjectsHeightmaps()
{
    GridF subsidenceFactor = EnvObject::scenario.computeSubsidence(initialHeightmap.getDimensions());
    subsidedHeightmap = initialHeightmap * subsidenceFactor;

    float absoluteWaterLevel = voxelGrid->getSizeZ() * voxelGrid->properties->waterLevel;

    GridF groundConstraintedHeights(subsidedHeightmap.getDimensions()); // Heightmaps from the ground
    GridF waterConstraintedHeights(subsidedHeightmap.getDimensions(), -100000.f); // Heightmaps from the water level
    GridF surfaceHeights(subsidedHeightmap.getDimensions());
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(obj->_patch)) {
            GridF grid = GridF(subsidedHeightmap.getDimensions(), 0.f);
            grid = grid.paste(obj->createHeightfield() * obj->computeGrowingState2(), patch->position.xy());
            if (flowErosionFactor != 0 && EnvObject::materials.count(toLower(stringFromMaterial(obj->material)))) {
                grid = grid.warpWith(EnvObject::flowfield * flowErosionFactor * EnvObject::materials[toLower(stringFromMaterial(obj->material))].waterTransport, 10);
            }
            if (obj->heightFrom == EnvObject::HeightmapFrom::SURFACE) {
                surfaceHeights = (surfaceHeights + grid * (isIn(obj->material, LayerBasedGrid::invisibleLayers) ? -1.f : 1.f)).max(-15.f);
            } else if (obj->heightFrom == EnvObject::HeightmapFrom::GROUND) {
                groundConstraintedHeights = groundConstraintedHeights.max(grid * subsidenceFactor, Vector3());
            } else if (obj->heightFrom == EnvObject::HeightmapFrom::WATER) {
                grid.iterateParallel([&] (size_t i) {
                    grid[i] = (std::abs(grid[i]) < 1e-4 ? -10000.f : grid[i]);
                });
                // std::cout << "Max height for " << obj->name << ": " << grid.max() << " while height = " << obj->height << "(grow = " <<  obj->computeGrowingState2() << ")" << std::endl;
                waterConstraintedHeights = waterConstraintedHeights.max((grid - (obj->height)) - (obj->name == "lagoon" || obj->name == "smalllagoon" ? 3.f : 1.f), Vector3()); // Not sure why I need to multiply by 2.0, but otherwise, maxHeight is heigher than obj->height...
            }
        }
    }
    // if (flowErosionFactor != 0) {
        // groundConstraintedHeights = groundConstraintedHeights.warpWith(EnvObject::flowfield * flowErosionFactor, 10);
        // waterConstraintedHeights = waterConstraintedHeights.warpWith(EnvObject::flowfield * flowErosionFactor, 10);
        // surfaceHeights = surfaceHeights.warpWith(EnvObject::flowfield * flowErosionFactor, 10);
    // }
    // Dirty, remove when you understand why lagoon get over the water...
    waterConstraintedHeights.iterateParallel([&](size_t i) {
        waterConstraintedHeights[i] = std::min(waterConstraintedHeights[i] + absoluteWaterLevel, absoluteWaterLevel - 1.f);
    });
    waterConstraintedHeights = waterConstraintedHeights.gaussianSmooth(1.f, true);

    bool modificationsAppliedToSurface = false;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (displayGrooves) {
            if (endsWith(toLower(obj->name), "reef")) {
                auto objAsEnvCurve = dynamic_cast<EnvCurve*>(obj);
                BSpline path = objAsEnvCurve->curve;
                float nbGrooves = path.length() / 10.f;
                float sigma = objAsEnvCurve->width;
                surfaceHeights.iterateParallel([&](const Vector3& pos) {
                    float closestT = path.estimateClosestTime(pos);
                    float closestGrooveStartT = float(int(closestT * nbGrooves)) / nbGrooves;
                    auto [closestPoint, direction, normal] = path.pointAndDerivativeAndSecondDerivative(closestT);
                    auto closestGrooveStartPoint = path.getPoint(closestGrooveStartT);
                    if (direction.norm2() == 0) return;
                    direction.normalize();
                    auto fakeNormal = direction.rotated90XY(); // (normal.norm2() > 0 ? normal.normalize() : direction.rotated90XY());
                    Vector3 newSpace = Vector3(pos - closestGrooveStartPoint).changeBasis(direction, fakeNormal, Vector3(0, 0, 1)); //.rotated(Vector3(0, 0, random_gen::generate_perlin(closestT * 500.f) * 0.2f));
                    float sizeX = 1.f/(nbGrooves * .5f), sizeY = 1.f/(sigma * 1.f);
                    // float initialDistance = std::clamp(1.f - (pos - closestPoint).norm() / sigma, 0.f, 1.f);
                    float grooves = std::max(0.f, 1.f - (sizeX * std::abs(newSpace.x - 1.f/sizeX) + std::pow(sizeY * newSpace.y, 2.f)));
                    // return std::max(grooves, initialDistance);
                    const Vector3& flow = EnvObject::flowfield(pos);
                    surfaceHeights(pos) += 2.f * grooves * std::max(abs(flow.dot(fakeNormal)), 0.f);
                });
                modificationsAppliedToSurface = true;
            }
        }
    }

    if (modificationsAppliedToSurface) {
        surfaceHeights = surfaceHeights.gaussianSmooth(1.f, true, true);
    }

    subsidedHeightmap = GridF::max(GridF::max(subsidedHeightmap, groundConstraintedHeights), waterConstraintedHeights).gaussianSmooth(2.f, true, true);
    subsidedHeightmap = (subsidedHeightmap.max(-15.f) + surfaceHeights).gaussianSmooth(1.f, true, true).max(-15.f);
    /*
    std::map<std::string, GridF> perTypeHeightmap;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(obj->_patch)) {
            if (perTypeHeightmap.count(obj->name) == 0) {
                perTypeHeightmap[obj->name] = GridF(subsidedHeightmap.getDimensions());
            }
            auto grid = obj->createHeightfield() * obj->computeGrowingState2();
            perTypeHeightmap[obj->name].max(grid, patch->position.xy());
        }
    }
    std::vector<std::pair<std::string, float>> supperposedGridsOrder;
    for (auto& [name, grid] : perTypeHeightmap) {
        supperposedGridsOrder.push_back({name, LayerBasedGrid::densityFromMaterial(EnvObject::availableObjects[name]->material)});
    }
    std::sort(supperposedGridsOrder.begin(), supperposedGridsOrder.end(), [&](const std::pair<std::string, float>& A, const std::pair<std::string, float>& B) {
        return A.second > B.second;
    });
    // for (auto& [name, grid] : perTypeHeightmap) {
    for (int i = 0; i < supperposedGridsOrder.size(); i++) {
        auto& name = supperposedGridsOrder[i].first;
        auto& density = supperposedGridsOrder[i].second;
        auto& grid = perTypeHeightmap[name];
        // std::cout << "For " << name << ", max is " << grid.max() << std::endl;
        subsidedHeightmap += grid * (isIn(EnvObject::availableObjects[name]->material, LayerBasedGrid::invisibleLayers) ? -1.f : 1.f);
    }
    */
}

void EnvObjsInterface::flowErosionSimulation()
{
    // this->heightmap->fromVoxelGrid(*this->voxelGrid);
    if (flowErosionFactor != 0) {
        // this->subsidedHeightmap = (this->initialHeightmap * EnvObject::scenario.computeSubsidence()).warpWith(EnvObject::flowfield * flowErosionFactor, 10);
        //this->subsidedHeightmap = (this->subsidedHeightmap).warpWith(EnvObject::flowfield * flowErosionFactor, 10);
    }
    this->heightmap->heights = subsidedHeightmap.gaussianSmooth(.5f, true) - 1.f;
}

void EnvObjsInterface::startNewObjectCreation()
{
    this->objectSkeletonCreation = BSpline();
    this->updateNewObjectMesh();
}

void EnvObjsInterface::addPointOnNewObjectCreation(const Vector3 &position, bool addPoint, float removeRadius)
{
    if (addPoint) {
        auto objectModel = EnvObject::availableObjects[getCurrentObjectName()];
        this->objectSkeletonCreation.points.push_back(position.xy());
        if (auto asPoint = dynamic_cast<EnvPoint*>(objectModel)) {
            this->endNewObjectCreation();
        } else {
            // Nothing to do
        }
    } else {
        for (int i = objectSkeletonCreation.size() - 1; i >= 0; i--) {
            if ((objectSkeletonCreation[i] - position).xy().norm2() < removeRadius * removeRadius) {
                objectSkeletonCreation.points.erase(objectSkeletonCreation.begin() + i);
            }
        }
    }
    this->updateNewObjectMesh();
}

void EnvObjsInterface::endNewObjectCreation()
{
    // Creation of an object
    auto objectModel = EnvObject::availableObjects[getCurrentObjectName()];

    if ((dynamic_cast<EnvCurve*>(objectModel) || dynamic_cast<EnvArea*>(objectModel)) && this->objectSkeletonCreation.empty())
        return;

    auto newObject = this->instantiateObjectUsingSpline(objectModel->name, this->objectSkeletonCreation);
    this->objectSkeletonCreation = BSpline();
    if (newObject) {
        newObject->createdManually = true;
    }
    if (!newObject) {
        this->log("Object not created");
        return;
    }
    newObject->recomputeEvaluationPoints();
    newObject->fittingFunction = EnvObject::parseFittingFunction(newObject->s_FittingFunction, newObject->name, true, newObject);
    newObject->fitnessFunction = EnvObject::parseFittingFunction(newObject->s_FitnessFunction, newObject->name, true, newObject);
    auto implicit = newObject->createImplicitPatch(subsidedHeightmap);
    this->implicitPatchesFromObjects[newObject] = implicit;
    if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
        this->implicitTerrain->addChild(this->rootPatch);

    if (implicit != nullptr) {
        rootPatch->addChild(implicit);
    }
    // Wait until the object is 100% grown:
    int maxIterations = 100;
    while (newObject->computeGrowingState() < 1.f) {
        this->updateEnvironmentFromEnvObjects(false, true);
        maxIterations--;
        if (maxIterations < 0) break;
    }
    this->currentSelections = {newObject};
    EnvObject::recomputeTerrainPropertiesForObject(newObject->name);
    this->updateEnvironmentFromEnvObjects(implicit != nullptr); // If implicit is null, don't update the map

    this->updateNewObjectMesh();
    this->updateSelectionMesh();

    Q_EMIT this->updated();
}

void EnvObjsInterface::startDraggingObject(const Vector3 &position, bool singleVertexMoved)
{
    if (singleVertexMoved) {
        draggingPoint = position.xy();
    } else {
        draggingFullObject = position.xy();
    }
    draggingHasBeenApplied = position.xy();
    draggingHasBeenApplied.setValid(false); // Keep position, but set to invalid
}

void EnvObjsInterface::moveDraggedObject(const Vector3 &position)
{
    if (draggingPoint.isValid()) {
        draggingHasBeenApplied.setValid(true);
        Vector3 translation = (position.xy() - draggingHasBeenApplied.xy());
        draggingHasBeenApplied = position.xy();
        float maxDistToPointSqr = 20.f * 20.f;

        for (auto currentSelection : currentSelections) {
            if (EnvPoint* point = dynamic_cast<EnvPoint*>(currentSelection)) {
                point->translate(translation);
            } else if (EnvCurve* curve = dynamic_cast<EnvCurve*>(currentSelection)) {
                auto newCurve = curve->curve;
                int pointIndexToMove = -1;
                float closestDistToPoint = std::numeric_limits<float>::max();

                for (int i = 0; i < curve->curve.size(); i++) {
                    float dist = (curve->curve[i] - position.xy()).norm2();
                    if (dist < maxDistToPointSqr && dist < closestDistToPoint) {
                        closestDistToPoint = dist;
                        pointIndexToMove = i;
                    }
                }

                if (pointIndexToMove > -1) {
                    newCurve[pointIndexToMove].translate(translation);
                }
                curve->updateCurve(newCurve);
            } else if (EnvArea* area = dynamic_cast<EnvArea*>(currentSelection)) {
                auto newCurve = area->curve;
                int pointIndexToMove = -1;
                float closestDistToPoint = std::numeric_limits<float>::max();

                for (int i = 0; i < area->curve.size(); i++) {
                    float dist = (area->curve[i] - position.xy()).norm2();
                    if (dist < maxDistToPointSqr && dist < closestDistToPoint) {
                        closestDistToPoint = dist;
                        pointIndexToMove = i;
                    }
                }

                if (pointIndexToMove > -1) {
                    newCurve[pointIndexToMove].translate(translation);
                }
                area->updateCurve(newCurve);
            }
        }
        // Also do it on the creation curve :

        auto newCurve = objectSkeletonCreation;
        int pointIndexToMove = -1;
        float closestDistToPoint = std::numeric_limits<float>::max();

        for (int i = 0; i < objectSkeletonCreation.size(); i++) {
            float dist = (objectSkeletonCreation[i] - position.xy()).norm2();
            if (dist < maxDistToPointSqr && dist < closestDistToPoint) {
                closestDistToPoint = dist;
                pointIndexToMove = i;
            }
        }

        if (pointIndexToMove > -1) {
            newCurve[pointIndexToMove].translate(translation);
        }
        objectSkeletonCreation = newCurve;

        this->updateNewObjectMesh();
        this->updateSelectionMesh();
    } else if (draggingFullObject.isValid()) {
        draggingHasBeenApplied.setValid(true);
        Vector3 translation = (position.xy() - draggingHasBeenApplied.xy());
        draggingHasBeenApplied = position.xy();

        for (auto currentSelection : currentSelections) {
            if (EnvPoint* point = dynamic_cast<EnvPoint*>(currentSelection)) {
                point->translate(translation);
            } else if (EnvCurve* curve = dynamic_cast<EnvCurve*>(currentSelection)) {
                curve->translate(translation);
            } else if (EnvArea* area = dynamic_cast<EnvArea*>(currentSelection)) {
                area->translate(translation);
            }
        }

        // Also do it for the creation curve
        objectSkeletonCreation.translate(translation);
    }
    this->updateNewObjectMesh();
    this->updateSelectionMesh();

    Q_EMIT this->updated();
}

void EnvObjsInterface::endDraggingObject(bool destroyObjects)
{
    if (draggingFullObject.isValid() && destroyObjects) {
        for (auto currentSelection : currentSelections) {
            this->destroyEnvObject(currentSelection);
        }
        this->updateEnvironmentFromEnvObjects(true, true);
    }
    if ((draggingPoint.isValid() || draggingFullObject.isValid()) && draggingHasBeenApplied.isValid()) {
        draggingPoint.setValid(false);
        draggingFullObject.setValid(false);
        if (!destroyObjects) {
            for (auto currentSelection : currentSelections) {
                currentSelection->age = 0.f;
                if (this->implicitPatchesFromObjects.count(currentSelection) != 0) {
                    currentSelection->geometryNeedsUpdate = true;
                    auto newPatch = currentSelection->createImplicitPatch(subsidedHeightmap * 0.f/*, dynamic_cast<ImplicitPrimitive*>(currentSelection->_patch)*/);
                    if (newPatch) {
                        *(this->implicitPatchesFromObjects[currentSelection]) = *newPatch;
                        // delete newPatch;
                    }
                }
                EnvObject::recomputeTerrainPropertiesForObject(currentSelection->name);
            }
        }
        this->materialSimulationStable = false;
        this->updateEnvironmentFromEnvObjects(true, true);

    }
    this->updateSelectionMesh();
    draggingPoint.setValid(false);
    draggingFullObject.setValid(false);
    draggingHasBeenApplied.setValid(false);

    Q_EMIT this->updated();
}

std::string EnvObjsInterface::getCurrentObjectName() const
{
    return objectCombobox->choices[objectCombobox->combobox()->currentIndex()].label;
}

void EnvObjsInterface::updateVectorFieldVisu()
{
    GridV3 velocities = EnvObject::flowfield.resize(Vector3(50, 50, 1));
    Mesh::createVectorField(velocities, this->voxelGrid->getDimensions(), &velocitiesMesh, 1.f, false, true);
}

StatsValues EnvObjsInterface::displayStatsForObjectCreation(std::string objectName, int nbSamples)
{
    std::vector<float> values(nbSamples);
    bool isPossible;
    GridF score = computeScoreMap(objectName, subsidedHeightmap.getDimensions(), isPossible) * focusedArea;
    if (isPossible) {
        // #pragma omp parallel for
        for (int i = 0; i < nbSamples; i++) {
            auto obj = fakeInstantiate(objectName, score);
            if (obj) {
                values[i] = obj->fitnessScoreAtCreation;
            }
        }
    }
    StatsValues stats = getStats(values);
    std::cout << objectName << ": " << stats.mean << " (" << stats.stdev << ")" << std::endl;
    return stats;

}

void EnvObjsInterface::fromGanUI()
{
    this->resetScene();
    EnvObject::readEnvObjectsFileContent(this->primitiveDefinitionFile.read());
    EnvObject::readEnvMaterialsFileContent(this->materialsDefinitionFile.read());

    std::string path = "Python_tests/test_island_heightmapfeatures/";
//    QString q_filename= QString::fromStdString(path + "_1.png");
    QString q_filename= QString::fromStdString(path + "1.png");  //QFileDialog::getOpenFileName(this, "Open feature map", QString::fromStdString(path), "*", nullptr);
    if (!q_filename.isEmpty()) {
        std::string file = q_filename.toStdString();
        GridV3 img = Image::readFromFile(file).colorImage;

        auto envObjects = CoralIslandGenerator::envObjsFromFeatureMap(img, voxelGrid->getDimensions());
        rootPatch->deleteAllChildren();
        for (auto& newObject : envObjects) {
            auto implicit = newObject->createImplicitPatch(subsidedHeightmap);
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (implicit != nullptr) {
                rootPatch->addChild(implicit);
            }
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);
        }

        implicitTerrain->updateCache();
        implicitTerrain->update();
        rootPatch->reevaluateAll();

        /*
        std::cout << "To voxels: " << showTime(timeIt([&]() {
            voxelGrid->fromImplicit(implicitTerrain.get(), 40);
        })) << std::endl;
        std::cout << "To heightmap: " << showTime(timeIt([&]() {
            heightmap->fromVoxelGrid(*voxelGrid.get());
        })) << std::endl;*/

        EnvObject::precomputeTerrainProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
        this->updateEnvironmentFromEnvObjects(true, true);
        displayProcessTime("Update object list... ", [&]() {
            updateObjectsList();
        });
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
