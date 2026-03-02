#include "EnvObjsInterface.h"

#include "GUIElements/InterfaceUtils.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "GUIElements/CommonInterface.h"
#include "TerrainModification/CoralIslandGenerator.h"
#include "DataStructure/Image.h"
#include "EnvObject/ExpressionParser.h"
#include "Utils/Voronoi.h"
#include "Interface/MeshInstanceAmplificationInterface.h"
#include "Interface/TerrainGenerationInterface.h"
#include "Utils/Delaunay.h"
#include "Utils/PSO.h"


#include "EnvObjGUI/EnvMaterialViewer.h"
#include "EnvObjGUI/FocusAreaViewer.h"
#include "EnvObjGUI/WaterFlowViewer.h"



EnvObjsInterface::EnvObjsInterface(QWidget *parent)
    : ActionInterface("envobjects", "Environmental Objects", "model", "Management of environmental objects generation", "envobjs_button.png", parent)
{
    this->scene = std::make_shared<EnvironmentalScene>();

    this->scene->readEnvMaterialsFile("EnvObjects/envMaterials.json");
    this->scene->readEnvObjectsFile("EnvObjects/primitives.json");

    primitiveDefinitionFile.onChange([&](const std::string& newDefinitions) { updateObjectsDefinitions(newDefinitions); });
    materialsDefinitionFile.onChange([&](const std::string& newDefinitions) { updateMaterialsDefinitions(newDefinitions); });
    transformationsFile.onChange([&](const std::string& newDefinitions) { updateMaterialsTransformationsDefinitions(newDefinitions); });
    scenarioFile.onChange([&](const std::string& newDefinitions) { updateScenarioDefinition(newDefinitions); });

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
    this->focusedArea = GridF(initialHeightmap.getDimensions(), .5f);
    this->userFlowField = GridV3(initialHeightmap.getDimensions());
    this->simulationFlowField = GridV3(initialHeightmap.getDimensions());
    this->scene->precomputeTerrainProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());

    QObject::connect(ImageViewer::get("Object Preview"), &ImageViewer::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* _event) {
        this->previewCurrentEnvObjectPlacement(clickPos);
    });

    /*QObject::connect(EnvMaterialViewer::get("Material"), &EnvMaterialViewer::movedOnImage, this, [&](const Vector3& clickPos, const Vector3& _prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;
        this->previewMaterialEdition(clickPos, leftPressed);
    });*/
    QObject::connect(EnvMaterialViewer::get("Material"), &EnvMaterialViewer::imagePainted, this, [&](const GridF& newDistrib) {
        this->scene->materials[this->currentMaterialEdited].currentState = newDistrib;
    });

    /*QObject::connect(FocusAreaViewer::get("Focus"), &ImageViewer::movedOnImage, this, [&](const Vector3& mousePos, const Vector3& prevPos, QMouseEvent* event) {
        bool leftPressed = event->buttons().testFlag(Qt::LeftButton);
        bool rightPressed = event->buttons().testFlag(Qt::RightButton);
        if (!leftPressed && !rightPressed) return;

        this->previewFocusAreaEdition(mousePos, leftPressed);
    });*/
    QObject::connect(FocusAreaViewer::get("Focus"), &FocusAreaViewer::imagePainted, this, [&](const GridF& newDistrib) {
        this->focusedArea = newDistrib;
    });

    QObject::connect(WaterFlowViewer::get("Flowfield"), &WaterFlowViewer::updated, this, [&]() {
        this->userKelvinlets = WaterFlowViewer::get("Flowfield")->kelvinletParams.kelvinlets;
        this->previewFlowEdition(Vector3::invalid, Vector3::invalid);
        Q_EMIT this->updated();
    });
}

void EnvObjsInterface::display(const Vector3 &camPos)
{
    if (!this->visible)
        return;

    if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
        this->implicitTerrain->addChild(this->rootPatch);

    if (this->waitAtEachFrame) {

        for (auto& obj : this->scene->instantiatedObjects) {
            obj->improvePositionning(1.f);
        }
        this->updateUntilStabilization();

        this->updateEnvironmentFromEnvObjects(true, false);
    }

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
        for (auto& obj : this->scene->instantiatedObjects) {
            obj->improvePositionning(1.f);
        }
        this->updateUntilStabilization();
        this->updateSelectionMesh();

        Q_EMIT this->updated();
        //this->saveScene("testEnvObjects.json");
    });
    CheckboxElement* waitAtEachFrameButton = new CheckboxElement("Auto wait", this->waitAtEachFrame);
//    ButtonElement* createFromGAN = new ButtonElement("From GAN", [&]() { this->fromGanUI(); });
    ButtonElement* createFromFile = new ButtonElement("From file", [&]() { this->loadScene("EnvObjects/testEnvObjects.json"); });
    TextEditElement* testingFittingFormula = new TextEditElement("", "Fitting func: ");
    // testingFittingFormula->setOnTextChange([&](const std::string& expression) { this->evaluateAndDisplayCustomFittingFormula(expression); });
    TextEditElement* testingFitnessFormula = new TextEditElement("", "Fitness func: ");
    // testingFitnessFormula->setOnTextChange([&](const std::string& expression) { this->evaluateAndDisplayCustomFitnessFormula(expression); });
    testingFitnessFormula->setOnTextChange([&](const std::string& expression) { this->testedFitnessFunction = expression; this->evaluateAndDisplayCustomFitnessAndFittingFormula(this->testedFitnessFunction, this->testedFittingFunction); });
    testingFittingFormula->setOnTextChange([&](const std::string& expression) { this->testedFittingFunction = expression; this->evaluateAndDisplayCustomFitnessAndFittingFormula(this->testedFitnessFunction, this->testedFittingFunction); });
    // ButtonElement* testPerformancesButton = new ButtonElement("Run test", [&]() { this->runPerformanceTest(); });
    ButtonElement* resetButton = new ButtonElement("Reset scene", [&]() { this->resetScene(); });
    CheckboxElement* addGroovesButton = new CheckboxElement("Spurs and grooves", displayGrooves);
    ButtonElement* saveForRendersButton = new ButtonElement("Save for render", [&]() { this->saveForRenders(); });

    LabelElement* label = new LabelElement("Objects: " + std::to_string(this->scene->instantiatedObjects.size()));

    objectsListWidget = new HierarchicalListUI;
    objectsListWidget->setSelectionMode(QAbstractItemView::SelectionMode::ExtendedSelection);
    updateObjectsList();
    objectsListWidget->setOnItemSelectionChanged([&]() { this->updateObjectsListSelection(); });


    std::vector<ComboboxLineElement<EnvObject*>*> objectsChoices;
    int selectionForCoral = 0;
    for (auto& [name, obj] : this->scene->availableObjects) {
        objectsChoices.push_back(new ComboboxLineElement<EnvObject*>(name, obj));
        if (toLower(name) == "coralpolyp") {
            selectionForCoral = objectsChoices.size() - 1;
        }
    }
    ButtonElement* showButton = new ButtonElement("Show", [&](){
        this->displayProbas(getCurrentObjectName());
    });
    ButtonElement* forceButton = new ButtonElement("Force", [&](){
        this->instantiateSpecific(getCurrentObjectName(), Vector3::invalid, GridF(), false, false); //true, true);
        updateObjectsList();
        Q_EMIT this->updated();
    });
    forceButton->setOnRepeat([&](){
        this->instantiateSpecific(getCurrentObjectName(), Vector3::invalid, GridF(), false, false); //true, true);
        updateObjectsList();
        Q_EMIT this->updated();
    });
    objectCombobox = new ComboboxElement("Objects", objectsChoices);
    objectCombobox->combobox()->setCurrentIndex(selectionForCoral);

    std::vector<UIElement*> materialsButtons;
    for (auto& [name, material] : this->scene->materials) {
        ButtonElement* showButton = new ButtonElement("Show " + toCapitalize(name), [&](){ this->displayMaterialDistrib(name); });
        materialsButtons.push_back(showButton);
    }

    ButtonElement* editFocusAreaButton = new ButtonElement("Edit focus", [&]() { this->manualModificationOfFocusArea(); });
    ButtonElement* editFlowfieldButton = new ButtonElement("Edit flowfield", [&]() { this->manualModificationOfFlowfield(); });
    ButtonElement* resetFlowfieldButton = new ButtonElement("Reset flow", [&]() { this->resetFlowfield(); });
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

    RadioButtonElement* grabKelvinlet = new RadioButtonElement("Grab");
    RadioButtonElement* scaleKelvinlet = new RadioButtonElement("Scale");
    RadioButtonElement* pinchKelvinlet = new RadioButtonElement("Pinch");
    RadioButtonElement* twistKelvinlet = new RadioButtonElement("Twist");
    grabKelvinlet->setOnChecked([&](bool checked) { if(checked) { this->KelvinletChoice = "grab"; } });
    scaleKelvinlet->setOnChecked([&](bool checked) { if(checked) { this->KelvinletChoice = "scale"; } });
    pinchKelvinlet->setOnChecked([&](bool checked) { if(checked) { this->KelvinletChoice = "pinch"; } });
    twistKelvinlet->setOnChecked([&](bool checked) { if(checked) { this->KelvinletChoice = "twist"; } });

    CheckboxElement* displayCurrentsButton = new CheckboxElement("Flow", this->displayFlow);

    ui->add({
             createHorizontalGroupUI({newObjectCreationBox, waitAtEachFrameButton}),
             createHorizontalGroupUI({spendTimeButton, nextStepButton, runButton}),
             createMultiColumnGroupUI(materialsButtons, 2),
             flowErosionSlider,
             objectCombobox,
             createMultiColumnGroupUI({showButton, forceButton}, 2),
             createHorizontalGroupUI({editFocusAreaButton, editFlowfieldButton, resetFlowfieldButton}),
             showElementsOnCanvasButton,
             objectsListWidget,
             createVerticalGroupUI({testingFitnessFormula, testingFittingFormula}),
             createHorizontalGroupUI({/*instantiaABCbutton, testPerformancesButton, */resetButton}),
             addGroovesButton,
             createVerticalGroupUI({grabKelvinlet, scaleKelvinlet, pinchKelvinlet, twistKelvinlet}),
             createHorizontalGroupUI({label, createFromFile, saveButton, displayCurrentsButton}),
             saveForRendersButton
    });

    return ui->get()->layout();
}

void EnvObjsInterface::createEnvObjectsFromImplicitTerrain()
{
    this->scene->removeAllObjects();
    auto tunnelsPatches = implicitTerrain->findAll(ImplicitPatch::ParametricTunnel);
    for (auto& tunnelPatch : tunnelsPatches) {
        auto asPrimitive = dynamic_cast<ImplicitPrimitive*>(tunnelPatch);
        if (asPrimitive && asPrimitive->material == WATER) {
            BSpline curve = asPrimitive->optionalCurve;
            for (auto& p : curve) {
                p = asPrimitive->getGlobalPositionOf(p);
            }
            EnvCurveInstance* passe = dynamic_cast<EnvCurveInstance*>(this->scene->instantiate("passe"));
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
            EnvAreaInstance* lagoon = dynamic_cast<EnvAreaInstance*>(this->scene->instantiate("lagoon"));
            lagoon->curve = curve;
        }
    }

    if (objectsListWidget != nullptr)
        updateObjectsList();
}

void EnvObjsInterface::setMaterialsDefinitionFile(const std::string& filename)
{
    this->materialsDefinitionFile.path = filename;
    this->scene->readEnvMaterialsFile(filename);
}

void EnvObjsInterface::setDefinitionFile(const std::string& filename)
{
    this->primitiveDefinitionFile.path = filename;
    this->scene->readEnvObjectsFile(filename);
}

void EnvObjsInterface::setTransformationsFile(const std::string& filename)
{
    this->transformationsFile.path = filename;
    this->scene->readEnvMaterialsTransformationsFile(filename);
}

void EnvObjsInterface::setScenarioFile(const std::string& filename)
{
    this->scenarioFile.path = filename;
    this->scene->readScenarioFile(filename);
}

void EnvObjsInterface::show()
{
    ActionInterface::show();
}

void EnvObjsInterface::hide()
{
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
    this->scene->recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());

    if (this->isVisible()) {
        if (ImageViewer::get("Object Preview")->isVisible()) {
            displayProbas(getCurrentObjectName());
        }
    }
}

void EnvObjsInterface::mouseClickedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    if (!this->isVisible()) return;
    if (!mouseInMap) return;

    if (event->button() == Qt::MouseButton::LeftButton) {

        bool moveSingleVertex = event->modifiers().testFlag(Qt::KeyboardModifier::ShiftModifier);
        bool moveWholeObject = event->modifiers().testFlag(Qt::KeyboardModifier::ControlModifier);

        if (moveSingleVertex || moveWholeObject) {
            this->startDraggingObject(mouseWorldPosition, moveSingleVertex);
        }
        else if (this->manuallyCreatingObject) {
            bool addingPoint = event->buttons().testFlag(Qt::LeftButton);
            this->addPointOnNewObjectCreation(mouseWorldPosition, addingPoint);
        }
    } else if (event->button() == Qt::MouseButton::RightButton) {
        if (event->modifiers().testFlag(Qt::KeyboardModifier::ShiftModifier)) {
            kelvinletDraggingPoint = mouseWorldPosition.xy();

            KelvinletPoint* _k;
            if (this->KelvinletChoice == "grab") {
                _k = new GrabKelvinlet();
            } else if (this->KelvinletChoice == "scale") {
                _k = new ScaleKelvinlet();
                _k->v = 0.1; // Don't keep the default value of 0.5 here!
            } else if (this->KelvinletChoice == "pinch") {
                _k = new PinchKelvinlet();
            } else if (this->KelvinletChoice == "twist") {
                _k = new TwistKelvinlet();
            }

            if (_k) {
                _k->pos = mouseWorldPosition.xy();
                this->userKelvinlets.push_back(_k);
            }
        }
    }
}

void EnvObjsInterface::mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model)
{
    if (!this->isVisible()) return;
    if (!mouseWorldPosition.isValid()) return;

    this->moveDraggedObject(mouseWorldPosition);


    if (this->kelvinletDraggingPoint.isValid()) {
        KelvinletPoint* _k = dynamic_cast<KelvinletPoint*>(this->userKelvinlets.back());
        Vector3 delta = mouseWorldPosition.xy() - _k->pos;

        if (_k != nullptr) {
            if (this->KelvinletChoice == "grab") {
                GrabKelvinlet* k = dynamic_cast<GrabKelvinlet*>(_k);
                k->force = delta;
            } else if (this->KelvinletChoice == "scale") {
                ScaleKelvinlet* k = dynamic_cast<ScaleKelvinlet*>(_k);
                // k->radialScale = 5.f;
                k->force = delta.norm();
            } else if (this->KelvinletChoice == "pinch") {
                PinchKelvinlet* k = dynamic_cast<PinchKelvinlet*>(_k);
                k->force = delta;
            } else if (this->KelvinletChoice == "twist") {
                TwistKelvinlet* k = dynamic_cast<TwistKelvinlet*>(_k);
                k->force = Vector3(0, 0, delta.norm() * sign(delta.x()));
                // k->radialScale = delta.norm();
            }

            this->scene->updateFlowfield(userFlowField + this->computeUserKelvinletField(), simulationFlowField);
            this->updateVectorFieldVisu();
            // Q_EMIT this->updated();
        }
    }
}

void EnvObjsInterface::mouseReleasedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model)
{
    if (!this->isVisible()) return;

    bool destroyObjects = !mouseInMap;

    this->endDraggingObject(destroyObjects);

    if (this->kelvinletDraggingPoint.isValid()) {
        this->kelvinletDraggingPoint = Vector3::invalid;

        this->updateVectorFieldVisu();
        Q_EMIT this->updated();
    }

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

GridF computeScoreMap(std::shared_ptr<EnvironmentalScene> scene, std::string objectName, const Vector3& dimensions, bool& possible, bool applyNormalization = false) {
    auto obj = scene->availableObjects[objectName]->instantiate();
    scene->recomputeFlow();
    GridF score = GridF(dimensions);
    score.iterateParallel([&](const Vector3i& pos) {
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


Vector3 bestPositionForInstantiationUniform(std::shared_ptr<EnvironmentalScene> scene, std::string objectName, const AABBox& bounds, const GridF& focusArea, const Vector3& focusedPosition = Vector3::invalid, int nbSamples = 20) {
    auto& func = scene->availableObjects[objectName]->fitnessFunction;
    std::vector<std::pair<Vector3, float>> evaluations(nbSamples);
    float minScoreThreshold = scene->availableObjects[objectName]->minScore;
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < nbSamples; i++) {
        float focus = 0.f;
        Vector3 pos;
        float score = 0.f;
        while (focus < random_gen::generate()) {
            pos = Vector3::random(bounds);
            focus = focusArea.at(pos) * (focusedPosition.isValid() ? 1.f / std::max(0.1f, std::log((pos - focusedPosition).norm2())) : 1.f);
        }
        score = func(pos);
        evaluations[i] = std::make_pair(pos, score > minScoreThreshold ? score : 0.f);
    }
    float totalScore = 0.f;
    for (int i = 0; i < nbSamples; i++) {
        totalScore += evaluations[i].second;
    }
    float cummul = 0.f;
    float target = random_gen::generate(0.f, totalScore);

    int iSample = 0;
    for (iSample = 0; iSample < nbSamples; iSample++) {
        cummul += evaluations[iSample].second;
        if (cummul > target) return evaluations[iSample].first;
    }
    return Vector3::invalid;
}

EnvObjectInstance* EnvObjsInterface::instantiateObjectAtBestPositionWithoutScoreMap(const std::string& objectName, Vector3 position, const Vector3& maxPos)
{
    bool verbose = false;
    Vector3 initialPosition = position;
    EnvObjectInstance* newObject = this->scene->instantiate(objectName);
    // std::cout << "Placing '" << newObject->name << "' at " << initialPosition << std::endl;

    if (!newObject->placeInTerrain(initialPosition)) {
        this->destroyEnvObject(newObject, false, false);
        delete newObject;
        return nullptr;
    }
    int nbOutside = 0;
    for (auto& p : newObject->evaluationPositions) {
        if (!Vector3::isInBox(p, Vector3(), maxPos)) {
            nbOutside++;
        }
        if (nbOutside > newObject->evaluationPositions.size() / 2) {
            log("Object is outside...", verbose);
            this->destroyEnvObject(newObject, false, false);
            return nullptr;
        }
    }
    this->log("Creation of obj at score = " + std::to_string(newObject->fitnessScoreAtCreation), verbose);
    return newObject;
}

EnvObjectInstance* EnvObjsInterface::instantiateObjectUsingSpline(const std::string& objectName, const BSpline &spline)
{
    EnvObjectInstance* newObject = this->scene->instantiate(objectName);
    newObject->snake.position = spline.getPoint(.5f);
    newObject->placeInTerrain(spline);
    this->log("Manual creation of obj at score = " + std::to_string(newObject->fitnessScoreAtCreation));
    return newObject;
}

EnvObjectInstance* EnvObjsInterface::instantiateSpecific(const std::string& _objectName, const Vector3 &targetPosition, const GridF &score, bool waitForFullyGrown, bool updateScreen, bool updateEnvironmentDirectly)
{
    EnvObjectInstance* result = nullptr;
    std::string objectName = _objectName;
    objectName = toLower(objectName);
    if (this->scene->availableObjects.count(objectName) == 0) {
        std::cerr << "No object '" << objectName << "' in database!" << std::endl;
        return result;
    }
    bool verbose = true;
    Vector3 position = bestPositionForInstantiationUniform(this->scene, objectName, AABBox(Vector3i::origin, this->heightmap->getDimensions()), this->focusedArea, targetPosition, 1000);

    if (position.isValid()) {
        EnvObjectInstance* newObject = instantiateObjectAtBestPositionWithoutScoreMap(objectName, position, this->heightmap->getDimensions());
        if (!newObject) {
            this->log("Object not created", verbose);
            return nullptr;
        }
        ImplicitPatch* implicit = newObject->createImplicitPatch(subsidedHeightmap);
        // newObject->getDefinition()->fittingFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FittingFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
        // newObject->getDefinition()->fitnessFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FitnessFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
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
            if (!isIn(newObject, this->scene->instantiatedObjects)) {
                return nullptr; // Object died in this process, stop this function now
            }
        }
        this->currentSelections = {newObject};
        if (updateEnvironmentDirectly) {
            this->scene->recomputeTerrainPropertiesForObject(objectName);
            this->updateEnvironmentFromEnvObjects(implicit != nullptr, updateScreen); // If implicit is null, don't update the map
        }
        result = newObject;
        this->materialSimulationStable = false; // We have to compute the simulation again
    } else {
        // std::cout << "Nope, impossible to instantiate..." << std::endl;
    }

    if (updateScreen)
        updateObjectsList();
    return result;
}

EnvObjectInstance* EnvObjsInterface::fakeInstantiate(const std::string& _objectName, const GridF &score)
{
    std::string objectName = _objectName;
    objectName = toLower(objectName);
    if (this->scene->availableObjects.count(objectName) == 0) {
        std::cerr << "No object '" << objectName << "' in database!" << std::endl;
        return nullptr;
    }
    Vector3 position = bestPositionForInstantiationUniform(this->scene, objectName, AABBox(Vector3i::origin, this->heightmap->getDimensions()), this->focusedArea, Vector3::invalid, 1000);
    if (!position.isValid()) {
        return nullptr;
    }
    EnvObjectInstance* newObject = instantiateObjectAtBestPositionWithoutScoreMap(objectName, position, this->heightmap->getDimensions());
    if (!newObject) {
        return nullptr;
    }
    // newObject->getDefinition()->fittingFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FittingFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
    // newObject->getDefinition()->fitnessFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FitnessFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
    this->destroyEnvObject(newObject, false, false);
    return newObject;
}

bool EnvObjsInterface::checkIfObjectShouldDie(EnvObjectInstance* obj, float limitFactorForDying)
{
    if (obj->createdManually) return false;
    return (obj->computeGrowingState2() <= obj->getDefinition()->minScore * limitFactorForDying);
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
    EnvObjectInstance* createdObject = nullptr;
    Scenario& scenario = this->scene->scenario;

    int nbObjectsCreated = 0;
    float t1 = timeIt([&]() {
        for (auto& nextObject : scenario.nextObjects()) {
            // bool possible;
            // GridF score;
            // displayProcessTime("Score map for " + nextObject.objectName + ": ", [&]() {
            // score = computeScoreMap(this->scene, nextObject.objectName, subsidedHeightmap.getDimensions(), possible) * focusedArea;
            // });
            // if (possible) {
                // displayProcessTime("Instantiation ", [&]() {
                    // createdObject = this->instantiateSpecific(nextObject.objectName, Vector3::invalid, score, false, false, false);
                // });
            // }
            createdObject = this->instantiateSpecific(nextObject.objectName, Vector3::invalid, GridF(), false, false, false);
            if (createdObject != nullptr)
                nbObjectsCreated++;
        }
        // this->updateEnvironmentFromEnvObjects(true, true, true);
    });
    float t2 = timeIt([&](){ updateEnvironmentFromEnvObjects(false, false, true); });
    this->scene->currentTime += scenario.dt;
    this->log("Step: " + std::to_string(this->scene->scenario.currentTime()) + " -- instantiation (" + std::to_string(nbObjectsCreated) + ") : " + showTime(t1) + " -- update : " + showTime(t2));
}

void EnvObjsInterface::runScenario()
{
    updateEnvironmentFromEnvObjects(true, false, false);
    Scenario& scenario = this->scene->scenario;
    this->forceScenarioInterruption = false;
    if (scenario.waterLevel >= 0) {
        (dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get()))->setWaterLevel(scenario.waterLevel);
    } else {
        scenario.waterLevel = heightmap->properties->waterLevel;
    }
    scenario.startTime = this->scene->currentTime;
    while (!scenario.finished() && !forceScenarioInterruption) {
        float time = scenario.currentTime();
        float waterLevel = scenario.computeWaterLevel();
        (dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get()))->setWaterLevel(waterLevel);
        runNextStep();
    }
    this->forceScenarioInterruption = true;
}

void EnvObjsInterface::updateEnvironmentFromEnvObjects(bool updateImplicitTerrain, bool emitUpdateSignal, bool killObjectsIfPossible)
{
    bool verbose = true;

    GridF subsidenceFactor = this->scene->scenario.computeSubsidence(initialHeightmap.getDimensions());
    subsidedHeightmap = initialHeightmap * subsidenceFactor;

    // std::vector<EnvObject*> immatureObjects;
    if (killObjectsIfPossible) {
        /*for (auto& obj : this->scene->instantiatedObjects) {
            if (obj->computeGrowingState() < 1.f) {
                materialSimulationStable = false;
                immatureObjects.push_back(obj);
            }
        }*/

        bool atLeastOneDeath = false;
        std::set<std::string> deadObjects;
        for (auto& obj : this->scene->instantiatedObjects) {
            // if (isIn(obj, immatureObjects)) continue;
            bool shouldDie = this->checkIfObjectShouldDie(obj, .1f);
            atLeastOneDeath |= shouldDie;
            if (shouldDie) {
                float startingScore = obj->fitnessScoreAtCreation;
                // float endingScore = obj->evaluate(obj->evaluationPosition);
                float endingScore = obj->evaluate();
                this->log(obj->getDefinition()->name + " went from " + std::to_string(startingScore) + " to " + std::to_string(endingScore) + " -> " + std::to_string(std::round(100.f * endingScore / startingScore)) + "%");
                this->destroyEnvObject(obj);
                deadObjects.insert(obj->getDefinition()->name);
                if (random_gen::generate() < .9)
                    this->instantiateSpecific(obj->getDefinition()->name, obj->evaluationPositions[0], GridF(), false, false, false);
            }
        }
        if (atLeastOneDeath) {
            updateImplicitTerrain = true;
            for (auto death : deadObjects) {
                this->scene->recomputeTerrainPropertiesForObject(death);
            }
            //updateEnvironmentFromEnvObjects(false, false, false);
        }
    }

    displayProcessTime("Get impacted... ", [&]() {
        this->scene->beImpactedByEvents();
    }, verbose);

    displayProcessTime("Apply effect (stabilization)... ", [&]() {
        this->scene->stabilizeMaterials(this->heightmap->getHeights());
    });

    if (!this->materialSimulationStable) { // If the simulation is stable, don't do anything
        /*
        displayProcessTime("Apply effects... ", [&]() {
                bool bigChangesInMaterials = this->scene->applyEffects(subsidedHeightmap, userFlowField + simulationFlowField + this->computeUserKelvinletField());
            //this->materialSimulationStable = !bigChangesInMaterials;
        }, true);
        */
        displayProcessTime("Recompute properties... ", [&]() {
            this->scene->recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
        }, verbose);
        // Get original flowfield, do not accumulate effects (for now).
        displayProcessTime("Get velocity... ", [&]() {
            if (!this->fluidSimulationIsStable) {
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->mainDirection = Vector3();
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->setObstacles(voxelGrid->getVoxelValues());
                dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->recomputeVelocities();
                this->simulationFlowField = dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->getVelocities(this->scene->flowfield.sizeX, this->scene->flowfield.sizeY, this->scene->flowfield.sizeZ);
                this->simulationFlowField *= .1f;
                this->fluidSimulationIsStable = true;

                this->simulationFlowField *= 0.f; // WARNING TEMP !!! TO REMOVE FAST!!! FOR DEBUG ONLY
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
            for (auto& [name, material] : this->scene->materials) {
                this->heightmap->heights += material.currentState * material.virtualHeight;
            }
        }
    }, verbose);
    this->scene->recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());

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

    /*
    for (auto& obj : immatureObjects) {
        if (obj->computeGrowingState() >= 1.f) {
            // Got mature during this process -> now let's save the fitting score
            // obj->fitnessScoreAtCreation = obj->evaluate(obj->evaluationPosition);
            obj->fitnessScoreAtCreation = obj->evaluate();
        }
    }
    */

    if (emitUpdateSignal) {
        Q_EMIT this->updated();
        updateObjectsList();
    }
}

void EnvObjsInterface::updateUntilStabilization()
{
    displayProcessTime("Stabilisation: ", [&]() {
        auto heights = heightmap->getHeights();
        this->scene->stabilizeMaterials(heights);
    });
}

void EnvObjsInterface::destroyEnvObject(EnvObjectInstance* object, bool applyDying, bool recomputeTerrainPropertiesForObject)
{
    if (!object) return;

    if (applyDying)
        object->die();
    for (size_t i = 0; i < this->scene->instantiatedObjects.size(); i++) {
        if (this->scene->instantiatedObjects[i] == object) {
            this->scene->instantiatedObjects.erase(this->scene->instantiatedObjects.begin() + i);
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
        this->scene->recomputeTerrainPropertiesForObject(object->getDefinition()->name);
}

void EnvObjsInterface::displayProbas(const std::string& objectName)
{
    focusAreaEditing = false;
    flowfieldEditing = false;
    previewingObjectInPlotter = true;
    // currentlyPreviewedObject = objectName;
    Vector3 dimensions = initialHeightmap.getDimensions();
    bool possible;
    GridF score = computeScoreMap(this->scene, objectName, dimensions, possible, false);
    if (!possible) {
        ImageViewer::get("Object Preview")->addImage(score * 0.f);
    } else {
        ImageViewer::get("Object Preview")->addImage(score);

        float smallestPositive = score.max();
        score.iterate([&](size_t i) {
            if (score[i] > 0.f && score[i] < smallestPositive)
                smallestPositive = score[i];
        });
        score.iterateParallel([&](size_t i) {
            score[i] = std::max(score[i], smallestPositive);
        });
    }
    ImageViewer::get("Object Preview")->show();
    dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(score);
    Q_EMIT updated();
}

void EnvObjsInterface::displayMaterialDistrib(const std::string& materialName)
{
    this->currentMaterialEdited = materialName;
    GridF distribution = this->scene->materials[materialName].currentState;
    EnvMaterialViewer::get("Material")->addImage(distribution);
    EnvMaterialViewer::get("Material")->show();
    dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(distribution);
    Q_EMIT updated();
}

void EnvObjsInterface::manualModificationOfFocusArea()
{
    this->focusAreaEditing = true;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = false;
    // FocusAreaViewer::get("Focus")->addImage(this->renderFocusArea());
    FocusAreaViewer::get("Focus")->addImage(this->focusedArea);
    FocusAreaViewer::get("Focus")->show();
}

void EnvObjsInterface::manualModificationOfFlowfield()
{
    // this->focusAreaEditing = false;
    // this->flowfieldEditing = true;
    // this->previewingObjectInPlotter = false;
    // this->scene->updateFlowfield(userFlowField + this->computeUserKelvinletField(), simulationFlowField);
    // WaterFlowViewer::get("Flowfield")->addImage(this->renderFlowfield());
    WaterFlowViewer::get("Flowfield")->addVectorField(this->scene->flowfield);
    WaterFlowViewer::get("Flowfield")->show();
}

void EnvObjsInterface::resetFlowfield()
{
    this->userFlowField.reset();
    for (int i = 0; i < this->userKelvinlets.size(); i++)
        delete this->userKelvinlets[i];
    this->userKelvinlets.resize(0);

    this->scene->updateFlowfield(GridV3(), simulationFlowField);
    this->addObjectsHeightmaps();
    this->flowErosionSimulation();
    this->updateVectorFieldVisu();

    WaterFlowViewer::get("Flowfield")->addVectorField(this->scene->flowfield);
    Q_EMIT this->updated();
}

void EnvObjsInterface::updateObjectsList()
{
    if (!this->isVisible()) return;
    if (!objectsListWidget) return;
    std::vector<int> currentSelectionsIDs;
    for (auto currentSelection : currentSelections) currentSelectionsIDs.push_back(currentSelection->ID);
    objectsListWidget->clear();
    auto list = this->scene->instantiatedObjects;

    for (auto& obj : list) {

        // float startingScore = obj->fitnessScoreAtCreation;
        // float endingScore = obj->evaluate(obj->evaluationPosition);

        std::string text = obj->getDefinition()->name;
        if (obj->createdManually)
            text += " [*]";
        else
            text += " (" + std::to_string(int(obj->computeGrowingState() * 100.f)) + "% -- " + std::to_string(int(100.f * obj->computeGrowingState2())) + "% -- " + std::to_string(obj->evaluate()) + "/" + std::to_string(obj->fitnessScoreAtCreation) + ")";
        objectsListWidget->addItem(new HierarchicalListWidgetItemBase(text, obj->ID, 0));
        if (isIn(obj->ID, currentSelectionsIDs))
            currentSelections.push_back(obj);
    }
    objectsListWidget->setCurrentItems(currentSelectionsIDs);
}

void EnvObjsInterface::updateObjectsListSelection(QListWidgetItem *__newSelectionItem)
{
    currentSelections.clear();

    for (auto newSelectionItem : objectsListWidget->selectedItems()) {
        auto newSelection = dynamic_cast<HierarchicalListWidgetItemBase*>(newSelectionItem);
        if (!newSelection) {
            continue; // Does it happen?
        }
        int objID = newSelection->ID;
        EnvObjectInstance* selection = nullptr;
        for (auto& obj : this->scene->instantiatedObjects) {
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
            evalPos.z() = subsidedHeightmap.interpolate(evalPos.x(), evalPos.y()) + offsetAbove;
            std::vector<Vector3> evalLines = {evalPos - Vector3(2, 2, 0), evalPos + Vector3(2, 2, 0), evalPos - Vector3(-2, 2, 0), evalPos + Vector3(-2, 2, 0)};
            lines.insert(lines.end(), evalLines.begin(), evalLines.end());
            std::vector<Vector3> evalColors = std::vector<Vector3>(evalLines.size(), Vector3(0.5, 0.5, 1));
            colors.insert(colors.end(), evalColors.begin(), evalColors.end());
        }

        if (auto asPoint = dynamic_cast<EnvPointInstance*>(currentSelection)) {
            continue; // Do not display the points!!
            selectionPos = asPoint->position;
            selectionPos.z() = subsidedHeightmap.interpolate(selectionPos.x(), selectionPos.y()) + offsetAbove;
            std::vector<Vector3> meshPoints = Mesh::getPointsForArrow(selectionPos + Vector3(0, 0, 20), selectionPos);
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(currentSelection)) {
            selectionPos = asCurve->curve.center();
            std::vector<Vector3> meshPoints;
            auto path = asCurve->curve.getPath(50);
            for (size_t i = 0; i < path.size() - 1; i++) {
                auto p1 = path[i];
                auto p2 = path[i + 1];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + offsetAbove);
                Vector3 p2leveled = p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x(), p2.y()) + offsetAbove);
                meshPoints.push_back(p1leveled);
                meshPoints.push_back(p2leveled);
            }
            for (int i = 0; i < asCurve->curve.size(); i++) {
                auto& p1 = asCurve->curve[i];
                auto& p2 = asCurve->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + offsetAbove);
                Vector3 perpendicular = (p2 - p1).rotate(0, 0, deg2rad(90)).normalized() * 1.f;
                meshPoints.push_back(p1leveled + perpendicular);
                meshPoints.push_back(p1leveled - perpendicular);
            }
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(currentSelection)) {
            selectionPos = asArea->curve.center();
            std::vector<Vector3> meshPoints;
            auto path = asArea->curve.getPath(20);
            for (size_t i = 0; i < path.size() - 1; i++) {
                auto p1 = path[i];
                auto p2 = path[i + 1];
                meshPoints.push_back(p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + 5.f));
                meshPoints.push_back(p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x(), p2.y()) + 5.f));
            }
            for (int i = 0; i < asArea->curve.size(); i++) {
                auto& p1 = asArea->curve[i];
                auto& p2 = asArea->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + offsetAbove);
                Vector3 perpendicular = (p2 - p1).rotate(0, 0, deg2rad(90)).normalized() * 1.f;
                meshPoints.push_back(p1leveled + perpendicular);
                meshPoints.push_back(p1leveled - perpendicular);
            }
            lines.insert(lines.end(), meshPoints.begin(), meshPoints.end());
            std::vector<Vector3> meshColors = std::vector<Vector3>(meshPoints.size(), Vector3(1, 0.5, 1));
            colors.insert(colors.end(), meshColors.begin(), meshColors.end());
        } else {
            std::cerr << "Object #" << currentSelection->ID << " (" << currentSelection->getDefinition()->name << ") could not be casted to Point, Curve or Area..." << std::endl;
//            return;
            continue;
        }
    }
    selectedObjectsMesh.colorsArray = colors;
    selectedObjectsMesh.fromArray(lines);

    Q_EMIT this->updated();
}

void EnvObjsInterface::updateNewObjectMesh()
{
    newObjectMesh.clear();

    std::vector<Vector3> lines;
    std::vector<Vector3> colors;
    float offsetAbove = 5.f;

    std::vector<Vector3> meshPoints;
    auto path = objectSkeletonCreation.getPath(50);
    for (size_t i = 0; i < path.size() - 1; i++) {
        auto p1 = path[i];
        auto p2 = path[i + 1];
        Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + offsetAbove);
        Vector3 p2leveled = p2 + Vector3(0, 0, subsidedHeightmap.interpolate(p2.x(), p2.y()) + offsetAbove);
        meshPoints.push_back(p1leveled);
        meshPoints.push_back(p2leveled);
    }
    for (int i = 0; i < objectSkeletonCreation.size(); i++) {
        auto& p1 = objectSkeletonCreation[i];
        auto& p2 = objectSkeletonCreation[std::abs(i - 1)];
        Vector3 p1leveled = p1 + Vector3(0, 0, subsidedHeightmap.interpolate(p1.x(), p1.y()) + offsetAbove);
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
        this->scene->readEnvObjectsFileContent(newDefinition);
        this->previousFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << primitiveDefinitionFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousFileContent != "")
            this->scene->readEnvObjectsFileContent(this->previousFileContent);
    }
}

void EnvObjsInterface::updateMaterialsDefinitions(const std::string &newDefinition)
{
    try {
        this->scene->readEnvMaterialsFileContent(newDefinition);
        this->previousMaterialsFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << materialsDefinitionFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousMaterialsFileContent != "")
            this->scene->readEnvMaterialsFileContent(this->previousMaterialsFileContent);
    }
}

void EnvObjsInterface::updateMaterialsTransformationsDefinitions(const std::string &newDefinition)
{
    try {
        this->scene->readEnvMaterialsTransformationsFileContent(newDefinition);
        this->previousMaterialsTransformationsFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << transformationsFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousMaterialsTransformationsFileContent != "")
            this->scene->readEnvMaterialsTransformationsFileContent(this->previousMaterialsTransformationsFileContent);
    }
}

void EnvObjsInterface::updateScenarioDefinition(const std::string &newDefinition)
{
    try {
        this->scene->readScenarioFileContent(newDefinition);
        this->previousScenarioFileContent = newDefinition;
    } catch (const nlohmann::detail::parse_error& exception) {
        std::cerr << "Error parsing " << scenarioFile.path << "... No change taken into account. Cause:\n" << exception.what() << std::endl;
        if (previousScenarioFileContent != "")
            this->scene->readScenarioFileContent(this->previousScenarioFileContent);
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFitnessFormula(const std::string& formula)
{
    EnvPoint fake;
    // fake.s_FittingFunction = formula;
    fake.s_FitnessFunction = formula;
    try {
        // fake.fittingFunction = EnvObject::parseFittingFunction(formula, "");
        fake.fitnessFunction = EnvObject::parseFittingFunction(formula, "", this->scene.get(), false);

        GridF eval(this->scene->flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3i& p) {
            eval(p) = fake.fitnessFunction(p);
        });
        ImageViewer::get("Fitness Function")->addImage(eval);
        ImageViewer::get("Fitness Function")->show();
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFittingFormula(const std::string& formula)
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = true;

    EnvPoint fake;
    fake.s_FittingFunction = formula;
    // fake.s_FitnessFunction = formula;
    try {
        fake.fittingFunction = EnvObject::parseFittingFunction(formula, "", this->scene.get(), false);
        // fake.fitnessFunction = EnvObject::parseFittingFunction(formula, "");

        GridF eval(this->scene->flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3i& p) {
            eval(p) = fake.fittingFunction(p);
        });
        ImageViewer::get("Fitting Function")->addImage(eval);
        ImageViewer::get("Fitting Function")->show();
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}

void EnvObjsInterface::evaluateAndDisplayCustomFitnessAndFittingFormula(const std::string& fitnessFuncFormula, std::string fittingFuncFormula)
{
    this->focusAreaEditing = false;
    this->flowfieldEditing = false;
    this->previewingObjectInPlotter = true;

    EnvPoint fake;
    fake.s_FittingFunction = (trim(fittingFuncFormula) == "" ? this->scene->availableObjects[getCurrentObjectName()]->s_FittingFunction : fittingFuncFormula);
    fake.s_FitnessFunction = (trim(fitnessFuncFormula) == "" ? this->scene->availableObjects[getCurrentObjectName()]->s_FitnessFunction : fitnessFuncFormula);
    try {
        fake.fittingFunction = EnvObject::parseFittingFunction(fake.s_FittingFunction, "", this->scene.get(), false);
        fake.fitnessFunction = EnvObject::parseFittingFunction(fake.s_FitnessFunction, "", this->scene.get(), false);

        GridV3 eval(this->scene->flowfield.getDimensions());
        eval.iterateParallel([&](const Vector3i& p) {
            eval(p).x() = fake.fitnessFunction(p);
            eval(p).y() = fake.fittingFunction(p);
        });
        ImageViewer::get("Object Preview")->addImage(eval);
        ImageViewer::get("Object Preview")->show();
        // dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->updateScalarFieldToDisplay(eval);
        Q_EMIT updated();
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    }
}


/*
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
}
*/

void EnvObjsInterface::resetScene()
{
    this->scene->reset();
    this->simulationFlowField.reset();
    this->userFlowField.reset();
    for (int i = 0; i < this->userKelvinlets.size(); i++)
        delete this->userKelvinlets[i];
    this->userKelvinlets.resize(0);
    this->materialSimulationStable = false; // We have to compute the simulation again
    for (auto& [obj, patch] : implicitPatchesFromObjects)
        delete patch;
    this->implicitPatchesFromObjects.clear();
//    this->rootPatch->deleteAllChildren();
    this->rootPatch->composables.clear();
    this->rootPatch->updateCache();
    this->currentSelections.clear();

    this->updateEnvironmentFromEnvObjects(true);

    this->focusedArea = GridF(initialHeightmap.getDimensions(), 1.f);

    Q_EMIT this->updated();
}

void EnvObjsInterface::loadScene(const std::string& filename)
{
    this->resetScene();
    nlohmann::json json = nlohmann::json::parse(std::ifstream(filename));

    std::vector<nlohmann::json> allObjects = json["objects"];
    std::vector<nlohmann::json> allMaterials = json["materials"];
    if (json.contains("initialflow")) {
        std::string flowStr = json["initialflow"];
        this->scene->initialFlowfield = json["initialflow"]; // loadGridV3(flowStr, false);
    }
    if (json.contains("userflow")) {
        this->userFlowField = json["userflow"]; //loadGridV3(json["userflow"], false);
    }

    if (json.contains("waterlevel")) {
        float waterLevel = json["waterlevel"];
        dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->setWaterLevel(waterLevel);
    }
    if (json.contains("heightmap")) {
        initialHeightmap = json["heightmap"]; // loadGridF(json["heightmap"], false);
    }

    for (auto mat : allMaterials) {
        this->scene->materials[mat["name"]] = mat; //.fromJSON(mat);
    }

    for (auto obj : allObjects) {
        std::string objectName = obj["name"];
        EnvObjectInstance* newObject = this->scene->instantiate(objectName);
        newObject->age = obj["age"];
        newObject->fitnessScoreAtCreation = obj["fitnessScoreAtCreation"];
        // newObject->evaluationPosition = obj["evaluationPosition"];
        /*std::vector<nlohmann::json> positions = obj["evaluationPositions"];
        for (auto position : positions) {

        }*/

        if (auto asPoint = dynamic_cast<EnvPointInstance*>(newObject)) {
            asPoint->position = obj["position"];
        } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(newObject)) {
            asCurve->curve = obj["curve"];
            newObject->createdManually = true;
        } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(newObject)) {
            asArea->curve = obj["curve"];
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
    for (auto& [objectName, obj] : this->scene->availableObjects) {
        this->scene->recomputeTerrainPropertiesForObject(objectName);
    }
    this->scene->recomputeFlowAndSandProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
    updateEnvironmentFromEnvObjects(true, false, false);
    for (int i = 0; i < 20; i++)
        updateEnvironmentFromEnvObjects(false, false, false);
    updateEnvironmentFromEnvObjects(false, true, false);
    updateObjectsList();
    updateSelectionMesh();

    for (auto& obj : this->scene->instantiatedObjects) {
        obj->premature = false;
        // newObject->createdManually = false;
    }
    updateEnvironmentFromEnvObjects(false, true, false);

    Q_EMIT this->updated();
}

void EnvObjsInterface::saveScene(const std::string& filename)
{
    nlohmann::json mainJson;
    std::vector<nlohmann::json> allObjects(this->scene->instantiatedObjects.size());
    std::vector<nlohmann::json> allMaterials(this->scene->materials.size());

    for (size_t i = 0; i < allObjects.size(); i++) {
        allObjects[i] = *this->scene->instantiatedObjects[i]; //->toJSON();
    }

    size_t i = 0;
    for (auto& [matName, material] : this->scene->materials) {
        allMaterials[i] = material; // .toJSON();
        i++;
    }

    mainJson["objects"] = allObjects;
    mainJson["materials"] = allMaterials;
    mainJson["initialflow"] = this->scene->initialFlowfield; // stringifyGridV3(this->scene->initialFlowfield, false);
    mainJson["userflow"] = this->userFlowField; // stringifyGridV3(this->userFlowField, false);
    mainJson["waterlevel"] = heightmap->properties->waterLevel;
    mainJson["heightmap"] = initialHeightmap; //stringifyGridF(initialHeightmap, false);
    std::ofstream out(filename);
    out << mainJson.dump(1, '\t');
    out.close();
}

void EnvObjsInterface::previewCurrentEnvObjectPlacement(const Vector3 &position)
{
    GridV3 dataV3 = ImageViewer::get("Object Preview")->dataModel->getImage();
    GridF fitnessScoreGrid(dataV3.getDimensions());
    GridF fittingScoreGrid(dataV3.getDimensions());
    dataV3.iterateParallel([&](size_t i) {
        fitnessScoreGrid[i] = dataV3[i].x();
        fittingScoreGrid[i] = dataV3[i].y();
    });

    auto obj = this->scene->availableObjects[getCurrentObjectName()]->instantiate();
    auto score = fittingScoreGrid;

    GridV3 result = GridV3(score.getDimensions());
    GridF resultAlpha = GridF(score.getDimensions(), 0.f);
    ShapeCurve isoline;
    if (auto objAsPoint = dynamic_cast<EnvPointInstance*>(obj)) {
        isoline = ShapeCurve::circle(objAsPoint->getDefinition()->radius, position, 20);
    } else if (auto objAsCurve = dynamic_cast<EnvCurveInstance*>(obj)) {
        BSpline initialCurve;
        SnakeSegmentation& s = obj->snake;
        float targetLength = objAsCurve->getDefinition()->length;
        if (objAsCurve->getDefinition()->curveFollow == EnvCurve::SKELETON) {
            Vector3 dir = gradientFromFieldFunction(obj->getDefinition()->fitnessFunction)(position).rotated90XY().normalize() * targetLength * .1f;
            initialCurve = BSpline({position - dir, position + dir}).resamplePoints(20);
        } else if (objAsCurve->getDefinition()->curveFollow == EnvCurve::ISOVALUE) {
            Vector3 dir = gradientFromFieldFunction(obj->getDefinition()->fitnessFunction)(position).rotated90XY().normalize() * targetLength * .1f;
            initialCurve = BSpline({position - dir, position + dir}).resamplePoints(20);
        } else if (objAsCurve->getDefinition()->curveFollow == EnvCurve::GRADIENTS) {
            Vector3 dir = gradientFromFieldFunction(obj->getDefinition()->fitnessFunction)(position).normalize() * targetLength * .1f;
            initialCurve = BSpline({position - dir, position + dir}).resamplePoints(20);
        }
        objAsCurve->updateCurve(initialCurve);
        s.position = position;

        int maxIterations = 5;
        for (int iteration = 0; iteration < maxIterations; iteration++) {
            objAsCurve->improvePositionning(5);
            initialCurve = objAsCurve->curve;
            BSpline display = initialCurve;
            display.points.push_back(display[0]);
            display.resamplePoints(100);
            for (size_t i = 0; i < display.size(); i++) {
                const auto& pos = display[i];
                result(pos) = colorPalette(float(iteration) / float(maxIterations - 1));
                resultAlpha(pos) = 1.f;
            }
        }
        isoline = initialCurve;
        isoline.closed = false;
    } else if (auto objAsArea = dynamic_cast<EnvAreaInstance*>(obj)) {
        ShapeCurve initialCurve;
        SnakeSegmentation& s = obj->snake;
        s.position = position;
        float fakeRadius = std::sqrt(s.params->targetArea * .5f / PI);

        ShapeCurve curve = ShapeCurve::circle(fakeRadius, position, 20);
        objAsArea->updateCurve(curve);
        int maxIterations = 10;
        for (int iteration = 0; iteration < maxIterations; iteration++) {
            obj->improvePositionning(5);
            initialCurve = objAsArea->curve;
            ShapeCurve display = initialCurve;
            display.points.push_back(display[0]);
            display.resamplePoints(100);
            if (iteration <= 9) {
                for (size_t i = 0; i < display.size(); i++) {
                    const auto& pos = display[i];
                    result(pos) = colorPalette(float(iteration) / float(maxIterations - 1));
                    resultAlpha(pos) = 1.f;
                }
            }
            if (iteration == 9) {
                for (const auto& v : s.randomGreenCoords) {
                    auto pos = computePointFromGreenCoordinates(v, ShapeCurve(s.contour));
                    result(pos).z() = 1;
                    resultAlpha(pos) = 1.f;
                }
            }
        }
        isoline = initialCurve;
    }

    if (isoline.closed) {
        dataV3.iterateParallel([&](const Vector3i& pos) {
            bool insideCurve = isoline.containsXY(pos, false);
            result(pos) += Vector3(.5f, .5f, .5f) * (insideCurve ? 1.f : 0.f);
            // resultAlpha(pos) = (insideCurve ? .5f : 0.f);
        });
    }
    int nbSamples = 500;
    auto path = isoline.getPath(nbSamples); // .resamplePoints(nbSamples).points;
    for (size_t i = 0; i < path.size(); i++) {
        result(path[i]) = Vector3(0, 0, 1.f); //colorPalette(float(i) / float(path.size() - 1));
        resultAlpha(path[i]) = 1.f;
    }
    for (size_t i = 0; i < isoline.size(); i++) {
        result(isoline[i]) = Vector3(1, 1, 1); //colorPalette(float(i) / float(path.size() - 1));
        resultAlpha(isoline[i]) = 1.f;
    }
    ImageViewer::get("Object Preview")->setOverlay(result, resultAlpha);
    ImageViewer::get("Object Preview")->show();
}

/*
void EnvObjsInterface::previewFocusAreaEdition(const Vector3 &mousePos, bool addingFocus)
{
    //        float velocity = (prevPos - mousePos).norm(); // Typically between 0.1 to 1.0
    auto brush = GridF::normalizedGaussian(30, 30, 1, 8.f) * (addingFocus ? 1.f : -1.f) * 8.f;
    this->focusedArea.add(brush, mousePos - brush.getDimensions().xy() * .5f);

    focusedArea.iterateParallel([&](size_t i) {
        focusedArea[i] = std::clamp(focusedArea[i], 0.f, 30.f);
    });
    FocusAreaViewer::get("Focus")->addImage(renderFocusArea());
    FocusAreaViewer::get("Focus")->show();
}
*/

void EnvObjsInterface::previewFlowEdition(const Vector3 &mousePos, const Vector3 &brushDir)
{
    displayProcessTime("updateFlowfield", [&]() {
        this->scene->updateFlowfield(userFlowField + this->computeUserKelvinletField(), simulationFlowField);
    });
    displayProcessTime("updateVectorFieldVisu", [&]() {
        this->updateVectorFieldVisu();
    });

    displayProcessTime("addObjectsHeightmaps", [&]() {
        this->addObjectsHeightmaps();
    });
    displayProcessTime("flowErosionSimulation", [&]() {
        this->flowErosionSimulation();
    });
}

void EnvObjsInterface::previewMaterialEdition(const Vector3 &position, bool addingMaterial)
{
    /*if (this->scene->materials.count(this->currentMaterialEdited) == 0) return;

    auto& materialContent = this->scene->materials[this->currentMaterialEdited].currentState;

    int radius = 10;
    float amount = 5.f;

    GridF mask = GridF::normalizedGaussian(2 * radius + 1, 2 * radius + 1, 1, float(radius / 4.f)) * amount;

    materialContent.add(mask, position - mask.getDimensions().xy() * .5f);

    EnvMaterialViewer::get("Material")->addImage(materialContent);
    EnvMaterialViewer::get("Material")->show();*/
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

    for (auto& obj : this->scene->instantiatedObjects) {
        TerrainTypes material = obj->getDefinition()->material;
        Vector3 col = materialToColor[material];
        if (auto asPoint = dynamic_cast<EnvPointInstance*>(obj)) {
            Vector3 pos = asPoint->position;
            for (int dx = -1; dx <= 2; dx++) {
                for (int dy = -1; dy <= 2; dy++) {
                    img(pos + Vector3(dx, dy)) = col;
                }
            }
        } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(obj)) {
            for (const auto& p : asCurve->curve.getPath(200)) {
                img(p) = col;
            }
        } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(obj)) {
            for (const auto& p : asArea->curve.getPath(200)) {
                img(p) = col;
            }
        }
    }

    ImageViewer::get("Topography")->addImage(img);
    ImageViewer::get("Topography")->show();
}

void EnvObjsInterface::addObjectsHeightmaps()
{
    float absoluteWaterLevel = voxelGrid->getSizeZ() * voxelGrid->properties->waterLevel;
    this->subsidedHeightmap = this->scene->getHeightmap(initialHeightmap, absoluteWaterLevel, flowErosionFactor, displayGrooves);
    /*
    GridF subsidenceFactor = this->scene->scenario.computeSubsidence(initialHeightmap.getDimensions());
    subsidedHeightmap = initialHeightmap * subsidenceFactor;

    float absoluteWaterLevel = voxelGrid->getSizeZ() * voxelGrid->properties->waterLevel;

    groundConstraintedHeights = GridF(subsidedHeightmap.getDimensions()); // Heightmaps from the ground
    waterConstraintedHeights = GridF(subsidedHeightmap.getDimensions(), -100000.f); // Heightmaps from the water level
    surfaceHeights = GridF(subsidedHeightmap.getDimensions());
    for (auto& obj : this->scene->instantiatedObjects) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(obj->_patch)) {
            GridF grid = GridF(subsidedHeightmap.getDimensions(), 0.f);
            grid = grid.paste(obj->createHeightfield() * obj->computeGrowingState2(), patch->position.xy());
            if (flowErosionFactor != 0 && this->scene->materials.count(toLower(stringFromMaterial(obj->material)))) {
                grid = grid.warpWith(this->scene->flowfield * flowErosionFactor * this->scene->materials[toLower(stringFromMaterial(obj->material))].waterTransport, 10);
            }
            if (obj->heightFrom == this->scene->HeightmapFrom::SURFACE) {
                surfaceHeights = (surfaceHeights + grid * (isIn(obj->material, LayerBasedGrid::invisibleLayers) ? -1.f : 1.f)).max(-15.f);
            } else if (obj->heightFrom == this->scene->HeightmapFrom::GROUND) {
                groundConstraintedHeights = groundConstraintedHeights.max(grid * subsidenceFactor, Vector3());
            } else if (obj->heightFrom == this->scene->HeightmapFrom::WATER) {
                grid.iterateParallel([&] (size_t i) {
                    grid[i] = (std::abs(grid[i]) < 1e-4 ? -10000.f : grid[i]);
                });
                // std::cout << "Max height for " << obj->name << ": " << grid.max() << " while height = " << obj->height << "(grow = " <<  obj->computeGrowingState2() << ")" << std::endl;
                waterConstraintedHeights = waterConstraintedHeights.max((grid - (obj->height)) - (obj->name == "lagoon" || obj->name == "smalllagoon" ? 3.f : 1.f), Vector3()); // Not sure why I need to multiply by 2.0, but otherwise, maxHeight is heigher than obj->height...
            }
        }
    }
    // if (flowErosionFactor != 0) {
        // groundConstraintedHeights = groundConstraintedHeights.warpWith(this->scene->flowfield * flowErosionFactor, 10);
        // waterConstraintedHeights = waterConstraintedHeights.warpWith(this->scene->flowfield * flowErosionFactor, 10);
        // surfaceHeights = surfaceHeights.warpWith(this->scene->flowfield * flowErosionFactor, 10);
    // }
    // Dirty, remove when you understand why lagoon get over the water...
    waterConstraintedHeights.iterateParallel([&](size_t i) {
        waterConstraintedHeights[i] = std::min(waterConstraintedHeights[i] + absoluteWaterLevel, absoluteWaterLevel - 1.f);
    });
    waterConstraintedHeights = waterConstraintedHeights.meanSmooth(3, 3, 1);
    // waterConstraintedHeights = waterConstraintedHeights.gaussianSmooth(1.f, true);

    bool modificationsAppliedToSurface = false;
    for (auto& obj : this->scene->instantiatedObjects) {
        if (displayGrooves) {
            if (endsWith(toLower(obj->name), "reef")) {
                auto objAsEnvCurve = dynamic_cast<EnvCurve*>(obj);
                BSpline path = objAsEnvCurve->curve;
                float nbGrooves = path.length() / 10.f;
                float sigma = objAsEnvCurve->width;
                surfaceHeights.iterateParallel([&](const Vector3i& pos) {
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
                    float grooves = std::max(0.f, 1.f - (sizeX * std::abs(newSpace.x() - 1.f/sizeX) + std::pow(sizeY * newSpace.y(), 2.f)));
                    // return std::max(grooves, initialDistance);
                    const Vector3& flow = this->scene->flowfield(pos);
                    surfaceHeights(pos) += 2.f * grooves * std::max(abs(flow.dot(fakeNormal)), 0.f);
                });
                modificationsAppliedToSurface = true;
            }
        }
    }

    if (modificationsAppliedToSurface) {
        surfaceHeights = surfaceHeights.meanSmooth(3, 3, 1);
        // surfaceHeights = surfaceHeights.gaussianSmooth(1.f, true, true);
    }

    subsidedHeightmap = GridF::max(GridF::max(subsidedHeightmap, groundConstraintedHeights), waterConstraintedHeights).meanSmooth(5, 5, 1);
    subsidedHeightmap = (subsidedHeightmap.max(-15.f) + surfaceHeights).meanSmooth(3, 3, 1).max(-15.f);
    // subsidedHeightmap = GridF::max(GridF::max(subsidedHeightmap, groundConstraintedHeights), waterConstraintedHeights).gaussianSmooth(2.f, true, true);
    // subsidedHeightmap = (subsidedHeightmap.max(-15.f) + surfaceHeights).gaussianSmooth(1.f, true, true).max(-15.f);
    */
}

void EnvObjsInterface::flowErosionSimulation()
{
    // this->heightmap->fromVoxelGrid(*this->voxelGrid);
    if (flowErosionFactor != 0) {
        // this->subsidedHeightmap = (this->initialHeightmap * this->scene->scenario.computeSubsidence()).warpWith(this->scene->flowfield * flowErosionFactor, 10);
        //this->subsidedHeightmap = (this->subsidedHeightmap).warpWith(this->scene->flowfield * flowErosionFactor, 10);
    }
    this->heightmap->heights = subsidedHeightmap - 1.f; // .gaussianSmooth(.5f, true) - 1.f;
}

void EnvObjsInterface::startNewObjectCreation()
{
    this->objectSkeletonCreation = BSpline();
    this->updateNewObjectMesh();
}

void EnvObjsInterface::addPointOnNewObjectCreation(const Vector3 &position, bool addPoint, float removeRadius)
{
    if (addPoint) {
        auto objectModel = this->scene->availableObjects[getCurrentObjectName()];
        this->objectSkeletonCreation.points.push_back(position.xy());
        if (objectModel->isPoint()) {
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
    auto objectModel = this->scene->availableObjects[getCurrentObjectName()];

    if ((objectModel->isCurve() || objectModel->isArea()) && this->objectSkeletonCreation.empty())
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
    // newObject->getDefinition()->fittingFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FittingFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
    // newObject->getDefinition()->fitnessFunction = EnvObject::parseFittingFunction(newObject->getDefinition()->s_FitnessFunction, newObject->getDefinition()->name, this->scene.get(), true, newObject);
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
    this->scene->recomputeTerrainPropertiesForObject(newObject->getDefinition()->name);
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
            if (currentSelection->getDefinition()->isPoint()) {
                currentSelection->translate(translation);
            } else if (EnvCurveInstance* curve = dynamic_cast<EnvCurveInstance*>(currentSelection)) {
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
            } else if (EnvAreaInstance* area = dynamic_cast<EnvAreaInstance*>(currentSelection)) {
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
            currentSelection->translate(translation);
        }

        // Also do it for the creation curve
        objectSkeletonCreation.translate(translation);
        this->updateNewObjectMesh();
        this->updateSelectionMesh();
    }

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
                this->scene->recomputeTerrainPropertiesForObject(currentSelection->getDefinition()->name);
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
    return objectCombobox->choices[objectCombobox->combobox()->currentIndex()]->label;
}

void EnvObjsInterface::updateVectorFieldVisu()
{
    GridV3 velocities = this->scene->flowfield.resize(Vector3(50, 50, 1));
    Mesh::createVectorField(velocities, this->voxelGrid->getDimensions(), &velocitiesMesh, 1.f, false, true);
}





void EnvObjsInterface::saveForRenders()
{
    Vector3 imageDimensions(256, 256, 1);
    Vector3 ratio = imageDimensions / subsidedHeightmap.getDimensions();
    float maxHeight = voxelGrid->getSizeZ();
    GridF depthMap = (heightmap->properties->waterLevel * maxHeight) - subsidedHeightmap; // .gaussianSmooth(1.f, true));;

    std::map<std::string, Vector3> featuresToColors = {
        {"abyss", Vector3(255,   0,   0)},
        {"reef", Vector3(  0,   0, 255)},
        {"lagoon", Vector3(  0, 255, 255)},
        {"coast", Vector3(  0, 255,   0)},
        {"island", Vector3(255, 255,   0)}
    };
    GridV3 labels(imageDimensions, featuresToColors["abyss"]);
    labels.iterateParallel([&](const Vector3i& _p) {
        Vector3 p = _p / ratio;
        float depth = depthMap(p);

        // Add colors for lagoons
        for (auto& obj : this->scene->instantiatedObjects) {
            if (toLower(obj->getDefinition()->name) != "lagoon" && toLower(obj->getDefinition()->name) != "smalllagoon") continue;
            EnvAreaInstance* asArea = dynamic_cast<EnvAreaInstance*>(obj);
            if (asArea->curve.containsXY(p, false)) {
                labels(_p) = featuresToColors["lagoon"];
            }
        }

        // Add colors for beaches
        for (auto& obj : this->scene->instantiatedObjects) {
            if (toLower(obj->getDefinition()->name) != "beach") continue;
            EnvAreaInstance* asArea = dynamic_cast<EnvAreaInstance*>(obj);
            if (asArea->curve.containsXY(p, false)) {
                labels(_p) = featuresToColors["coast"];
            }
        }


        // Add colors for reefs
        for (auto& obj : this->scene->instantiatedObjects) {
            if (toLower(obj->getDefinition()->name) != "reef" && toLower(obj->getDefinition()->name) != "greatreef") continue;
            EnvCurveInstance* asCurve= dynamic_cast<EnvCurveInstance*>(obj);
            if (asCurve->curve.estimateDistanceFrom(p) < asCurve->getDefinition()->width * .5f) {
                labels(_p) = featuresToColors["reef"];
            }
        }

        // Add colors for islands
        for (auto& obj : this->scene->instantiatedObjects) {
            if (toLower(obj->getDefinition()->name) != "island") continue;
            EnvAreaInstance* asArea = dynamic_cast<EnvAreaInstance*>(obj);
            if (asArea->curve.containsXY(p, false)) {
                if (depth > 10)
                    labels(_p) = featuresToColors["lagoon"];
                else if (depthMap(p) > 0)
                    labels(_p) = featuresToColors["coast"];
                else
                    labels(_p) = featuresToColors["island"];
            }
        }
    });

    auto flowToRGB = [&](const GridV3& flow) -> GridV3 {
        GridV3 rgb = flow;
        Vector3 mid(0.5, 0.5, 0.5), mini(-1, -1, -1), maxi(1, 1, 1), flowScale(0.2, 0.2, 0.2);
        rgb.iterateParallel([&](size_t i) {
            rgb[i] = mid + (Vector3::min(Vector3::max(rgb[i], mini), maxi) * flowScale);
        });
        return rgb;
    };

    auto heightmapScaling = [&](const GridF& height) -> GridF {
        auto newHeights = (height / maxHeight);
        newHeights.iterateParallel([&](size_t i) {
            newHeights[i] = std::clamp(newHeights[i], 0.f, 1.f);
        });
        return newHeights;
    };

    std::map<std::string, GridF> objectsScores;
    // Create all maps
    for (auto& [name, obj] : this->scene->availableObjects) {
        objectsScores[name] = GridF(heightmap->getDimensions());
    }
    // Fill each maps
    for (auto& obj : this->scene->instantiatedObjects) {
        auto& grid = objectsScores[obj->getDefinition()->name];
        float evaluationScore = obj->evaluate();
        if (auto asPoint = dynamic_cast<EnvPointInstance*>(obj)) {
            grid(asPoint->position) = std::max(grid(asPoint->position), evaluationScore);
        } else if (auto asCurve = dynamic_cast<EnvCurveInstance*>(obj)) {
            for (auto& p : asCurve->curve.getPath(500)){
                grid(p) = std::max(grid(p), evaluationScore);
            }
        } else if (auto asArea = dynamic_cast<EnvAreaInstance*>(obj)) {
            grid.iterateParallel([&](const Vector3i& p) {
                if (asArea->curve.containsXY(p)) {
                    grid(p) = std::max(grid(p), evaluationScore);
                }
            });
        }
    }
    // Add blur + dilatation to maps
    for (auto& [name, g] : objectsScores) {
        GridF grid = g;
        g = grid.dilate(true, 3.f);
        grid = g.gaussianSmooth(2.f, true);
        g.iterateParallel([&](size_t i) {
            g[i] = std::max(g[i], grid[i]);
        });
        g.normalize();
    }

    auto saveEnvObjsToJSON = [&](const std::string& filename) {
        auto meshInterface = std::static_pointer_cast<MeshInstanceAmplificationInterface>(this->findOtherInterface("meshinstance"));
        nlohmann::json json;
        for (auto& meshType : meshInterface->meshesOptions) {
            json["assets"][meshType.name] = meshType.folderPath;
            if (meshType.positions.empty()) continue;

            bool ok = true;
            auto scoreMap = computeScoreMap(this->scene, meshType.name, this->scene->flowfield.getDimensions(), ok);
            json["instances"][meshType.name] = meshType.currentInstancesToJSON(this->scene, scoreMap);
        }

        std::ofstream out(filename);
        out << json.dump(1, '\t');
        out.close();
    };


    // Make folder "EnvObjRendering/{time}/"
    time_t now = std::time(0);
    tm *gmtm = std::gmtime(&now);
    char s_time[80];
    std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
    std::string folder = "EnvObjRendering/" + std::string(s_time) + "/";
    makedir(folder);

    Image(labels).writeToFile(folder + "/input_label.png");
    Image(heightmapScaling(initialHeightmap)).writeToFile(folder + "heightmap_initial.png");
    Image(heightmapScaling(subsidedHeightmap)).writeToFile(folder + "heightmap_subsided.png");
    Image(heightmapScaling(groundConstraintedHeights)).writeToFile(folder + "heightmap_ground-constraint.png");
    Image(heightmapScaling(waterConstraintedHeights)).writeToFile(folder + "heightmap_water-constraint.png");
    Image(heightmapScaling(surfaceHeights + maxHeight * 0.5)).writeToFile(folder + "heightmap_surface-constraint.png");
    Image(flowToRGB(userFlowField)).writeToFile(folder + "flowfield_user.png");
    Image(flowToRGB(simulationFlowField)).writeToFile(folder + "flowfield_simu.png");
    Image(flowToRGB(this->scene->flowfield)).writeToFile(folder + "flowfield_total.png");
    saveEnvObjsToJSON(folder + "terrain_saved.json");

    for (auto& [name, g] : objectsScores) {
        Image(g).writeToFile(folder + "obj_score_" + name + ".png");
    }

    for (auto& [name, material] : this->scene->materials) {
        Image(material.currentState.normalized()).writeToFile(folder + "material_" + name + ".png");
    }
    Image(this->scene->allScalarProperties["current.vel"]).writeToFile(folder + "material_current.png");

    log("Exported to '" + folder + "'.");
}





/*
StatsValues EnvObjsInterface::displayStatsForObjectCreation(const std::string& objectName, int nbSamples)
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
*/

GridV3 EnvObjsInterface::computeUserKelvinletField() const
{
    GridV3 flow(this->userFlowField.getDimensions());
    flow.iterateParallel([&](const Vector3i& p) {
        Vector3 v;
        for (auto& kelvinlet : this->userKelvinlets) {
            v += kelvinlet->evaluate(p);
        }
        flow(p) = v;
    });
    return flow;
}

void EnvObjsInterface::fromGanUI()
{
    this->resetScene();
    this->scene->readEnvObjectsFileContent(this->primitiveDefinitionFile.read());
    this->scene->readEnvMaterialsFileContent(this->materialsDefinitionFile.read());

    std::string path = "Python_tests/test_island_heightmapfeatures/";
//    QString q_filename= QString::fromStdString(path + "_1.png");
    QString q_filename= QString::fromStdString(path + "1.png");  //QFileDialog::getOpenFileName(this, "Open feature map", QString::fromStdString(path), "*", nullptr);
    if (!q_filename.isEmpty()) {
        std::string file = q_filename.toStdString();
        GridV3 img = Image::readFromFile(file).colorImage;

        auto envObjects = CoralIslandGenerator::envObjsFromFeatureMap(img, voxelGrid->getDimensions(), this->scene);
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

        this->scene->precomputeTerrainProperties(subsidedHeightmap, heightmap->properties->waterLevel, voxelGrid->getSizeZ());
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
                                                                Vector3::invalid,
                                                                Vector3::invalid,
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
