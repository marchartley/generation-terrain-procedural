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
    primitiveDefinitionFile.onChange([&](std::string newDefinitions) { updateObjectsDefinitions(newDefinitions); });
    materialsDefinitionFile.onChange([&](std::string newDefinitions) { updateMaterialsDefinitions(newDefinitions); });
    transformationsFile.onChange([&](std::string newDefinitions) { updateMaterialsTransformationsDefinitions(newDefinitions); });

    QTimer* hotreloadTimer = new QTimer(this);
    hotreloadTimer->setInterval(500);
    QObject::connect(hotreloadTimer, &QTimer::timeout, this, [&]() {
        materialsDefinitionFile.check();
        primitiveDefinitionFile.check();
        transformationsFile.check();
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
    this->objectsMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader), true, GL_LINES);
    this->objectsMesh.useIndices = false;

    this->rootPatch = new ImplicitNaryOperator;
    this->implicitTerrain->addChild(rootPatch);
//    createEnvObjectsFromImplicitTerrain();
//    recomputeErosionValues();
    EnvObject::precomputeTerrainProperties(*heightmap);

    this->focusedArea = GridF(heightmap->heights.getDimensions(), 1.f);


    QObject::connect(Plotter::get(), &Plotter::clickedOnImage, this, [&](const Vector3& clickPos, Vector3 value) {
        if (!this->isVisible()) return;
        if (focusAreaEditing) return;
        GridV3 dataV3 = Plotter::get()->displayedImage;
        GridF data(dataV3.getDimensions());
        dataV3.iterateParallel([&](size_t i) {
            data[i] = dataV3[i].x;
        });
        GridV3 gradients = data.gradient();
        data.raiseErrorOnBadCoord = false;
        data.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        gradients.raiseErrorOnBadCoord = false;
        gradients.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        auto obj = EnvObject::availableObjects[currentlyPreviewedObject];
        GridV3 result(data.getDimensions());
        ShapeCurve isoline;
        if (auto asArea = dynamic_cast<EnvArea*>(obj)) {
            isoline = this->computeNewObjectsShapeAtPosition(clickPos, gradients, data, 1000.f, 1.f).close();
            result.iterateParallel([&](const Vector3& pos) {
                if (isoline.containsXY(pos)) {
                    result(pos) = Vector3(.5f, .5f, .5f);
                }
            });
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(obj)) {
            isoline = this->computeNewObjectsCurveAtPosition(clickPos, gradients, data, asCurve->length, asCurve->width, asCurve->curveFollow == EnvCurve::ISOVALUE);
        } else if (auto asPoint = dynamic_cast<EnvPoint*>(obj)) {
            isoline = ShapeCurve::circle(asPoint->radius, clickPos, 20);
            result.iterateParallel([&](const Vector3& pos) {
                if (isoline.containsXY(pos)) {
                    result(pos) = Vector3(.5f, .5f, .5f);
                }
            });
        }

        result += dataV3;


        int nbSamples = 500;
        auto path = isoline.resamplePoints(nbSamples).points;
        for (size_t i = 0; i < path.size(); i++) {
            result(path[i]) = colorPalette(float(i) / float(path.size() - 1)); //Vector3(1.f, 0, 0);
        }
        Plotter::get()->addImage(result);
        Plotter::get()->show();
        Plotter::get()->addImage(dataV3);
    });

    QObject::connect(Plotter::get(), &Plotter::movedOnImage, this, [&](const Vector3& mousePos, bool pressed, const Vector3& prevPos, QMouseEvent* event) {
        if (!this->isVisible()) return;
        if (!this->focusAreaEditing) return;

        if (!pressed) return;

        auto brush = GridF::normalizedGaussian(20, 20, 1, 8.f) * (event->modifiers().testFlag(Qt::ShiftModifier) ? -1.f : 1.f);
        this->focusedArea.add(brush, mousePos - brush.getDimensions().xy() * .5f);
        Plotter::get()->addImage(focusedArea);
        Plotter::get()->show();
    });
}

void EnvObjsInterface::display(const Vector3 &camPos)
{
    if (!this->visible)
        return;

    if (this->waitAtEachFrame)
        this->updateEnvironmentFromEnvObjects(true, false);

    velocitiesMesh.shader->setVector("color", std::vector<float>{.2f, .2f, .8f, .5f});
    velocitiesMesh.display(GL_LINES, 5.f);

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
    ButtonElement* spendTimeButton = new ButtonElement("Wait", [&]() {
        this->updateEnvironmentFromEnvObjects(/* meh... I don't know if I should update the terrain or not */);
        this->updateSelectionMesh();
        this->saveScene("testEnvObjects.json");
    });
    CheckboxElement* waitAtEachFrameButton = new CheckboxElement("Auto wait", this->waitAtEachFrame);
//    ButtonElement* createFromGAN = new ButtonElement("From GAN", [&]() { this->fromGanUI(); });
    ButtonElement* createFromFile = new ButtonElement("From file", [&]() { this->loadScene("testEnvObjects.json"); });
    TextEditElement* testingFormula = new TextEditElement("", "Try: ");
    testingFormula->setOnTextChange([&](std::string expression) { this->evaluateAndDisplayCustomCostFormula(expression); });
    ButtonElement* testPerformancesButton = new ButtonElement("Run test", [&]() { this->runPerformanceTest(); });
    ButtonElement* resetButton = new ButtonElement("Reset scene", [&]() { this->resetScene(); });

    ButtonElement* instantiaABCbutton = new ButtonElement("Instantiate ABC", [&]() {
        this->instantiateSpecific("CoralPolypFlatA", false);
        this->instantiateSpecific("CoralPolypFlatB", false);
        this->instantiateSpecific("CoralPolypFlatC", false);
    });

    QLabel* label = new QLabel(QString::fromStdString("Objects: " + std::to_string(EnvObject::instantiatedObjects.size())));

    objectsListWidget = new HierarchicalListWidget;
    objectsListWidget->setSelectionMode(QAbstractItemView::SelectionMode::ExtendedSelection);
    updateObjectsList();
//    QObject::connect(objectsListWidget, &HierarchicalListWidget::clicked, this, [&](QModelIndex item) { qDebug() << item; }); //&EnvObjsInterface::updateObjectsListSelection);
//    QObject::connect(objectsListWidget, &HierarchicalListWidget::currentItemChanged, this, [&](QListWidgetItem* current, QListWidgetItem* previous) { this->updateObjectsListSelection(current); });
    QObject::connect(objectsListWidget, &HierarchicalListWidget::itemSelectionChanged, this, [&]() { this->updateObjectsListSelection(); });


//    std::vector<QWidget*> probaButtons;
    std::vector<ComboboxLineElement> objectsChoices;
    for (auto& [name, obj] : EnvObject::availableObjects) {
//        ButtonElement* showButton = new ButtonElement("Show " + obj->name, [&](){ this->displayProbas(name); });
//        ButtonElement* forceButton = new ButtonElement("Force", [&](){ this->instantiateSpecific(name); });
//        probaButtons.push_back(showButton->get());
//        probaButtons.push_back(forceButton->get());
        objectsChoices.push_back(ComboboxLineElement{name, 0});
    }
    ButtonElement* showButton = new ButtonElement("Show", [&](){ this->displayProbas(objectCombobox->choices[objectCombobox->combobox()->currentIndex()].label); });
    ButtonElement* forceButton = new ButtonElement("Force", [&](){ this->instantiateSpecific(objectCombobox->choices[objectCombobox->combobox()->currentIndex()].label); });
    objectCombobox = new ComboboxElement("Objects", objectsChoices);
//    objectCombobox->setOnSelectionChanged([=](int newSelectionIndex) {
//        std::cout << newSelectionIndex << " " << objectCombobox->choices[newSelectionIndex].label << std::endl;
//    });

    std::vector<QWidget*> materialsButtons;
    for (auto& [name, material] : EnvObject::materials) {
        ButtonElement* showButton = new ButtonElement("Show " + toCapitalize(name), [&](){ this->displayMaterialDistrib(name); });
        materialsButtons.push_back(showButton->get());
    }

    ButtonElement* editFocusAreaButton = new ButtonElement("Edit focus", [&]() { this->manualModificationOfFocusArea(); });

    layout->addWidget(spendTimeButton->get());
    layout->addWidget(waitAtEachFrameButton->get());
//    layout->addWidget(showDepositionButton->get());
//    layout->addWidget(showPolypButton->get());
    layout->addWidget(createMultiColumnGroup(materialsButtons, 2));
//    layout->addWidget(showFlowfieldButton->get());
//    layout->addWidget(createFromGAN->get());
//    layout->addWidget(instantiateButton->get());
//    layout->addWidget(createMultiColumnGroup(probaButtons, 2));
    layout->addWidget(objectCombobox->get());
    layout->addWidget(createMultiColumnGroup({showButton->get(), forceButton->get()}, 2));
    layout->addWidget(editFocusAreaButton->get());
    layout->addWidget(objectsListWidget);
//    layout->addWidget(recomputeErosionButton->get());
    layout->addWidget(testingFormula->get());
    layout->addWidget(createHorizontalGroupUI({instantiaABCbutton, testPerformancesButton, resetButton})->get());
    layout->addWidget(createHorizontalGroup({label, createFromFile->get()}));

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

void EnvObjsInterface::show()
{
    ActionInterface::show();
}

void EnvObjsInterface::hide()
{
//    updateObjectsListSelection(nullptr); // Hide single object's data
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
    if (event->modifiers().testFlag(Qt::KeyboardModifier::ShiftModifier)) {
        draggingPoint = mouseWorldPosition.xy();
    } else if (event->modifiers().testFlag(Qt::KeyboardModifier::ControlModifier)) {
        draggingFullObject = mouseWorldPosition.xy();
    }
    draggingHasBeenApplied = mouseWorldPosition.xy();
    draggingHasBeenApplied.setValid(false); // Keep position, but set to invalid
}

void EnvObjsInterface::mouseMovedOnMapEvent(const Vector3& mouseWorldPosition, TerrainModel* model)
{
    if (!this->isVisible()) return;
    if (!mouseWorldPosition.isValid()) return;

    if (draggingPoint.isValid()) {
        draggingHasBeenApplied.setValid(true);
        Vector3 translation = (mouseWorldPosition.xy() - draggingHasBeenApplied.xy());
        draggingHasBeenApplied = mouseWorldPosition.xy();

        for (auto currentSelection : currentSelections) {
            float maxDistToPointSqr = 10.f * 10.f;
            if (EnvPoint* point = dynamic_cast<EnvPoint*>(currentSelection)) {
                point->translate(translation);
            } else if (EnvCurve* curve = dynamic_cast<EnvCurve*>(currentSelection)) {
                auto newCurve = curve->curve;
                int pointIndexToMove = -1;
                float closestDistToPoint = std::numeric_limits<float>::max();

                for (int i = 0; i < curve->curve.size(); i++) {
                    float dist = (curve->curve[i] - mouseWorldPosition.xy()).norm2();
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
                    float dist = (area->curve[i] - mouseWorldPosition.xy()).norm2();
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
        this->updateSelectionMesh();
    } else if (draggingFullObject.isValid()) {
        draggingHasBeenApplied.setValid(true);
        Vector3 translation = (mouseWorldPosition.xy() - draggingHasBeenApplied.xy());
        draggingHasBeenApplied = mouseWorldPosition.xy();

        for (auto currentSelection : currentSelections) {
            if (EnvPoint* point = dynamic_cast<EnvPoint*>(currentSelection)) {
                point->translate(translation);
            } else if (EnvCurve* curve = dynamic_cast<EnvCurve*>(currentSelection)) {
                curve->translate(translation);
            } else if (EnvArea* area = dynamic_cast<EnvArea*>(currentSelection)) {
                area->translate(translation);
            }
        }
        this->updateSelectionMesh();
    }
}

void EnvObjsInterface::mouseReleasedOnMapEvent(const Vector3& mouseWorldPosition, bool mouseInMap, QMouseEvent* event, TerrainModel* model)
{
    if (!this->isVisible()) return;

    if (draggingFullObject.isValid() && !mouseInMap) {
        for (auto currentSelection : currentSelections) {
            this->destroyEnvObject(currentSelection);
        }
        this->updateEnvironmentFromEnvObjects(true, true);
    }
    if ((draggingPoint.isValid() || draggingFullObject.isValid()) && draggingHasBeenApplied.isValid()) {
        draggingPoint.setValid(false);
        draggingFullObject.setValid(false);
        if (mouseInMap) {
            for (auto currentSelection : currentSelections) {
                currentSelection->age = 0.f;
                if (this->implicitPatchesFromObjects.count(currentSelection) != 0) {
                    auto newPatch = currentSelection->createImplicitPatch(heightmap->heights);
                    if (newPatch) {
                        *(this->implicitPatchesFromObjects[currentSelection]) = *newPatch;
                        delete newPatch;
                    }
                }
                EnvObject::recomputeTerrainPropertiesForObject(*heightmap, currentSelection->name);
            }
        }
        this->materialSimulationStable = false;
        this->updateEnvironmentFromEnvObjects(true, true);

        this->updateSelectionMesh();
    }
    draggingPoint.setValid(false);
    draggingFullObject.setValid(false);
    draggingHasBeenApplied.setValid(false);
}

GridF computeScoreMap(std::string objectName, const Vector3& dimensions, bool& possible, bool applyNormalization = false) {
    auto obj = EnvObject::availableObjects[objectName];
    GridF score = GridF(dimensions);
    score.iterateParallel([&](const Vector3& pos) {
        score(pos) = std::max(obj->evaluate(pos), 0.f);
    });
    if (abs(score.min() - score.max()) > 1e-5) {
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
    float scoreToObtain = random_gen::generate(score.sum());
    float bestScore = 10000;
    Vector3 bestPos = Vector3::invalid();

    int nbIterationsToFindBestPos = 0;
    while (scoreToObtain > 0) {
        Vector3 p = Vector3::random(Vector3(), score.getDimensions().xy());
        /*for (int iter = 0; iter < 50; iter++) {
            p += gradients(p).normalize();
        }*/
        float scoreObtained = score(p);
        scoreToObtain -= scoreObtained;
        bestPos = p;
        nbIterationsToFindBestPos++;
    }
    return bestPos;
}

EnvObject* EnvObjsInterface::instantiateObjectAtBestPosition(std::string objectName, const Vector3& position, const GridF& score) {
    EnvObject* newObject = EnvObject::instantiate(objectName);

    auto objAsPoint = dynamic_cast<EnvPoint*>(newObject);
    auto objAsCurve = dynamic_cast<EnvCurve*>(newObject);
    auto objAsArea = dynamic_cast<EnvArea*>(newObject);

    if (objAsPoint) {
        // Nothing to do...
    } else if (objAsCurve) {
        GridV3 gradients = score.gaussianSmooth(2.f).gradient();
        gradients.raiseErrorOnBadCoord = false;
        gradients.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        BSpline initialCurve = this->computeNewObjectsCurveAtPosition(position, gradients, score, objAsCurve->length, objAsCurve->width, objAsCurve->curveFollow == EnvCurve::ISOVALUE);
        BSpline curve = initialCurve;
        curve.translate(-position);
        objAsCurve->curve = curve.resamplePoints(10);
    } else if (objAsArea) {
        if (toLower(objectName) == "island") {
            float width = objAsArea->width;
            objAsArea->curve = ShapeCurve::circle(width, Vector3(), 10);
            for (auto& p : objAsArea->curve) {
                p *= random_gen::generate(0.3f, 1.f);
            }
        } else {
            GridV3 gradients = score/*.gaussianSmooth(2.f)*/.gradient();
            gradients.raiseErrorOnBadCoord = false;
            gradients.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
            ShapeCurve initialCurve = this->computeNewObjectsShapeAtPosition(position, gradients, score, 1000.f, 1.f).close();
            ShapeCurve curve = initialCurve.resamplePoints();
            curve.translate(-position);
            objAsArea->curve = curve.resamplePoints(10);
        }
    }

    newObject->translate(position.xy());
    newObject->fittingScoreAtCreation = newObject->evaluate(newObject->evaluationPosition);
    std::cout << "Creation of obj at score = " << newObject->fittingScoreAtCreation << std::endl;
    return newObject;
}

void EnvObjsInterface::instantiateObject(bool waitForFullyGrown)
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
            score.raiseErrorOnBadCoord = false;
            score.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        }

        std::vector<std::string> possibleObjects;
        for (auto& [name, isPossible] : possible)
            if (isPossible)
                possibleObjects.push_back(name);

        if (possibleObjects.size() > 0) {
            this->materialSimulationStable = false; // We have to compute the simulation again
            std::shuffle(possibleObjects.begin(), possibleObjects.end(), random_gen::random_generator);
            std::string name = possibleObjects[0];
            auto& score = scores[name];

            Vector3 bestPos = bestPositionForInstantiation(name, score * focusedArea);
            EnvObject* newObject = instantiateObjectAtBestPosition(name, bestPos, score);
            auto implicit = newObject->createImplicitPatch(heightmap->heights);
//            dynamic_cast<ImplicitPrimitive*>(implicit)->position.z = heightmap->getHeight(bestPos);
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);

            if (implicit != nullptr) {
                rootPatch->addChild(implicit);
            }
            std::cout << "Instantiating " << name << " at position " << bestPos << std::endl;
            // Wait until the object is 100% grown:
            int maxIterations = 100;
            while (waitForFullyGrown && newObject->computeGrowingState() < 1.f) {
                this->updateEnvironmentFromEnvObjects(false, true);
                maxIterations--;
                if (maxIterations < 0) break;
            }
            EnvObject::recomputeTerrainPropertiesForObject(*heightmap, name);
            EnvObject::recomputeFlowAndSandProperties(*heightmap);
            updateEnvironmentFromEnvObjects(implicit != nullptr, true);
        } else {
            std::cout << "No object to instantiate..." << std::endl;
        }
    });
    updateObjectsList();
}

void EnvObjsInterface::instantiateSpecific(std::string objectName, bool waitForFullyGrown)
{
    objectName = toLower(objectName);
    if (EnvObject::availableObjects.count(objectName) == 0) {
        std::cerr << "No object " << objectName << " in database!" << std::endl;
        return;
    }
    bool verbose = false;
    displayProcessTime("Instantiate new " + objectName + " (forced)... ", [&]() {
        GridF heights = heightmap->getHeights();

        bool possible;
        GridF score = computeScoreMap(objectName, heights.getDimensions(), possible) * focusedArea;
        score.raiseErrorOnBadCoord = false;
        score.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;

        if (possible) {
            this->materialSimulationStable = false; // We have to compute the simulation again
            Vector3 bestPos = bestPositionForInstantiation(objectName, score);
            EnvObject* newObject = instantiateObjectAtBestPosition(objectName, bestPos, score);
            auto implicit = newObject->createImplicitPatch(heightmap->heights);
            this->implicitPatchesFromObjects[newObject] = implicit;
            if (!isIn((ImplicitPatch*)this->rootPatch, this->implicitTerrain->composables))
                this->implicitTerrain->addChild(this->rootPatch);

            if (implicit != nullptr) {
                rootPatch->addChild(implicit);
            }
            // Wait until the object is 100% grown:
            int maxIterations = 100;
            while (waitForFullyGrown && newObject->computeGrowingState() < 1.f) {
                this->updateEnvironmentFromEnvObjects(false, true);
                maxIterations--;
                if (maxIterations < 0) break;
            }
            this->currentSelections = {newObject};
            EnvObject::recomputeTerrainPropertiesForObject(*heightmap, objectName);
            this->updateEnvironmentFromEnvObjects(implicit != nullptr); // If implicit is null, don't update the map
        } else {
            std::cout << "Nope, impossible to instantiate..." << std::endl;
        }
    }, verbose);
    updateObjectsList();
    updateSelectionMesh();
}

bool EnvObjsInterface::checkIfObjectShouldDie(EnvObject *obj, float limitFactorForDying)
{
    Vector3 evaluationPosition = obj->evaluationPosition;
    float score = obj->evaluate(evaluationPosition);

    return score <= obj->fittingScoreAtCreation * limitFactorForDying;
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

void EnvObjsInterface::updateEnvironmentFromEnvObjects(bool updateImplicitTerrain, bool emitUpdateSignal, bool killObjectsIfPossible)
{
    bool verbose = false;


    std::vector<EnvObject*> immatureObjects;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (obj->computeGrowingState() < 1.f) {
            materialSimulationStable = false;
            immatureObjects.push_back(obj);
        }
    }

    if (killObjectsIfPossible) {
        for (auto& obj : EnvObject::instantiatedObjects) {
            if (isIn(obj, immatureObjects)) continue;
            bool shouldDie = this->checkIfObjectShouldDie(obj, .001f);
            if (shouldDie) {
                float startingScore = obj->fittingScoreAtCreation;
                float endingScore = obj->evaluate(obj->evaluationPosition);
                std::cout << "Object went from " << startingScore << " to " << endingScore << " -> " << std::round(100.f * endingScore / startingScore) << "%" << std::endl;
                this->destroyEnvObject(obj);
            }
        }
    }

    displayProcessTime("Get impacted... ", [&]() {
        EnvObject::beImpactedByEvents();
    }, verbose);

    if (!this->materialSimulationStable) { // If the simulation is stable, don't do anything
        // Get original flowfield, do not accumulate effects (for now).
        displayProcessTime("Get velocity... ", [&]() {
            EnvObject::flowfield = dynamic_cast<WarpedFluidSimulation*>(GlobalTerrainProperties::get()->simulations[WARP])->getVelocities(EnvObject::flowfield.sizeX, EnvObject::flowfield.sizeY, EnvObject::flowfield.sizeZ);
        }, verbose);
        displayProcessTime("Apply effects... ", [&]() {
            bool bigChangesInMaterials = EnvObject::applyEffects(heightmap->heights);
            this->materialSimulationStable = !bigChangesInMaterials;
        }, verbose);
        displayProcessTime("Recompute properties... ", [&]() {
            EnvObject::recomputeFlowAndSandProperties(*heightmap);
        }, verbose);
    }

    if (updateImplicitTerrain) {
        for (auto& [obj, implicit] : this->implicitPatchesFromObjects) {
            obj->createImplicitPatch(heightmap->heights, dynamic_cast<ImplicitPrimitive*>(implicit));
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
    }
    for (auto& obj : immatureObjects) {
        if (obj->computeGrowingState() >= 1.f) {
            // Got mature during this process -> now let's save the fitting score
            obj->fittingScoreAtCreation = obj->evaluate(obj->evaluationPosition);
        }
    }
    if (emitUpdateSignal) {
        Q_EMIT this->updated();
        updateObjectsList();
    }
}

void EnvObjsInterface::onlyUpdateFlowAndSandFromEnvObjects()
{
    EnvObject::applyEffects(heightmap->heights);
}

void EnvObjsInterface::destroyEnvObject(EnvObject *object)
{
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
}

void EnvObjsInterface::displayProbas(std::string objectName)
{
    focusAreaEditing = false;
    currentlyPreviewedObject = objectName;
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

void EnvObjsInterface::displayMaterialDistrib(std::string materialName)
{
    GridF distribution = EnvObject::materials[materialName].currentState;
    Plotter::get()->addImage(distribution);
    Plotter::get()->show();
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

void EnvObjsInterface::displayFlowfieldAsImage()
{
    GridV3 flow = EnvObject::flowfield;
    Plotter::get()->addImage(flow);
    Plotter::get()->show();
}

void EnvObjsInterface::manualModificationOfFocusArea()
{
    this->focusAreaEditing = true;
    Plotter::get()->addImage(this->focusedArea);
    Plotter::get()->show();
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

        float startingScore = obj->fittingScoreAtCreation;
        float endingScore = obj->evaluate(obj->evaluationPosition);

        objectsListWidget->addItem(new HierarchicalListWidgetItem(obj->name + " (" + std::to_string(int(obj->computeGrowingState() * 100.f)) + "% -- " + std::to_string(int(100.f * endingScore / startingScore)) + "%)", obj->ID, 0));
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
        this->checkIfObjectShouldDie(selection);
    }

    this->updateSelectionMesh();
}

void EnvObjsInterface::updateSelectionMesh()
{
    if (currentSelections.empty()) {
        objectsMesh.clear();
        velocitiesMesh.clear();
        return;
    }

    this->velocitiesMesh.fromArray(std::vector<float>{});
    this->objectsMesh.fromArray(std::vector<float>{});
    Vector3 selectionPos;
    std::vector<Vector3> lines;
    std::vector<Vector3> colors;
    float offsetAbove = 5.f;
    for (auto currentSelection : currentSelections) {
        Vector3 evalPos = currentSelection->evaluationPosition;
        evalPos.z = heightmap->getHeight(evalPos.x, evalPos.y) + offsetAbove;
        std::vector<Vector3> evalLines = {evalPos - Vector3(2, 2, 0), evalPos + Vector3(2, 2, 0), evalPos - Vector3(-2, 2, 0), evalPos + Vector3(-2, 2, 0)};
        lines.insert(lines.end(), evalLines.begin(), evalLines.end());
        std::vector<Vector3> evalColors = std::vector<Vector3>(evalLines.size(), Vector3(0.5, 0.5, 1));
        colors.insert(colors.end(), evalColors.begin(), evalColors.end());

        if (auto asPoint = dynamic_cast<EnvPoint*>(currentSelection)) {
            selectionPos = asPoint->position;
            selectionPos.z = heightmap->getHeight(selectionPos.x, selectionPos.y) + offsetAbove;
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
                Vector3 p1leveled = p1 + Vector3(0, 0, heightmap->getHeight(p1.x, p1.y) + offsetAbove);
                Vector3 p2leveled = p2 + Vector3(0, 0, heightmap->getHeight(p2.x, p2.y) + offsetAbove);
                meshPoints.push_back(p1leveled);
                meshPoints.push_back(p2leveled);
            }
            for (int i = 0; i < asCurve->curve.size(); i++) {
                auto& p1 = asCurve->curve[i];
                auto& p2 = asCurve->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, heightmap->getHeight(p1.x, p1.y) + offsetAbove);
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
                meshPoints.push_back(p1 + Vector3(0, 0, heightmap->getHeight(p1.x, p1.y) + 5.f));
                meshPoints.push_back(p2 + Vector3(0, 0, heightmap->getHeight(p2.x, p2.y) + 5.f));
            }
            for (int i = 0; i < asArea->curve.size(); i++) {
                auto& p1 = asArea->curve[i];
                auto& p2 = asArea->curve[std::abs(i - 1)];
                Vector3 p1leveled = p1 + Vector3(0, 0, heightmap->getHeight(p1.x, p1.y) + offsetAbove);
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
    objectsMesh.colorsArray = colors;
    objectsMesh.fromArray(lines);
    Q_EMIT this->updated();
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
    Vector3 dir(1, 0, 0);
    bool didAFullCircle = false;
    float totalDistance = 0.f;
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
    }
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

ShapeCurve EnvObjsInterface::computeNewObjectsShapeAtPosition(const Vector3 &seedPosition, const GridV3& gradients, const GridF& score, float directionLength, float widthMaxLength)
{
    BSpline isoline = followIsovalue(score, gradients, seedPosition, directionLength);
    return isoline;
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
    this->materialSimulationStable = false; // We have to compute the simulation again
    for (auto& [obj, patch] : implicitPatchesFromObjects)
        delete patch;
    this->implicitPatchesFromObjects.clear();
//    this->rootPatch->deleteAllChildren();
    this->rootPatch->composables.clear();
    this->rootPatch->updateCache();

    this->updateEnvironmentFromEnvObjects(true);

    this->focusedArea = GridF(heightmap->heights.getDimensions(), 1.f);

    Q_EMIT this->updated();
}

void EnvObjsInterface::loadScene(std::string filename)
{
    this->resetScene();
    nlohmann::json json = nlohmann::json::parse(std::ifstream(filename));

    std::vector<nlohmann::json> allObjects = json["objects"];
    std::vector<nlohmann::json> allMaterials = json["materials"];

    for (auto mat : allMaterials) {
        EnvObject::materials[mat["name"]].fromJSON(mat);
    }

    for (auto obj : allObjects) {
        std::string objectName = obj["name"];
        EnvObject* newObject = EnvObject::instantiate(objectName);
        newObject->age = obj["age"];
        newObject->fittingScoreAtCreation = obj["fittingScoreAtCreation"];
        newObject->evaluationPosition = json_to_vec3(obj["evaluationPosition"]);

        if (auto asPoint = dynamic_cast<EnvPoint*>(newObject)) {
            asPoint->position = json_to_vec3(obj["position"]);
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(newObject)) {
            asCurve->curve = json_to_bspline(obj["curve"]);
        } else if (auto asArea = dynamic_cast<EnvArea*>(newObject)) {
            asArea->curve = json_to_bspline(obj["curve"]);
        }
        auto implicit = newObject->createImplicitPatch(heightmap->heights);
        this->implicitPatchesFromObjects[newObject] = implicit;
        if (implicit != nullptr) {
            rootPatch->addChild(implicit);
        }
        this->currentSelections = {};
        EnvObject::recomputeTerrainPropertiesForObject(*heightmap, objectName);
    }
    EnvObject::recomputeFlowAndSandProperties(*heightmap);
    updateEnvironmentFromEnvObjects(true, true, false);
    updateObjectsList();
    updateSelectionMesh();
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

    std::ofstream out(filename);
    out << mainJson.dump(1, '\t');
    out.close();
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
            auto implicit = newObject->createImplicitPatch(heightmap->heights);
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

        EnvObject::precomputeTerrainProperties(*heightmap);
        this->updateEnvironmentFromEnvObjects(true, true);
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
