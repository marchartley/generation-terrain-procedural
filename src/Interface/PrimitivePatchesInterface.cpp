#include "PrimitivePatchesInterface.h"
#include "Interface/InterfaceUtils.h"
//#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"

PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("PrimitivePatchesInterface", parent)
{
    previewMesh.cullFace = false;

    this->mainPatch = new ImplicitPrimitive;
    this->mainPatch->name = "Identity";
//    dynamic_cast<ImplicitPrimitive*>(this->mainPatch)->material = WATER;


//    this->mainFilename = "saved_maps/implicit_patches.json";
    this->mainFilename = "saved_maps/trench.json";

    QTimer* hotreloadTimer = new QTimer(this);
    hotreloadTimer->setInterval(500);
    QObject::connect(hotreloadTimer, &QTimer::timeout, this, &PrimitivePatchesInterface::hotReloadFile);
    hotreloadTimer->start();
}

void PrimitivePatchesInterface::display()
{
    if (this->isVisible()) {
        if (previewMesh.shader != nullptr) {
            this->setSelectedShape(this->currentShapeSelected);
            previewMesh.shader->setVector("color", std::vector<float>({0.f, 0.4f, 0.8f, 0.5f}));
            previewMesh.display();
        }
        if (!patchAABBoxMesh.vertexArray.empty() && patchAABBoxMesh.shader != nullptr) {
            GLint polygonMode;
            glGetIntegerv(GL_POLYGON_MODE, &polygonMode);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            patchAABBoxMesh.shader->setVector("color", std::vector<float>({0.f, 0.8f, 0.4f, 1.0f}));
            patchAABBoxMesh.display(GL_LINES, 3.f);
            glPolygonMode(GL_FRONT_AND_BACK, polygonMode);
        }
        this->primitiveControlPoint->display();

        if (this->debugMeshDisplayed) {
//            this->debuggingVoxelsMesh.shader->setVector("color", std::vector<float> {0.f, 0.f, 0.f, 1.f});
            displayDebuggingVoxels();
        }
    }
}

void PrimitivePatchesInterface::reloadShaders()
{
}

void PrimitivePatchesInterface::replay(nlohmann::json action)
{

}

void PrimitivePatchesInterface::affectTerrains(std::shared_ptr<Grid> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid);

    // Waiting for OpenGL to be available to create the shaders...
    std::string shaderPath = "src/Shaders/";
    std::string vertexFilename = shaderPath + "MarchingCubes.vert";
    std::string geometryFilename = shaderPath + "MarchingCubes.geom";
    std::string fragmentFilename = shaderPath + "no_shader.frag";

    debuggingVoxelsMesh.shader = std::make_shared<Shader>(vertexFilename, fragmentFilename, geometryFilename);
    debuggingVoxelsMesh.useIndices = false;
    debuggingVoxelsMesh.cullFace = true;

    this->loadPatchesFromFile(this->mainFilename);
}

QLayout *PrimitivePatchesInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    QRadioButton* createSphereButton = new QRadioButton("Sphere");
    QRadioButton* createBlockButton = new QRadioButton("Block");
    QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
    QRadioButton* createRockButton = new QRadioButton("Rock");
    QRadioButton* createMountainButton = new QRadioButton("Mountain");
    QRadioButton* createDuneButton = new QRadioButton("Dune");
    QRadioButton* createBasinButton = new QRadioButton("Basin");
    QRadioButton* createCaveButton = new QRadioButton("Cave");
    QRadioButton* createArchButton = new QRadioButton("Arch");
    QRadioButton* createNoise2DButton = new QRadioButton("Noise");
    QRadioButton* createFromFileButton = new QRadioButton("...");

//    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");

    QRadioButton* abovePosButton = new QRadioButton("Above");
    QRadioButton* insideTopPosButton = new QRadioButton("Inside top");
    QRadioButton* insideBottomPosButton = new QRadioButton("Inside bottom");
    QRadioButton* fixedPosButton = new QRadioButton("Fixed");

    FancySlider* blendingFactorSlider = new FancySlider(Qt::Orientation::Horizontal, 2.f, 10.f, .1f);

    FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
    FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
    FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
    FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 20.f, .1f);

    QPushButton* addNoiseButton = new QPushButton("Add noise");

    QPushButton* resetButton = new QPushButton("Reset");

    QLabel* selectedFilenameLabel = new QLabel((this->mainFilename == "" ? "No file selected" : QString::fromStdString(this->mainFilename).split("/").back().split("\\").back()));
    QPushButton* fileSelectionButton = new QPushButton("...");
    QCheckBox* enableHotreloadButton = new QCheckBox("Hot reloading ?");

    QLabel* patchesCounterLabel = new QLabel(QString::fromStdString("Prim: " + std::to_string(nbPrimitives) + ", Unary: " + std::to_string(nbUnaryOperators) + ", Binary: " + std::to_string(nbBinaryOperators)));

    QRadioButton* airDensityCheckbox = new QRadioButton("Air");
    QRadioButton* waterDensityCheckbox = new QRadioButton("Water");
    QRadioButton* coralDensityCheckbox = new QRadioButton("Coral");
    QRadioButton* sandDensityCheckbox = new QRadioButton("Sand");
    QRadioButton* dirtDensityCheckbox = new QRadioButton("Dirt");
    QRadioButton* rockDensityCheckbox = new QRadioButton("Rock");
    QRadioButton* bedrockDensityCheckbox = new QRadioButton("Bedrock");

    QCheckBox* applyIntersectionButton = new QCheckBox("Intersection");

    primitiveSelectionGui = new HierarchicalListWidget(this);

    layout->addWidget(createMultiColumnGroup({
                                                    createSphereButton,
                                                    createBlockButton,
                                                    createGaussianButton,
                                                    createRockButton,
                                                    createMountainButton,
                                                    createDuneButton,
                                                    createBasinButton,
                                                    createCaveButton,
                                                    createArchButton,
                                                    createNoise2DButton,
                                                    createFromFileButton
                                          }));
//    layout->addWidget(createSliderGroup("Density", densitySlider));
    layout->addWidget(createMultiColumnGroup({
                                                airDensityCheckbox,
                                                waterDensityCheckbox,
                                                coralDensityCheckbox,
                                                sandDensityCheckbox,
                                                dirtDensityCheckbox,
                                                rockDensityCheckbox,
                                                bedrockDensityCheckbox
                                            }, 3));
    layout->addWidget(createSliderGroup("Blend factor", blendingFactorSlider));
    layout->addWidget(createHorizontalGroup({createVerticalGroup({
                                              stackingButton,
                                              blendingButton,
                                              replacingButton,
                                          }),
                                             createVerticalGroup({
                                                 abovePosButton,
                                                 insideTopPosButton,
                                                 insideBottomPosButton,
                                                 fixedPosButton
                                             })
                                            }));
    layout->addWidget(createMultipleSliderGroup({
                                              {"Width / radius", widthSlider},
                                                    {"Depth", depthSlider},
                                                    {"Height", heightSlider},
                                                    {"Sigma", sigmaSlider}
                                          }));
    layout->addWidget(applyIntersectionButton);
    layout->addWidget(addNoiseButton);
    layout->addWidget(resetButton);
    layout->addWidget(primitiveSelectionGui);

    layout->addWidget(createHorizontalGroup({selectedFilenameLabel, fileSelectionButton}));
    layout->addWidget(enableHotreloadButton);
    layout->addWidget(patchesCounterLabel);

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Sphere); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Block); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Gaussian); });
    QObject::connect(createRockButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Rock); });
    QObject::connect(createMountainButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Mountain); });
    QObject::connect(createDuneButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Dune); });
    QObject::connect(createBasinButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Basin); });
    QObject::connect(createCaveButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Cave); });
    QObject::connect(createArchButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Arch); });
    QObject::connect(createNoise2DButton, &QRadioButton::toggled, this, [=]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Noise2D); });
    QObject::connect(createFromFileButton, &QRadioButton::toggled, this, [=](bool checked) {
        if (checked) {
            this->openFileForNewPatch();
            createSphereButton->setChecked(true);
        }
    });
//    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [=](float newDensity) { this->selectedDensity = newDensity; });

    QObject::connect(stackingButton, &QRadioButton::toggled, this, [=]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [=]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [=]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::REPLACE); });

    QObject::connect(abovePosButton, &QRadioButton::toggled, this, [=]() { this->setCurrentPositioning(ImplicitPatch::PositionalLabel::ABOVE); });
    QObject::connect(insideTopPosButton, &QRadioButton::toggled, this, [=]() { this->setCurrentPositioning(ImplicitPatch::PositionalLabel::INSIDE_TOP); });
    QObject::connect(insideBottomPosButton, &QRadioButton::toggled, this, [=]() { this->setCurrentPositioning(ImplicitPatch::PositionalLabel::INSIDE_BOTTOM); });
    QObject::connect(fixedPosButton, &QRadioButton::toggled, this, [=]() { this->setCurrentPositioning(ImplicitPatch::PositionalLabel::FIXED_POS); });


    QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { this->selectedWidth = newVal; });
    QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { this->selectedHeight = newVal; });
    QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { this->selectedDepth = newVal; });
    QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { this->selectedSigma = newVal; });
    QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { this->selectedBlendingFactor = newVal; });

    QObject::connect(airDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = AIR; });
    QObject::connect(waterDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = WATER; });
    QObject::connect(coralDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = CORAL; });
    QObject::connect(sandDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = SAND; });
    QObject::connect(dirtDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = DIRT; });
    QObject::connect(rockDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = ROCK; });
    QObject::connect(bedrockDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->selectedTerrainType = BEDROCK; });

    QObject::connect(applyIntersectionButton, &QCheckBox::toggled, this, [=](bool checked) { this->applyIntersection = checked; });

    QObject::connect(addNoiseButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::addNoiseOnSelectedPatch);

    QObject::connect(resetButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::resetPatch);

    QObject::connect(fileSelectionButton, &QPushButton::pressed, this, [this, selectedFilenameLabel]() {
        std::string path = "saved_maps/";
        QString fileSelection = QFileDialog::getSaveFileName(this, "Saving file", QString::fromStdString(path), "*.json", nullptr, QFileDialog::DontConfirmOverwrite);
        if (!fileSelection.isEmpty()) {
            this->mainFilename = fileSelection.toStdString();
            selectedFilenameLabel->setText(QString::fromStdString(this->mainFilename).split("/").back().split("\\").back());
            if (QFileInfo::exists(QString::fromStdString(this->mainFilename))) {
                this->loadPatchesFromFile(this->mainFilename);
            } else {
                this->savePatchesAsFile(this->mainFilename);
            }
        }
    });
    QObject::connect(enableHotreloadButton, &QCheckBox::toggled, this, [=](bool checked) {
        this->enableHotReloading = checked;
    });

    createSphereButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Sphere);
    createBlockButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Block);
    createGaussianButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Gaussian);
    createRockButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Rock);
    createMountainButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Mountain);
    createDuneButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Dune);
    createBasinButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Basin);
    createCaveButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Cave);
    createArchButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Arch);
    createNoise2DButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Noise2D);

    stackingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::STACK);
    blendingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::BLEND);
    replacingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::REPLACE);

    abovePosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::ABOVE);
    insideTopPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::INSIDE_TOP);
    insideBottomPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM);
    fixedPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::FIXED_POS);

    blendingFactorSlider->setfValue(this->selectedBlendingFactor);

//    densitySlider->setfValue(this->selectedDensity);
//    TerrainTypes currentMat = LayerBasedGrid::materialFromDensity(this->selectedDensity);
    airDensityCheckbox->setChecked(this->selectedTerrainType == AIR);
    waterDensityCheckbox->setChecked(this->selectedTerrainType == WATER);
    coralDensityCheckbox->setChecked(this->selectedTerrainType == CORAL);
    sandDensityCheckbox->setChecked(this->selectedTerrainType == SAND);
    dirtDensityCheckbox->setChecked(this->selectedTerrainType == DIRT);
    rockDensityCheckbox->setChecked(this->selectedTerrainType == ROCK);
    bedrockDensityCheckbox->setChecked(this->selectedTerrainType == BEDROCK);

    widthSlider->setfValue(this->selectedWidth);
    heightSlider->setfValue(this->selectedHeight);
    depthSlider->setfValue(this->selectedDepth);
    sigmaSlider->setfValue(this->selectedSigma);

    applyIntersectionButton->setChecked(this->applyIntersection);

    enableHotreloadButton->setChecked(this->enableHotReloading);

    this->updatePrimitiveList();

    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::currentItemChanged, this, &PrimitivePatchesInterface::updateSelectedPrimitiveItem);
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemDoubleClicked, this, &PrimitivePatchesInterface::openPrimitiveModificationDialog);
//    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemChangedHierarchy, this, &PrimitivePatchesInterface::modifyPrimitiveHierarchy);


    // This initialization should be in the constructor I guess...
    this->primitiveControlPoint = std::make_unique<ControlPoint>(Vector3(), 5.f);
    this->primitiveControlPoint->allowAllAxisRotations(false);
    this->primitiveControlPoint->allowAllAxisTranslation(true);
    this->primitiveControlPoint->hide();

    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::modified, this, [=](){
        // Get patch being manipulated
        if (this->currentlyManipulatedPatch != nullptr) {
            // Display modified AABBox
            auto AABBox = currentlyManipulatedPatch->getBBox();
            Vector3 minPos = AABBox.first;
            Vector3 maxPos = AABBox.second;
            Vector3 dim = maxPos - minPos;
            std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
            for (auto& p : box) {
                p = (p * dim) + /*minPos + */primitiveControlPoint->getPosition(); //p * currentlyManipulatedPatch->getDimensions() + (primitiveControlPoint->getPosition() - currentlyManipulatedPatch->getDimensions().xy() * .5f);//p * currentlyManipulatedPatch->getDimensions() + (primitiveControlPoint->getPosition() - (currentlyManipulatedPatch->getDimensions() * .5f).xy());
            }
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
        }
    });
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::translationApplied, this, [=](Vector3 translation){
//        QObject::blockSignals(true);
        // Get patch being manipulated
        if (this->currentlyManipulatedPatch != nullptr) {
            primitiveControlPoint->blockSignals(true);
            primitiveSelectionGui->blockSignals(true);

            ImplicitUnaryOperator* manipulatedAsUnary = dynamic_cast<ImplicitUnaryOperator*>(this->currentlyManipulatedPatch);
            if (manipulatedAsUnary != nullptr) { // We are updating an unary operator
                manipulatedAsUnary->translate(translation); // Just update it
            } else { // Otherwise, create a new Unary operator
                ImplicitUnaryOperator* translate = new ImplicitUnaryOperator;
                translate->composableA = this->currentlyManipulatedPatch;
                translate->translate(translation);
                translate->name = "Translation";
                if (currentlyManipulatedPatch == this->mainPatch) {
                    this->mainPatch = translate;
                } else {
                    ImplicitOperator* parent = (ImplicitOperator*)this->naiveApproachToGetParent(currentlyManipulatedPatch);
                    if (currentlyManipulatedPatch == parent->composableA) {
                        parent->composableA = translate;
                    } else {
                        parent->composableB = translate;
                    }
                }
                this->storedPatches.push_back(translate);
                this->updatePrimitiveList();
                this->primitiveSelectionGui->setCurrentItem(translate->index);
            }
            this->updateMapWithCurrentPatch();

            primitiveControlPoint->blockSignals(false);
            primitiveSelectionGui->blockSignals(false);
        }
    });
    /*
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::rotationApplied, this, [=](Vector3 rotation){
        primitiveControlPoint->blockSignals(true);
        primitiveSelectionGui->blockSignals(true);
//        QObject::blockSignals(true);
        // Get patch being manipulated
        if (this->currentlyManipulatedPatch != nullptr) {
            std::cout << "Modifying the terrain..." << std::endl;
//            Vector3 initialControlPosition = currentlyManipulatedPatch->getBBox().first;
//            Vector3 translation = primitiveControlPoint->getPosition() - initialControlPosition;
            // When released, recreate patch + update terrain
            //currentlyManipulatedPatch->position = primitiveControlPoint->getPosition() - (currentlyManipulatedPatch->getDimensions() * .5f).xy();
            ImplicitUnaryOperator* rotate = new ImplicitUnaryOperator;
            rotate->composableA = this->currentlyManipulatedPatch;
            rotate->rotate(rotation.x, rotation.y, rotation.z);
            rotate->name = "Rotation";
            if (currentlyManipulatedPatch == this->mainPatch) {
                this->mainPatch = rotate; // TODO : Not right, but for tests
            } else {
                ImplicitOperator* parent = (ImplicitOperator*)this->naiveApproachToGetParent(currentlyManipulatedPatch);
                if (currentlyManipulatedPatch == parent->composableA) {
                    parent->composableA = rotate;
                } else {
                    parent->composableB = rotate;
                }
            }
            this->updateMapWithCurrentPatch();
            this->storedPatches.push_back(rotate);
            this->updatePrimitiveList();
        }
        primitiveControlPoint->blockSignals(false);
        primitiveSelectionGui->blockSignals(false);
    });
    */


    return layout;
}

void PrimitivePatchesInterface::show()
{
    this->previewMesh.show();
    this->debuggingVoxelsMesh.show();
    this->patchAABBoxMesh.show();
    ActionInterface::show();
}

void PrimitivePatchesInterface::hide()
{
    this->previewMesh.hide();
    this->debuggingVoxelsMesh.hide();
    this->patchAABBoxMesh.hide();
    ActionInterface::hide();
}

void PrimitivePatchesInterface::mouseMovedOnMap(Vector3 newPosition)
{
    if (this->isVisible()) {
        this->currentPos = newPosition;
        this->setSelectedShape(this->currentShapeSelected, newPosition);
    }
}

void PrimitivePatchesInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent *event)
{
    if (this->isVisible() && mouseInMap) {

        if (this->primitiveControlPoint->grabsMouse())
            return; // Don't create patches if we try to displace another patch
        ImplicitPatch* patch = this->createPatchFromParameters(mousePosInMap - this->functionSize.xy() * .5f);
        ImplicitPatch* previousMain;
        if (this->currentlySelectedPatch != nullptr) {
            previousMain = this->currentlySelectedPatch;
        } else {
            previousMain = mainPatch;
        }
        ImplicitPatch* _parent = this->naiveApproachToGetParent(previousMain);
        ImplicitPatch* operation = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;

        if (_parent == nullptr) {
            this->mainPatch = operation;
        } else {
            ImplicitOperator* parent = dynamic_cast<ImplicitOperator*>(_parent);
            if (previousMain == parent->composableA)
                parent->composableA = operation;
            else
                parent->composableB = operation;
        }
//        ImplicitPatch* newOperation = this->createOperationPatchFromParameters(previousMain, patch);
//        this->mainPatch = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;
        this->updateMapWithCurrentPatch();

        this->storedPatches.push_back(patch);
        this->storedPatches.push_back(operation);
        this->storedPatches.push_back(mainPatch);
        this->updatePrimitiveList();

        if (this->mainFilename != "")
            this->savePatchesAsFile(this->mainFilename);
    }
}

void PrimitivePatchesInterface::setSelectedShape(ImplicitPatch::PredefinedShapes newShape, Vector3 newPosition)
{
    newPosition = this->currentPos;
    this->currentShapeSelected = newShape;

    auto vertices = CubeMesh::cubesVertices;
    this->functionSize = Vector3(selectedWidth, selectedDepth, selectedHeight);

    if (this->desiredPatchFromFile != nullptr) {
        this->functionSize = desiredPatchFromFile->getDimensions();
    }
    for (auto& vert : vertices) {
        vert *= functionSize;
        vert += newPosition;
        vert -= functionSize.xy() * .5f;
    }
    this->previewMesh.fromArray(vertices);
    this->previewMesh.update();
}

void PrimitivePatchesInterface::setCurrentOperation(ImplicitPatch::CompositionFunction newOperation)
{
    this->currentOperation = newOperation;
}

void PrimitivePatchesInterface::setCurrentPositioning(ImplicitPatch::PositionalLabel positioning)
{
    this->currentPositioning = positioning;
}

void PrimitivePatchesInterface::resetPatch()
{
    this->storedPatches.clear();
    delete this->mainPatch;
    this->mainPatch = new ImplicitPrimitive;
    this->mainPatch->name = "Identity";
    this->updateMapWithCurrentPatch();
    this->storedPatches.push_back(mainPatch);
    this->updatePrimitiveList();
//    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    //    voxelGrid->fromIsoData();
}

void PrimitivePatchesInterface::updateMapWithCurrentPatch()
{
    this->layerGrid->layers = this->layerGrid->previousState;
//    this->mainPatch->cleanCache();
    this->layerGrid->add(this->mainPatch/*, SAND, false*/);
    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    voxelGrid->fromIsoData();
    heightmap->fromLayerGrid(*layerGrid);

    Q_EMIT this->terrainUpdated();
}

void PrimitivePatchesInterface::setSelectedWidth(float newVal) {
    this->selectedWidth = newVal;
    this->setSelectedShape(this->currentShapeSelected); // Just update the mesh and "functionSize"
}
void PrimitivePatchesInterface::setSelectedHeight(float newVal) {
    this->selectedHeight = newVal;
    this->setSelectedShape(this->currentShapeSelected); // Just update the mesh and "functionSize"
}
void PrimitivePatchesInterface::setSelectedDepth(float newVal) {
    this->selectedDepth = newVal;
    this->setSelectedShape(this->currentShapeSelected); // Just update the mesh and "functionSize"
}
void PrimitivePatchesInterface::setSelectedSigma(float newVal) {
    this->selectedSigma = newVal;
    this->setSelectedShape(this->currentShapeSelected); // Just update the mesh and "functionSize"
}
void PrimitivePatchesInterface::setSelectedBlendingFactor(float newVal) {
    this->selectedBlendingFactor = newVal;
    this->setSelectedShape(this->currentShapeSelected); // Just update the mesh and "functionSize"
}

void PrimitivePatchesInterface::addNoiseOnSelectedPatch()
{
    ImplicitPatch* selectedPatch;
    if (this->primitiveSelectionGui->selectedItems().isEmpty()) {
        selectedPatch = this->mainPatch;
    } else {
        auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(this->primitiveSelectionGui->selectedItems()[0]);
        int selectedPatchID = selectedPatchItem->ID;
        selectedPatch = this->findPrimitiveById(selectedPatchID);
    }

    ImplicitUnaryOperator* noisePatch = new ImplicitUnaryOperator;
    noisePatch->composableA = selectedPatch;
    if (selectedPatch == this->mainPatch) {
        this->mainPatch = noisePatch;
    } else {
        ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
        ImplicitOperator* parent = dynamic_cast<ImplicitOperator*>(_parent);
        if (parent) {
            if (selectedPatch == parent->composableA)
                parent->composableA = noisePatch;
            else
                parent->composableB = noisePatch;
        }
    }
    noisePatch->noiseFunction = [=](Vector3 pos) -> float {
        FastNoiseLite noise;
    //        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float period = 20.f;
        float noiseOffset = 10.f; // sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseVal = noise.GetNoise(pos.x * period + noiseOffset, pos.y * period + noiseOffset, pos.z * period + noiseOffset) * 1.f;
        return noiseVal;
    };
    this->storedPatches.push_back(noisePatch);
    this->updateMapWithCurrentPatch();
    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::savePatchesAsFile(std::string filename)
{
    this->lastTimeFileHasBeenModified = QDateTime::currentDateTime().addDays(1); // Just to avoid the reloading during the save
    nlohmann::json content = this->mainPatch->toJson();
    std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc);
    out << nlohmann::json({ {ImplicitPatch::json_identifier, content} }).dump(1, '\t');
    out.flush();
    this->lastTimeFileHasBeenModified = QDateTime::currentDateTime();
}

void PrimitivePatchesInterface::loadPatchesFromFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    nlohmann::json json_content = nlohmann::json::parse(content);
    this->lastTimeFileHasBeenModified = QFileInfo(QString::fromStdString(filename)).lastModified();
    if (json_content.contains(ImplicitPatch::json_identifier)) {
        this->mainPatch = ImplicitPatch::fromJson(json_content[ImplicitPatch::json_identifier]);
        this->updateMapWithCurrentPatch();
        this->updatePrimitiveList();
    } else {
        std::cerr << "No patches defined in file " << filename << " (no '" << ImplicitPatch::json_identifier << "' found in file)" << std::endl;
    }
}

void PrimitivePatchesInterface::hotReloadFile()
{
    if (this->isVisible() && this->enableHotReloading) {
        if (this->mainFilename != "") {
            QFileInfo infos(QString::fromStdString(this->mainFilename));
            auto lastModifTime = infos.lastModified();
            if (lastModifTime > this->lastTimeFileHasBeenModified || !this->lastTimeFileHasBeenModified.isValid()) {
                this->loadPatchesFromFile(this->mainFilename);
            }
        }
    }
}

void PrimitivePatchesInterface::updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous) {
    bool newSelectionIsExisting = false;
    this->currentlySelectedPatch = nullptr;
    if (current) {
        auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(current);
        int selectedPatchID = selectedPatchItem->ID;
        ImplicitPatch* selectedPatch = this->findPrimitiveById(selectedPatchID);
        this->currentlySelectedPatch = selectedPatch;
        this->currentlyManipulatedPatch = selectedPatch;
        if (selectedPatch != nullptr) {
            auto patchAABBox = selectedPatch->getBBox();
            auto patchSupportedAABBox = selectedPatch->getSupportBBox();
            Vector3 patchDimensions = patchAABBox.second - patchAABBox.first;
            Vector3 patchSupportedDimensions = patchSupportedAABBox.second - patchSupportedAABBox.first;
            Vector3 controlPosition = patchAABBox.first;
            this->primitiveControlPoint->setPosition(controlPosition /*+ (selectedPatch->getDimensions() * .5f).xy()*/);
            std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
            Vector3 minPos = patchAABBox.first; //selectedPatch->getMinPosition();
            Vector3 maxPos = patchAABBox.second; //selectedPatch->getMaxPosition();
            Vector3 minimalDimensions = maxPos - minPos;
            for (auto& p : box) {
                p = p * minimalDimensions + minPos; //p * selectedPatch->getDimensions() + selectedPatch->position;
            }
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
            this->primitiveControlPoint->show();
            newSelectionIsExisting = true;

            Vector3 resolution = Vector3(20, 20, 20);
            Vector3 ratio = (patchSupportedDimensions) / (resolution /*+ Vector3(1, 1, 1)*/);
            debuggingVoxels = Matrix3<float>(resolution);
            debuggingVoxels.raiseErrorOnBadCoord = false;
            for (int x = 0; x < resolution.x; x++) {
                for (int y = 0; y < resolution.y; y++) {
                    for (int z = 0; z < resolution.z; z++) {
                        Vector3 pos(x, y, z);
                        debuggingVoxels.at(pos) = selectedPatch->evaluate((pos) * ratio + patchSupportedAABBox.first);
                    }
                }
            }
            this->debugMeshDisplayed = true;
        }
    }
    if (!newSelectionIsExisting) {
        this->patchAABBoxMesh.fromArray(std::vector<Vector3>{});
        this->patchAABBoxMesh.update();
        this->primitiveControlPoint->hide();
        this->debugMeshDisplayed = false;
    }
    Q_EMIT updated();
}

void PrimitivePatchesInterface::openPrimitiveModificationDialog(QListWidgetItem *item) {
//        int selectedIndex = this->biomeSelectionGui->currentRow();
    auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(item);
    int selectedPatchID = selectedPatchItem->ID;

    ImplicitPatch* selectedPatch = this->findPrimitiveById(selectedPatchID);
    for (const auto& p : storedPatches) {
        if (p->index == selectedPatchID) {
            selectedPatch = p;
            break;
        }
    }
    if (selectedPatch) {
        std::cout << selectedPatch->toString() << std::endl;
        PatchReplacementDialog* dialog = new PatchReplacementDialog(this, selectedPatch);
        int result = dialog->exec();
//            this->debugMeshDisplayed = false;
        if (result > 0) {
            this->updateMapWithCurrentPatch();
        } else {
            // Action cancelled
        }
        delete dialog;
    }
    updatePrimitiveList();
}

void PrimitivePatchesInterface::modifyPrimitiveHierarchy(int ID_to_move, int relatedID, HIERARCHY_TYPE relation, QDropEvent *event) {
    /*
    ImplicitPatch* toMove = this->findPrimitiveById(ID_to_move);
//        bool isLeaf = !toMove->isOperation();
    bool createCopy = event != nullptr && event->keyboardModifiers().testFlag(Qt::KeyboardModifier::ControlModifier);
    // If Ctrl is held, create a copy of the instance
    if (createCopy) {
        toMove = new ImplicitPatch(toMove->clone());
        this->storedPatches.push_back(toMove);
    }
    ImplicitPatch* related = this->findPrimitiveById(relatedID);

    if (relation == HIERARCHY_TYPE::SIBLING) {
        // Critical case : trying to make sibling with the root, just consider it as an error and make it a child
        if (related != this->mainPatch)
            related = this->naiveApproachToGetParent(related);
        relation = HIERARCHY_TYPE::CHILD;
    } else if (relation == HIERARCHY_TYPE::PARENT) {
        ImplicitPatch* tmp = toMove;
        toMove = related;
        related = tmp;
        relation = HIERARCHY_TYPE::CHILD;
    }
    ImplicitPatch* previousParent = this->naiveApproachToGetParent(toMove);
    ImplicitPatch* newGrandparent = this->naiveApproachToGetParent(related);

    std::cout << "Patch " << toMove->toString() << " (prev. child of " << previousParent->toString() << ") is now the child of " << related->toString() << std::endl;

    // Add an identity instead of "toMove"
    ImplicitPatch* identity = new ImplicitPatch;
    identity->name = "Identity";
    if (previousParent->composableA == toMove)
        previousParent->composableA = identity;
    else
        previousParent->composableB = identity;

    // Add a node to join "toMove" with the new position in tree
    ImplicitPatch* newNodeToAdd = new ImplicitPatch(related->clone());
    newNodeToAdd->composableA = toMove;
    newNodeToAdd->composableB = related;

    if (newGrandparent == nullptr){
        // No grandparent = "related" is the main patch
        newGrandparent = new ImplicitPatch(related->clone());
        this->mainPatch = newGrandparent;
    }

    if (newGrandparent->composableA == related)
        newGrandparent->composableA = newNodeToAdd;
    else
        newGrandparent->composableB = newNodeToAdd;

    this->updatePrimitiveList();
    */
}

void PrimitivePatchesInterface::openFileForNewPatch()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Open a composition"), QString::fromStdString("saved_maps/"));
    if (!q_filename.isEmpty()) {
        if (QFileInfo(q_filename).absoluteFilePath() != QFileInfo(QString::fromStdString(this->mainFilename)).absoluteFilePath()) {
            nlohmann::json content = nlohmann::json::parse(std::ifstream(q_filename.toStdString()));
            if (content.contains(ImplicitPatch::json_identifier)) {
                this->desiredPatchFromFile = ImplicitPatch::fromJson(content[ImplicitPatch::json_identifier]);
                this->desiredPatchFromFile->used_json_filename = q_filename.toStdString();
                return;
            }
        } else {
            std::cerr << "Cannot include this file in the here : current file would cause infinite recursivity" << std::endl;
        }
    }
    this->desiredPatchFromFile = nullptr;
}

ImplicitPatch* PrimitivePatchesInterface::createPatchFromParameters(Vector3 position, ImplicitPatch *replacedPatch)
{
    if (desiredPatchFromFile != nullptr) {
        ImplicitPatch* patch = desiredPatchFromFile;
        desiredPatchFromFile = nullptr;
//        return patch;
        // Set the center of the patch at mouse position
        auto [BBoxMin, BBoxMax] = patch->getBBox();
        Vector3 center = (BBoxMin + BBoxMax) * .5f;
        Vector3 translation = position.xy() - center.xy();
        // Create a translation operation to move the new patch
        ImplicitUnaryOperator* translate = new ImplicitUnaryOperator;
        translate->translate(translation);
        translate->composableA = patch;
        // Return the translation operation
        translate->updateCache();
        return translate;
    } else {
        ImplicitPrimitive* patch = (ImplicitPrimitive*) ImplicitPatch::createPredefinedShape(currentShapeSelected, functionSize, selectedSigma);
        patch->position = position.xy();
        patch->setDimensions(functionSize);
        patch->predefinedShape = this->currentShapeSelected;
    //    TerrainMaterial m;
    //    m.density = selectedDensity;
        patch->material = this->selectedTerrainType;
        patch->name = toCapitalize(stringFromPredefinedShape(currentShapeSelected));

        return patch;
    }
    /*
    ImplicitPatch* patch;
    if (replacedPatch != nullptr)
        patch = replacedPatch;
    else
        patch = new ImplicitPatch;

    Vector3 pos = position - this->functionSize.xy() * .5f;
    Vector3 useless;
    Vector3 size = this->functionSize;

    // Just to init positions and dimensions :
    *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createBlockFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
    patch->shape = this->currentShapeSelected;
    patch->densityValue = this->selectedDensity;
    patch->sigmaValue = this->selectedSigma;
    patch->name = stringFromPredefinedShape(patch->shape); //(currentShapeSelected == ImplicitPatch::PredefinedShapes::BLOCK ? "Block" : currentShapeSelected == ImplicitPatch::PredefinedShapes::GAUSSIAN ? "Gaussian" : "Sphere");
    patch->defineFunctionsBasedOnPredefinedShape();
    return patch;
    */
}

ImplicitPatch *PrimitivePatchesInterface::createOperationPatchFromParameters(ImplicitPatch *composableA, ImplicitPatch *composableB, ImplicitPatch *replacedPatch)
{
    ImplicitOperator* operation = new ImplicitOperator;
    operation->composableA = composableA;
    operation->composableB = composableB;
    operation->composeFunction = currentOperation;
    operation->positionalB = currentPositioning;
    operation->blendingFactor = selectedBlendingFactor;
    operation->withIntersectionOnB = applyIntersection;

    operation->name = toCapitalize(stringFromCompositionOperation(currentOperation));
    operation->updateCache();
    return operation;
    /*
    ImplicitPatch* operation = nullptr;

    if (currentOperation == ImplicitPatch::CompositionFunction::STACK) {
        operation = new ImplicitPatch(ImplicitPatch::createStack(composableA, composableB, currentPositioning));
    } else if (currentOperation == ImplicitPatch::CompositionFunction::BLEND) {
        operation = new ImplicitPatch(ImplicitPatch::createBlending(composableA, composableB, currentPositioning));
        operation->blendingFactor = this->selectedBlendingFactor;
    } else if (currentOperation == ImplicitPatch::CompositionFunction::REPLACE) {
        operation = new ImplicitPatch(ImplicitPatch::createReplacement(composableA, composableB, currentPositioning));
    }
    operation->name = stringFromCompositionOperation(operation->compositionOperation);
    return operation;
    */
}

ImplicitPatch* PrimitivePatchesInterface::findPrimitiveById(int ID)
{
    ImplicitPatch* patch = nullptr;
    for (const auto& p : storedPatches) {
        if (p->index == ID) {
            patch = p;
            break;
        }
    }
    return patch;
}

ImplicitPatch* PrimitivePatchesInterface::naiveApproachToGetParent(ImplicitPatch* child)
{
    std::vector<ImplicitPatch*> waitingPatches = {this->mainPatch};
    std::set<ImplicitPatch*> visitedPatches;
    while (!waitingPatches.empty()) {
        auto current = waitingPatches.back();
        waitingPatches.pop_back();
        if (visitedPatches.find(current) != visitedPatches.end())
            continue;
        visitedPatches.insert(current);
        // Check children if they exist
        ImplicitOperator* currentAsOperator = dynamic_cast<ImplicitOperator*>(current);
        if (currentAsOperator == nullptr)
            continue;
        if (currentAsOperator->composableA) {
            if (currentAsOperator->composableA == child)
                return current;
//            if (visitedPatches.find(current->composableA) == visitedPatches.end())
            waitingPatches.push_back(currentAsOperator->composableA);
//            visitedPatches.insert(current->composableA);
        }
        if (currentAsOperator->composableB) {
            if (currentAsOperator->composableB == child)
                return current;
//            if (visitedPatches.find(current->composableB) == visitedPatches.end())
            waitingPatches.push_back(currentAsOperator->composableB);
//            visitedPatches.insert(current->composableB);
        }
    }
    return nullptr;
}


void PrimitivePatchesInterface::updateFunctionSize()
{

}

void PrimitivePatchesInterface::updatePrimitiveList()
{
    if (this->isVisible()) {
        primitiveSelectionGui->clear();
        this->storedPatches.clear();
    //    this->cleanPatch(this->mainPatch);
        bool patchesCanBeLocked = this->mainFilename != "";
        std::vector<std::tuple<ImplicitPatch*, int, bool>> waitingPatches = { {this->mainPatch, 0, false} };

        this->nbPrimitives = 0;
        this->nbBinaryOperators = 0;
        this->nbUnaryOperators = 0;

        while (!waitingPatches.empty()) {
            auto [current, level, locked] = waitingPatches.back();
            bool childrenAreLocked = locked || (current->used_json_filename != "" && patchesCanBeLocked && current->used_json_filename != this->mainFilename);
            waitingPatches.pop_back();
            primitiveSelectionGui->addItem(new HierarchicalListWidgetItem((locked ? "*" : "") + current->toString(), current->index, level));
            this->storedPatches.push_back(current);
            ImplicitOperator* currentAsOperator = dynamic_cast<ImplicitOperator*>(current);
            if (currentAsOperator) {
                if (currentAsOperator->composableA) {
                    waitingPatches.push_back({currentAsOperator->composableA, level+1, childrenAreLocked});
                }
                if (currentAsOperator->composableB) {
                    waitingPatches.push_back({currentAsOperator->composableB, level+1, childrenAreLocked});
                }
            }

            // Counting patches
            if (dynamic_cast<ImplicitPrimitive*>(current) != nullptr)
                nbPrimitives ++;
            else if (dynamic_cast<ImplicitUnaryOperator*>(current) != nullptr)
                nbUnaryOperators ++;
            else
                nbBinaryOperators ++;
        }
        if (nbPrimitives > 1) // If there is more than just the Identity
            this->savePatchesAsFile(this->mainFilename);
    }
}

void PrimitivePatchesInterface::cleanPatch(ImplicitPatch *_patch)
{
    ImplicitOperator* patch = dynamic_cast<ImplicitOperator*>(_patch);
    if (!patch) // Not a composition ? Nothing to do
        return;

    // Check if it only contains Identity patches. In this case, transform it into an Identity function
    if (patch->composableA && patch->composableB) {
        if (patch->composableA->name == "Identity" && patch->composableB->name == "Identity") {
            patch->name = "Identity";
            patch->composeFunction = ImplicitPatch::CompositionFunction::NONE;
//            patch->dimension = Vector3(0, 0, 0);
            patch->composableA = nullptr;
            patch->composableB = nullptr;
        }
    }

    if (patch->composableA && patch->composableB) {
        // If it's a composite, clean children before, so we can use their clean data
        cleanPatch(patch->composableA);
        cleanPatch(patch->composableB);

        // Recompute the position
        /*if (patch->composableA->isIdentity()) { // If child A is useless, just use child B's data
            patch->position = patch->composableB->position;
            patch->dimension = patch->composableB->dimension;
        } else if (patch->composableB->isIdentity()) { // If child B is useless, just use child A's data
            patch->position = patch->composableA->position;
            patch->dimension = patch->composableA->dimension;
        } else { // Normal case, use child A and B's data
            patch->position = Vector3::min(patch->composableA->position, patch->composableB->position);
            patch->dimension = Vector3::max(patch->composableA->position + patch->composableA->dimension, patch->composableB->position + patch->composableB->dimension) - patch->position;
        }*/
    } else {
        // Nothing to do (?)
    }
}

void PrimitivePatchesInterface::displayDebuggingVoxels()
{
    std::vector<float> isoValues = {0.9f, 0.5f, 0.2f};
//    std::vector<float> isoValues = {0.9f, 0.8f, 0.65f, 0.5f, 0.35f, 0.2f, 0.1f};

    std::vector<Vector3> positions(debuggingVoxels.size());
    for (size_t i = 0; i < positions.size(); i++) {
        positions[i] = debuggingVoxels.getCoordAsVector3(i);
    }
    debuggingVoxelsMesh.fromArray(positions);
    debuggingVoxelsMesh.update();
    GlobalsGL::f()->glBindVertexArray(debuggingVoxelsMesh.vao);
    debuggingVoxelsMesh.shader->setTexture3D("dataFieldTex", 0, debuggingVoxels + .5f);
    debuggingVoxelsMesh.shader->setInt("dataFieldTex", 0);
    debuggingVoxelsMesh.shader->setInt("edgeTableTex", 1);
    debuggingVoxelsMesh.shader->setInt("triTableTex", 2);
    debuggingVoxelsMesh.shader->setFloat("isolevel", 0.f);
    debuggingVoxelsMesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    debuggingVoxelsMesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
    debuggingVoxelsMesh.shader->setBool("useMarchingCubes", true);
    //Edge Table texture//
    //This texture store the 256 different configurations of a marching cube.
    //This is a table accessed with a bitfield of the 8 cube edges states
    //(edge cut by isosurface or totally in or out).
    //(cf. MarchingCubes.cpp)

    GLuint edgeTableTex, triTableTex;
    GlobalsGL::f()->glGenTextures(1, &edgeTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, edgeTableTex);
    //Integer textures must use nearest filtering mode

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 256, 1, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::cubeEdges));

    //Triangle Table texture//
    //This texture store the vertex index list for
    //generating the triangles of each configurations.
    //(cf. MarchingCubes.cpp)

    glGenTextures(1, &triTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE2);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, triTableTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 16, 256, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::triangleTable));


    GlobalsGL::f()->glActiveTexture(GL_TEXTURE0);

    // Scale the debugging mesh to the right size
    Vector3 scaleToDisplayPatch = currentlyManipulatedPatch->getSupportDimensions() / debuggingVoxels.getDimensions(); //currentlyManipulatedPatch->getDimensions() / debuggingVoxels.getDimensions();
    debuggingVoxelsMesh.shader->setVector("scale", scaleToDisplayPatch);

    // Translate the debugging mesh to the right position
    Vector3 positionToDisplayPatch = currentlyManipulatedPatch->getSupportBBox().first; // currentlyManipulatedPatch->getBBox().first; //currentlyManipulatedPatch->position - currentlyManipulatedPatch->getDimensions() - Vector3(1.f, 1.f, 1.f) * scaleToDisplayPatch;
    debuggingVoxelsMesh.shader->setFloat("offsetX", positionToDisplayPatch.x);
    debuggingVoxelsMesh.shader->setFloat("offsetY", positionToDisplayPatch.y);
    debuggingVoxelsMesh.shader->setFloat("offsetZ", positionToDisplayPatch.z);

    // Ignore parameters to hide some voxels
    debuggingVoxelsMesh.shader->setVector("min_vertice_positions", Vector3::min());
    debuggingVoxelsMesh.shader->setVector("max_vertice_positions", Vector3::max());
    debuggingVoxelsMesh.shader->setFloat("min_isolevel", -1000.f);
    debuggingVoxelsMesh.shader->setFloat("max_isolevel", 1000.f);

    for (size_t i = 0; i < isoValues.size(); i++) {
        float iso = isoValues[i];
        Vector3 color = HSVtoRGB(i / float(isoValues.size()), 1.f, 1.f);
        debuggingVoxelsMesh.shader->setVector("color", std::vector<float> {color.x, color.y, color.z, .3f});
        debuggingVoxelsMesh.shader->setFloat("isolevel", iso);

        // display the mesh
        debuggingVoxelsMesh.display(GL_POINTS);
    }
}









PatchReplacementDialog::PatchReplacementDialog(PrimitivePatchesInterface* caller, ImplicitPatch *patchToModify)
    : QDialog(), caller(caller), patch(patchToModify)
{
    bool patchIsOperation = (dynamic_cast<ImplicitOperator*>(patch) != nullptr);

    QVBoxLayout* layout = new QVBoxLayout();

    if (patchIsOperation) {
        ImplicitOperator* patchAsOperator = dynamic_cast<ImplicitOperator*>(patch);
        QRadioButton* stackingButton = new QRadioButton("Stack");
        QRadioButton* blendingButton = new QRadioButton("Blend");
        QRadioButton* replacingButton = new QRadioButton("Replace");

        QRadioButton* abovePosButton = new QRadioButton("Above");
        QRadioButton* insideTopPosButton = new QRadioButton("Inside top");
        QRadioButton* insideBottomPosButton = new QRadioButton("Inside bottom");
        QRadioButton* fixedPosButton = new QRadioButton("Fixed");

        FancySlider* blendingFactorSlider = new FancySlider(Qt::Orientation::Horizontal, 0.1f, 10.f, .1f);
        QCheckBox* applyIntersectionButton = new QCheckBox("Intersection");

        layout->addWidget(createSliderGroup("Blend factor", blendingFactorSlider));
        layout->addWidget(createHorizontalGroup({createVerticalGroup({
                                                  stackingButton,
                                                  blendingButton,
                                                  replacingButton,
                                              }),
                                                 createVerticalGroup({
                                                     abovePosButton,
                                                     insideTopPosButton,
                                                     insideBottomPosButton,
                                                     fixedPosButton
                                                 })
                                                }));
        layout->addWidget(applyIntersectionButton);

        stackingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::STACK);
        blendingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::BLEND);
        replacingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::REPLACE);

        abovePosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::ABOVE);
        insideTopPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::INSIDE_TOP);
        insideBottomPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM);
        fixedPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::FIXED_POS);

        blendingFactorSlider->setfValue(patchAsOperator->blendingFactor);
        applyIntersectionButton->setChecked(patchAsOperator->withIntersectionOnB);

        QObject::connect(stackingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::STACK; });
        QObject::connect(blendingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::BLEND; });
        QObject::connect(replacingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::REPLACE; });

        QObject::connect(abovePosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::ABOVE; });
        QObject::connect(insideTopPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::INSIDE_TOP; });
        QObject::connect(insideBottomPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::INSIDE_BOTTOM; });
        QObject::connect(fixedPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::FIXED_POS; });
        QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [=](float newVal) {
            patchAsOperator->blendingFactor = newVal;
        });
        QObject::connect(applyIntersectionButton, &QCheckBox::toggled, this, [=](bool checked) { patchAsOperator->withIntersectionOnB = checked; });


    } else {
        ImplicitPrimitive* patchAsPrimitive = dynamic_cast<ImplicitPrimitive*>(patch);
        QRadioButton* createSphereButton = new QRadioButton("Sphere");
        QRadioButton* createBlockButton = new QRadioButton("Block");
        QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
        QRadioButton* createRockButton = new QRadioButton("Rock");
        QRadioButton* createMountainButton = new QRadioButton("Mountain");
        QRadioButton* createDuneButton = new QRadioButton("Dune");
        QRadioButton* createBasinButton = new QRadioButton("Basin");
        QRadioButton* createCaveButton = new QRadioButton("Cave");
        QRadioButton* createArchButton = new QRadioButton("Arch");
        QRadioButton* createNoise2DButton = new QRadioButton("Noise");

//        FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

        QRadioButton* airDensityCheckbox = new QRadioButton("Air");
        QRadioButton* waterDensityCheckbox = new QRadioButton("Water");
        QRadioButton* coralDensityCheckbox = new QRadioButton("Coral");
        QRadioButton* sandDensityCheckbox = new QRadioButton("Sand");
        QRadioButton* dirtDensityCheckbox = new QRadioButton("Dirt");
        QRadioButton* rockDensityCheckbox = new QRadioButton("Rock");
        QRadioButton* bedrockDensityCheckbox = new QRadioButton("Bedrock");

        FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 20.f, .1f);

        layout->addWidget(createMultiColumnGroup({
                                                        createSphereButton,
                                                        createBlockButton,
                                                        createGaussianButton,
                                                        createRockButton,
                                                        createMountainButton,
                                                        createDuneButton,
                                                        createBasinButton,
                                                        createCaveButton,
                                                        createArchButton,
                                                        createNoise2DButton
                                              }));
//        layout->addWidget(createSliderGroup("Density", densitySlider));
        layout->addWidget(createMultiColumnGroup({
                                                    airDensityCheckbox,
                                                    waterDensityCheckbox,
                                                    coralDensityCheckbox,
                                                    sandDensityCheckbox,
                                                    dirtDensityCheckbox,
                                                    rockDensityCheckbox,
                                                    bedrockDensityCheckbox
                                                }, 3));
        layout->addWidget(createMultipleSliderGroup({
                                                  {"Width / radius", widthSlider},
                                                        {"Depth", depthSlider},
                                                        {"Height", heightSlider},
                                                        {"Sigma", sigmaSlider}
                                              }));

        createSphereButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Sphere);
        createBlockButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Block);
        createGaussianButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Gaussian);
        createRockButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Rock);
        createMountainButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Mountain);
        createDuneButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Dune);
        createBasinButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Basin);
        createCaveButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Cave);
        createArchButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Arch);
        createNoise2DButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Noise2D);

//        densitySlider->setfValue(patchAsPrimitive->material.density);
//        TerrainTypes currentMat = LayerBasedGrid::materialFromDensity(patchAsPrimitive->material.density);
        airDensityCheckbox->setChecked(patchAsPrimitive->material == AIR);
        waterDensityCheckbox->setChecked(patchAsPrimitive->material == WATER);
        coralDensityCheckbox->setChecked(patchAsPrimitive->material == CORAL);
        sandDensityCheckbox->setChecked(patchAsPrimitive->material == SAND);
        dirtDensityCheckbox->setChecked(patchAsPrimitive->material == DIRT);
        rockDensityCheckbox->setChecked(patchAsPrimitive->material == ROCK);
        bedrockDensityCheckbox->setChecked(patchAsPrimitive->material == BEDROCK);

        widthSlider->setfValue(patchAsPrimitive->dimensions.x);
        heightSlider->setfValue(patchAsPrimitive->dimensions.z);
        depthSlider->setfValue(patchAsPrimitive->dimensions.y);
        sigmaSlider->setfValue(patchAsPrimitive->parametersProvided[0]);

        QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Sphere; });
        QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Block; });
        QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Gaussian; });
        QObject::connect(createRockButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Rock; });
        QObject::connect(createMountainButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Mountain; });
        QObject::connect(createDuneButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Dune; });
        QObject::connect(createBasinButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Basin; });
        QObject::connect(createCaveButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Cave; });
        QObject::connect(createArchButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Arch; });
        QObject::connect(createNoise2DButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Noise2D; });
//        QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [=](float newDensity) { patchAsPrimitive->material.density = newDensity; });
        QObject::connect(airDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = AIR; });
        QObject::connect(waterDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = WATER; });
        QObject::connect(coralDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = CORAL; });
        QObject::connect(sandDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = SAND; });
        QObject::connect(dirtDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = DIRT; });
        QObject::connect(rockDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = ROCK; });
        QObject::connect(bedrockDensityCheckbox, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = BEDROCK; });


        QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.x = newVal; });
        QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.z = newVal; });
        QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.y = newVal; });
        QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->parametersProvided[0] = newVal; });

    }



    cancelButton = new QPushButton("Annuler", this);
    validButton = new QPushButton("Confirmer", this);

    layout->addWidget(createHorizontalGroup({cancelButton, validButton}));
    QObject::connect(cancelButton, &QPushButton::pressed, this, &PatchReplacementDialog::cancel);
    QObject::connect(validButton, &QPushButton::pressed, this, &PatchReplacementDialog::confirm);
    setLayout(layout);
    setSizeGripEnabled(true);
}

void PatchReplacementDialog::open()
{
    QDialog::open();
}

void PatchReplacementDialog::cancel()
{
    setResult(-1);
//    this->close();
    this->done(-1);
}

void PatchReplacementDialog::confirm()
{
    patch->update();
    /*ImplicitPatch* replacement;
    if (patch->isOperation()) {
        replacement = caller->createOperationPatchFromParameters(patch->composableA, patch->composableB);
    } else {
        replacement = caller->createPatchFromParameters(this->patch->position + (this->patch->getDimensions() * .5f).xy()); // Position at the center of current patch
    }
    *this->patch = *replacement;
    delete replacement;*/
    setResult(1);
//    this->close();
    this->done(1);
}
