#include "PrimitivePatchesInterface.h"
#include "Graphics/DisplayGraphics.h"
#include "Interface/InterfaceUtils.h"
//#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"
//#include "Utils/stb_image.h"

PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("PrimitivePatchesInterface", parent)
{
    previewMesh.cullFace = false;

//    this->implicitTerrain = new ImplicitPrimitive;
//    this->implicitTerrain->name = "Identity";
//    dynamic_cast<ImplicitPrimitive*>(this->implicitTerrain)->material = WATER;


//    this->setMainFilename("saved_maps/implicit_patches.json");
    this->setMainFilename("saved_maps/trench.json");
//    this->setMainFilename("saved_maps/water_currents2.json");
//    this->setMainFilename("saved_maps/cave.json");
//    this->setMainFilename("saved_maps/rock_in_sand.json");

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

        if (parametricCurveMesh.shader != nullptr) {
            parametricCurveMesh.shader->setVector("color", std::vector<float> {1.f, 1.f, 1.f, 1.f});
            this->displayParametricCurve();
        }
    }
}

void PrimitivePatchesInterface::reloadShaders()
{
}

void PrimitivePatchesInterface::replay(nlohmann::json action)
{

}

void PrimitivePatchesInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    // Waiting for OpenGL to be available to create the shaders...
    std::string shaderPath = "src/Shaders/";
    std::string vertexFilename = shaderPath + "MarchingCubes.vert";
    std::string geometryFilename = shaderPath + "MarchingCubes.geom";
    std::string fragmentFilename = shaderPath + "no_shader.frag";

    std::string noShaderVertexFilename = shaderPath + "no_shader.vert";
    std::string noShaderFragmentFilename = shaderPath + "no_shader.frag";

    debuggingVoxelsMesh.shader = std::make_shared<Shader>(vertexFilename, fragmentFilename, geometryFilename);
    debuggingVoxelsMesh.useIndices = false;
    debuggingVoxelsMesh.cullFace = true;

    parametricCurveMesh = Mesh(std::make_shared<Shader>(noShaderVertexFilename, noShaderFragmentFilename), true, GL_LINES);
    parametricCurveMesh.useIndices = false;

    this->loadTransformationRules();
//    this->loadPatchesFromFile(this->mainFilename);
}

QLayout *PrimitivePatchesInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    QRadioButton* createSphereButton = new QRadioButton("Sphere");
    QRadioButton* createBlockButton = new QRadioButton("Block");
    QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
    QRadioButton* createCylinderButton = new QRadioButton("Cylinder");
    QRadioButton* createRockButton = new QRadioButton("Rock");
    QRadioButton* createMountainButton = new QRadioButton("Mountain");
    QRadioButton* createDuneButton = new QRadioButton("Dune");
    QRadioButton* createBasinButton = new QRadioButton("Basin");
    QRadioButton* createCaveButton = new QRadioButton("Cave");
    QRadioButton* createArchButton = new QRadioButton("Arch");
    QRadioButton* createNoise2DButton = new QRadioButton("Noise");
    QRadioButton* createMountainChainButton = new QRadioButton("Mountains");
    QRadioButton* createPolygonButton = new QRadioButton("Polygon");
    QRadioButton* createTunnelButton = new QRadioButton("Tunnel");
    QRadioButton* createFromFileButton = new QRadioButton("...");

//    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");
    QRadioButton* oneSideBlendButton = new QRadioButton("1-blend");

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
    QPushButton* addDistoButton = new QPushButton("Distord");
    QPushButton* addSpreadButton = new QPushButton("Spread");

    QPushButton* resetButton = new QPushButton("Reset");

    QLabel* selectedFilenameLabel = new QLabel((this->mainFilename == "" ? "No file selected" : QString::fromStdString(this->mainFilename).split("/").back().split("\\").back()));
    QPushButton* fileSelectionButton = new QPushButton("...");
    QCheckBox* enableHotreloadButton = new QCheckBox("Hot reloading ?");

    QLabel* patchesCounterLabel = new QLabel(QString::fromStdString("Prim: " + std::to_string(nbPrimitives) + ", Unary: " + std::to_string(nbUnaryOperators) + ", Binary: " + std::to_string(nbBinaryOperators) + ", N-ary: " + std::to_string(nbNaryOperators)));

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
                                                    createCylinderButton,
                                                    createRockButton,
                                                    createMountainButton,
                                                    createDuneButton,
                                                    createBasinButton,
                                                    createCaveButton,
                                                    createArchButton,
                                                    createNoise2DButton,
                                                    createMountainChainButton,
                                                    createPolygonButton,
                                                    createTunnelButton,
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
                                              oneSideBlendButton
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
    layout->addWidget(createHorizontalGroup({addNoiseButton, addDistoButton, addSpreadButton}));
    layout->addWidget(resetButton);
    layout->addWidget(primitiveSelectionGui);

    layout->addWidget(createHorizontalGroup({selectedFilenameLabel, fileSelectionButton}));
    layout->addWidget(enableHotreloadButton);
    layout->addWidget(patchesCounterLabel);

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Sphere); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Block); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Gaussian); });
    QObject::connect(createCylinderButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Cylinder); });
    QObject::connect(createRockButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Rock); });
    QObject::connect(createMountainButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Mountain); });
    QObject::connect(createDuneButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Dune); });
    QObject::connect(createBasinButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Basin); });
    QObject::connect(createCaveButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Cave); });
    QObject::connect(createArchButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Arch); });
    QObject::connect(createNoise2DButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Noise2D); });
    QObject::connect(createMountainChainButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::MountainChain); });
    QObject::connect(createPolygonButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Polygon); });
    QObject::connect(createTunnelButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::ParametricTunnel); });
    QObject::connect(createFromFileButton, &QRadioButton::toggled, this, [=](bool checked) {
        if (checked) {
            this->setSelectedShape(ImplicitPatch::PredefinedShapes::ImplicitHeightmap);
            this->openFileForNewPatch();
//            createSphereButton->setChecked(true);
        }
    });
//    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [=](float newDensity) { this->selectedDensity = newDensity; });

    QObject::connect(stackingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setCurrentOperation(ImplicitPatch::CompositionFunction::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setCurrentOperation(ImplicitPatch::CompositionFunction::REPLACE); });
    QObject::connect(oneSideBlendButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setCurrentOperation(ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND); });

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
    QObject::connect(addDistoButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::addDistortionOnSelectedPatch);
    QObject::connect(addSpreadButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::addSpreadOnSelectedPatch);

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
    createCylinderButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Cylinder);
    createRockButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Rock);
    createMountainButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Mountain);
    createDuneButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Dune);
    createBasinButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Basin);
    createCaveButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Cave);
    createArchButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Arch);
    createMountainChainButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::MountainChain);
    createPolygonButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Polygon);
    createTunnelButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::ParametricTunnel);
    createNoise2DButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Noise2D);

    stackingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::STACK);
    blendingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::BLEND);
    replacingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::REPLACE);
    oneSideBlendButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND);

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
    this->primitiveControlPoint->allowAllAxisRotations(true);
    this->primitiveControlPoint->allowAllAxisTranslation(true);
    this->primitiveControlPoint->hide();

    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::modified, this, [=](){
        // Get patch being manipulated
        if (this->selectedPatch() != nullptr) {
            // Display modified AABBox
            auto AABBox = this->selectedPatch()->getBBox();
            Vector3 minPos = AABBox.min();
            Vector3 maxPos = AABBox.max();
            Vector3 dim = maxPos - minPos;
            std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
            for (auto& p : box) {
                p = (p * dim) + /*minPos + */primitiveControlPoint->getPosition(); //p * currentlySelectedPatch->getDimensions() + (primitiveControlPoint->getPosition() - currentlySelectedPatch->getDimensions().xy() * .5f);//p * currentlySelectedPatch->getDimensions() + (primitiveControlPoint->getPosition() - (currentlySelectedPatch->getDimensions() * .5f).xy());
            }
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
        }
    });
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::translationApplied, this, [=](Vector3 translation){
//        QObject::blockSignals(true);
        // Get patch being manipulated
        if (this->selectedPatch() != nullptr) {
            primitiveControlPoint->blockSignals(true);
            primitiveSelectionGui->blockSignals(true);

            ImplicitUnaryOperator* manipulatedAsUnary = dynamic_cast<ImplicitUnaryOperator*>(this->selectedPatch());
            if (manipulatedAsUnary == nullptr) {
                // It's not an unary yet, but the parent might be
                manipulatedAsUnary = dynamic_cast<ImplicitUnaryOperator*>(this->naiveApproachToGetParent(this->selectedPatch()));
            }

            if (manipulatedAsUnary != nullptr) { // We are updating an unary operator
                manipulatedAsUnary->translate(translation); // Just update it
            } else { // Otherwise, create a new Unary operator
                ImplicitUnaryOperator* translate = new ImplicitUnaryOperator;
                translate->composableA = this->selectedPatch();
                translate->translate(translation);
                translate->name = "Translation";
                if (this->selectedPatch() == this->implicitTerrain.get()) {
                    this->implicitTerrain->composables = {translate};
                } else {
                    ImplicitOperator* parentAsOperator = dynamic_cast<ImplicitOperator*>(this->naiveApproachToGetParent(this->selectedPatch()));
                    if (parentAsOperator) {
                        if (this->selectedPatch() == parentAsOperator->composableA) {
                            parentAsOperator->composableA = translate;
                        } else {
                            parentAsOperator->composableB = translate;
                        }
                    }
                    ImplicitNaryOperator* parentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(this->naiveApproachToGetParent(this->selectedPatch()));
                    if (parentAsNaryOperator) {
                        for (auto& c : parentAsNaryOperator->composables)
                            if (this->selectedPatch() == c)
                                c = translate;
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
        if (this->currentlySelectedPatch != nullptr) {
            std::cout << "Modifying the terrain..." << std::endl;
//            Vector3 initialControlPosition = currentlySelectedPatch->getBBox().first;
//            Vector3 translation = primitiveControlPoint->getPosition() - initialControlPosition;
            // When released, recreate patch + update terrain
            //currentlySelectedPatch->position = primitiveControlPoint->getPosition() - (currentlySelectedPatch->getDimensions() * .5f).xy();
            ImplicitUnaryOperator* rotate = new ImplicitUnaryOperator;
            rotate->composableA = this->currentlySelectedPatch;
            rotate->rotate(rotation.x, rotation.y, rotation.z);
            rotate->name = "Rotation";
            if (currentlySelectedPatch == this->implicitTerrain) {
                this->implicitTerrain = rotate; // TODO : Not right, but for tests
            } else {
                ImplicitOperator* parent = (ImplicitOperator*)this->naiveApproachToGetParent(currentlySelectedPatch);
                if (currentlySelectedPatch == parent->composableA) {
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

void PrimitivePatchesInterface::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key::Key_Q && (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel)) {
        this->createPatchWithOperation(Vector3());
        this->parametricCurve.points.clear();
        this->parametricCurveMesh.clear();
    }
    if (e->key() == Qt::Key::Key_W) {
        this->displayPatchesTree();
    }
    if (this->selectedPatch() == nullptr) return;
    ImplicitPrimitive* patch = dynamic_cast<ImplicitPrimitive*>(this->selectedPatch());
    if (patch == nullptr) return;
    if (e->key() == Qt::Key::Key_Plus) {
        Vector3 dimensions = patch->getDimensions();
        dimensions.z += 1.f;
        patch->setDimensions(dimensions);
        patch->update();
        this->updateMapWithCurrentPatch();
    } else if (e->key() == Qt::Key::Key_Minus) {
        Vector3 dimensions = patch->getDimensions();
        dimensions.z -= 1.f;
        patch->setDimensions(dimensions);
        patch->update();
        this->updateMapWithCurrentPatch();
    }
}

ImplicitPatch *PrimitivePatchesInterface::selectedPatch()
{
    if (this->currentlySelectedPatch != nullptr && this->currentlySelectedPatch != this->implicitTerrain.get()) {
        return this->currentlySelectedPatch;
    } else if (!this->implicitTerrain->composables.empty()) {
        return this->implicitTerrain->composables[0];
    } else {
        return nullptr;
    }
}

void PrimitivePatchesInterface::mouseMovedOnMapEvent(Vector3 newPosition, TerrainModel* model)
{
    if (this->isVisible()) {
        this->currentPos = newPosition;
        this->setSelectedShape(this->currentShapeSelected, newPosition);
    }
}

void PrimitivePatchesInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel* model)
{
    if (this->isVisible() && mouseInMap) {

        if (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel) {
            this->addParametricPoint(model->getTerrainPos(mousePosInMap));
            return;
        }

        if (this->primitiveControlPoint->grabsMouse())
            return; // Don't create patches if we try to displace another patch
        this->createPatchWithOperation(model->getTerrainPos(mousePosInMap));
    }
}

void PrimitivePatchesInterface::createPatchWithOperation(Vector3 pos)
{
    ImplicitPatch* patch = this->createPatchFromParameters(pos - this->functionSize.xy() * .5f);
    ImplicitPatch* previousMain = this->selectedPatch();
    ImplicitPatch* _parent = this->naiveApproachToGetParent(previousMain);
    if (!previousMain) {
        if (this->implicitTerrain->composables.empty()) {
            previousMain = ImplicitPatch::createIdentity();
        } else {
            previousMain = this->implicitTerrain->composables.front();
        }
    }
    ImplicitPatch* operation = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;

    if (_parent == nullptr) {
        this->implicitTerrain->composables = {operation};
    } else {
        ImplicitOperator* parentAsOperator = dynamic_cast<ImplicitOperator*>(_parent);
        if (parentAsOperator) {
            if (previousMain == parentAsOperator->composableA)
                parentAsOperator->composableA = operation;
            else
                parentAsOperator->composableB = operation;
        }
        ImplicitNaryOperator* parentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(_parent);
        if (parentAsNaryOperator) {
            int iComposable = 0;
            for (iComposable = 0; iComposable < parentAsNaryOperator->composables.size(); iComposable++) {
                if (previousMain == parentAsNaryOperator->composables[iComposable])
                    break;
            }
            parentAsNaryOperator->composables[iComposable] = operation;

        }

    }
//        ImplicitPatch* newOperation = this->createOperationPatchFromParameters(previousMain, patch);
//        this->implicitTerrain = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;
    this->updateMapWithCurrentPatch();

    this->storedPatches.push_back(patch);
    this->storedPatches.push_back(operation);
    this->storedPatches.push_back(implicitTerrain.get());
    this->updatePrimitiveList();

    if (this->mainFilename != "")
        this->savePatchesAsFile(this->mainFilename);

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
    if (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::ParametricTunnel) {
        this->functionSize = Vector3(selectedSigma, selectedSigma, selectedSigma);
    }
    for (auto& vert : vertices) {
        vert *= functionSize;
        vert += newPosition;
        vert -= functionSize.xy() * .5f;
    }

//    ImplicitPrimitive* previewPatch = ImplicitPatch::createPredefinedShape(newShape, this->functionSize, this->selectedSigma);
//    previewPatch->material = this->selectedTerrainType;
//    Vector3 resolution = Vector3(20, 20, 20);
//    Vector3 ratio = (patchSupportedDimensions) / (resolution /*+ Vector3(1, 1, 1)*/);
//    debuggingVoxelsPosition = newPosition - functionSize * Vector3(1.f, 1.f, 0.f);
//    debuggingVoxelsScale = Vector3(2.f, 2.f, 2.f);
//    debuggingVoxels = previewPatch->getVoxelized(Vector3(false), debuggingVoxelsScale);
//    debuggingVoxels = debuggingVoxels.resize(debuggingVoxelsScale);
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
//    delete this->implicitTerrain;
    ImplicitPatch::currentMaxIndex = -1;
    this->implicitTerrain->composables.clear(); // = new ImplicitPrimitive;
    this->implicitTerrain->name = "Identity";
    this->updateMapWithCurrentPatch();
    this->storedPatches.push_back(implicitTerrain.get());
    this->updatePrimitiveList();
//    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    //    voxelGrid->fromIsoData();
}

void PrimitivePatchesInterface::updateMapWithCurrentPatch()
{
    this->implicitTerrain->updateCache();
//    this->layerGrid->layers = this->layerGrid->previousState;
    this->layerGrid->reset();
//    this->implicitTerrain->cleanCache();
    this->layerGrid->add(this->implicitTerrain.get()/*, SAND, false*/);
    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
//    voxelGrid->fromCachedData();
    heightmap->fromLayerGrid(*layerGrid);
    this->savePatchesAsFile(this->mainFilename);

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
    ImplicitPatch* selectedPatch = this->selectedPatch();
    ImplicitUnaryOperator* asOp = dynamic_cast<ImplicitUnaryOperator*>(selectedPatch);
    ImplicitUnaryOperator* noisePatch;

    if (asOp == nullptr) {
        asOp = dynamic_cast<ImplicitUnaryOperator*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (asOp != nullptr) {
        noisePatch = asOp;
    } else {
        noisePatch = new ImplicitUnaryOperator;
        noisePatch->composableA = selectedPatch;
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->composables = {noisePatch};
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitOperator* parentAsOperator = dynamic_cast<ImplicitOperator*>(_parent);
            if (parentAsOperator) {
                if (selectedPatch == parentAsOperator->composableA)
                    parentAsOperator->composableA = noisePatch;
                else
                    parentAsOperator->composableB = noisePatch;
            }
        }
    }
    float amplitude = this->selectedSigma * .01f;
    noisePatch->addRandomNoise(amplitude);
    noisePatch->name = "Noise";
    this->storedPatches.push_back(noisePatch);
    this->updateMapWithCurrentPatch();
    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::addDistortionOnSelectedPatch()
{
    ImplicitPatch* selectedPatch = this->selectedPatch();
    ImplicitUnaryOperator* asOp = dynamic_cast<ImplicitUnaryOperator*>(selectedPatch);
    ImplicitUnaryOperator* distortionPatch;

    if (asOp == nullptr) {
        asOp = dynamic_cast<ImplicitUnaryOperator*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (asOp != nullptr) {
        distortionPatch = asOp;
    } else {
        distortionPatch = new ImplicitUnaryOperator;
        distortionPatch->composableA = selectedPatch;
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->composables = {distortionPatch};
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitOperator* parent = dynamic_cast<ImplicitOperator*>(_parent);
            if (parent) {
                if (selectedPatch == parent->composableA)
                    parent->composableA = distortionPatch;
                else
                    parent->composableB = distortionPatch;
            }
        }
    }
    Matrix3<Vector3> disto = Matrix3<Vector3>({
//                                                  { Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), /*Vector3(0.0, 0.0, 0), */Vector3(0.0, -.5, 0)/*, Vector3(0.0, 0.5, 0) */},
//                                                  { Vector3(0.0, 0.0, 0), Vector3(0.0, -.2, 0), Vector3(0.0, -.0, 0), /*Vector3(0.0, -.3, 0), */Vector3(0.0, -.5, 0)/*, Vector3(0.0, 0.5, 0) */},
//                                                  { Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), /*Vector3(0.0, 0.0, 0), */Vector3(0.0, 0.0, 0)/*, Vector3(0.0, 0.0, 0) */},
//                                                  { Vector3(0.0, 0.0, 0), Vector3(0.0, 0.2, 0), Vector3(0.0, 0.0, 0), /*Vector3(0.0, 0.3, 0), */Vector3(0.0, 0.5, 0)/*, Vector3(0.0, -.5, 0) */},
//                                                  { Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), Vector3(0.0, 0.0, 0), /*Vector3(0.0, 0.0, 0), */Vector3(0.0, 0.5, 0)/*, Vector3(0.0, -.5, 0) */}
                                                  { Vector3(0.2, -.5, 0), Vector3(0.0, -.5, 0), Vector3(0.0, -.5, 0), Vector3(0.0, -.5, 0) },
                                                  { Vector3(0.2, 0.5, 0), Vector3(0.0, 0.5, 0), Vector3(0.0, 0.5, 0), Vector3(0.0, 0.5, 0) }
    }).flip(false, true, false);
//    std::cout << disto.displayValues() << std::endl;
    float amplitude = this->selectedSigma;
//    distortionPatch->addRandomWrap(amplitude);
    distortionPatch->addWrapFunction(disto);
    distortionPatch->name = "Wrapping";
    this->storedPatches.push_back(distortionPatch);
    this->updateMapWithCurrentPatch();
    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::addSpreadOnSelectedPatch()
{
    ImplicitPatch* selectedPatch = this->selectedPatch();
    ImplicitUnaryOperator* asOp = dynamic_cast<ImplicitUnaryOperator*>(selectedPatch);
    ImplicitUnaryOperator* distortionPatch;

    if (asOp == nullptr) {
        asOp = dynamic_cast<ImplicitUnaryOperator*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (asOp != nullptr) {
        distortionPatch = asOp;
    } else {
        distortionPatch = new ImplicitUnaryOperator;
        distortionPatch->composableA = selectedPatch;
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->composables = {distortionPatch};
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitOperator* parent = dynamic_cast<ImplicitOperator*>(_parent);
            if (parent) {
                if (selectedPatch == parent->composableA)
                    parent->composableA = distortionPatch;
                else
                    parent->composableB = distortionPatch;
            }
        }
    }
    float amplitude = this->selectedSigma / 10.f;
    distortionPatch->spread(amplitude);
    distortionPatch->name = "Spread";
    this->storedPatches.push_back(distortionPatch);
    this->updateMapWithCurrentPatch();
    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::savePatchesAsFile(std::string filename)
{
    this->lastTimeFileHasBeenModified = QDateTime::currentDateTime().addDays(1); // Just to avoid the reloading during the save
    nlohmann::json content = this->implicitTerrain->toJson();
    std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc);
    out << nlohmann::json({ {ImplicitPatch::json_identifier, content} }).dump(1, '\t');
    out.flush();
    this->lastTimeFileHasBeenModified = QDateTime::currentDateTime();
}

void PrimitivePatchesInterface::loadPatchesFromFile(std::string filename)
{
    this->loadTransformationRules();
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    nlohmann::json json_content = nlohmann::json::parse(content);
    this->lastTimeFileHasBeenModified = QFileInfo(QString::fromStdString(filename)).lastModified();
    if (json_content.contains(ImplicitPatch::json_identifier)) {
        this->implicitTerrain->composables = {ImplicitPatch::fromJson(json_content[ImplicitPatch::json_identifier])};
        this->updateMapWithCurrentPatch();
        this->updatePrimitiveList();
    } else {
        std::cerr << "No patches defined in file " << filename << " (no '" << ImplicitPatch::json_identifier << "' found in file)" << std::endl;
    }
}

void PrimitivePatchesInterface::hotReloadFile()
{
    if (/*this->isVisible() && */this->enableHotReloading) {
        bool needReload = false;
        if (this->mainFilename != "") {
            QFileInfo infos(QString::fromStdString(this->mainFilename));
            auto lastModifTime = infos.lastModified();
            if (lastModifTime > this->lastTimeFileHasBeenModified || !this->lastTimeFileHasBeenModified.isValid()) {
                needReload = true;
            }
        }
        if (this->rulesFilename != "") {
            QFileInfo infos(QString::fromStdString(this->rulesFilename));
            auto lastModifTime = infos.lastModified();
            if (lastModifTime > this->lastTimeFileHasBeenModified || !this->lastTimeFileHasBeenModified.isValid()) {
                needReload = true;
            }
        }
        for (auto file : allSubfiles) {
            QFileInfo infos(QString::fromStdString(file));
            auto lastModifTime = infos.lastModified();
            if (lastModifTime > this->lastTimeFileHasBeenModified || !this->lastTimeFileHasBeenModified.isValid()) {
                needReload = true;
                break;
            }
        }
        if (needReload)
            this->loadPatchesFromFile(this->mainFilename);
    }
}

void PrimitivePatchesInterface::updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous) {
    bool newSelectionIsExisting = false;
    this->currentlySelectedPatch = nullptr;
    if (current) {
        auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(current);
        int selectedPatchID = selectedPatchItem->ID;
        this->currentlySelectedPatch = this->findPrimitiveById(selectedPatchID);
        if (currentlySelectedPatch != nullptr) {
            auto patchAABBox = currentlySelectedPatch->getBBox();
            auto patchSupportedAABBox = currentlySelectedPatch->getSupportBBox();
            Vector3 patchDimensions = patchAABBox.dimensions();
            Vector3 patchSupportedDimensions = patchSupportedAABBox.dimensions();
            Vector3 controlPosition = patchAABBox.min();
            this->primitiveControlPoint->setPosition(controlPosition /*+ (selectedPatch->getDimensions() * .5f).xy()*/);
            std::vector<Vector3> box = CubeMesh::createTriangles(patchAABBox);
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
            this->primitiveControlPoint->show();
            newSelectionIsExisting = true;

            Vector3 resolution = Vector3(20, 20, 20);
            Vector3 ratio = (patchSupportedDimensions) / (resolution /*+ Vector3(1, 1, 1)*/);
            this->debuggingVoxelsPosition = this->selectedPatch()->getSupportBBox().min();
            this->debuggingVoxelsScale = ratio;
            debuggingVoxels = Matrix3<float>(resolution);
            debuggingVoxels.raiseErrorOnBadCoord = false;
            for (int x = 0; x < resolution.x; x++) {
                for (int y = 0; y < resolution.y; y++) {
                    for (int z = 0; z < resolution.z; z++) {
                        Vector3 pos(x, y, z);
                        debuggingVoxels.at(pos) = currentlySelectedPatch->evaluate((pos) * ratio + patchSupportedAABBox.min());
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
        if (related != this->implicitTerrain)
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
        this->implicitTerrain = newGrandparent;
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
    std::string filename = q_filename.toStdString();
    if (!q_filename.isEmpty()) {
        if (QFileInfo(q_filename).absoluteFilePath() != QFileInfo(QString::fromStdString(this->mainFilename)).absoluteFilePath()) {
            this->desiredPatchFilename = q_filename.toStdString();
//            if (q_filename.endsWith(".json")) {
//                nlohmann::json content = nlohmann::json::parse(std::ifstream(q_filename.toStdString()));
//                if (content.contains(ImplicitPatch::json_identifier)) {
//                    this->desiredPatchFromFile = ImplicitPatch::fromJson(content[ImplicitPatch::json_identifier]);
//                    this->desiredPatchFromFile->used_json_filename = q_filename.toStdString();
//                    return;
//                }
//            } else if (q_filename.endsWith(".png")) {
//                this->desiredPatchFromFile = ImplicitPrimitive::fromHeightmap(filename, Vector3(this->selectedWidth, this->selectedDepth, this->selectedHeight));
//                this->desiredPatchFromFile->used_json_filename = q_filename.toStdString();
//                return;
//            }
        } else {
            std::cerr << "Cannot include this file in the here : current file would cause infinite recursivity" << std::endl;
        }
    }
    this->desiredPatchFromFile = nullptr;
}

void PrimitivePatchesInterface::setMainFilename(std::string newMainFilename)
{
    this->mainFilename = newMainFilename;
}

void PrimitivePatchesInterface::findAllSubfiles()
{
    this->allSubfiles.clear();
    std::vector<ImplicitPatch*> pendingPatches = {this->implicitTerrain.get()};
    while (!pendingPatches.empty()) {
        ImplicitPatch* current = pendingPatches.back();
        pendingPatches.pop_back();

        if (current->used_json_filename != "")
            allSubfiles.push_back(current->used_json_filename);

        ImplicitOperator* asOp = dynamic_cast<ImplicitOperator*>(current);
        if (asOp != nullptr) {
            if (asOp->composableA != nullptr)
                pendingPatches.push_back(asOp->composableA);
            if (asOp->composableB != nullptr)
                pendingPatches.push_back(asOp->composableB);
        }
    }
}

void PrimitivePatchesInterface::loadTransformationRules()
{
    std::cout << "Reading and interpreting file " << this->rulesFilename << "..." << std::endl;
    std::vector<std::pair<std::map<TerrainTypes, float>, std::map<TerrainTypes, float>>> rules;
    std::ifstream file(this->rulesFilename);
    std::string sline;
//    bool takesInputs = true;
    while (std::getline(file, sline)) {
        if (sline.empty() || sline[0] == '#') continue; // Comments with "#"
        std::istringstream line(sline);
        std::map<TerrainTypes, float> inputs, outputs;
        std::string value, word, operation;
        // Get inputs
        while (operation != "=") {
            line >> value;
            line >> word;
            line >> operation;
            inputs[materialFromString(word)] = std::stof(value);
        }
        // Get outputs
        while (true) {
            line >> value;
            line >> word;
            outputs[materialFromString(word)] = std::stof(value);
            if (!(line >> operation)) break;
        }
        rules.push_back({inputs, outputs});
    }
    this->layerGrid->transformationRules = rules;
}

void PrimitivePatchesInterface::addParametricPoint(Vector3 point)
{
    this->parametricCurve.points.push_back(point);
    auto path = parametricCurve.getPath(100);
    std::vector<Vector3> allPointsDoubled;
    for (size_t i = 0; i < path.size() - 1; i++) {
        allPointsDoubled.push_back(path[i] + Vector3(0, 0, 2.f));
        allPointsDoubled.push_back(path[i + 1] + Vector3(0, 0, 2.f));
    }
    parametricCurveMesh.fromArray(allPointsDoubled);
}

void PrimitivePatchesInterface::displayParametricCurve()
{
    parametricCurveMesh.display(GL_LINES, 5.f);
}

void PrimitivePatchesInterface::displayPatchesTree()
{
    Plotter* plt = Plotter::getInstance();

    std::vector<std::pair<ImplicitPatch*, std::vector<int>>> allPatchesAndDirections;
    std::vector<std::tuple<ImplicitPatch*, std::vector<int>, int>> queueWithDirections = {{implicitTerrain.get(), {}, -1}};
    std::vector<Vector3> positions;
    std::vector<std::string> labels;
    std::vector<std::pair<int, int>> edges;

    int maxDepth = 0;
    while (!queueWithDirections.empty()) {
        auto [current, directions, parentIndex] = queueWithDirections.front();
        queueWithDirections.erase(queueWithDirections.begin());
        allPatchesAndDirections.push_back({current, directions});
        int currentIndex = int(allPatchesAndDirections.size()) - 1;
        edges.push_back({parentIndex, currentIndex});
        std::string label = current->name;
        ImplicitPrimitive* asPrim = dynamic_cast<ImplicitPrimitive*>(current);
        if (asPrim != nullptr)
            label += "\n(" + toCapitalize(stringFromMaterial(asPrim->material)) + ")";
        labels.push_back(label);
        maxDepth = std::max(maxDepth, int(directions.size()));

        ImplicitOperator* asOp = dynamic_cast<ImplicitOperator*>(current);
        if (asOp != nullptr) {
            if (asOp->composableA != nullptr) {
                auto directionsForA = directions;
                if (asOp->composableB == nullptr) { // A is unique child
                    directionsForA.push_back(0);
                } else {
                    directionsForA.push_back(-1);
                }
                queueWithDirections.push_back({asOp->composableA, directionsForA, currentIndex});
            }
            if (asOp->composableB != nullptr) {
                auto directionsForB = directions;
                directionsForB.push_back(1);
                queueWithDirections.push_back({asOp->composableB, directionsForB, currentIndex});
            }
        }
    }

    for (const auto& [patch, directions] : allPatchesAndDirections) {
        float x = 0.f;
        float y = 0.f;
        for (size_t depth = 0; depth < directions.size(); depth++) {
            x += float(directions[depth]) * std::pow(0.9f, depth);
            y -= 1.f;
        }
        positions.push_back(Vector3(x, y));
    }
    for (auto& [e0, e1] : edges) {
        plt->addPlot({positions[e0], positions[e1]});
    }
    plt->addScatter(positions, "", labels);
    plt->draw();
    plt->show();
//    QObject::connect(plt, &Plotter::finished, this, [=](int result) {
//        std::cout << "Closed with result " << result << std::endl;
//        plt->close();
//    });
}

ImplicitPatch* PrimitivePatchesInterface::createPatchFromParameters(Vector3 position, ImplicitPatch *replacedPatch)
{
//    if (desiredPatchFromFile != nullptr) {
    if (this->currentShapeSelected == ImplicitPatch::ImplicitHeightmap) {
        ImplicitPatch* patch; // = desiredPatchFromFile;

        if (endsWith(desiredPatchFilename, ".json")) {
            // JSON process
            nlohmann::json content = nlohmann::json::parse(std::ifstream(desiredPatchFilename));
            patch = ImplicitPatch::fromJson(content[ImplicitPatch::json_identifier]);
        } else {
            // Consider it's an image
            ImplicitPrimitive* patchAsPrimitive = ImplicitPrimitive::fromHeightmap(desiredPatchFilename, Vector3(this->selectedWidth, this->selectedDepth, this->selectedHeight));
            patchAsPrimitive->material = this->selectedTerrainType;
            patch = patchAsPrimitive;
        }
        patch->used_json_filename = desiredPatchFilename;

        desiredPatchFromFile = nullptr;
//        return patch;
        // Set the center of the patch at mouse position
        auto [BBoxMin, BBoxMax] = patch->getBBox();
//        Vector3 center = (BBoxMin + BBoxMax) * .5f;
        Vector3 translation = position;
        // Create a translation operation to move the new patch
        ImplicitUnaryOperator* translate = new ImplicitUnaryOperator;
        translate->translate(translation);
        translate->composableA = patch;
        // Return the translation operation
        translate->updateCache();
        return translate;
    } else {
        Vector3 funcSize = this->functionSize * layerGrid->scaling;
        if (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel) {
            auto AABBox = this->parametricCurve.AABBox();
            position = std::get<0>(AABBox) - funcSize * .5f;
            funcSize = (std::get<1>(AABBox) + funcSize * .5f) - position;
        }
        BSpline translatedSpline = this->parametricCurve;
        for (auto& p : translatedSpline.points)
            p -= position;
        ImplicitPrimitive* primitive = (ImplicitPrimitive*) ImplicitPatch::createPredefinedShape(currentShapeSelected, funcSize, selectedSigma, translatedSpline);
        primitive->position = position.xy();
        primitive->setDimensions(funcSize);
        primitive->predefinedShape = this->currentShapeSelected;
        primitive->parametersProvided = {this->selectedSigma};

        primitive->material = this->selectedTerrainType;
        primitive->name = toCapitalize(stringFromPredefinedShape(currentShapeSelected));
        if (this->currentPositioning == ImplicitPatch::PositionalLabel::FIXED_POS && position.z != 0.f) {
            ImplicitUnaryOperator* translate = new ImplicitUnaryOperator;
            translate->translate(Vector3(0.f, 0.f, position.z));
            translate->composableA = primitive;
            // Return the translation operation
            translate->updateCache();
            return translate;
        }
        return primitive;
    }
}

ImplicitPatch *PrimitivePatchesInterface::createOperationPatchFromParameters(ImplicitPatch *composableA, ImplicitPatch *composableB, ImplicitOperator *replacedPatch)
{
    ImplicitOperator* operation;
    if (replacedPatch == nullptr)
        operation = new ImplicitOperator;
    else
        operation = replacedPatch;

    operation->composableA = composableA;
    operation->composableB = composableB;
    operation->composeFunction = currentOperation;
    operation->positionalB = currentPositioning;
    operation->blendingFactor = selectedBlendingFactor;
    operation->withIntersectionOnB = applyIntersection;

    operation->name = toCapitalize(stringFromCompositionOperation(currentOperation));
    operation->updateCache();
    return operation;
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
    std::vector<ImplicitPatch*> waitingPatches = {this->implicitTerrain.get()};
    std::set<ImplicitPatch*> visitedPatches;
    while (!waitingPatches.empty()) {
        auto current = waitingPatches.back();
        waitingPatches.pop_back();
        if (visitedPatches.find(current) != visitedPatches.end())
            continue;
        visitedPatches.insert(current);
        // Check children if they exist
        ImplicitOperator* currentAsOperator = dynamic_cast<ImplicitOperator*>(current);
        if (currentAsOperator) {
            if (currentAsOperator->composableA) {
                if (currentAsOperator->composableA == child)
                    return current;
                waitingPatches.push_back(currentAsOperator->composableA);
            }
            if (currentAsOperator->composableB) {
                if (currentAsOperator->composableB == child)
                    return current;
                waitingPatches.push_back(currentAsOperator->composableB);
            }
        }
        ImplicitNaryOperator* currentAsNary = dynamic_cast<ImplicitNaryOperator*>(current);
        if (currentAsNary) {
            for (auto& c : currentAsNary->composables) {
                if (c == child)
                    return current;
                waitingPatches.push_back(c);
            }
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
    //    this->cleanPatch(this->implicitTerrain);
        bool patchesCanBeLocked = this->mainFilename != "";
        std::vector<std::tuple<ImplicitPatch*, int, bool>> waitingPatches = { {this->implicitTerrain.get(), 0, false} };

        this->nbPrimitives = 0;
        this->nbBinaryOperators = 0;
        this->nbUnaryOperators = 0;
        this->nbNaryOperators = 0;
        int currentIndex = 0;

        while (!waitingPatches.empty()) {
            auto [current, level, locked] = waitingPatches.back();
            current->index = currentIndex++;
            bool childrenAreLocked = locked || (current->used_json_filename != "" && patchesCanBeLocked && current->used_json_filename != this->mainFilename);
            waitingPatches.pop_back();
            primitiveSelectionGui->addItem(new HierarchicalListWidgetItem((locked ? "*" : "") + current->toString(), current->index, level));
            this->storedPatches.push_back(current);
            ImplicitOperator* currentAsOperator = dynamic_cast<ImplicitOperator*>(current);
            ImplicitNaryOperator* currentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(current);
            if (currentAsOperator) {
                if (currentAsOperator->composableA) {
                    waitingPatches.push_back({currentAsOperator->composableA, level+1, childrenAreLocked});
                }
                if (currentAsOperator->composableB) {
                    waitingPatches.push_back({currentAsOperator->composableB, level+1, childrenAreLocked});
                }
            } else if (currentAsNaryOperator) {
                for (const auto& c : currentAsNaryOperator->composables)
                    waitingPatches.push_back({c, level+1, childrenAreLocked});
            }

            // Counting patches
            if (dynamic_cast<ImplicitPrimitive*>(current) != nullptr)
                nbPrimitives ++;
            else if (dynamic_cast<ImplicitUnaryOperator*>(current) != nullptr)
                nbUnaryOperators ++;
            else if (currentAsOperator)
                nbBinaryOperators ++;
            else
                nbNaryOperators ++;
        }
//        if (nbPrimitives > 1) // If there is more than just the Identity
//            this->savePatchesAsFile(this->mainFilename);
        ImplicitPatch::currentMaxIndex = currentIndex;
        this->findAllSubfiles();
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
    auto suppDim = this->selectedPatch()->getSupportDimensions();
    Vector3 scaleToDisplayPatch = suppDim / debuggingVoxels.getDimensions(); // Vector3(1.f, 1.f, 1.f); // 1.f / this->debuggingVoxelsScale; //currentlySelectedPatch->getDimensions() / debuggingVoxels.getDimensions();
    debuggingVoxelsMesh.shader->setVector("scale", scaleToDisplayPatch);

    // Translate the debugging mesh to the right position
    Vector3 positionToDisplayPatch = this->selectedPatch()->getSupportBBox().min(); // this->debuggingVoxelsPosition * debuggingVoxels.getDimensions() * .5f; // currentlySelectedPatch->getBBox().first; //currentlySelectedPatch->position - currentlySelectedPatch->getDimensions() - Vector3(1.f, 1.f, 1.f) * scaleToDisplayPatch;
    debuggingVoxelsMesh.shader->setFloat("offsetX", positionToDisplayPatch.x);
    debuggingVoxelsMesh.shader->setFloat("offsetY", positionToDisplayPatch.y);
    debuggingVoxelsMesh.shader->setFloat("offsetZ", positionToDisplayPatch.z);
    std::cout << positionToDisplayPatch << " " << scaleToDisplayPatch << std::endl;

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
        QRadioButton* oneSideBlendButton = new QRadioButton("1-blend");

        QRadioButton* abovePosButton = new QRadioButton("Above");
        QRadioButton* insideTopPosButton = new QRadioButton("Inside top");
        QRadioButton* insideBottomPosButton = new QRadioButton("Inside bottom");
        QRadioButton* fixedPosButton = new QRadioButton("Fixed");

        FancySlider* blendingFactorSlider = new FancySlider(Qt::Orientation::Horizontal, 0.1f, 10.f, .1f);
        QCheckBox* applyIntersectionButton = new QCheckBox("Intersection");
        QPushButton* swapButton = new QPushButton("Swap A and B");

        layout->addWidget(createSliderGroup("Blend factor", blendingFactorSlider));
        layout->addWidget(createHorizontalGroup({createVerticalGroup({
                                                  stackingButton,
                                                  blendingButton,
                                                  replacingButton,
                                                  oneSideBlendButton,
                                              }),
                                                 createVerticalGroup({
                                                     abovePosButton,
                                                     insideTopPosButton,
                                                     insideBottomPosButton,
                                                     fixedPosButton
                                                 })
                                                }));
        layout->addWidget(createHorizontalGroup({ applyIntersectionButton, swapButton }));

        stackingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::STACK);
        blendingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::BLEND);
        replacingButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::REPLACE);
        oneSideBlendButton->setChecked(patchAsOperator->composeFunction == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND);

        abovePosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::ABOVE);
        insideTopPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::INSIDE_TOP);
        insideBottomPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM);
        fixedPosButton->setChecked(patchAsOperator->positionalB == ImplicitPatch::PositionalLabel::FIXED_POS);

        blendingFactorSlider->setfValue(patchAsOperator->blendingFactor);
        applyIntersectionButton->setChecked(patchAsOperator->withIntersectionOnB);

        QObject::connect(stackingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::STACK; });
        QObject::connect(blendingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::BLEND; });
        QObject::connect(replacingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::REPLACE; });
        QObject::connect(oneSideBlendButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->composeFunction = ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND; });

        QObject::connect(abovePosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::ABOVE; });
        QObject::connect(insideTopPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::INSIDE_TOP; });
        QObject::connect(insideBottomPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::INSIDE_BOTTOM; });
        QObject::connect(fixedPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsOperator->positionalB = ImplicitPatch::PositionalLabel::FIXED_POS; });
        QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [=](float newVal) {
            patchAsOperator->blendingFactor = newVal;
        });
        QObject::connect(applyIntersectionButton, &QCheckBox::toggled, this, [=](bool checked) { patchAsOperator->withIntersectionOnB = checked; });
        QObject::connect(swapButton, &QPushButton::pressed, this, [=]() { patchAsOperator->swapAB(); });


    } else {
        ImplicitPrimitive* patchAsPrimitive = dynamic_cast<ImplicitPrimitive*>(patch);
        QRadioButton* createSphereButton = new QRadioButton("Sphere");
        QRadioButton* createBlockButton = new QRadioButton("Block");
        QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
        QRadioButton* createCylinderButton = new QRadioButton("Cylinder");
        QRadioButton* createRockButton = new QRadioButton("Rock");
        QRadioButton* createMountainButton = new QRadioButton("Mountain");
        QRadioButton* createDuneButton = new QRadioButton("Dune");
        QRadioButton* createBasinButton = new QRadioButton("Basin");
        QRadioButton* createCaveButton = new QRadioButton("Cave");
        QRadioButton* createArchButton = new QRadioButton("Arch");
        QRadioButton* createNoise2DButton = new QRadioButton("Noise");
        QRadioButton* createFromFileButton = new QRadioButton("From file");

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
                                                        createCylinderButton,
                                                        createRockButton,
                                                        createMountainButton,
                                                        createDuneButton,
                                                        createBasinButton,
                                                        createCaveButton,
                                                        createArchButton,
                                                        createNoise2DButton,
                                                        createFromFileButton
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
        createCylinderButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Cylinder);
        createRockButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Rock);
        createMountainButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Mountain);
        createDuneButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Dune);
        createBasinButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Basin);
        createCaveButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Cave);
        createArchButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Arch);
        createNoise2DButton->setChecked(patchAsPrimitive->predefinedShape == ImplicitPatch::PredefinedShapes::Noise2D);
        createFromFileButton->setChecked(patchAsPrimitive->used_json_filename != "");

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

        QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Sphere; });
        QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Block; });
        QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Gaussian; });
        QObject::connect(createCylinderButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Cylinder; });
        QObject::connect(createRockButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Rock; });
        QObject::connect(createMountainButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Mountain; });
        QObject::connect(createDuneButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Dune; });
        QObject::connect(createBasinButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Basin; });
        QObject::connect(createCaveButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Cave; });
        QObject::connect(createArchButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Arch; });
        QObject::connect(createNoise2DButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Noise2D; });
        QObject::connect(createFromFileButton, &QRadioButton::toggled, this, [=](bool checked) {
            if (checked) {
                patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::ImplicitHeightmap;
                QString q_filename = QFileDialog::getOpenFileName(this, QString("Open a composition"), QString::fromStdString("saved_maps/"));
                std::string filename = q_filename.toStdString();
                if (!q_filename.isEmpty()) {
                    patchAsPrimitive->used_json_filename = q_filename.toStdString();
                }
            }
        });

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
    this->done(-1);
}

void PatchReplacementDialog::confirm()
{
    patch->update();
    setResult(1);
    this->done(1);
}
