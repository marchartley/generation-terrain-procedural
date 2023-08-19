#include "PrimitivePatchesInterface.h"
#include "Graphics/DisplayGraphics.h"
#include "Interface/InterfaceUtils.h"
//#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"
//#include "Utils/stb_image.h"
#include "Utils/ShapeCurve.h"

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

void PrimitivePatchesInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

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
    QRadioButton* createRippleButton = new QRadioButton("Ripple");

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
    QPushButton* addRipplesButton = new QPushButton("Ripples");
    QPushButton* deformFromFlowButton = new QPushButton("Deform flow");

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

    QPushButton* createStructureButton = new QPushButton("Auto gen.");

    primitiveSelectionGui = new HierarchicalListWidget(this);

    layout->addWidget(createMultiColumnGroup({
                                                    createSphereButton,
                                                    createBlockButton,
//                                                    createGaussianButton,
                                                    createCylinderButton,
                                                    createRockButton,
                                                    createMountainButton,
//                                                    createDuneButton,
//                                                    createBasinButton,
                                                    createCaveButton,
                                                    createArchButton,
                                                    createNoise2DButton,
                                                    createMountainChainButton,
                                                    createPolygonButton,
                                                    createTunnelButton,
                                                    createFromFileButton,
//                                                    createRippleButton
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
    layout->addWidget(createHorizontalGroup({applyIntersectionButton, createStructureButton}));
    layout->addWidget(createHorizontalGroup({addNoiseButton, addDistoButton, addSpreadButton, addRipplesButton, deformFromFlowButton}));
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
        }
    });
    QObject::connect(createRippleButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) this->setSelectedShape(ImplicitPatch::PredefinedShapes::Ripple); });
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
    QObject::connect(addRipplesButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::rippleScene);
    QObject::connect(deformFromFlowButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::deformationFromFlow);

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
    QObject::connect(createStructureButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::structureAutoGeneration);

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
    createRippleButton->setChecked(this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Ripple);

    stackingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::STACK);
    blendingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::BLEND);
    replacingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::REPLACE);
    oneSideBlendButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND);

    abovePosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::ABOVE);
    insideTopPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::INSIDE_TOP);
    insideBottomPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM);
    fixedPosButton->setChecked(this->currentPositioning == ImplicitPatch::PositionalLabel::FIXED_POS);

    blendingFactorSlider->setfValue(this->selectedBlendingFactor);

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

//    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::currentItemChanged, this, &PrimitivePatchesInterface::updateSelectedPrimitiveItem);
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemClicked, this, [=](QListWidgetItem* item) { updateSelectedPrimitiveItem(item); });
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemDoubleClicked, this, &PrimitivePatchesInterface::openPrimitiveModificationDialog);
//    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemChangedHierarchy, this, &PrimitivePatchesInterface::modifyPrimitiveHierarchy);


    // This initialization should be in the constructor I guess...
    this->primitiveControlPoint = std::make_unique<ControlPoint>(Vector3(), 5.f);
    this->primitiveControlPoint->allowAllAxisRotations(true);
    this->primitiveControlPoint->allowAllAxisTranslation(true);
    this->primitiveControlPoint->hide();

    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::modified, this, &PrimitivePatchesInterface::moveDebugBoxWithControlPoint);
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::translationApplied, this, &PrimitivePatchesInterface::translatePatch);
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::rotationApplied, this, &PrimitivePatchesInterface::rotatePatch);
    /*
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::rotationApplied, this, [=](const Vector3& rotation){
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

void PrimitivePatchesInterface::mouseMovedOnMapEvent(const Vector3& newPosition, TerrainModel* model)
{
    if (this->isVisible()) {
        this->currentPos = newPosition;
        this->setSelectedShape(this->currentShapeSelected);
    }
}

void PrimitivePatchesInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel* model)
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

void PrimitivePatchesInterface::createPatchWithOperation(const Vector3& pos)
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
        this->implicitTerrain->deleteAllChildren();
        this->implicitTerrain->addChild(operation);
    } else {
        ImplicitBinaryOperator* parentAsBinary = dynamic_cast<ImplicitBinaryOperator*>(_parent);
        if (parentAsBinary) {
            if (previousMain == parentAsBinary->composableA())
                parentAsBinary->addChild(operation, 0);
            else
                parentAsBinary->addChild(operation, 1);
        }
        ImplicitNaryOperator* parentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(_parent);
        if (parentAsNaryOperator) {
            int iComposable = 0;
            for (iComposable = 0; iComposable < parentAsNaryOperator->composables.size(); iComposable++) {
                if (previousMain == parentAsNaryOperator->composables[iComposable])
                    break;
            }
            parentAsNaryOperator->addChild(operation, iComposable);

        }

    }
//        ImplicitPatch* newOperation = this->createOperationPatchFromParameters(previousMain, patch);
//        this->implicitTerrain = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;
    this->updateMapWithCurrentPatch();

    this->storedPatches.push_back(patch);
    this->storedPatches.push_back(operation);
    this->storedPatches.push_back(implicitTerrain.get());
    this->updatePrimitiveList();

//    if (this->mainFilename != "")
//        this->savePatchesAsFile(this->mainFilename);

}

void PrimitivePatchesInterface::setSelectedShape(ImplicitPatch::PredefinedShapes newShape)
{
    Vector3 newPosition = this->currentPos;
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
    this->implicitTerrain->deleteAllChildren(); // = new ImplicitPrimitive;
    this->implicitTerrain->name = "Identity";
    this->updateMapWithCurrentPatch();
    this->storedPatches.push_back(implicitTerrain.get());
    this->updatePrimitiveList();
//    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    //    voxelGrid->fromIsoData();

//    this->addParametricPoint(Vector3(45, 0, 0));
//    this->addParametricPoint(Vector3(35, 0, 0));
//    this->addParametricPoint(Vector3(25, 0, 0));
//    this->currentShapeSelected = ImplicitPatch::ParametricTunnel;
//    this->selectedSigma = 5.f;
//    this->createPatchFromParameters(Vector3());
//    this->updateMapWithCurrentPatch();
//    this->storedPatches.push_back(implicitTerrain.get());
//    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::updateMapWithCurrentPatch()
{
    this->implicitTerrain->updateCache();
    this->implicitTerrain->augment();
//    this->layerGrid->layers = this->layerGrid->previousState;
    this->layerGrid->reset();
//    this->implicitTerrain->cleanCache();
    this->layerGrid->add(this->implicitTerrain.get()/*, SAND, false*/);
    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    voxelGrid->smoothVoxels();
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
    ImplicitNoise* AsBin = dynamic_cast<ImplicitNoise*>(selectedPatch);
    ImplicitNoise* noisePatch;

    if (AsBin == nullptr) {
        AsBin = dynamic_cast<ImplicitNoise*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (AsBin != nullptr) {
        noisePatch = AsBin;
    } else {
        noisePatch = new ImplicitNoise;
        noisePatch->addChild(selectedPatch);
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->deleteAllChildren();
            this->implicitTerrain->addChild(noisePatch);
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitBinaryOperator* parentAsBinary = dynamic_cast<ImplicitBinaryOperator*>(_parent);
            if (parentAsBinary) {
                if (selectedPatch == parentAsBinary->composableA())
                    parentAsBinary->addChild(noisePatch, 0);
                else
                    parentAsBinary->addChild(noisePatch, 1);
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
    ImplicitWraping* AsBin = dynamic_cast<ImplicitWraping*>(selectedPatch);
    ImplicitWraping* distortionPatch;

    if (AsBin == nullptr) {
        AsBin = dynamic_cast<ImplicitWraping*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (AsBin != nullptr) {
        distortionPatch = AsBin;
    } else {
        distortionPatch = new ImplicitWraping;
        distortionPatch->addChild(selectedPatch);
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->deleteAllChildren();
            this->implicitTerrain->addChild(distortionPatch);
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitBinaryOperator* parent = dynamic_cast<ImplicitBinaryOperator*>(_parent);
            if (parent) {
                if (selectedPatch == parent->composableA())
                    parent->addChild(distortionPatch, 0);
                else
                    parent->addChild(distortionPatch, 1);
            }
        }
    }
    GridV3 disto = GridV3({
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
    ImplicitSpread* AsBin = dynamic_cast<ImplicitSpread*>(selectedPatch);
    ImplicitSpread* distortionPatch;

    if (AsBin == nullptr) {
        AsBin = dynamic_cast<ImplicitSpread*>(this->naiveApproachToGetParent(selectedPatch));
    }
    if (AsBin != nullptr) {
        distortionPatch = AsBin;
    } else {
        distortionPatch = new ImplicitSpread;
        distortionPatch->addChild(selectedPatch);
        if (selectedPatch == this->implicitTerrain.get()) {
            this->implicitTerrain->deleteAllChildren();
            this->implicitTerrain->addChild(distortionPatch);
        } else {
            ImplicitPatch* _parent = this->naiveApproachToGetParent(selectedPatch);
            ImplicitBinaryOperator* parent = dynamic_cast<ImplicitBinaryOperator*>(_parent);
            if (parent) {
                if (selectedPatch == parent->composableA())
                    parent->addChild(distortionPatch, 0);
                else
                    parent->addChild(distortionPatch, 1);
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

void PrimitivePatchesInterface::rippleScene()
{
    GridI terrainSurface = this->voxelGrid->getVoxelValues().binarize();
    GridF waterSpeed(terrainSurface.getDimensions());
    GridV3 flowfield = this->voxelGrid->getFlowfield();
    GridV3 normals = terrainSurface.gradient();
    for(auto& n : normals)
        n.normalize();
    terrainSurface.raiseErrorOnBadCoord = false;
    terrainSurface.defaultValueOnBadCoord = 0;

    std::vector<Vector3> allSurfaceVoxels;
    for (int x = 0; x < terrainSurface.sizeX; x++) {
        for (int y = 0; y < terrainSurface.sizeY; y++) {
            for (int z = 0; z < terrainSurface.sizeZ; z++) {
                Vector3 pos(x, y, z);
                if (terrainSurface.at(pos) == 0) continue;
                bool isSurface = (terrainSurface.at(pos + Vector3(0, 0, 1)) == 0);
                if (isSurface) {
                    terrainSurface.at(pos) = 1;
                    waterSpeed.at(pos) = std::abs(1 - flowfield.at(pos).dot(normals.at(pos)));
                    if (waterSpeed.at(pos) > 0.01)
                        allSurfaceVoxels.push_back(pos);
                } else {
                    terrainSurface.at(pos) = 0;
                }
            }
        }
    }
    ImplicitNaryOperator* allRipples = new ImplicitNaryOperator;

    std::shuffle(allSurfaceVoxels.begin(), allSurfaceVoxels.end(), random_gen::random_generator);
    for (int i = 0; i < std::min(100, int(allSurfaceVoxels.size())); i++) {
        Vector3 pos = allSurfaceVoxels[i];
        Vector3 waterFlow = flowfield.at(pos);
        float waterSpeed = waterFlow.norm();
        float rippleWidth = 20 * waterSpeed;
        float rippleDepth = 40; // * waterSpeed;
        Vector3 rippleDims = Vector3(rippleWidth, rippleDepth, 4);
        ImplicitPrimitive* ripple = ImplicitPatch::createPredefinedShape(ImplicitPatch::PredefinedShapes::Ripple, rippleDims, 0.f);
        ripple->position = (pos - rippleDims.xy() * .5f) - Vector3(0, 0, 3);
        ripple->material = TerrainTypes::DIRT;

//        UnaryOpRotate rotationTransfo(Vector3(0, 0, waterFlow.toEulerAngles().z), Vector3(rippleWidth * .5f, rippleDepth * .5f, 0));
        ImplicitRotation* rotation = new ImplicitRotation;
        rotation->addChild(ripple);
        rotation->rotate(Vector3(0, 0, waterFlow.toEulerAngles().z));
        allRipples->addChild(rotation);
//        std::cout << pos << std::endl;
    }
    /*for (int x = 0; x < terrainSurface.sizeX; x++) {
        for (int y = 0; y < terrainSurface.sizeY; y++) {
            for (int z = 0; z < terrainSurface.sizeZ; z++) {
                Vector3 pos(x, y, z);
                if (terrainSurface.at(pos) == 0)
                    continue;

                ImplicitPrimitive* ripple = ImplicitPatch::createPredefinedShape(ImplicitPatch::PredefinedShapes::Ripple, Vector3(3, 3, 5), 0.f);
                ripple->position = pos;
                allRipples->composables.push_back(ripple);
                std::cout << pos << std::endl;
            }
        }
    }*/
    std::cout << "Added " << allRipples->composables.size() << " elements" << std::endl;




    this->implicitTerrain->addChild(allRipples);
    this->storedPatches.push_back(allRipples);
    this->updateMapWithCurrentPatch();
    this->updatePrimitiveList();
}

void PrimitivePatchesInterface::deformationFromFlow()
{
    std::vector<ImplicitPatch*> queue = {this->implicitTerrain.get()};
    while (!queue.empty()) {
        auto patch = queue[0];
        queue.erase(queue.begin());
        ImplicitWraping* wrap = new ImplicitWraping;
        AABBox bbox = patch->getSupportBBox();
        if (bbox.dimensions().norm2() < 50*50*50) {
            wrap->addWrapFunction(voxelGrid->getFlowfield().subset(bbox.min(), bbox.max()).resize(10, 10, 10) * 0.1f);
            wrap->addChild(patch);
            ImplicitPatch* parent = this->naiveApproachToGetParent(patch);
            ImplicitNaryOperator* asNary = dynamic_cast<ImplicitNaryOperator*>(parent);
            ImplicitUnaryOperator* asUnary = dynamic_cast<ImplicitUnaryOperator*>(parent);
            ImplicitBinaryOperator* asBinary = dynamic_cast<ImplicitBinaryOperator*>(parent);
            if (asUnary) {
                asUnary->addChild(wrap);
            } else if (asBinary) {
                if (patch == asBinary->composableA())
                    asBinary->addChild(wrap, 0);
                else
                    asBinary->addChild(wrap, 1);
            } else if (asNary) {
                for (size_t i = 0; i < asNary->composables.size(); i++)
                    if (asNary->composables[i] == patch)
                        asNary->addChild(wrap, i);
            }
        } else {
            ImplicitNaryOperator* asNary = dynamic_cast<ImplicitNaryOperator*>(patch);
            ImplicitUnaryOperator* asUnary = dynamic_cast<ImplicitUnaryOperator*>(patch);
            ImplicitBinaryOperator* asBinary = dynamic_cast<ImplicitBinaryOperator*>(patch);
            if (asUnary) {
                queue.push_back(asUnary->composableA());
            } else if (asBinary) {
                queue.push_back(asBinary->composableA());
                queue.push_back(asBinary->composableB());
            } else if (asNary) {
                for (size_t i = 0; i < asNary->composables.size(); i++)
                    queue.push_back(asNary->composables[i]);
            }
        }
//        this->storedPatches.push_back(wrap);
    }
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
        *this->implicitTerrain = *dynamic_cast<ImplicitNaryOperator*>(ImplicitNaryOperator::fromJson(json_content[ImplicitPatch::json_identifier]));
//        this->implicitTerrain->composables = {ImplicitPatch::fromJson(json_content[ImplicitPatch::json_identifier])};
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

void PrimitivePatchesInterface::autoGeneratePasse(const GridF &terrainSurface, const GridV3 &waterFlow, const GridV3 &surfaceNormals, const std::vector<Vector3> &availablePositions)
{
    auto lagoons = this->implicitTerrain->findAll(ImplicitPatch::Polygon);
    auto barriers = this->implicitTerrain->findAll(ImplicitPatch::MountainChain);
    ImplicitPrimitive* barrier = dynamic_cast<ImplicitPrimitive*>(barriers[0]);
    BSpline transformedBarrierCurve;
    for (size_t i = 0; i < barrier->optionalCurve.points.size(); i++) {
        Vector3 firstPointInCurve = barrier->optionalCurve.points[i];
        Vector3 transformedPoint = barrier->getGlobalPositionOf(firstPointInCurve);
        transformedBarrierCurve.points.push_back(transformedPoint);
    }

    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    ShapeCurve transformedLagoonArea;
    for (size_t i = 0; i < lagoon->optionalCurve.points.size(); i++) {
        Vector3 firstPointInCurve = lagoon->optionalCurve.points[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.points.push_back(transformedPoint.xy());
    }

    GridF sqrDistToBarrier(voxelGrid->getDimensions(), -1.f);
    GridV3 barrierCurveNormals(voxelGrid->getDimensions());

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];


        if (!transformedLagoonArea.contains(pos.xy())) {
            float closestTime = transformedBarrierCurve.estimateClosestTime(pos);
            Vector3 closestPoint = transformedBarrierCurve.getPoint(closestTime);
            sqrDistToBarrier.at(pos) = (pos - closestPoint).norm2(); //transformedBarrierCurve.estimateSqrDistanceFrom(pos);
            barrierCurveNormals.at(pos) = transformedBarrierCurve.getNormal(closestTime);
        } else {
            sqrDistToBarrier.at(pos) = 100000.f;
            barrierCurveNormals.at(pos) = Vector3(false);
        }
    }

    Vector3 bestFitPos;
    float bestFitScore = std::numeric_limits<float>::max();

    float distToBarrier;
    Vector3 surfaceNormal;
    Vector3 waterVel;

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        float score = 0;

        if (sqrDistToBarrier.at(pos) > 15*15) score += 100;
        if (waterFlow.at(pos).norm2() < 0.01) score += 200;
//        else score -= waterFlow.at(pos).normalized().dot(surfaceNormals.at(pos)) * 100.f;
        else score -= waterFlow.at(pos).normalized().dot(barrierCurveNormals.at(pos)) * 100.f;

        if (score < bestFitScore) {
            bestFitScore = score;
            bestFitPos = pos;

            distToBarrier = std::sqrt(sqrDistToBarrier.at(pos));
            waterVel = waterFlow.at(pos);
            surfaceNormal = surfaceNormals.at(pos);
        }
    }
    std::cout << "Adding at " << bestFitPos << " because score = " << bestFitScore << std::endl;
    std::cout << "Distance to a barrier : " << distToBarrier <<
                 "\nWater velocity : " << waterVel <<
                 "\nSurface normal : " << surfaceNormal <<
                 "\nScore = - vel . normal * 100.0 = " << -waterVel.normalized().dot(surfaceNormal) * 100.f
              << std::endl;

    Vector3 startingPosition = bestFitPos;
    Vector3 passDirection = barrierCurveNormals.at(bestFitPos).xy().normalize(); //surfaceNormals.at(startingPosition).xy().normalize();
    Vector3 endingPosition = startingPosition + passDirection * 40.f;

    Vector3 diff = endingPosition - startingPosition;

    int nbPoints = 4;
    std::vector<Vector3> passPositions;
    for (int i = 0; i < nbPoints + 2; i++) {
        Vector3 originalPoint = startingPosition + i * diff / (nbPoints + 1.f) + Vector3::random().xy() * diff.norm() / (2.f * nbPoints);
        originalPoint.z = this->implicitTerrain->getHeight(originalPoint);
        passPositions.push_back(originalPoint);
    }

    float width = diff.norm() * .1f;
    AABBox dimensions(passPositions);

//    dimensions.mini -= Vector3(width, width, width);
//    dimensions.maxi += Vector3(width, width, width);
    Vector3 start = dimensions.min();
    for (auto& p : passPositions)
        p = p - dimensions.min(); // + Vector3(width, width, width);
    ImplicitPrimitive* pass = ImplicitPrimitive::createPredefinedShape(ImplicitPatch::ParametricTunnel, dimensions.dimensions() + Vector3(width, width, width), width, passPositions);
    pass->material = WATER;
    pass->position = start; //Vector3::min(startingPosition, endingPosition); // - dimensions.dimensions() * .5f;
    pass->name = "Pass (autogen)";

    this->implicitTerrain->addChild(pass);
}

void PrimitivePatchesInterface::autoGenerateDelta(const GridF &terrainSurface, const GridV3 &waterFlow, const GridV3 &surfaceNormals, const std::vector<Vector3> &availablePositions)
{
    auto lagoons = this->implicitTerrain->findAll(ImplicitPatch::Polygon);
    auto passes = this->implicitTerrain->findAll(ImplicitPatch::ParametricTunnel);
    ImplicitPrimitive* passe = dynamic_cast<ImplicitPrimitive*>(passes[0]);
    BSpline transformedPasseCurve;
    for (size_t i = 0; i < passe->optionalCurve.points.size(); i++) {
        Vector3 firstPointInCurve = passe->optionalCurve.points[i];
        Vector3 transformedPoint = passe->getGlobalPositionOf(firstPointInCurve);
        transformedPasseCurve.points.push_back(transformedPoint);
    }
    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    ShapeCurve transformedLagoonArea;
    for (size_t i = 0; i < lagoon->optionalCurve.points.size(); i++) {
        Vector3 firstPointInCurve = lagoon->optionalCurve.points[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.points.push_back(transformedPoint.xy());
    }

    GridF sqrDistToPasse(voxelGrid->getDimensions(), -1.f);
    GridV3 passeDirections(voxelGrid->getDimensions());

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        Vector3 passStarting = transformedPasseCurve.points.front();
        Vector3 passEnding = transformedPasseCurve.points.back();

        if (transformedLagoonArea.contains(pos.xy())) {
//            if ((pos - passStarting).norm2() < (pos - passEnding).norm2()) {
//                sqrDistToPasse.at(pos) = (pos - passStarting).norm2();
//                passeDirections.at(pos) = -transformedPasseCurve.getDirection(0.f);
//            } else {
                sqrDistToPasse.at(pos) = (pos - passEnding).norm2();
                passeDirections.at(pos) = transformedPasseCurve.getDirection(1.f);
//            }
        } else {
            sqrDistToPasse.at(pos) = 10000.f;
            passeDirections.at(pos) = Vector3(false);
        }
    }

    Vector3 bestFitPos;
    float bestFitScore = std::numeric_limits<float>::max();

    float distToPasse;
    Vector3 surfaceNormal;
    Vector3 waterVel;

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        float score = 0;

        if (!passeDirections.at(pos).isValid()) score = 1000;
        else score = sqrDistToPasse.at(pos);

        if (score < bestFitScore) {
            bestFitScore = score;
            bestFitPos = pos;

            distToPasse = std::sqrt(sqrDistToPasse.at(pos));
            waterVel = waterFlow.at(pos);
            surfaceNormal = surfaceNormals.at(pos);
        }
    }
    std::cout << "Adding at " << bestFitPos << " because score = " << bestFitScore << std::endl;
    std::cout << "Distance to a passe : " << distToPasse <<
                 "\nScore = dist = " << sqrDistToPasse.at(bestFitPos)
              << std::endl;

    Vector3 startingPosition = bestFitPos;

    int nbRivers = 2;
    int nbPoints = 4;
    ImplicitNaryOperator* deltaContainer = new ImplicitNaryOperator;
    deltaContainer->name = "Delta (autogen)";
    float riverLength = (2.f / float(nbRivers)) * passe->optionalCurve.length() * .2f;
    for (int iRiver = 0; iRiver < nbRivers; iRiver++) {
        Vector3 passDirection = passeDirections.at(startingPosition).xy().normalize();
        passDirection.rotate(0, 0, interpolation::inv_linear(float(iRiver) / float(nbRivers - 1), -M_PI*.25f, M_PI*.25f));
        Vector3 endingPosition = startingPosition + passDirection * riverLength;

        Vector3 diff = endingPosition - startingPosition;

        std::vector<Vector3> riverPositions;
        for (int i = 0; i < nbPoints + 2; i++) {
            Vector3 originalPoint = startingPosition + i * diff / (nbPoints + 1.f) + Vector3::random().xy() * diff.norm() / (2.f * nbPoints);
            originalPoint.z = this->implicitTerrain->getHeight(originalPoint);
            riverPositions.push_back(originalPoint);
        }

        float width = diff.norm() * .5f;
        AABBox dimensions(riverPositions);

        Vector3 start = dimensions.min();
        for (auto& p : riverPositions)
            p = p - dimensions.min();
        ImplicitPrimitive* river = ImplicitPrimitive::createPredefinedShape(ImplicitPatch::ParametricTunnel, dimensions.dimensions() + Vector3(width, width, width), width * .3f, riverPositions);
        river->material = WATER;
        river->position = start;
        deltaContainer->addChild(river);
    }

    this->implicitTerrain->addChild(deltaContainer);
}

void PrimitivePatchesInterface::autoGenerateMotu(const GridF &terrainSurface, const GridV3 &waterFlow, const GridV3 &surfaceNormals, const std::vector<Vector3> &availablePositions)
{
    auto passes = this->implicitTerrain->findAll(ImplicitPatch::ParametricTunnel);
    std::vector<BSpline> transformedPassesCurves(passes.size());
    for (size_t i = 0; i < passes.size(); i++) {
        ImplicitPrimitive* passe = dynamic_cast<ImplicitPrimitive*>(passes[i]);
        BSpline transformedPasseCurve;
        for (size_t i = 0; i < passe->optionalCurve.points.size(); i++) {
            Vector3 firstPointInCurve = passe->optionalCurve.points[i];
            Vector3 transformedPoint = passe->getGlobalPositionOf(firstPointInCurve);
            transformedPasseCurve.points.push_back(transformedPoint);
        }
        transformedPassesCurves[i] = transformedPasseCurve;
    }

    auto lagoons = this->implicitTerrain->findAll(ImplicitPatch::Polygon);
    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    ShapeCurve transformedLagoonArea;
    for (size_t i = 0; i < lagoon->optionalCurve.points.size(); i++) {
        Vector3 firstPointInCurve = lagoon->optionalCurve.points[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.points.push_back(transformedPoint.xy());
    }

    GridF sqrDistTo3Passes(voxelGrid->getDimensions(), 100000);
    Matrix3<std::vector<float>> sqrDistToEachPasses(voxelGrid->getDimensions(), std::vector<float>(passes.size(), 1000000));
    Matrix3<std::vector<int>> indicesToClosestPasses(voxelGrid->getDimensions(), std::vector<int>(passes.size()));
    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        auto& distToEach = sqrDistToEachPasses.at(pos);
        auto& indices = indicesToClosestPasses.at(pos);

        for (size_t iPasse = 0; iPasse < transformedPassesCurves.size(); iPasse++) {
            distToEach[iPasse] = transformedPassesCurves[iPasse].estimateSqrDistanceFrom(pos);
            indices[iPasse] = iPasse;
        }
//        std::sort(sqrDistToEachPasses.at(pos).begin(), sqrDistToEachPasses.at(pos).end());
        std::sort(indices.begin(), indices.end(), [&](int a, int b) { return distToEach[a] < distToEach[b]; });
        sqrDistTo3Passes.at(pos) = distToEach[indices[0]]; // + distToEach[indices[1]] + distToEach[indices[2]];
    }

    auto reefs = this->implicitTerrain->findAll(ImplicitPatch::MountainChain);
    std::vector<BSpline> transformedReefCurves(passes.size());
    for (size_t i = 0; i < reefs.size(); i++) {
        ImplicitPrimitive* reef = dynamic_cast<ImplicitPrimitive*>(reefs[i]);
        BSpline transformedReefCurve;
        for (size_t i = 0; i < reef->optionalCurve.points.size(); i++) {
            Vector3 firstPointInCurve = reef->optionalCurve.points[i];
            Vector3 transformedPoint = reef->getGlobalPositionOf(firstPointInCurve);
            transformedReefCurve.points.push_back(transformedPoint);
        }
        transformedReefCurves[i] = transformedReefCurve;
    }

    GridF sqrDistToReef(voxelGrid->getDimensions(), 100000);
    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];

        for (size_t iReef = 0; iReef < transformedReefCurves.size(); iReef++) {
            sqrDistToReef.at(pos) = (transformedLagoonArea.contains(pos) ? std::min(sqrDistToReef.at(pos), transformedReefCurves[iReef].estimateSqrDistanceFrom(pos)) : 100000.f);
        }
    }

    Vector3 bestFitPos;
    float bestFitScore = std::numeric_limits<float>::max();

    float distToPasses = -1;
    Vector3 surfaceNormal = Vector3(false);
    Vector3 waterVel = Vector3(false);

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        float score = 0;

//        if (!passeDirections.at(pos).isValid()) score = 1000;
//        else score = sqrDistToPasse.at(pos);
//        score = sqrDistTo3Passes.at(pos);
        score = std::abs(sqrDistTo3Passes.at(pos) - 10 * 10) + sqrDistToReef.at(pos);

        if (score < bestFitScore) {
            bestFitScore = score;
            bestFitPos = pos;

            distToPasses = std::sqrt(sqrDistTo3Passes.at(pos));
            waterVel = waterFlow.at(pos);
            surfaceNormal = surfaceNormals.at(pos);
        }
    }
    std::cout << "Adding at " << bestFitPos << " because score = " << bestFitScore << std::endl;
    std::cout << "Distance to a passe : " << distToPasses <<
                 "\nScore = dist to pass + dist to reef = " << sqrDistTo3Passes.at(bestFitPos) << " + " << sqrDistToReef.at(bestFitPos)
              << std::endl;

//    Vector3 avgPassesDirections;
//    auto& distToEach = sqrDistToEachPasses.at(bestFitPos);
//    auto& indices = indicesToClosestPasses.at(bestFitPos);
//    for (int i = 0; i < 3; i++) {
//        avgPassesDirections += transformedPassesCurves[indices[i]].getDirection(transformedPassesCurves[indices[i]].estimateClosestTime(bestFitPos));
//    }
//    avgPassesDirections /= 3.f;
    Vector3 avgPassesDirections = Vector3(-1, 1, 0).normalized(); // transformedReefCurves[0].getDirection(transformedReefCurves[0].estimateClosestTime(bestFitPos));
    Vector3 startingPosition = bestFitPos + avgPassesDirections * 5.f;

    ImplicitPrimitive* motu = ImplicitPrimitive::createPredefinedShape(ImplicitPatch::Dune, Vector3(30, 30, 10), 0.f);
    motu->position = startingPosition - motu->dimensions.xy() * .5f;
    motu->material = SAND;
    motu->name = "Motu (autogen)";

    this->implicitTerrain->addChild(motu);
}

void PrimitivePatchesInterface::autoGenerateLongShore(const GridF &terrainSurface, const GridV3 &waterFlow, const GridV3 &surfaceNormals, const std::vector<Vector3> &availablePositions)
{

}

void PrimitivePatchesInterface::structureAutoGeneration()
{    
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
    }

    this->updatePrimitiveList();
    this->updateMapWithCurrentPatch();
}

void PrimitivePatchesInterface::updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous) {
    bool newSelectionIsExisting = false;
    if (current) {
        auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(current);
        int selectedPatchID = selectedPatchItem->ID;
        auto newlySelectedPatch = this->findPrimitiveById(selectedPatchID);
        if (newlySelectedPatch != nullptr && newlySelectedPatch != currentlySelectedPatch) {
            this->currentlySelectedPatch = newlySelectedPatch; //this->findPrimitiveById(selectedPatchID);
            auto patchAABBox = currentlySelectedPatch->getBBox();
            auto patchSupportedAABBox = currentlySelectedPatch->getBBox();
//            auto patchSupportedAABBox = currentlySelectedPatch->getSupportBBox();
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
            this->debuggingVoxelsPosition = this->selectedPatch()->getBBox().min();
//            this->debuggingVoxelsPosition = this->selectedPatch()->getSupportBBox().min();
            this->debuggingVoxelsScale = ratio;
            debuggingVoxels = GridF(resolution);
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
            auto hierarchy = naiveApproachToGetAllParents(currentlySelectedPatch);
            for (auto& node : hierarchy) {
                std::cout << node->name << "\n";
            }
            std::cout << std::endl;
        }
    }
    if (!newSelectionIsExisting) {
        this->patchAABBoxMesh.fromArray(std::vector<Vector3>{});
        this->patchAABBoxMesh.update();
        this->primitiveControlPoint->hide();
        this->debugMeshDisplayed = false;
        this->currentlySelectedPatch = nullptr;
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

        ImplicitBinaryOperator* AsBin = dynamic_cast<ImplicitBinaryOperator*>(current);
        if (AsBin != nullptr) {
            if (AsBin->composableA() != nullptr)
                pendingPatches.push_back(AsBin->composableA());
            if (AsBin->composableB() != nullptr)
                pendingPatches.push_back(AsBin->composableB());
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

void PrimitivePatchesInterface::addParametricPoint(const Vector3& point)
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

        ImplicitBinaryOperator* AsBin = dynamic_cast<ImplicitBinaryOperator*>(current);
        if (AsBin != nullptr) {
            if (AsBin->composableA() != nullptr) {
                auto directionsForA = directions;
                if (AsBin->composableB() == nullptr) { // A is unique child
                    directionsForA.push_back(0);
                } else {
                    directionsForA.push_back(-1);
                }
                queueWithDirections.push_back({AsBin->composableA(), directionsForA, currentIndex});
            }
            if (AsBin->composableB() != nullptr) {
                auto directionsForB = directions;
                directionsForB.push_back(1);
                queueWithDirections.push_back({AsBin->composableB(), directionsForB, currentIndex});
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

void PrimitivePatchesInterface::moveDebugBoxWithControlPoint()
{
    // Get patch being manipulated
    if (this->selectedPatch() != nullptr) {
        // Display modified AABBox
        auto BBox = this->selectedPatch()->getBBox();
        Vector3 minPos = BBox.min();
        Vector3 maxPos = BBox.max();
        Vector3 dim = maxPos - minPos;
        std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
        for (auto& p : box) {
            p = (p * dim) + /*minPos + */primitiveControlPoint->getPosition(); //p * currentlySelectedPatch->getDimensions() + (primitiveControlPoint->getPosition() - currentlySelectedPatch->getDimensions().xy() * .5f);//p * currentlySelectedPatch->getDimensions() + (primitiveControlPoint->getPosition() - (currentlySelectedPatch->getDimensions() * .5f).xy());
            p = (p - BBox.center()).rotate(-primitiveControlPoint->getRotation()) + BBox.center();
        }
        this->patchAABBoxMesh.fromArray(box);
        this->patchAABBoxMesh.update();
    }
}

void PrimitivePatchesInterface::translatePatch(const Vector3& translation)
{
//        QObject::blockSignals(true);
    // Get patch being manipulated
    if (this->selectedPatch() != nullptr) {
        primitiveControlPoint->blockSignals(true);
        primitiveSelectionGui->blockSignals(true);

        ImplicitTranslation* manipulatedAsUnary = dynamic_cast<ImplicitTranslation*>(this->selectedPatch());
        if (manipulatedAsUnary == nullptr) {
            // It's not an unary yet, but the parent might be
            manipulatedAsUnary = dynamic_cast<ImplicitTranslation*>(this->naiveApproachToGetParent(this->selectedPatch()));
        }

        if (manipulatedAsUnary != nullptr) { // We are updating an unary operator
            manipulatedAsUnary->translate(translation); // Just update it
        } else { // Otherwise, create a new Unary operator
            ImplicitTranslation* translate = new ImplicitTranslation;
            if (this->selectedPatch() == this->implicitTerrain.get()) {
                this->implicitTerrain->deleteAllChildren();
                this->implicitTerrain->addChild(translate);
            } else {
                ImplicitNaryOperator* parentAsNaryOperator = this->selectedPatch()->getParent(); //dynamic_cast<ImplicitNaryOperator*>(this->naiveApproachToGetParent(this->selectedPatch()));
                if (parentAsNaryOperator) {
                    for (size_t i = 0; i < parentAsNaryOperator->composables.size(); i++) {
                        if (parentAsNaryOperator->composables[i] == selectedPatch()) {
                            parentAsNaryOperator->addChild(translate, i);
                            break;
                        }
                    }
                }
            }
            translate->addChild(this->selectedPatch());
            translate->translate(translation);
            this->storedPatches.push_back(translate);
            this->updatePrimitiveList();
            this->primitiveSelectionGui->setCurrentItem(translate->index);


        }
        this->updateMapWithCurrentPatch();

        primitiveControlPoint->blockSignals(false);
        primitiveSelectionGui->blockSignals(false);
    }
}

void PrimitivePatchesInterface::rotatePatch(const Vector3& rotation)
{
//        QObject::blockSignals(true);
    // Get patch being manipulated
    if (this->selectedPatch() != nullptr) {
        primitiveControlPoint->blockSignals(true);
        primitiveSelectionGui->blockSignals(true);

        ImplicitRotation* manipulatedAsUnary = dynamic_cast<ImplicitRotation*>(this->selectedPatch());
        if (manipulatedAsUnary == nullptr) {
            // It's not an unary yet, but the parent might be
            manipulatedAsUnary = dynamic_cast<ImplicitRotation*>(this->naiveApproachToGetParent(this->selectedPatch()));
        }

        if (manipulatedAsUnary != nullptr) { // We are updating an unary operator
            manipulatedAsUnary->rotate(rotation); // Just update it
        } else { // Otherwise, create a new Unary operator
            ImplicitRotation* translate = new ImplicitRotation;
            translate->addChild(this->selectedPatch());
            translate->rotate(rotation);
            if (this->selectedPatch() == this->implicitTerrain.get()) {
                this->implicitTerrain->deleteAllChildren();
                this->implicitTerrain->addChild(translate);
            } else {
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
}

ImplicitPatch* PrimitivePatchesInterface::createPatchFromParameters(const Vector3& position, ImplicitPatch *replacedPatch)
{
    Vector3 finalPos = position;
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
        Vector3 translation = finalPos;
        // Create a translation operation to move the new patch
        ImplicitTranslation* translate = new ImplicitTranslation;
        translate->translate(translation);
        translate->addChild(patch);
        // Return the translation operation
        translate->updateCache();
        return translate;
    } else {
        Vector3 funcSize = this->functionSize * layerGrid->scaling;
        if (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel) {
            auto AABBox = this->parametricCurve.AABBox();
            finalPos = std::get<0>(AABBox) - funcSize * .5f;
            funcSize = (std::get<1>(AABBox) + funcSize * .5f) - finalPos;
        }
        BSpline translatedSpline = this->parametricCurve;
        for (auto& p : translatedSpline.points)
            p -= finalPos;
        ImplicitPrimitive* primitive = (ImplicitPrimitive*) ImplicitPatch::createPredefinedShape(currentShapeSelected, funcSize, selectedSigma, translatedSpline);
        primitive->position = finalPos.xy();
        primitive->setDimensions(funcSize);
        primitive->predefinedShape = this->currentShapeSelected;
        primitive->parametersProvided = {this->selectedSigma};

        primitive->material = this->selectedTerrainType;
        primitive->name = toCapitalize(stringFromPredefinedShape(currentShapeSelected));
        if (this->currentPositioning == ImplicitPatch::PositionalLabel::FIXED_POS && finalPos.z != 0.f) {
            ImplicitTranslation* translate = new ImplicitTranslation;
            translate->translate(Vector3(0.f, 0.f, finalPos.z));
            translate->addChild(primitive);
            // Return the translation operation
            translate->updateCache();
            return translate;
        }
        return primitive;
    }
}

ImplicitPatch *PrimitivePatchesInterface::createOperationPatchFromParameters(ImplicitPatch *composableA, ImplicitPatch *composableB, ImplicitBinaryOperator *replacedPatch)
{
    ImplicitBinaryOperator* operation;
    if (replacedPatch == nullptr)
        operation = new ImplicitBinaryOperator;
    else
        operation = replacedPatch;

    operation->addChild(composableA, 0);
    operation->addChild(composableB, 1);
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

std::vector<ImplicitPatch *> PrimitivePatchesInterface::naiveApproachToGetAllParents(ImplicitPatch *child)
{
    std::vector<ImplicitPatch*> parents = {child};
    auto parent = naiveApproachToGetParent(child);
    while(parent != nullptr) {
        parents.push_back(parent);
        parent = naiveApproachToGetParent(parent);
    }
    return parents;
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
            if (current != this->implicitTerrain.get() && current->parent == nullptr)
                    std::cout << current->toString() << " has no parent" << std::endl;
            bool childrenAreLocked = locked || (current->used_json_filename != "" && patchesCanBeLocked && current->used_json_filename != this->mainFilename);
            waitingPatches.pop_back();
            primitiveSelectionGui->addItem(new HierarchicalListWidgetItem((locked ? "*" : "") + current->toString(), current->index, level));
            this->storedPatches.push_back(current);
            ImplicitBinaryOperator* currentAsBinary = dynamic_cast<ImplicitBinaryOperator*>(current);
            ImplicitUnaryOperator* currentAsUnary = dynamic_cast<ImplicitUnaryOperator*>(current);
            ImplicitNaryOperator* currentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(current);
            if (currentAsNaryOperator) {
                for (const auto& c : currentAsNaryOperator->composables)
                    waitingPatches.push_back({c, level+1, childrenAreLocked});
            }

            // Counting patches
            if (dynamic_cast<ImplicitPrimitive*>(current) != nullptr)
                nbPrimitives ++;
            else if (dynamic_cast<ImplicitUnaryOperator*>(current) != nullptr)
                nbUnaryOperators ++;
            else if (currentAsBinary)
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
    auto Nary = dynamic_cast<ImplicitNaryOperator*>(_patch);
    Nary->deleteAllChildren();
    /*
    ImplicitBinaryOperator* patch = dynamic_cast<ImplicitBinaryOperator*>(_patch);
    if (!patch) // Not a composition ? Nothing to do
        return;

    // Check if it only contains Identity patches. In this case, transform it into an Identity function
    if (patch->composableA() && patch->composableB()) {
        if (patch->composableA()->name == "Identity" && patch->composableB()->name == "Identity") {
            patch->name = "Identity";
            patch->composeFunction = ImplicitPatch::CompositionFunction::NONE;
//            patch->dimension = Vector3(0, 0, 0);
            patch->composableA() = nullptr;
            patch->composableB() = nullptr;
        }
    }

    if (patch->composableA() && patch->composableB()) {
        // If it's a composite, clean children before, so we can use their clean data
        cleanPatch(patch->composableA());
        cleanPatch(patch->composableB());
    } else {
        // Nothing to do (?)
    }*/
}

void PrimitivePatchesInterface::displayDebuggingVoxels()
{
    if (!this->selectedPatch())
        return;
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
//    auto suppDim = this->selectedPatch()->getSupportDimensions();
    auto suppDim = this->selectedPatch()->getDimensions();
    Vector3 scaleToDisplayPatch = suppDim / debuggingVoxels.getDimensions(); // Vector3(1.f, 1.f, 1.f); // 1.f / this->debuggingVoxelsScale; //currentlySelectedPatch->getDimensions() / debuggingVoxels.getDimensions();
    debuggingVoxelsMesh.shader->setVector("scale", scaleToDisplayPatch);

    // Translate the debugging mesh to the right position
    //    Vector3 positionToDisplayPatch = this->selectedPatch()->getSupportBBox().min(); // this->debuggingVoxelsPosition * debuggingVoxels.getDimensions() * .5f; // currentlySelectedPatch->getBBox().first; //currentlySelectedPatch->position - currentlySelectedPatch->getDimensions() - Vector3(1.f, 1.f, 1.f) * scaleToDisplayPatch;
    Vector3 positionToDisplayPatch = this->selectedPatch()->getBBox().min(); // this->debuggingVoxelsPosition * debuggingVoxels.getDimensions() * .5f; // currentlySelectedPatch->getBBox().first; //currentlySelectedPatch->position - currentlySelectedPatch->getDimensions() - Vector3(1.f, 1.f, 1.f) * scaleToDisplayPatch;
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
    bool patchIsOperation = (dynamic_cast<ImplicitBinaryOperator*>(patch) != nullptr);

    QVBoxLayout* layout = new QVBoxLayout();

    if (patchIsOperation) {
        ImplicitBinaryOperator* patchAsBinary = dynamic_cast<ImplicitBinaryOperator*>(patch);
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

        stackingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::STACK);
        blendingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::BLEND);
        replacingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::REPLACE);
        oneSideBlendButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND);

        abovePosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::ABOVE);
        insideTopPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::INSIDE_TOP);
        insideBottomPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM);
        fixedPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::FIXED_POS);

        blendingFactorSlider->setfValue(patchAsBinary->blendingFactor);
        applyIntersectionButton->setChecked(patchAsBinary->withIntersectionOnB);

        QObject::connect(stackingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::STACK; });
        QObject::connect(blendingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::BLEND; });
        QObject::connect(replacingButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::REPLACE; });
        QObject::connect(oneSideBlendButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND; });

        QObject::connect(abovePosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::ABOVE; });
        QObject::connect(insideTopPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::INSIDE_TOP; });
        QObject::connect(insideBottomPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::INSIDE_BOTTOM; });
        QObject::connect(fixedPosButton, &QRadioButton::toggled, this, [=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::FIXED_POS; });
        QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [=](float newVal) {
            patchAsBinary->blendingFactor = newVal;
        });
        QObject::connect(applyIntersectionButton, &QCheckBox::toggled, this, [=](bool checked) { patchAsBinary->withIntersectionOnB = checked; });
        QObject::connect(swapButton, &QPushButton::pressed, this, [=]() { patchAsBinary->swapAB(); });


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
                QString q_filename = QDir(QDir::currentPath()).relativeFilePath(QFileDialog::getOpenFileName(this, QString("Open a composition"), QString::fromStdString("saved_maps/")));
                std::string filename = q_filename.toStdString();
                if (!q_filename.isEmpty()) {
                    patchAsPrimitive->used_json_filename = filename;
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
