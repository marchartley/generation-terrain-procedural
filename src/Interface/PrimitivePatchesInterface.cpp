#include "PrimitivePatchesInterface.h"
#include "Interface/InterfaceUtils.h"
//#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"

PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("PrimitivePatchesInterface", parent)
{
    previewMesh.cullFace = false;

    this->mainPatch = new ImplicitPatch;
    this->mainPatch->name = "Identity";


#ifdef linux
    this->mainFilename = "saved_maps/implicit_patches.json";
#else
    this->mainFilename = "saved_maps/implicit_patches.json";
#endif

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

void PrimitivePatchesInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    ActionInterface::affectVoxelGrid(voxelGrid);

    // Waiting for OpenGL to be available to create the shaders...
#ifdef linux
        std::string shaderPath = "src/Shaders/";
#else
        std::string shaderPath = "src/Shaders/";
#endif
        std::string vertexFilename = shaderPath + "MarchingCubes.vert";
        std::string geometryFilename = shaderPath + "MarchingCubes.geom";
        std::string fragmentFilename = shaderPath + "no_shader.frag";

        debuggingVoxelsMesh.shader = std::make_shared<Shader>(vertexFilename, fragmentFilename, geometryFilename);
        debuggingVoxelsMesh.useIndices = false;
        debuggingVoxelsMesh.cullFace = true;
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
    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");
    QRadioButton* negStackingButton = new QRadioButton("Negative stack");
    QRadioButton* stackInWaterButton = new QRadioButton("Stack in water");
    FancySlider* blendingFactorSlider = new FancySlider(Qt::Orientation::Horizontal, 0.1f, 10.f, .1f);

    FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);

    QPushButton* resetButton = new QPushButton("Reset");

    QLabel* selectedFilenameLabel = new QLabel((this->mainFilename == "" ? "No file selected" : QString::fromStdString(this->mainFilename).split("/").back().split("\\").back()));
    QPushButton* fileSelectionButton = new QPushButton("...");
    QCheckBox* enableHotreloadButton = new QCheckBox("Hot reloading ?");

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
                                                    createArchButton
                                          }));
    layout->addWidget(createSliderGroup("Density", densitySlider));
    layout->addWidget(createSliderGroup("Blend factor", blendingFactorSlider));
    layout->addWidget(createVerticalGroup({
                                              stackingButton,
                                              blendingButton,
                                              replacingButton,
                                              negStackingButton,
                                              stackInWaterButton
                                          }));
    layout->addWidget(createMultipleSliderGroup({
                                              {"Width / radius", widthSlider},
                                                    {"Depth", depthSlider},
                                                    {"Height", heightSlider},
                                                    {"Sigma", sigmaSlider}
                                          }));
    layout->addWidget(resetButton);
    layout->addWidget(primitiveSelectionGui);

    layout->addWidget(createHorizontalGroup({selectedFilenameLabel, fileSelectionButton}));
    layout->addWidget(enableHotreloadButton);

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Sphere); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Block); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Gaussian); });
    QObject::connect(createRockButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Rock); });
    QObject::connect(createMountainButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Mountain); });
    QObject::connect(createDuneButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Dune); });
    QObject::connect(createBasinButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Basin); });
    QObject::connect(createCaveButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Cave); });
    QObject::connect(createArchButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(ImplicitPatch::PredefinedShapes::Arch); });
    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [&](float newDensity) { this->selectedDensity = newDensity; });

    QObject::connect(stackingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::REPLACE); });
    QObject::connect(negStackingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::NEG_STACKING); });
    QObject::connect(stackInWaterButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK_IN_WATER); });

    QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedWidth = newVal; });
    QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedHeight = newVal; });
    QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedDepth = newVal; });
    QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedSigma = newVal; });
    QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedBlendingFactor = newVal; });

    QObject::connect(resetButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::resetPatch);

    QObject::connect(fileSelectionButton, &QPushButton::pressed, this, [this, selectedFilenameLabel]() {
        #ifdef linux
            std::string path = "saved_maps/";
        #else
            std::string path = "saved_maps/";
        #endif
        QString fileSelection = QFileDialog::getSaveFileName(this, "Saving file", QString::fromStdString(path), "*.json", nullptr, QFileDialog::DontConfirmOverwrite);
        if (!fileSelection.isEmpty()) {
            this->mainFilename = fileSelection.toStdString();
            selectedFilenameLabel->setText(QString::fromStdString(this->mainFilename).split("/").back().split("\\").back());
            this->savePatchesAsFile(this->mainFilename);
        }
    });
    QObject::connect(enableHotreloadButton, &QCheckBox::toggled, this, [&](bool checked) {
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

    stackingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::STACK);
    blendingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::BLEND);
    replacingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::REPLACE);
    negStackingButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::NEG_STACKING);
    stackInWaterButton->setChecked(this->currentOperation == ImplicitPatch::CompositionFunction::STACK_IN_WATER);
    blendingFactorSlider->setfValue(this->selectedBlendingFactor);

    densitySlider->setfValue(this->selectedDensity);
    widthSlider->setfValue(this->selectedWidth);
    heightSlider->setfValue(this->selectedHeight);
    depthSlider->setfValue(this->selectedDepth);
    sigmaSlider->setfValue(this->selectedSigma);

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
            std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
            for (auto& p : box) {
                p = p * currentlyManipulatedPatch->getDimensions() + (primitiveControlPoint->getPosition() - (currentlyManipulatedPatch->getDimensions() * .5f).xy());
            }
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
        }
    });
    QObject::connect(this->primitiveControlPoint.get(), &ControlPoint::afterModified, this, [=](){
        // Get patch being manipulated
        if (this->currentlyManipulatedPatch != nullptr) {
            // When released, recreate patch + update terrain
            currentlyManipulatedPatch->position = primitiveControlPoint->getPosition() - (currentlyManipulatedPatch->getDimensions() * .5f).xy();
            this->updateMapWithCurrentPatch();
        }
    });



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
        ImplicitPatch* patch = this->createPatchFromParameters(mousePosInMap);
        ImplicitPatch* previousMain = mainPatch;
//        ImplicitPatch* newOperation = this->createOperationPatchFromParameters(previousMain, patch);
        this->mainPatch = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;
        this->updateMapWithCurrentPatch();

        this->storedPatches.push_back(patch);
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
    for (auto& vert : vertices) {
        vert *= functionSize;
        vert += newPosition;
        vert -= functionSize.xy() * .5f;
    }
    this->previewMesh.fromArray(vertices);
    this->previewMesh.update();

    /*
    if (this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Sphere) {
        this->currentSelectionFunction = this->sphereFunction(selectedWidth * .5f);
    }
    if (this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Block) {
        this->currentSelectionFunction = this->blockFunction(selectedHeight);
    }
    if (this->currentShapeSelected == ImplicitPatch::PredefinedShapes::Gaussian) {
        this->currentSelectionFunction = this->gaussianFunction(selectedSigma, selectedHeight);
    }*/
}

void PrimitivePatchesInterface::setCurrentOperation(ImplicitPatch::CompositionFunction newOperation)
{
    this->currentOperation = newOperation;
}

void PrimitivePatchesInterface::resetPatch()
{
    this->storedPatches.clear();
    delete this->mainPatch;
    this->mainPatch = new ImplicitPatch;
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
    this->mainPatch->cleanCache();
    this->layerGrid->add(this->mainPatch, SAND, false);
    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    voxelGrid->fromIsoData();
    heightmap->fromLayerGrid(*layerGrid);
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

void PrimitivePatchesInterface::savePatchesAsFile(std::string filename)
{
    nlohmann::json content = this->mainPatch->toJson();
    std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc);
    out << nlohmann::json({ {"implicitPatches", content} }).dump(1, '\t');
    out.flush();
}

void PrimitivePatchesInterface::loadPatchesFromFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    nlohmann::json json_content = nlohmann::json::parse(content);
    this->lastTimeFileHasBeenModified = QFileInfo(QString::fromStdString(filename)).lastModified();
    if (json_content.contains("implicitPatches")) {
        this->mainPatch = ImplicitPatch::fromJson(json_content["implicitPatches"]);
        this->updateMapWithCurrentPatch();
        this->updatePrimitiveList();
    } else {
        std::cerr << "No patches defined in file " << filename << std::endl;
    }
}

void PrimitivePatchesInterface::hotReloadFile()
{
    if (this->isVisible() && this->enableHotReloading) {
        if (this->mainFilename != "") {
            QFileInfo infos(QString::fromStdString(this->mainFilename));
            auto lastModifTime = infos.lastModified();
            if (lastModifTime != this->lastTimeFileHasBeenModified) {
                this->loadPatchesFromFile(this->mainFilename);
            }
        }
    }
}

void PrimitivePatchesInterface::updateSelectedPrimitiveItem(QListWidgetItem *current, QListWidgetItem *previous) {
    bool newSelectionIsExisting = false;
    if (current) {
        auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem*>(current);
        int selectedPatchID = selectedPatchItem->ID;
        ImplicitPatch* selectedPatch = this->findPrimitiveById(selectedPatchID);
        this->currentlyManipulatedPatch = selectedPatch;
        if (selectedPatch != nullptr) {
            this->primitiveControlPoint->setPosition(selectedPatch->position + (selectedPatch->getDimensions() * .5f).xy());
            std::vector<Vector3> box = CubeMesh::cubesEdgesVertices;
            for (auto& p : box) {
                p = p * selectedPatch->getDimensions() + selectedPatch->position;
            }
            this->patchAABBoxMesh.fromArray(box);
            this->patchAABBoxMesh.update();
            this->primitiveControlPoint->show();
            newSelectionIsExisting = true;

            Vector3 resolution = selectedPatch->getDimensions(); // Vector3(20, 20, 20);
            Vector3 ratio = selectedPatch->getDimensions() * 3.f / (resolution + Vector3(1, 1, 1));
            debuggingVoxels = Matrix3<float>(resolution);
            debuggingVoxels.raiseErrorOnBadCoord = false;
            for (int x = 0; x < resolution.x; x++) {
                for (int y = 0; y < resolution.y; y++) {
                    for (int z = 0; z < resolution.z; z++) {
                        Vector3 pos(x, y, z);
                        debuggingVoxels.at(pos) = selectedPatch->evaluate((pos) * ratio/* + selectedPatch->position*/  - selectedPatch->getDimensions());
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
}

ImplicitPatch* PrimitivePatchesInterface::createPatchFromParameters(Vector3 position, ImplicitPatch *replacedPatch)
{
//    if (this->currentOperation == PATCH_OPERATION::STACK) // Remove the z offset when stacking... Should be done in the Blending behavior?
//        position.z = 0.f;


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


    /*
    switch (this->currentShapeSelected) {
    case ImplicitPatch::PredefinedShapes::Sphere:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createSphereFunction(selectedSigma, selectedWidth * .5f, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Block:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createBlockFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Gaussian:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createGaussianFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Rock:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createRockFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Mountain:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createMountainFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Dune:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createDuneFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Basin:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createBasinFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Cave:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createCaveFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::Arch:
        *patch = ImplicitPatch(pos, useless, size, ImplicitPatch::createArchFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    case ImplicitPatch::PredefinedShapes::None:
        break; // Code should not come here
    }*/
    patch->shape = this->currentShapeSelected;
    patch->densityValue = this->selectedDensity;
    patch->sigmaValue = this->selectedSigma;
    patch->name = stringFromPredefinedShape(patch->shape); //(currentShapeSelected == ImplicitPatch::PredefinedShapes::BLOCK ? "Block" : currentShapeSelected == ImplicitPatch::PredefinedShapes::GAUSSIAN ? "Gaussian" : "Sphere");
    patch->defineFunctionsBasedOnPredefinedShape();
    return patch;
}

ImplicitPatch *PrimitivePatchesInterface::createOperationPatchFromParameters(ImplicitPatch *composableA, ImplicitPatch *composableB, ImplicitPatch *replacedPatch)
{
    ImplicitPatch* operation = nullptr;/*
    if (replacedPatch != nullptr)
        operation = replacedPatch;
    else
        operation = ne;*/

    if (currentOperation == ImplicitPatch::CompositionFunction::STACK) {
        operation = new ImplicitPatch(ImplicitPatch::createStack(composableA, composableB));
    } else if (currentOperation == ImplicitPatch::CompositionFunction::BLEND) {
        operation = new ImplicitPatch(ImplicitPatch::createBlending(composableA, composableB));
        operation->blendingFactor = this->selectedBlendingFactor;
    } else if (currentOperation == ImplicitPatch::CompositionFunction::REPLACE) {
        operation = new ImplicitPatch(ImplicitPatch::createReplacement(composableA, composableB));
    } else if (currentOperation == ImplicitPatch::CompositionFunction::NEG_STACKING) {
        operation = new ImplicitPatch(ImplicitPatch::createNegStacking(composableA, composableB));
    } else if (currentOperation == ImplicitPatch::CompositionFunction::STACK_IN_WATER) {
        operation = new ImplicitPatch(ImplicitPatch::createStackInWater(composableA, composableB));
    }
    operation->name = stringFromCompositionOperation(operation->compositionOperation);
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
    std::vector<ImplicitPatch*> waitingPatches = {this->mainPatch};
    std::set<ImplicitPatch*> visitedPatches;
    while (!waitingPatches.empty()) {
        auto current = waitingPatches.back();
        waitingPatches.pop_back();
        if (visitedPatches.find(current) != visitedPatches.end())
            continue;
        visitedPatches.insert(current);
        if (current->composableA) {
            if (current->composableA == child)
                return current;
//            if (visitedPatches.find(current->composableA) == visitedPatches.end())
            waitingPatches.push_back(current->composableA);
//            visitedPatches.insert(current->composableA);
        }
        if (current->composableB) {
            if (current->composableB == child)
                return current;
//            if (visitedPatches.find(current->composableB) == visitedPatches.end())
            waitingPatches.push_back(current->composableB);
//            visitedPatches.insert(current->composableB);
        }
    }
    return nullptr;
}

std::function<float (Vector3)> PrimitivePatchesInterface::blockFunction(float height)
{
    std::function levelFunction = [height](Vector3 pos) {
        return height - pos.z;
    };
    return levelFunction;
}

std::function<float (Vector3)> PrimitivePatchesInterface::sphereFunction(float radius)
{
    std::function sphere = [radius](Vector3 pos) {
        Vector3 center(radius, radius, radius);
        float value = (pos - radius).norm2() - radius * radius;
        return -value;
    };
    return sphere;
}

std::function<float (Vector3)> PrimitivePatchesInterface::gaussianFunction(float sigma, float height)
{
    std::function gauss = [sigma, height, this](Vector3 pos) {
        return (normalizedGaussian(this->functionSize.xy(), pos.xy(), sigma) * height) - pos.z;
    };
    return gauss;
}

void PrimitivePatchesInterface::updateFunctionSize()
{

}

void PrimitivePatchesInterface::updatePrimitiveList()
{
    primitiveSelectionGui->clear();
    this->storedPatches.clear();
    this->cleanPatch(this->mainPatch);
    bool patchesCanBeLocked = this->mainFilename != "";
    std::vector<std::tuple<ImplicitPatch*, int, bool>> waitingPatches = { {this->mainPatch, 0, false} };
    while (!waitingPatches.empty()) {
        auto [current, level, locked] = waitingPatches.back();
        bool childrenAreLocked = locked || (current->used_json_filename != "" && patchesCanBeLocked && current->used_json_filename != this->mainFilename);
        waitingPatches.pop_back();
        primitiveSelectionGui->addItem(new HierarchicalListWidgetItem((locked ? "*" : "") + current->toString(), current->index, level));
        this->storedPatches.push_back(current);
        if (current->composableA) {
            waitingPatches.push_back({current->composableA, level+1, childrenAreLocked});
        }
        if (current->composableB) {
            waitingPatches.push_back({current->composableB, level+1, childrenAreLocked});
        }
    }
}

void PrimitivePatchesInterface::cleanPatch(ImplicitPatch *patch)
{
    // Check if it only contains Identity patches. In this case, transform it into an Identity function
    if (patch->composableA && patch->composableB) {
        if (patch->composableA->name == "Identity" && patch->composableB->name == "Identity") {
            patch->name = "Identity";
            patch->compositionOperation = ImplicitPatch::CompositionFunction::NONE;
            patch->dimension = Vector3(0, 0, 0);
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
    Vector3 scaleToDisplayPatch = (currentlyManipulatedPatch->getDimensions() * 3.f) / debuggingVoxels.getDimensions();
    debuggingVoxelsMesh.shader->setVector("scale", scaleToDisplayPatch);

    // Translate the debugging mesh to the right position
    Vector3 positionToDisplayPatch = currentlyManipulatedPatch->position - currentlyManipulatedPatch->getDimensions() - Vector3(1.f, 1.f, 1.f) * scaleToDisplayPatch;
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
    bool patchIsOperation = patch->isOperation();


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
    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");
    QRadioButton* negStackingButton = new QRadioButton("Negative stack");
    QRadioButton* stackInWaterButton = new QRadioButton("Stack in water");
    FancySlider* blendingFactorSlider = new FancySlider(Qt::Orientation::Horizontal, 0.1f, 10.f, .1f);

    FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);

    if (!patchIsOperation) {
        layout->addWidget(createMultiColumnGroup({
                                                        createSphereButton,
                                                        createBlockButton,
                                                        createGaussianButton,
                                                        createRockButton,
                                                        createMountainButton,
                                                        createDuneButton,
                                                        createBasinButton,
                                                        createCaveButton,
                                                        createArchButton
                                              }));
        layout->addWidget(createSliderGroup("Density", densitySlider));
        layout->addWidget(createMultipleSliderGroup({
                                                  {"Width / radius", widthSlider},
                                                        {"Depth", depthSlider},
                                                        {"Height", heightSlider},
                                                        {"Sigma", sigmaSlider}
                                              }));
    } else {
        layout->addWidget(createVerticalGroup({
                                                  stackingButton,
                                                  blendingButton,
                                                  replacingButton,
                                                  negStackingButton,
                                                  stackInWaterButton
                                              }));
        layout->addWidget(createSliderGroup("Blend factor", blendingFactorSlider));
    }

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Sphere); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Block); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Gaussian); });
    QObject::connect(createRockButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Rock); });
    QObject::connect(createMountainButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Mountain); });
    QObject::connect(createDuneButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Dune); });
    QObject::connect(createBasinButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Basin); });
    QObject::connect(createCaveButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Cave); });
    QObject::connect(createArchButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(ImplicitPatch::PredefinedShapes::Arch); });
    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [=](float newDensity) { caller->selectedDensity = newDensity; });

    QObject::connect(stackingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(ImplicitPatch::CompositionFunction::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(ImplicitPatch::CompositionFunction::REPLACE); });
    QObject::connect(negStackingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(ImplicitPatch::CompositionFunction::NEG_STACKING); });
    QObject::connect(stackInWaterButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(ImplicitPatch::CompositionFunction::STACK_IN_WATER); });

    QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { caller->selectedWidth = newVal; });
    QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { caller->selectedHeight = newVal; });
    QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { caller->selectedDepth = newVal; });
    QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { caller->selectedSigma = newVal; });
    QObject::connect(blendingFactorSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { caller->selectedBlendingFactor = newVal; });

    createSphereButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Sphere);
    createBlockButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Block);
    createGaussianButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Gaussian);
    createRockButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Rock);
    createMountainButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Mountain);
    createDuneButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Dune);
    createBasinButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Basin);
    createCaveButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Cave);
    createArchButton->setChecked(caller->currentShapeSelected == ImplicitPatch::PredefinedShapes::Arch);

    stackingButton->setChecked(caller->currentOperation == ImplicitPatch::CompositionFunction::STACK);
    blendingButton->setChecked(caller->currentOperation == ImplicitPatch::CompositionFunction::BLEND);
    replacingButton->setChecked(caller->currentOperation == ImplicitPatch::CompositionFunction::REPLACE);
    negStackingButton->setChecked(caller->currentOperation == ImplicitPatch::CompositionFunction::NEG_STACKING);
    stackInWaterButton->setChecked(caller->currentOperation == ImplicitPatch::CompositionFunction::STACK_IN_WATER);
    blendingFactorSlider->setfValue(caller->selectedBlendingFactor);

    densitySlider->setfValue(caller->selectedDensity);
    widthSlider->setfValue(caller->selectedWidth);
    heightSlider->setfValue(caller->selectedHeight);
    depthSlider->setfValue(caller->selectedDepth);
    sigmaSlider->setfValue(caller->selectedSigma);

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
    ImplicitPatch* replacement;
    if (patch->isOperation()) {
        replacement = caller->createOperationPatchFromParameters(patch->composableA, patch->composableB);
    } else {
        replacement = caller->createPatchFromParameters(this->patch->position + (this->patch->getDimensions() * .5f).xy()); // Position at the center of current patch
    }
    *this->patch = *replacement;
    delete replacement;
    setResult(1);
//    this->close();
    this->done(1);
}
