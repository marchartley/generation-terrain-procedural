#include "PrimitivePatchesInterface.h"
#include <QFileDialog>
#include "GUIElements/ImageViewer.h"
#include "GUIElements/InterfaceUtils.h"
//#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"
//#include "Utils/stb_image.h"
#include "Curves/Contour.h"
#include "Graphics/MarchingCubes.h"
#include "Interface/EnvObjsInterface.h"
#include "TerrainGen/ImplicitTree.h"
#include <QSaveFile>
#include <QSignalBlocker>
#include <memory>

namespace {
void wrapImplicitChild(ImplicitPatch* child, ImplicitUnaryOperator* wrapper)
{
    if (!child || !child->parent) throw std::invalid_argument("Wrapping requires an attached implicit child");
    auto* parent = child->parent;
    auto slot = std::find(parent->composables.begin(), parent->composables.end(), child);
    if (slot == parent->composables.end()) throw std::logic_error("Missing implicit child slot");
    parent->addChild(wrapper, int(slot - parent->composables.begin()));
    wrapper->addChild(child);
}
}


PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("primitivepatches", "Implicit modeling", "model", "Feature primitives", "feature_primitive_button.png", parent)
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
    GlobalsGL::ContextScope contextScope;
    if (!this->isVisible())
        return;

    if (constructionMode && previewMesh.shader != nullptr) {
        previewMesh.shader->setVector("color", std::vector<float>({0.f, 0.4f, 0.8f, 0.5f}));
        previewMesh.display();
    }
    if (!patchAABBoxMesh.vertexArray.empty() && patchAABBoxMesh.shader != nullptr) {
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        patchAABBoxMesh.shader->setVector("color", std::vector<float>({0.f, 0.8f, 0.4f, 1.0f}));
        patchAABBoxMesh.display(GL_LINES, 3.f);
        glPolygonMode(GL_FRONT, polygonMode[0]);
        glPolygonMode(GL_BACK, polygonMode[1]);
    }
    if (this->primitiveControlPoint) this->primitiveControlPoint->display();

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

    debuggingVoxelsMesh.shader = std::make_shared<Shader>(noShaderVertexFilename, noShaderFragmentFilename);
    patchAABBoxMesh.shader = debuggingVoxelsMesh.shader;
    previewMesh.shader = debuggingVoxelsMesh.shader;
    debuggingVoxelsMesh.useIndices = false;
    debuggingVoxelsMesh.cullFace = true;

    parametricCurveMesh = Mesh(std::make_shared<Shader>(noShaderVertexFilename, noShaderFragmentFilename), true, GL_LINES);
    parametricCurveMesh.useIndices = false;

    this->loadTransformationRules();
//    this->loadPatchesFromFile(this->mainFilename);
}

InterfaceUI* PrimitivePatchesInterface::createGUI()
{
    auto UI = new InterfaceUI();

    auto constructionModeButton = new CheckboxElement("Construction mode", constructionMode);

    auto shapeCombobox = new ComboboxElement("Shape");
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Sphere", ImplicitPatch::Sphere));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Block", ImplicitPatch::Block));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Gaussian", ImplicitPatch::Gaussian));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Cylinder", ImplicitPatch::Cylinder));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Rock", ImplicitPatch::Rock));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Mountain", ImplicitPatch::Mountain));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Dune", ImplicitPatch::Dune));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Basin", ImplicitPatch::Basin));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Cave", ImplicitPatch::Cave));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Arch", ImplicitPatch::Arch));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Noise", ImplicitPatch::Noise2D));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Mountains", ImplicitPatch::MountainChain));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Polygon", ImplicitPatch::Polygon));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Tunnel", ImplicitPatch::ParametricTunnel));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Ripple", ImplicitPatch::Ripple));
    shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("...", ImplicitPatch::ImplicitHeightmap));
    for (size_t i = 0; i < shapeCombobox->choices.size(); ++i)
        if (static_cast<ComboboxLineElement<ImplicitPatch::PredefinedShapes>*>(shapeCombobox->choices[i])->value == currentShapeSelected)
            shapeCombobox->combobox()->setCurrentIndex(int(i));
    shapeCombobox->setOnSelectionChanged([=](int) {
        auto shape = shapeCombobox->getSelection<ImplicitPatch::PredefinedShapes>();
        this->setSelectedShape(shape);
        if (shape == ImplicitPatch::ImplicitHeightmap) {
            this->openFileForNewPatch();
        }
    });

    auto stackingButton = new RadioButtonElement("Stack", ImplicitPatch::CompositionFunction::STACK, currentOperation);
    auto blendingButton = new RadioButtonElement("Blend", ImplicitPatch::CompositionFunction::BLEND, currentOperation);
    auto replacingButton = new RadioButtonElement("Replace", ImplicitPatch::CompositionFunction::REPLACE, currentOperation);
    auto oneSideBlendButton = new RadioButtonElement("1-blend", ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND, currentOperation);
    auto differenceButton = new RadioButtonElement("A minus B", ImplicitPatch::DIFFERENCE, currentOperation);
    auto intersectionButton = new RadioButtonElement("A intersect B", ImplicitPatch::INTERSECTION, currentOperation);

    auto abovePosButton = new RadioButtonElement("Above", ImplicitPatch::PositionalLabel::ABOVE, currentPositioning);
    auto insideTopPosButton = new RadioButtonElement("Inside top", ImplicitPatch::PositionalLabel::INSIDE_TOP, currentPositioning);
    auto insideBottomPosButton = new RadioButtonElement("Inside bottom", ImplicitPatch::PositionalLabel::INSIDE_BOTTOM, currentPositioning);
    auto fixedPosButton = new RadioButtonElement("Fixed", ImplicitPatch::PositionalLabel::FIXED_POS, currentPositioning);

    auto blendingFactorSlider = new SliderElement("Blend factor", 2.f, 10.f, .1f, selectedBlendingFactor);

    auto widthSlider = new SliderElement("Width / radius", 0.f, 50.f, .1f, selectedWidth);
    auto depthSlider = new SliderElement("Depth", 0.f, 50.f, .1f, selectedDepth);
    auto heightSlider = new SliderElement("Height", 0.f, 50.f, .1f, selectedHeight);
    auto sigmaSlider = new SliderElement("Sigma", 0.f, 20.f, .1f, selectedSigma);
    widthSlider->setOnValueChanged([this](float value) { setSelectedWidth(value); Q_EMIT updated(); });
    depthSlider->setOnValueChanged([this](float value) { setSelectedDepth(value); Q_EMIT updated(); });
    heightSlider->setOnValueChanged([this](float value) { setSelectedHeight(value); Q_EMIT updated(); });
    sigmaSlider->setOnValueChanged([this](float value) { setSelectedSigma(value); Q_EMIT updated(); });

    auto addNoiseButton = new ButtonElement("Add noise", [=]() { this->addNoiseOnSelectedPatch(); });
    auto addDistoButton = new ButtonElement("Distord", [=]() { this->addDistortionOnSelectedPatch(); });
    auto addSpreadButton = new ButtonElement("Spread", [=]() { this->addSpreadOnSelectedPatch(); });
    auto addRipplesButton = new ButtonElement("Ripples", [=]() { this->rippleScene(); });
    auto deformFromFlowButton = new ButtonElement("Deform flow", [=]() { this->deformationFromFlow(); });

    auto resetButton = new ButtonElement("Reset", [=]() { this->resetPatch(); });

    auto selectedFilenameLabel = new LabelElement((this->mainFilename == "" ? "No file selected" : split(split(this->mainFilename, "/").back(), "\\").back()));
    auto fileSelectionButton = new ButtonElement("...", [this, selectedFilenameLabel]() {
        std::string path = "saved_maps/";
        QString fileSelection = QFileDialog::getSaveFileName(this, "Saving file", QString::fromStdString(path), "*.json", nullptr, QFileDialog::DontConfirmOverwrite);
        if (!fileSelection.isEmpty()) {
            this->mainFilename = fileSelection.toStdString();
            selectedFilenameLabel->setText(split(split(this->mainFilename, "/").back(), "\\").back());
            if (QFileInfo::exists(QString::fromStdString(this->mainFilename))) {
                this->loadPatchesFromFile(this->mainFilename);
            } else {
                this->savePatchesAsFile(this->mainFilename);
            }
        }
    });
    auto enableHotreloadButton = new CheckboxElement("Hot reloading ?", enableHotReloading);

    auto patchesCounterLabel = new LabelElement("Prim: " + std::to_string(nbPrimitives) + ", Unary: " + std::to_string(nbUnaryOperators) + ", Binary: " + std::to_string(nbBinaryOperators) + ", N-ary: " + std::to_string(nbNaryOperators));

    auto densityCombobox = new ComboboxElement("Density");
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Air", TerrainTypes::AIR));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Water", TerrainTypes::WATER));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Coral", TerrainTypes::CORAL));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Sand", TerrainTypes::SAND));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Dirt", TerrainTypes::DIRT));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Rock", TerrainTypes::ROCK));
    densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Bedrock", TerrainTypes::BEDROCK));
    for (size_t i = 0; i < densityCombobox->choices.size(); ++i)
        if (static_cast<ComboboxLineElement<TerrainTypes>*>(densityCombobox->choices[i])->value == selectedTerrainType)
            densityCombobox->combobox()->setCurrentIndex(int(i));
    densityCombobox->setOnSelectionChanged([=](int) { this->selectedTerrainType = densityCombobox->getSelection<TerrainTypes>(); });

    auto applyIntersectionButton = new CheckboxElement("Always intersect", applyIntersection);

    auto createStructureButton = new ButtonElement("Auto gen.", [=]() { this->structureAutoGeneration(); });

    primitiveSelectionGui = new HierarchicalListWidget;


    UI->add({
        constructionModeButton,
        shapeCombobox,
        densityCombobox,
        blendingFactorSlider,

        stackingButton,
        blendingButton,
        replacingButton,
        oneSideBlendButton,
        differenceButton,
        intersectionButton,

        abovePosButton,
        insideTopPosButton,
        insideBottomPosButton,
        fixedPosButton,

        widthSlider,
        depthSlider,
        heightSlider,
        sigmaSlider,

        createHorizontalGroup({applyIntersectionButton, createStructureButton}),
        createHorizontalGroup({addNoiseButton, addDistoButton, addSpreadButton, addRipplesButton, deformFromFlowButton}),
        resetButton,
        new UIElement(primitiveSelectionGui),

        createHorizontalGroup({selectedFilenameLabel, fileSelectionButton}),
        enableHotreloadButton,
        patchesCounterLabel
    });

    this->updatePrimitiveList();

//    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::currentItemChanged, this, &PrimitivePatchesInterface::updateSelectedPrimitiveItem);
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemClicked, this, [=](QListWidgetItem* item) { updateSelectedPrimitiveItem(item); });
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemDoubleClicked, this, &PrimitivePatchesInterface::openPrimitiveModificationDialog);
//    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemChangedHierarchy, this, &PrimitivePatchesInterface::modifyPrimitiveHierarchy);


    // This initialization should be in the constructor I guess...
    this->primitiveControlPoint = std::make_shared<ControlPoint>(Vector3(), 5.f);
    this->primitiveControlPoint->allowAllAxisRotations(true);
    this->primitiveControlPoint->allowAllAxisTranslation(true);
    this->primitiveControlPoint->hide();

    // QObject::connect(this->primitiveControlPoint.get(), &ControlPoint3D::pointModified, this, &PrimitivePatchesInterface::moveDebugBoxWithControlPoint);
    primitiveControlPoint->setOnPointModified([=]() { moveDebugBoxWithControlPoint(); });
    // QObject::connect(this->primitiveControlPoint.get(), &ControlPoint3D::translationApplied, this, &PrimitivePatchesInterface::translatePatch);
    primitiveControlPoint->setOnPointTranslated([=]() { translatePatch(primitiveControlPoint->state().transform.position - primitiveControlPoint->state().drag.startWorldPosition); });
    // QObject::connect(this->primitiveControlPoint.get(), &ControlPoint3D::rotationApplied, this, &PrimitivePatchesInterface::rotatePatch);
    primitiveControlPoint->setOnPointRotated([=]() { rotatePatch(primitiveControlPoint->state().transform.rotation.toVector3()); });
    return UI;
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
    if (!this->isVisible()) return;
    if (constructionMode && e->key() == Qt::Key::Key_Q && (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel)) {
        this->createPatchWithOperation(Vector3());
        this->parametricCurve.reset();
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
        dimensions.z() += 1.f;
        patch->setDimensions(dimensions);
        patch->update();
        this->updateMapWithCurrentPatch();
    } else if (e->key() == Qt::Key::Key_Minus) {
        Vector3 dimensions = patch->getDimensions();
        dimensions.z() -= 1.f;
        patch->setDimensions(dimensions);
        patch->update();
        this->updateMapWithCurrentPatch();
    }
}

void PrimitivePatchesInterface::clearSelection()
{
    currentlySelectedPatch = nullptr;
    debugMeshDisplayed = false;
    debuggingVoxelsMesh.clear();
    patchAABBoxMesh.clear();
    if (primitiveControlPoint) primitiveControlPoint->hide();
    storedPatches.clear();
    if (primitiveSelectionGui) primitiveSelectionGui->clear();
}

ImplicitPatch *PrimitivePatchesInterface::selectedPatch()
{
    if (currentlySelectedPatch && currentlySelectedPatch != implicitTerrain.get()) return currentlySelectedPatch;
    for (auto* child : implicitTerrain->composables) if (child) return child;
    return nullptr;
}

void PrimitivePatchesInterface::mouseMovedOnMapEvent(const Vector3& newPosition, TerrainModel* model)
{
    if (this->isVisible() && newPosition.isValid()) {
        this->currentPos = model ? model->getTerrainPos(newPosition) : newPosition;
        this->setSelectedShape(this->currentShapeSelected);
    }
}

void PrimitivePatchesInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel* model)
{
    if (this->isVisible() && mouseInMap) {

        if (constructionMode && (this->currentShapeSelected == ImplicitPatch::MountainChain || this->currentShapeSelected == ImplicitPatch::Polygon || this->currentShapeSelected == ImplicitPatch::ParametricTunnel)) {
            this->addParametricPoint(model->getTerrainPos(mousePosInMap));
            return;
        }

        if (this->primitiveControlPoint->grabsMouse())
            return; // Don't create patches if we try to displace another patch
        if (constructionMode)
            this->createPatchWithOperation(model->getTerrainPos(mousePosInMap));
    }
}

void PrimitivePatchesInterface::createPatchWithOperation(const Vector3& pos)
{
    setSelectedShape(currentShapeSelected);
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
    auto* originalParent = dynamic_cast<ImplicitNaryOperator*>(_parent);
    int originalSlot = -1;
    if (originalParent) {
        auto slot = std::find(originalParent->composables.begin(), originalParent->composables.end(), previousMain);
        if (slot == originalParent->composables.end()) throw std::logic_error("Missing selected implicit child");
        originalSlot = int(slot - originalParent->composables.begin());
    }
    ImplicitPatch* operation = this->createOperationPatchFromParameters(previousMain, patch); // *newOperation;

    if (_parent == nullptr) {
        this->implicitTerrain->deleteAllChildren();
        this->implicitTerrain->addChild(operation);
    } else {
        originalParent->addChild(operation, originalSlot);
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

    if (vertices != previewMesh.vertexArray) {
        this->previewMesh.fromArray(vertices);
        this->previewMesh.update();
    }
}

/*
void PrimitivePatchesInterface::setCurrentOperation(ImplicitPatch::CompositionFunction newOperation)
{
    this->currentOperation = newOperation;
}

void PrimitivePatchesInterface::setCurrentPositioning(ImplicitPatch::PositionalLabel positioning)
{
    this->currentPositioning = positioning;
}
*/

void PrimitivePatchesInterface::resetPatch()
{
    clearSelection();
    if (auto environment = std::dynamic_pointer_cast<EnvObjsInterface>(findOtherInterface("envobjects"))) environment->detachImplicitRoot();
//    delete this->implicitTerrain;
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
    // Implicit editing renders the procedural source directly. Other
    // representations can still be explicitly converted through the terrain UI.
    if (!viewer || viewer->mapMode != IMPLICIT_MODE) {
        this->layerGrid->fromImplicit(this->implicitTerrain.get());
        voxelGrid->fromImplicit(implicitTerrain.get());
    }
    refreshSelectedPatchPreview();
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
        if (!selectedPatch) return;
        noisePatch = new ImplicitNoise;
        wrapImplicitChild(selectedPatch, noisePatch);
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
        if (!selectedPatch) return;
        distortionPatch = new ImplicitWraping;
        wrapImplicitChild(selectedPatch, distortionPatch);
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
        if (!selectedPatch) return;
        distortionPatch = new ImplicitSpread;
        wrapImplicitChild(selectedPatch, distortionPatch);
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
        rotation->rotate(Vector3(0, 0, waterFlow.toEulerAngles().z()));
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

void PrimitivePatchesInterface::savePatchesAsFile(const std::string& filename)
{
    const auto payload = nlohmann::json({{ImplicitPatch::json_identifier, implicitTerrain->toJson()}}).dump(1, '\t');
    QSaveFile file(QString::fromStdString(filename));
    if (!file.open(QIODevice::WriteOnly) || file.write(payload.data(), qint64(payload.size())) != qint64(payload.size()) || !file.commit())
        throw std::runtime_error("Cannot atomically save implicit terrain: " + filename);
    lastTimeFileHasBeenModified = QFileInfo(QString::fromStdString(filename)).lastModified();
}

void PrimitivePatchesInterface::loadPatchesFromFile(const std::string& filename)
{
    this->loadTransformationRules();
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    nlohmann::json json_content = nlohmann::json::parse(content);
    this->lastTimeFileHasBeenModified = QFileInfo(QString::fromStdString(filename)).lastModified();
    if (json_content.contains(ImplicitPatch::json_identifier)) {
        std::unique_ptr<ImplicitPatch> decoded(ImplicitPatch::fromJson(json_content.at(ImplicitPatch::json_identifier)));
        // Prepare a fully typed copy before changing scene aliases.
        ImplicitNaryOperator replacement;
        replaceImplicitRoot(replacement, *decoded);
        clearSelection();
        if (auto environment = std::dynamic_pointer_cast<EnvObjsInterface>(findOtherInterface("envobjects")))
            environment->detachImplicitRoot();
        replaceImplicitRoot(*implicitTerrain, replacement);
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
    CatmullRomSpline transformedBarrierCurve;
    const auto barrierPoints = barrier->optionalCurve.getPath();
    for (size_t i = 0; i < barrierPoints.size(); i++) {
        Vector3 firstPointInCurve = barrierPoints[i];
        Vector3 transformedPoint = barrier->getGlobalPositionOf(firstPointInCurve);
        transformedBarrierCurve.addPoint(transformedPoint);
    }

    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    Contour transformedLagoonArea;
    const auto lagoonPoints = lagoon->optionalCurve.getPath();
    for (size_t i = 0; i < lagoonPoints.size(); i++) {
        Vector3 firstPointInCurve = lagoonPoints[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.curve->addPoint(transformedPoint.xy());
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
            barrierCurveNormals.at(pos) = Vector3::invalid;
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
        originalPoint.z() = this->implicitTerrain->getHeight(originalPoint);
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
    CatmullRomSpline transformedPasseCurve;
    const auto passePoints = passe->optionalCurve.getPath();
    for (size_t i = 0; i < passePoints.size(); i++) {
        Vector3 firstPointInCurve = passePoints[i];
        Vector3 transformedPoint = passe->getGlobalPositionOf(firstPointInCurve);
        transformedPasseCurve.addPoint(transformedPoint);
    }
    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    Contour transformedLagoonArea;
    const auto lagoonPoints = lagoon->optionalCurve.getPath();
    for (size_t i = 0; i < lagoonPoints.size(); i++) {
        Vector3 firstPointInCurve = lagoonPoints[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.curve->addPoint(transformedPoint.xy());
    }

    GridF sqrDistToPasse(voxelGrid->getDimensions(), -1.f);
    GridV3 passeDirections(voxelGrid->getDimensions());

    for (size_t i = 0; i < availablePositions.size(); i++) {
        const Vector3& pos = availablePositions[i];
        Vector3 passStarting = transformedPasseCurve.front();
        Vector3 passEnding = transformedPasseCurve.back();

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
            passeDirections.at(pos) = Vector3::invalid;
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
            originalPoint.z() = this->implicitTerrain->getHeight(originalPoint);
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
    std::vector<CatmullRomSpline> transformedPassesCurves(passes.size());
    for (size_t i = 0; i < passes.size(); i++) {
        ImplicitPrimitive* passe = dynamic_cast<ImplicitPrimitive*>(passes[i]);
        CatmullRomSpline transformedPasseCurve;
        const auto passePoints = passe->optionalCurve.getPath();
        for (size_t i = 0; i < passePoints.size(); i++) {
            Vector3 firstPointInCurve = passePoints[i];
            Vector3 transformedPoint = passe->getGlobalPositionOf(firstPointInCurve);
            transformedPasseCurve.addPoint(transformedPoint);
        }
        transformedPassesCurves[i] = transformedPasseCurve;
    }

    auto lagoons = this->implicitTerrain->findAll(ImplicitPatch::Polygon);
    ImplicitPrimitive* lagoon = dynamic_cast<ImplicitPrimitive*>(lagoons[0]);
    Contour transformedLagoonArea;
    const auto lagoonPoints = lagoon->optionalCurve.getPath();
    for (size_t i = 0; i < lagoonPoints.size(); i++) {
        Vector3 firstPointInCurve = lagoonPoints[i];
        Vector3 transformedPoint = lagoon->getGlobalPositionOf(firstPointInCurve);
        transformedLagoonArea.curve->addPoint(transformedPoint.xy());
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
    std::vector<CatmullRomSpline> transformedReefCurves(passes.size());
    for (size_t i = 0; i < reefs.size(); i++) {
        ImplicitPrimitive* reef = dynamic_cast<ImplicitPrimitive*>(reefs[i]);
        CatmullRomSpline transformedReefCurve;
        const auto reefPoints = reef->optionalCurve.getPath();
        for (size_t i = 0; i < reefPoints.size(); i++) {
            Vector3 firstPointInCurve = reefPoints[i];
            Vector3 transformedPoint = reef->getGlobalPositionOf(firstPointInCurve);
            transformedReefCurve.addPoint(transformedPoint);
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
    Vector3 surfaceNormal = Vector3::invalid;
    Vector3 waterVel = Vector3::invalid;

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
    Q_UNUSED(previous);
    auto* item = dynamic_cast<HierarchicalListWidgetItem<ImplicitPatch*>*>(current);
    auto* patch = item ? findPrimitiveById(item->ID) : nullptr;
    if (patch == currentlySelectedPatch) return;
    currentlySelectedPatch = patch;
    refreshSelectedPatchPreview();
    Q_EMIT updated();
}

void PrimitivePatchesInterface::refreshSelectedPatchPreview()
{
    debugMeshDisplayed = false;
    debuggingVoxelsMesh.clear();
    patchAABBoxMesh.clear();
    if (!currentlySelectedPatch) {
        if (primitiveControlPoint) primitiveControlPoint->hide();
        return;
    }
    const auto bounds = currentlySelectedPatch->getBBox();
    auto boxVertices = CubeMesh::createTriangles(bounds);
    for (auto& vertex : boxVertices)
        vertex = currentlySelectedPatch->getGlobalPositionOf(vertex - bounds.min());
    const AABBox sceneBounds(boxVertices);
    if (primitiveControlPoint) {
        primitiveControlPoint->setPosition(sceneBounds.min());
        primitiveControlPoint->show();
    }
    patchAABBoxMesh.fromArray(boxVertices);
    const auto support = currentlySelectedPatch->getSupportBBox();
    // Keep preview cost bounded, and leave air samples around bounded solids.
    const float spacing = std::max(.25f, support.dimensions().maxComp() / 24.f);
    const Vector3 padding(spacing * 2.f, spacing * 2.f, spacing * 2.f);
    auto mesh = currentlySelectedPatch->getGeometryInBounds(
        AABBox(support.min() - padding, support.max() + padding), Vector3(1, 1, 1) / spacing);
    for (auto& vertex : mesh.vertexArray)
        vertex = currentlySelectedPatch->getGlobalPositionOf(vertex - bounds.min());
    mesh.vertexArrayFloat = Vector3::toArray(mesh.vertexArray);
    mesh.computeNormals();
    mesh.shader = debuggingVoxelsMesh.shader;
    mesh.cullFace = false;
    debuggingVoxelsMesh = std::move(mesh);
    debugMeshDisplayed = true;
}

void PrimitivePatchesInterface::openPrimitiveModificationDialog(QListWidgetItem *item) {
    auto selectedPatchItem = dynamic_cast<HierarchicalListWidgetItem<ImplicitPatch*>*>(item);
    int selectedPatchID = selectedPatchItem->ID;

    ImplicitPatch* selectedPatch = this->findPrimitiveById(selectedPatchID);
    for (const auto& p : storedPatches) {
        if (p->index == selectedPatchID) {
            selectedPatch = p;
            break;
        }
    }
    if (selectedPatch && selectedPatch->parent &&
        (dynamic_cast<ImplicitPrimitive*>(selectedPatch) || dynamic_cast<ImplicitBinaryOperator*>(selectedPatch))) {
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

void PrimitivePatchesInterface::setMainFilename(const std::string& newMainFilename)
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

        if (auto* group = dynamic_cast<ImplicitNaryOperator*>(current))
            for (auto* child : group->composables) if (child) pendingPatches.push_back(child);
    }
}

void PrimitivePatchesInterface::loadTransformationRules()
{
    displayProcessTime("Reading and interpreting file " + this->rulesFilename + "... ", [&]() {
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
    });
}

void PrimitivePatchesInterface::addParametricPoint(const Vector3& point)
{
    this->parametricCurve.addPoint(point);
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
    if (!implicitTerrain) return;
    auto& plot = ImageViewer::get("implicit-tree");
    plot.dataModel->reset();

    // Operators share the owning child API, including unary and fixed-arity
    // operators with vacant slots. A root has no incoming edge.
    std::vector<ImplicitPatch*> nodes{implicitTerrain.get()};
    std::vector<int> parents{-1};
    std::vector<size_t> depths{0};
    std::vector<size_t> levelCounts{1};
    for (size_t cursor = 0; cursor < nodes.size(); ++cursor) {
        auto* group = dynamic_cast<ImplicitNaryOperator*>(nodes[cursor]);
        if (!group) continue;
        const size_t depth = depths[cursor] + 1;
        if (levelCounts.size() <= depth) levelCounts.resize(depth + 1);
        for (auto* child : group->composables) {
            if (!child) continue;
            nodes.push_back(child);
            parents.push_back(int(cursor));
            depths.push_back(depth);
            ++levelCounts[depth];
        }
    }

    std::vector<Vector3> positions;
    std::vector<std::string> labels;
    std::vector<size_t> levelOffsets(levelCounts.size(), 0);
    for (size_t i = 0; i < nodes.size(); ++i) {
        const size_t depth = depths[i];
        positions.emplace_back(float(levelOffsets[depth]++) - float(levelCounts[depth] - 1) * .5f, -float(depth));
        std::string label = nodes[i]->name;
        if (auto* primitive = dynamic_cast<ImplicitPrimitive*>(nodes[i]))
            label += "\n(" + toCapitalize(stringFromMaterial(primitive->material)) + ")";
        labels.push_back(label);
    }
    for (size_t i = 1; i < nodes.size(); ++i)
        plot.addPlot({positions[size_t(parents[i])], positions[i]});
    plot.addScatter(positions, "", labels);
    plot.draw();
    plot.show();
}

void PrimitivePatchesInterface::moveDebugBoxWithControlPoint()
{
    if (!currentlySelectedPatch || !primitiveControlPoint || patchAABBoxMesh.vertexArray.empty()) return;
    const auto bounds = AABBox(patchAABBoxMesh.vertexArray);
    patchAABBoxMesh.translate(primitiveControlPoint->getPosition() - bounds.min());
}

void PrimitivePatchesInterface::translatePatch(const Vector3& translation)
{
    auto* selected = selectedPatch();
    if (!selected) return;
    QSignalBlocker blocker(primitiveSelectionGui);
    auto* operation = dynamic_cast<ImplicitTranslation*>(selected);
    if (!operation) operation = dynamic_cast<ImplicitTranslation*>(selected->parent);
    if (!operation) {
        operation = new ImplicitTranslation;
        wrapImplicitChild(selected, operation);
    }
    operation->translate(translation);
    currentlySelectedPatch = operation;
    updateMapWithCurrentPatch();
    updatePrimitiveList();
    primitiveSelectionGui->setCurrentItem(operation->index);
}

void PrimitivePatchesInterface::rotatePatch(const Vector3& rotation)
{
    auto* selected = selectedPatch();
    if (!selected) return;
    QSignalBlocker blocker(primitiveSelectionGui);
    auto* operation = dynamic_cast<ImplicitRotation*>(selected);
    if (!operation) operation = dynamic_cast<ImplicitRotation*>(selected->parent);
    if (!operation) {
        operation = new ImplicitRotation;
        wrapImplicitChild(selected, operation);
    }
    operation->rotate(rotation);
    currentlySelectedPatch = operation;
    updateMapWithCurrentPatch();
    updatePrimitiveList();
    primitiveSelectionGui->setCurrentItem(operation->index);
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
        // auto [BBoxMin, BBoxMax] = patch->getBBox();
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
            auto AABBox = this->parametricCurve.getAABBox();
            finalPos = AABBox.min() - funcSize * .5f;
            funcSize = (AABBox.max() + funcSize * .5f) - finalPos;
        }
        CatmullRomSpline translatedSpline = this->parametricCurve;
        for (auto& p : translatedSpline)
            p -= finalPos;
        ImplicitPrimitive* primitive = (ImplicitPrimitive*) ImplicitPatch::createPredefinedShape(currentShapeSelected, funcSize, selectedSigma, translatedSpline);
        primitive->position = finalPos.xy();
        primitive->setDimensions(funcSize);
        primitive->predefinedShape = this->currentShapeSelected;
        primitive->parametersProvided = {this->selectedSigma};

        primitive->material = this->selectedTerrainType;
        primitive->name = toCapitalize(stringFromPredefinedShape(currentShapeSelected));
        if (this->currentPositioning == ImplicitPatch::PositionalLabel::FIXED_POS && finalPos.z() != 0.f) {
            ImplicitTranslation* translate = new ImplicitTranslation;
            translate->translate(Vector3(0.f, 0.f, finalPos.z()));
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
    if (!child) return nullptr;
    if (child->parent != nullptr)
        return child->parent;

    // Not sure if the rest is useful, it's an old code...
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
    if (this->isVisible() && primitiveSelectionGui) {
        primitiveSelectionGui->clear();
        this->storedPatches.clear();
        bool patchesCanBeLocked = this->mainFilename != "";
        std::vector<std::tuple<ImplicitPatch*, int, bool>> waitingPatches = { {this->implicitTerrain.get(), 0, false} };

        this->nbPrimitives = 0;
        this->nbBinaryOperators = 0;
        this->nbUnaryOperators = 0;
        this->nbNaryOperators = 0;

        while (!waitingPatches.empty()) {
            auto [current, level, locked] = waitingPatches.back();
            if (current != this->implicitTerrain.get() && current->parent == nullptr)
                    std::cout << current->toString() << " has no parent" << std::endl;
            bool childrenAreLocked = locked || (current->used_json_filename != "" && patchesCanBeLocked && current->used_json_filename != this->mainFilename);
            waitingPatches.pop_back();
            primitiveSelectionGui->addItem(new HierarchicalListWidgetItem<ImplicitPatch*>((locked ? "*" : "") + current->toString(), current, current->index, level));
            this->storedPatches.push_back(current);
            ImplicitBinaryOperator* currentAsBinary = dynamic_cast<ImplicitBinaryOperator*>(current);
            ImplicitUnaryOperator* currentAsUnary = dynamic_cast<ImplicitUnaryOperator*>(current);
            ImplicitNaryOperator* currentAsNaryOperator = dynamic_cast<ImplicitNaryOperator*>(current);
            if (currentAsNaryOperator) {
                for (const auto& c : currentAsNaryOperator->composables)
                    if (c) waitingPatches.push_back({c, level+1, childrenAreLocked});
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
        this->findAllSubfiles();
        if (currentlySelectedPatch) primitiveSelectionGui->setCurrentItem(currentlySelectedPatch->index);
    }
}


void PrimitivePatchesInterface::displayDebuggingVoxels()
{
    if (!debugMeshDisplayed || !debuggingVoxelsMesh.shader) return;
    debuggingVoxelsMesh.shader->setVector("color", Vector4(0.f, .6f, .8f, .5f));
    debuggingVoxelsMesh.display(GL_TRIANGLES);
}


PatchReplacementDialog::PatchReplacementDialog(PrimitivePatchesInterface* caller, ImplicitPatch *patchToModify)
    : QDialog(), caller(caller), patch(patchToModify)
{
    if (!patchToModify) throw std::invalid_argument("Missing implicit primitive to edit");
    originalPatch = patchToModify;
    originalIndex = patchToModify->index;
    draftPatch.reset(patchToModify->copy());
    patch = draftPatch.get();
    bool patchIsOperation = (dynamic_cast<ImplicitBinaryOperator*>(patch) != nullptr);
    if (auto* binary = dynamic_cast<ImplicitBinaryOperator*>(patch)) firstDraftChild = binary->composableA();

    auto UI = new InterfaceUI();

    if (patchIsOperation) {
        ImplicitBinaryOperator* patchAsBinary = dynamic_cast<ImplicitBinaryOperator*>(patch);
        auto stackingButton = new RadioButtonElement("Stack", ImplicitPatch::CompositionFunction::STACK, patchAsBinary->composeFunction);
        auto blendingButton = new RadioButtonElement("Blend", ImplicitPatch::CompositionFunction::BLEND, patchAsBinary->composeFunction);
        auto replacingButton = new RadioButtonElement("Replace", ImplicitPatch::CompositionFunction::REPLACE, patchAsBinary->composeFunction);
        auto oneSideBlendButton = new RadioButtonElement("1-blend", ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND, patchAsBinary->composeFunction);
        auto differenceButton = new RadioButtonElement("A minus B", ImplicitPatch::DIFFERENCE, patchAsBinary->composeFunction);
        auto intersectionButton = new RadioButtonElement("A intersect B", ImplicitPatch::INTERSECTION, patchAsBinary->composeFunction);

        auto abovePosButton = new RadioButtonElement("Above", ImplicitPatch::PositionalLabel::ABOVE, patchAsBinary->positionalB);
        auto insideTopPosButton = new RadioButtonElement("Inside top", ImplicitPatch::PositionalLabel::INSIDE_TOP, patchAsBinary->positionalB);
        auto insideBottomPosButton = new RadioButtonElement("Inside bottom", ImplicitPatch::PositionalLabel::INSIDE_BOTTOM, patchAsBinary->positionalB);
        auto fixedPosButton = new RadioButtonElement("Fixed", ImplicitPatch::PositionalLabel::FIXED_POS, patchAsBinary->positionalB);

        auto blendingFactorSlider = new SliderElement("Blend factor", 0.1f, 10.f, .1f, patchAsBinary->blendingFactor);
        auto applyIntersectionButton = new CheckboxElement("Always intersect", patchAsBinary->withIntersectionOnB);
        auto swapButton = new ButtonElement("Swap A and B", [=]() { patchAsBinary->swapAB(); });

        /*
        stackingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::STACK).setOnChecked([=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::STACK; });
        blendingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::BLEND).setOnChecked([=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::BLEND; });
        replacingButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::REPLACE).setOnChecked([=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::REPLACE; });
        oneSideBlendButton->setChecked(patchAsBinary->composeFunction == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND).setOnChecked([=](bool checked) { if (checked) patchAsBinary->composeFunction = ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND; });


        abovePosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::ABOVE)
            .setOnChecked([=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::ABOVE; });
        insideTopPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::INSIDE_TOP)
            .setOnChecked([=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::INSIDE_TOP; });
        insideBottomPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM)
            .setOnChecked([=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::INSIDE_BOTTOM; });
        fixedPosButton->setChecked(patchAsBinary->positionalB == ImplicitPatch::PositionalLabel::FIXED_POS)
            .setOnChecked([=](bool checked) { if (checked) patchAsBinary->positionalB = ImplicitPatch::PositionalLabel::FIXED_POS; });
        */

        UI->add({
            blendingFactorSlider,

            stackingButton,
            blendingButton,
            replacingButton,
            oneSideBlendButton,
            differenceButton,
            intersectionButton,

            abovePosButton,
            insideTopPosButton,
            insideBottomPosButton,
            fixedPosButton,

            createHorizontalGroup({ applyIntersectionButton, swapButton })
        });

    } else {
        ImplicitPrimitive* patchAsPrimitive = dynamic_cast<ImplicitPrimitive*>(patch);


        auto shapeCombobox = new ComboboxElement("Shape");
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Sphere", ImplicitPatch::Sphere), patchAsPrimitive->predefinedShape == ImplicitPatch::Sphere);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Block", ImplicitPatch::Block), patchAsPrimitive->predefinedShape == ImplicitPatch::Block);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Gaussian", ImplicitPatch::Gaussian), patchAsPrimitive->predefinedShape == ImplicitPatch::Gaussian);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Cylinder", ImplicitPatch::Cylinder), patchAsPrimitive->predefinedShape == ImplicitPatch::Cylinder);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Rock", ImplicitPatch::Rock), patchAsPrimitive->predefinedShape == ImplicitPatch::Rock);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Mountain", ImplicitPatch::Mountain), patchAsPrimitive->predefinedShape == ImplicitPatch::Mountain);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Dune", ImplicitPatch::Dune), patchAsPrimitive->predefinedShape == ImplicitPatch::Dune);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Basin", ImplicitPatch::Basin), patchAsPrimitive->predefinedShape == ImplicitPatch::Basin);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Cave", ImplicitPatch::Cave), patchAsPrimitive->predefinedShape == ImplicitPatch::Cave);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Arch", ImplicitPatch::Arch), patchAsPrimitive->predefinedShape == ImplicitPatch::Arch);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Noise", ImplicitPatch::Noise2D), patchAsPrimitive->predefinedShape == ImplicitPatch::Noise2D);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Mountains", ImplicitPatch::MountainChain), patchAsPrimitive->predefinedShape == ImplicitPatch::MountainChain);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Polygon", ImplicitPatch::Polygon), patchAsPrimitive->predefinedShape == ImplicitPatch::Polygon);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Tunnel", ImplicitPatch::ParametricTunnel), patchAsPrimitive->predefinedShape == ImplicitPatch::ParametricTunnel);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("Ripple", ImplicitPatch::Ripple), patchAsPrimitive->predefinedShape == ImplicitPatch::Ripple);
        shapeCombobox->addChoice(new ComboboxLineElement<ImplicitPatch::PredefinedShapes>("...", ImplicitPatch::ImplicitHeightmap), patchAsPrimitive->predefinedShape == ImplicitPatch::ImplicitHeightmap);
        shapeCombobox->setOnSelectionChanged([=](int) {
            auto shape = shapeCombobox->getSelection<ImplicitPatch::PredefinedShapes>();
            // this->setSelectedShape(shape);
            patchAsPrimitive->predefinedShape = shape;
            if (shape == ImplicitPatch::ImplicitHeightmap) {
                QString q_filename = QDir(QDir::currentPath()).relativeFilePath(QFileDialog::getOpenFileName(this, QString("Open a composition"), QString::fromStdString("saved_maps/")));
                std::string filename = q_filename.toStdString();
                if (!q_filename.isEmpty()) {
                    patchAsPrimitive->used_json_filename = filename;
                }
            }
        });

        auto widthSlider = new SliderElement("Width / radius", 0.f, 50.f, .1f, patchAsPrimitive->dimensions.x());
        auto depthSlider = new SliderElement("Depth", 0.f, 50.f, .1f, patchAsPrimitive->dimensions.y());
        auto heightSlider = new SliderElement("Height", 0.f, 50.f, .1f, patchAsPrimitive->dimensions.z());
        if (patchAsPrimitive->parametersProvided.empty()) patchAsPrimitive->parametersProvided.push_back(1.f);
        auto sigmaSlider = new SliderElement("Sigma", 0.f, 20.f, .1f, patchAsPrimitive->parametersProvided[0]);

        auto densityCombobox = new ComboboxElement("Density");
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Air", TerrainTypes::AIR), patchAsPrimitive->material == TerrainTypes::AIR);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Water", TerrainTypes::WATER), patchAsPrimitive->material == TerrainTypes::WATER);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Coral", TerrainTypes::CORAL), patchAsPrimitive->material == TerrainTypes::CORAL);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Sand", TerrainTypes::SAND), patchAsPrimitive->material == TerrainTypes::SAND);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Dirt", TerrainTypes::DIRT), patchAsPrimitive->material == TerrainTypes::DIRT);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Rock", TerrainTypes::ROCK), patchAsPrimitive->material == TerrainTypes::ROCK);
        densityCombobox->addChoice(new ComboboxLineElement<TerrainTypes>("Bedrock", TerrainTypes::BEDROCK), patchAsPrimitive->material == TerrainTypes::BEDROCK);
        densityCombobox->setOnSelectionChanged([=](int) {
            patchAsPrimitive->material = densityCombobox->getSelection<TerrainTypes>();
            patchAsPrimitive->sampledMaterialColumns = ImplicitPrimitive::MaterialColumns();
        });


        UI->add({
            shapeCombobox,
            densityCombobox,

            widthSlider,
            depthSlider,
            heightSlider,
            sigmaSlider
        });

        /*
        auto createSphereButton = new RadioButtonElement("Sphere"));
        auto createBlockButton = new RadioButtonElement("Block"));
        auto createGaussianButton = new RadioButtonElement("Gaussian"));
        auto createCylinderButton = new RadioButtonElement("Cylinder"));
        auto createRockButton = new RadioButtonElement("Rock"));
        auto createMountainButton = new RadioButtonElement("Mountain"));
        auto createDuneButton = new RadioButtonElement("Dune"));
        auto createBasinButton = new RadioButtonElement("Basin"));
        auto createCaveButton = new RadioButtonElement("Cave"));
        auto createArchButton = new RadioButtonElement("Arch"));
        auto createNoise2DButton = new RadioButtonElement("Noise"));
        auto createFromFileButton = new RadioButtonElement("From file"));

        auto airDensityCheckbox = new RadioButtonElement("Air"));
        auto waterDensityCheckbox = new RadioButtonElement("Water"));
        auto coralDensityCheckbox = new RadioButtonElement("Coral"));
        auto sandDensityCheckbox = new RadioButtonElement("Sand"));
        auto dirtDensityCheckbox = new RadioButtonElement("Dirt"));
        auto rockDensityCheckbox = new RadioButtonElement("Rock"));
        auto bedrockDensityCheckbox = new RadioButtonElement("Bedrock"));

        auto widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        auto depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        auto heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 50.f, .1f);
        auto sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 20.f, .1f);

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

        widthSlider->setfValue(patchAsPrimitive->dimensions.x());
        heightSlider->setfValue(patchAsPrimitive->dimensions.z());
        depthSlider->setfValue(patchAsPrimitive->dimensions.y());
        sigmaSlider->setfValue(patchAsPrimitive->parametersProvided[0]);

        QObject::connect(createSphereButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Sphere; });
        QObject::connect(createBlockButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Block; });
        QObject::connect(createGaussianButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Gaussian; });
        QObject::connect(createCylinderButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Cylinder; });
        QObject::connect(createRockButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Rock; });
        QObject::connect(createMountainButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Mountain; });
        QObject::connect(createDuneButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Dune; });
        QObject::connect(createBasinButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Basin; });
        QObject::connect(createCaveButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Cave; });
        QObject::connect(createArchButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Arch; });
        QObject::connect(createNoise2DButton, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->predefinedShape = ImplicitPatch::PredefinedShapes::Noise2D; });
        QObject::connect(createFromFileButton, &RadioButtonElement::toggled, this, [=](bool checked) {
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
        QObject::connect(airDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = AIR; });
        QObject::connect(waterDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = WATER; });
        QObject::connect(coralDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = CORAL; });
        QObject::connect(sandDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = SAND; });
        QObject::connect(dirtDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = DIRT; });
        QObject::connect(rockDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = ROCK; });
        QObject::connect(bedrockDensityCheckbox, &RadioButtonElement::toggled, this, [=](bool checked) { if (checked) patchAsPrimitive->material = BEDROCK; });


        QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.x() = newVal; });
        QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.z() = newVal; });
        QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->dimensions.y() = newVal; });
        QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [=](float newVal) { patchAsPrimitive->parametersProvided[0] = newVal; });
        */

    }



    auto cancelButton = new ButtonElement("Annuler", [=]() { this->cancel(); });
    auto validButton = new ButtonElement("Confirmer", [=]() { this->confirm(); });
    auto error = new LabelElement("");
    validationError = error->label();
    validationError->setWordWrap(true);
    validationError->hide();

    UI->add({error, createHorizontalGroup({cancelButton, validButton})});
    setLayout(UI->get()->layout());
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
    try {
        validationError->clear();
        validationError->hide();
        if (!originalPatch) throw std::logic_error("This edit has already finished");
        if (caller) {
            // Hot reload can destroy the original while this dialog's event loop
            // remains active. Compare identities against live nodes before touching
            // the saved pointer; IDs also distinguish an address reused by new data.
            bool found = false;
            std::vector<ImplicitPatch*> pending;
            if (caller->implicitTerrain) pending.push_back(caller->implicitTerrain.get());
            while (!pending.empty()) {
                auto* node = pending.back();
                pending.pop_back();
                if (node == originalPatch && node->index == originalIndex) { found = true; break; }
                if (auto* group = dynamic_cast<ImplicitNaryOperator*>(node))
                    for (auto* child : group->composables) if (child) pending.push_back(child);
            }
            if (!found) throw std::logic_error("This primitive was removed or reloaded. Close this editor and select it again.");
        }
        patch->update();
        if (!originalPatch || !originalPatch->parent) throw std::logic_error("Edited implicit patch is detached");
        auto* parent = originalPatch->parent;
        if (std::find(parent->composables.begin(), parent->composables.end(), originalPatch) == parent->composables.end())
            throw std::logic_error("Edited implicit patch is missing from its parent");
        // Scene objects retain aliases to primitive nodes and descendants. Commit
        // parameters in place so editing a parent never replaces their identities.
        if (auto* target = dynamic_cast<ImplicitPrimitive*>(originalPatch)) {
            auto* source = dynamic_cast<ImplicitPrimitive*>(patch);
            if (!source || typeid(*source) != typeid(*target)) throw std::invalid_argument("Edited primitive type changed");
            const int index = target->index;
            *target = *source;
            target->parent = parent;
            target->index = index;
            target->update();
        } else if (auto* target = dynamic_cast<ImplicitBinaryOperator*>(originalPatch)) {
            auto* source = dynamic_cast<ImplicitBinaryOperator*>(patch);
            if (!source) throw std::invalid_argument("Edited operation type changed");
            if (source->composableA() != firstDraftChild) target->swapAB();
            target->composeFunction = source->composeFunction;
            target->positionalB = source->positionalB;
            target->blendingFactor = source->blendingFactor;
            target->withIntersectionOnB = source->withIntersectionOnB;
            target->name = source->name;
            target->mirrored = source->mirrored;
            target->optionalCurve = source->optionalCurve;
        } else throw std::invalid_argument("Unsupported implicit parameter editor node");
        originalPatch->used_json_filename.clear();
        originalPatch->referenceProperties = nlohmann::json::object();
        for (ImplicitPatch* node = originalPatch; node; node = node->parent) node->_cached = false;
        if (caller) caller->clearSelection();
        originalPatch = nullptr;
        setResult(1);
        done(1);
    } catch (const std::exception& error) {
        validationError->setText(QString::fromUtf8(error.what()));
        validationError->show();
    }

}
