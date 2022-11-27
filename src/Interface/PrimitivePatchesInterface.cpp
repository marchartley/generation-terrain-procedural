#include "PrimitivePatchesInterface.h"
#include "Interface/InterfaceUtils.h"
#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"

PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("PrimitivePatchesInterface", parent)
{
    previewMesh.cullFace = false;
}

void PrimitivePatchesInterface::display()
{
    if (this->isVisible()) {
        if (previewMesh.shader != nullptr) {
            this->setSelectedShape(this->currentShapeSelected);
            previewMesh.shader->setVector("color", std::vector<float>({0.f, 0.4f, 0.8f, 0.5f}));
            previewMesh.display();
        }
    }
}

void PrimitivePatchesInterface::reloadShaders()
{

}

void PrimitivePatchesInterface::replay(nlohmann::json action)
{

}

QLayout *PrimitivePatchesInterface::createGUI()
{

    QVBoxLayout* layout = new QVBoxLayout();

    QRadioButton* createSphereButton = new QRadioButton("Sphere");
    QRadioButton* createBlockButton = new QRadioButton("Block");
    QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");

    FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);

    QPushButton* resetButton = new QPushButton("Reset");

    layout->addWidget(createVerticalGroup({
                                                    createSphereButton,
                                                    createBlockButton,
                                                    createGaussianButton
                                                }));
    layout->addWidget(createSliderGroup("Density", densitySlider));
    layout->addWidget(createVerticalGroup({
                                                    stackingButton,
                                                    blendingButton,
                                                    replacingButton
                                                }));
    layout->addWidget(createMultipleSliderGroup({
                                              {"Width / radius", widthSlider},
                                                    {"Depth", depthSlider},
                                                    {"Height", heightSlider},
                                                    {"Sigma", sigmaSlider}
                                          }));
    layout->addWidget(resetButton);

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(PRIMITIVE_SHAPE::SPHERE); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(PRIMITIVE_SHAPE::BLOCK); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [&]() { this->setSelectedShape(PRIMITIVE_SHAPE::GAUSSIAN); });
    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [&](float newDensity) { this->selectedDensity = newDensity; });
    QObject::connect(stackingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(PATCH_OPERATION::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(PATCH_OPERATION::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [&]() { this->setCurrentOperation(PATCH_OPERATION::REPLACE); });

    QObject::connect(widthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedWidth = newVal; });
    QObject::connect(heightSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedHeight = newVal; });
    QObject::connect(depthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedDepth = newVal; });
    QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->selectedSigma = newVal; });

    QObject::connect(resetButton, &QPushButton::pressed, this, &PrimitivePatchesInterface::resetPatch);

    createSphereButton->setChecked(this->currentShapeSelected == PRIMITIVE_SHAPE::SPHERE);
    createBlockButton->setChecked(this->currentShapeSelected == PRIMITIVE_SHAPE::BLOCK);
    createGaussianButton->setChecked(this->currentShapeSelected == PRIMITIVE_SHAPE::GAUSSIAN);

    stackingButton->setChecked(this->currentOperation == PATCH_OPERATION::STACK);
    blendingButton->setChecked(this->currentOperation == PATCH_OPERATION::BLEND);
    replacingButton->setChecked(this->currentOperation == PATCH_OPERATION::REPLACE);

    densitySlider->setfValue(this->selectedDensity);
    widthSlider->setfValue(this->selectedWidth);
    heightSlider->setfValue(this->selectedHeight);
    depthSlider->setfValue(this->selectedDepth);
    sigmaSlider->setfValue(this->selectedSigma);

    return layout;
}

void PrimitivePatchesInterface::show()
{
    this->previewMesh.show();
    ActionInterface::show();
}

void PrimitivePatchesInterface::hide()
{
    this->previewMesh.hide();
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
        this->layerGrid->layers = this->layerGrid->previousState;

        Patch3D* patch = new Patch3D(mousePosInMap - (this->functionSize.xy() * .5f), Vector3(), this->functionSize, this->currentSelectionFunction, this->selectedDensity);
        Patch3D thisResult;

        if (storedPatches.size() > 0) {
            Patch3D* previousPatch = storedPatches.back();
            if (currentOperation == PATCH_OPERATION::STACK)
                thisResult = Patch3D::stack(previousPatch, patch);
            if (currentOperation == PATCH_OPERATION::BLEND)
                thisResult = Patch3D::blend(previousPatch, patch);
            if (currentOperation == PATCH_OPERATION::REPLACE)
                thisResult = Patch3D::replace(previousPatch, patch);

            currentPatch = new Patch3D(thisResult);
        } else {
            currentPatch = patch;
            thisResult = *patch;
        }
        storedPatches.push_back(currentPatch);

        this->layerGrid->add(*currentPatch, TerrainTypes::AIR, false);
        auto newVoxels = layerGrid->voxelize(this->voxelGrid->getSizeZ(), 3.f);
        voxelGrid->applyModification(newVoxels - voxelGrid->getVoxelValues());
        heightmap->fromLayerGrid(*this->layerGrid);
    }
}

void PrimitivePatchesInterface::setSelectedShape(PRIMITIVE_SHAPE newShape, Vector3 newPosition)
{
    newPosition = this->currentPos;
    this->currentShapeSelected = newShape;

    if (this->currentShapeSelected == PRIMITIVE_SHAPE::SPHERE) {
        this->functionSize = Vector3(selectedWidth, selectedWidth, selectedWidth);
        auto vertices = CubeMesh::cubesVertices;
        for (auto& vert : vertices) {
            vert *= functionSize;
            vert += newPosition;
            vert -= functionSize.xy() * .5f;
        }
        this->previewMesh.fromArray(vertices);
        this->previewMesh.update();

        this->currentSelectionFunction = this->sphereFunction(selectedWidth * .5f);
    }
    if (this->currentShapeSelected == PRIMITIVE_SHAPE::BLOCK) {
        this->functionSize = Vector3(selectedWidth, selectedDepth, selectedHeight);
        auto vertices = CubeMesh::cubesVertices;
        for (auto& vert : vertices) {
            vert *= functionSize;
            vert += newPosition;
            vert -= functionSize.xy() * .5f;
        }
        this->previewMesh.fromArray(vertices);
        this->previewMesh.update();

        this->currentSelectionFunction = this->blockFunction(selectedHeight);
    }
    if (this->currentShapeSelected == PRIMITIVE_SHAPE::GAUSSIAN) {
        this->functionSize = Vector3(selectedWidth, selectedDepth, selectedHeight);
        auto vertices = CubeMesh::cubesVertices;
        for (auto& vert : vertices) {
                vert *= functionSize;
                vert += newPosition;
                vert -= functionSize.xy() * .5f;
            }
        this->previewMesh.fromArray(vertices);
        this->previewMesh.update();

        this->currentSelectionFunction = this->gaussianFunction(selectedSigma, selectedHeight);
    }
}

void PrimitivePatchesInterface::setCurrentOperation(PATCH_OPERATION newOperation)
{
    this->currentOperation = newOperation;
}

void PrimitivePatchesInterface::resetPatch()
{
    this->storedPatches.clear();
    this->layerGrid->layers = this->layerGrid->previousState;
    voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
    voxelGrid->fromIsoData();
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
