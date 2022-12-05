#include "PrimitivePatchesInterface.h"
#include "Interface/InterfaceUtils.h"
#include "Graphics/Sphere.h"
#include "Graphics/CubeMesh.h"

PrimitivePatchesInterface::PrimitivePatchesInterface(QWidget *parent)
    : ActionInterface("PrimitivePatchesInterface", parent)
{
    previewMesh.cullFace = false;

    this->mainPatch = new ImplicitPatch;
    this->mainPatch->name = "Identity";
}

void PrimitivePatchesInterface::display()
{
    if (this->isVisible()) {
        if (previewMesh.shader != nullptr) {
            this->setSelectedShape(this->currentShapeSelected);
            previewMesh.shader->setVector("color", std::vector<float>({0.f, 0.4f, 0.8f, 0.5f}));
            previewMesh.display();
//            if (patchAABBoxMesh.shader == nullptr) {
//                patchAABBoxMesh.shader = previewMesh.shader;
//            }
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

    primitiveSelectionGui = new HierarchicalListWidget(this);

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
    layout->addWidget(primitiveSelectionGui);

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

    this->updatePrimitiveList();

    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::currentItemChanged, this, [&](QListWidgetItem *current, QListWidgetItem *previous) {
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
            }
        }
        if (!newSelectionIsExisting) {
            this->patchAABBoxMesh.fromArray(std::vector<Vector3>{});
            this->patchAABBoxMesh.update();
            this->primitiveControlPoint->hide();
        }
        Q_EMIT updated();
    });
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemDoubleClicked, this, [&](QListWidgetItem* item) -> void {
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
            if (result > 0) { // TODO : find why "0" is returned
                this->updateMapWithCurrentPatch();
            } else {
                this->updateMapWithCurrentPatch();
            }
        }
        updatePrimitiveList();
    });
    QObject::connect(primitiveSelectionGui, &HierarchicalListWidget::itemChangedHierarchy, this, [&] (int ID_to_move, int relatedID, HIERARCHY_TYPE relation, QDropEvent* event) {
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
        /*
        if (!createCopy) {
            // Remove the instance from the parent's tree (and save the index to possibly plug another biome to the parent)
            auto placeToInsert = previousParent->instances.erase(std::find(previousParent->instances.begin(), previousParent->instances.end(), toMove));
            // If we moved a biome in a lower level of his hierarchy...
            if (toMove->getPathToChild(related).size() > 1) {
                // Get the first child that leads to the target and change his parent
                auto child = toMove->getPathToChild(related)[1];
                child->parent = previousParent;
                previousParent->instances.insert(placeToInsert, child);
            }
        }

//        toMove->parent = related;
//        related->instances.push_back(toMove);
        related->addInstance(toMove);

        previousParent->updateSubInstances();
//        this->rootBiome = this->rootBiome->clone(rootBiome->area, rootBiome->position);
        this->updatePrimitiveList();
        */
    });


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

ImplicitPatch* PrimitivePatchesInterface::createPatchFromParameters(Vector3 position, ImplicitPatch *replacedPatch)
{
//    if (this->currentOperation == PATCH_OPERATION::STACK) // Remove the z offset when stacking... Should be done in the Blending behavior?
//        position.z = 0.f;


    ImplicitPatch* patch;
    if (replacedPatch != nullptr)
        patch = replacedPatch;
    else
        patch = new ImplicitPatch;

    switch (this->currentShapeSelected) {
    case SPHERE:
        *patch = ImplicitPatch(position - this->functionSize.xy() * .5f, Vector3(), this->functionSize, ImplicitPatch::createSphereFunction(selectedWidth * .5f));
        break;
    case BLOCK:
        *patch = ImplicitPatch(position - this->functionSize.xy() * .5f, Vector3(), this->functionSize, ImplicitPatch::createBlockFunction(selectedWidth, selectedDepth, selectedHeight));
        break;
    case GAUSSIAN:
        *patch = ImplicitPatch(position - this->functionSize.xy() * .5f, Vector3(), this->functionSize, ImplicitPatch::createGaussianFunction(selectedSigma, selectedWidth, selectedDepth, selectedHeight));
        break;
    }
    patch->densityValue = this->selectedDensity;
    patch->name = (currentShapeSelected == PRIMITIVE_SHAPE::BLOCK ? "Block" : currentShapeSelected == PRIMITIVE_SHAPE::GAUSSIAN ? "Gaussian" : "Sphere");
    return patch;
}

ImplicitPatch *PrimitivePatchesInterface::createOperationPatchFromParameters(ImplicitPatch *composableA, ImplicitPatch *composableB, ImplicitPatch *replacedPatch)
{
    ImplicitPatch* operation = nullptr;

    if (currentOperation == PATCH_OPERATION::STACK) {
        operation = new ImplicitPatch(ImplicitPatch::createStack(composableA, composableB));
        operation->name = "Stacking";
    } else if (currentOperation == PATCH_OPERATION::BLEND) {
        operation = new ImplicitPatch(ImplicitPatch::createBlending(composableA, composableB));
        operation->name = "Blending";
    } else if (currentOperation == PATCH_OPERATION::REPLACE) {
        operation = new ImplicitPatch(ImplicitPatch::createReplacement(composableA, composableB));
        operation->name = "Replacement";
    }
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
    std::vector<std::tuple<ImplicitPatch*, int>> waitingPatches = { {this->mainPatch, 0} };
    while (!waitingPatches.empty()) {
        auto [current, level] = waitingPatches.back();
        waitingPatches.pop_back();
        primitiveSelectionGui->addItem(new HierarchicalListWidgetItem(current->toString(), current->index, level));
        this->storedPatches.push_back(current);
        if (current->composableA)
            waitingPatches.push_back({current->composableA, level+1});
        if (current->composableB)
            waitingPatches.push_back({current->composableB, level+1});
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









PatchReplacementDialog::PatchReplacementDialog(PrimitivePatchesInterface* caller, ImplicitPatch *patchToModify)
    : QDialog(), caller(caller), patch(patchToModify)
{
    QVBoxLayout * layout = new QVBoxLayout(this);

    QRadioButton* createSphereButton = new QRadioButton("Sphere");
    QRadioButton* createBlockButton = new QRadioButton("Block");
    QRadioButton* createGaussianButton = new QRadioButton("Gaussian");
//    QRadioButton* createIdentityButton = new QRadioButton("Identity");
    FancySlider* densitySlider = new FancySlider(Qt::Orientation::Horizontal, -5.f, 5.f, .1f);

    QRadioButton* stackingButton = new QRadioButton("Stack");
    QRadioButton* blendingButton = new QRadioButton("Blend");
    QRadioButton* replacingButton = new QRadioButton("Replace");
//    QRadioButton* noneButton = new QRadioButton("None");

    FancySlider* widthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* depthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* heightSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);
    FancySlider* sigmaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 25.f, .1f);

    cancelButton = new QPushButton("Annuler", this);
    validButton = new QPushButton("Confirmer", this);

//    QPushButton* resetButton = new QPushButton("Reset");

//    primitiveSelectionGui = new HierarchicalListWidget(this);

    bool patchIsOperation = patch->isOperation();
    if (patchIsOperation) {
        layout->addWidget(createVerticalGroup({
                                                        stackingButton,
                                                        blendingButton,
                                                        replacingButton
                                                    }));
    } else {
        layout->addWidget(createVerticalGroup({
                                                        createSphereButton,
                                                        createBlockButton,
                                                        createGaussianButton
                                                    }));
        layout->addWidget(createSliderGroup("Density", densitySlider));
        layout->addWidget(createMultipleSliderGroup({
                                                  {"Width / radius", widthSlider},
                                                        {"Depth", depthSlider},
                                                        {"Height", heightSlider},
                                                        {"Sigma", sigmaSlider}
                                              }));
    }

    layout->addWidget(cancelButton);
    layout->addWidget(validButton);

    createSphereButton->setChecked(caller->currentShapeSelected == PRIMITIVE_SHAPE::SPHERE);
    createBlockButton->setChecked(caller->currentShapeSelected == PRIMITIVE_SHAPE::BLOCK);
    createGaussianButton->setChecked(caller->currentShapeSelected == PRIMITIVE_SHAPE::GAUSSIAN);

    stackingButton->setChecked(caller->currentOperation == PATCH_OPERATION::STACK);
    blendingButton->setChecked(caller->currentOperation == PATCH_OPERATION::BLEND);
    replacingButton->setChecked(caller->currentOperation == PATCH_OPERATION::REPLACE);

    densitySlider->setfValue(patch->densityValue);
    widthSlider->setfValue(patch->dimension.x);
    heightSlider->setfValue(patch->dimension.z);
    depthSlider->setfValue(patch->dimension.y);
    sigmaSlider->setfValue(1.f);

    QObject::connect(cancelButton, &QPushButton::pressed, this, &PatchReplacementDialog::cancel);
    QObject::connect(validButton, &QPushButton::pressed, this, &PatchReplacementDialog::confirm);

    QObject::connect(createSphereButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(PRIMITIVE_SHAPE::SPHERE); });
    QObject::connect(createBlockButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(PRIMITIVE_SHAPE::BLOCK); });
    QObject::connect(createGaussianButton, &QRadioButton::toggled, this, [=]() { caller->setSelectedShape(PRIMITIVE_SHAPE::GAUSSIAN); });
    QObject::connect(densitySlider, &FancySlider::floatValueChanged, this, [=](float newDensity) { caller->selectedDensity = newDensity; });
    QObject::connect(stackingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(PATCH_OPERATION::STACK); });
    QObject::connect(blendingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(PATCH_OPERATION::BLEND); });
    QObject::connect(replacingButton, &QRadioButton::toggled, this, [=]() { caller->setCurrentOperation(PATCH_OPERATION::REPLACE); });

    QObject::connect(widthSlider, &FancySlider::floatValueChanged, this->caller, &PrimitivePatchesInterface::setSelectedWidth);
    QObject::connect(heightSlider, &FancySlider::floatValueChanged, this->caller, &PrimitivePatchesInterface::setSelectedHeight);
    QObject::connect(depthSlider, &FancySlider::floatValueChanged, this->caller, &PrimitivePatchesInterface::setSelectedDepth);
    QObject::connect(sigmaSlider, &FancySlider::floatValueChanged, this->caller, &PrimitivePatchesInterface::setSelectedSigma);

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
    this->close();
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
    this->close();
}
