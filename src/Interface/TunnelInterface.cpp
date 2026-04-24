#include "TunnelInterface.h"

#include "GUIElements/InterfaceUtils.h"
#include "Curves/BSpline.h"
#include "TerrainModification/UnderwaterErosion.h"

TunnelInterface::TunnelInterface(QWidget *parent)
    : ActionInterface("tunnel", "Tunnel generation", "digging", "Tunnels creation", "tunnel_button.png", parent),
      startingShape(KarstHolePredefinedShapes::TUBE), endingShape(KarstHolePredefinedShapes::TUBE)
{

}

/*void TunnelInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;

//    this->currentTunnelPoints = {Vector3(0, 0, 0), Vector3(10, 0, 0)};
//    createTunnel();
    //    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );
}*/

void TunnelInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    for (auto& ctrl : this->controlPoints) {
        ctrl->display();
    }
    if (controlPoints.size() > 1) {
        if (this->tunnelPreview.shader != nullptr) {
            this->tunnelPreview.shader->setVector("color", std::vector<float>({0.1f, 0.2f, 0.7f, 0.6f}));
        }
        this->tunnelPreview.reorderTriangles(camPos);
        this->tunnelPreview.display(); //GL_LINES, 5.f);
    }
}

void TunnelInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        bool removingMatter = parameters.at("removing_matter").get<bool>();
        KarstHolePredefinedShapes startingShape = parameters.at("starting_shape").get<KarstHolePredefinedShapes>();
        KarstHolePredefinedShapes endingShape = parameters.at("ending_shape").get<KarstHolePredefinedShapes>();
        float width = parameters.at("width").get<float>() * random_gen::generate(0.1f, 2.f);
        float height = parameters.at("height").get<float>() * random_gen::generate(0.1f, 2.f);
        float erosionStrength = parameters.at("erosion_strength").get<float>() * random_gen::generate(0.1f, 2.f);
        BSpline path = parameters.at("path");

        UnderwaterErosion erod(this->voxelGrid.get(), 0, erosionStrength, 0);
        KarstHole hole(path, width, height, startingShape, endingShape);
        erod.CreateTunnel(hole, !removingMatter);
    }
}

void TunnelInterface::hide()
{
    for (auto& ctrl : this->controlPoints)
        ctrl->hide();
    this->tunnelPreview.hide();
    CustomInteractiveObject::hide();
}

void TunnelInterface::mouseClickedOnMapEvent(const Vector3& mousePosInWorld, bool mouseInMap, QMouseEvent* event, TerrainModel* model)
{
    if (this->isVisible() && mouseInMap && event->button() == Qt::MouseButton::LeftButton)
        this->addCurvesControlPoint(model->getTerrainPos(mousePosInWorld));
}

void TunnelInterface::show()
{
    for (auto& ctrl : this->controlPoints)
        ctrl->show();
    this->tunnelPreview.show();
    CustomInteractiveObject::show();
}

InterfaceUI* TunnelInterface::createGUI()
{
    auto UI = new InterfaceUI();

    auto tunnelClearControlPointButton = new ButtonElement("Tout retirer", [&](){this->clearTunnelPoints(); /*computeTunnelPreview();*/ });
    auto tunnelWidthSlider = new SliderElement("Width", 1, 30, 1, tunnelWidth, [=](float) { computeTunnelPreview(); });
    auto tunnelHeightSlider = new SliderElement("Height", 1, 30, 1, tunnelHeight, [=](float) { computeTunnelPreview();}, UIElement::VERTICAL);
    auto tunnelStrengthSlider = new SliderElement("", 0.0f, 3.0f, 0.1f, erosionStrength, [=](float) { computeTunnelPreview(); });
    auto tunnelCreateMatter = new ButtonElement("Arche", [&]() { this->createTunnel(false); });
    auto tunnelRemoveMatter = new ButtonElement("Tunnel", [&]() { this->createTunnel(true); });
//    CheckboxElement* tunnelDisplayButton = std::make_shared<CheckboxElement>("Afficher");

    this->shapes = {
        new ComboboxLineElement<KarstHolePredefinedShapes>("Tube", ":/tunnels/src/assets/tunnels_icons/tunnel_type_tube.png", TUBE),
        new ComboboxLineElement<KarstHolePredefinedShapes>("Soluble bed", ":/tunnels/src/assets/tunnels_icons/tunnel_type_soluble_bed.png", SOLUBLE_BED),
        new ComboboxLineElement<KarstHolePredefinedShapes>("Keyhole", ":/tunnels/src/assets/tunnels_icons/tunnel_type_keyhole.png", KEYHOLE),
        new ComboboxLineElement<KarstHolePredefinedShapes>("Canyon", ":/tunnels/src/assets/tunnels_icons/tunnel_type_canyon.png", CANYON),
        new ComboboxLineElement<KarstHolePredefinedShapes>("Crack", ":/tunnels/src/assets/tunnels_icons/tunnel_type_fracture.png", CRACK),
        new ComboboxLineElement<KarstHolePredefinedShapes>("Flat", ":/tunnels/src/assets/tunnels_icons/tunnel_type_fracture_flat.png", STAR)
    };

    auto startingShapeCombobox = new ComboboxElement("Inlet", shapes, startingShapeIndex);
    auto endingShapeCombobox = new ComboboxElement("Outlet", shapes, endingShapeIndex);


    UI->add({
        tunnelCreateMatter,
        tunnelRemoveMatter,
        tunnelClearControlPointButton,
        tunnelWidthSlider,
        tunnelHeightSlider,
        startingShapeCombobox,
        endingShapeCombobox
    });

    startingShapeCombobox->setOnSelectionChanged([=](int) { this->startingShape = startingShapeCombobox->getSelection<KarstHolePredefinedShapes>(); this->computeTunnelPreview(); });
    endingShapeCombobox->setOnSelectionChanged([=](int) { this->endingShape = endingShapeCombobox->getSelection<KarstHolePredefinedShapes>(); this->computeTunnelPreview(); });

    return UI;
}

void TunnelInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";

    this->tunnelPreview = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->tunnelPreview.useIndices = false;
}


void TunnelInterface::addCurvesControlPoint(const Vector3& pos, bool justUpdatePath)
{
    if (!justUpdatePath)
    {
        bool addTheNewPoint = true;
        for (auto& controls : this->controlPoints) {
            if (controls->isManipulated()) {
                addTheNewPoint = false;
                break;
            }
        }
        if (addTheNewPoint) {
            this->controlPoints.push_back(std::make_shared<ControlPoint>(pos, 2.5f));
            std::shared_ptr<ControlPoint>& newCtrl = this->controlPoints.back();
            newCtrl->allowAllAxisTranslation(true);
            newCtrl->setDisplayOnTop(true);
            // QObject::connect(newCtrl.get(), &ControlPoint3D::pointModified,
                             // this, [=](){
            newCtrl->setOnPointModified([=]() {
                this->addCurvesControlPoint(Vector3(), true);
            });
        }
    }
    this->currentTunnelPoints.clear();
//    bool atLeastOnePointIsManipulated = false;
    for (auto& control : this->controlPoints) {
        this->currentTunnelPoints.push_back(control->getPosition());
        if (control->isManipulated() && this->voxelGrid->contains(control->getPosition())) {
//            atLeastOnePointIsManipulated = true;
            Q_EMIT this->needToClipView(
                        control->getFluidTranslation(),
                        control->getPosition(),
                        true
                        );
            // QObject::connect(control.get(), &ControlPoint3D::pointReleased,
                             // this, [&]() -> void { Q_EMIT this->needToClipView(Vector3(), Vector3(), false); });
            control->setOnPointReleased([=]() { Q_EMIT this->needToClipView(Vector3(), Vector3(), false); });
        }
    }
    this->computeTunnelPreview();

    Q_EMIT updated();
}

void TunnelInterface::updateStartingShape()
{
    startingShape = shapes[startingShapeIndex]->value; //(this->startingShapeCombobox->currentData().toInt());
    this->computeTunnelPreview();
}

void TunnelInterface::updateEndingShape()
{
    endingShape = shapes[endingShapeIndex]->value; //(this->endingShapeCombobox->currentData().toInt());
    this->computeTunnelPreview();
}

void TunnelInterface::clearTunnelPoints()
{
    this->currentTunnelPoints.clear();
    this->controlPoints.clear();
    this->tunnelPreview.clear();

    Q_EMIT updated();
}
void TunnelInterface::createTunnel(bool removingMatter)
{
    if (this->currentTunnelPoints.empty()) return;

    UnderwaterErosion erod(this->voxelGrid.get(), 0, erosionStrength, 0);
    BSpline path(this->currentTunnelPoints);
    KarstHole hole(path, this->tunnelWidth, this->tunnelHeight, startingShape, endingShape);
    this->tunnelPreview.fromArray(erod.CreateTunnel(hole, !removingMatter, true));
    this->currentTunnelPoints.clear();
    this->controlPoints.clear();
    this->tunnelPreview.update();

    this->addTerrainAction(nlohmann::json({
                                              {"removing_matter", removingMatter},
                                              {"starting_shape", startingShape},
                                              {"ending_shape", endingShape},
                                              {"height", tunnelHeight},
                                              {"width", tunnelWidth},
                                              {"erosion_strength", erosionStrength},
                                              {"path", path}
                                          }));
    Q_EMIT tunnelCreated(hole);
    Q_EMIT updated();
}
void TunnelInterface::createCrack(bool removingMatter)
{/*
    if (this->currentTunnelPoints.size() < 2) return;

    UnderwaterErosion erod(this->voxelGrid, 0, erosionStrength, 0);
    BSpline path(this->currentTunnelPoints);
    KarstHole hole(path, this->tunnelWidth, this->tunnelHeight, CRACK, CRACK);
    this->tunnelPreview.fromArray(erod.CreateTunnel(hole, !removingMatter, true));
    this->currentTunnelPoints.clear();
    this->controlPoints.clear();
    this->tunnelPreview.update();

    Q_EMIT updated();*/
}

void TunnelInterface::computeTunnelPreview() {
    if (this->currentTunnelPoints.size() > 1) {
        BSpline path(this->currentTunnelPoints);
        KarstHole previewHole(path, this->tunnelWidth, this->tunnelHeight, startingShape, endingShape);
        std::vector<std::vector<Vector3>> vertices = previewHole.generateMesh();
        std::vector<Vector3> meshVertices;
        for (const auto& triangle : vertices) {
            meshVertices.push_back(triangle[0]);
            meshVertices.push_back(triangle[1]);
            meshVertices.push_back(triangle[2]);
        }
        this->tunnelPreview.fromArray(meshVertices);
    } else {
        this->tunnelPreview.clear();
    }
    this->tunnelPreview.update();
    Q_EMIT this->updated();
}

void TunnelInterface::wheelEvent(QWheelEvent* event)
{
    if (this->isHidden()) return;
    if (event->modifiers().testFlag(Qt::ControlModifier)) {
        this->setTunnelWidth(this->tunnelWidth - event->angleDelta().y() / 2);
        Q_EMIT this->updated();
    } else if (event->modifiers().testFlag(Qt::ShiftModifier)) {
        this->setTunnelHeight(this->tunnelHeight - event->angleDelta().y() / 2);
        Q_EMIT this->updated();
    }
    CustomInteractiveObject::wheelEvent(event);
}


void TunnelInterface::setTunnelWidth(int newSize)
{
    if (newSize < 0) newSize = 0;
    this->tunnelWidth = newSize;
//    this->tunnelWidthSlider->setValue(newSize);
    computeTunnelPreview();
}
void TunnelInterface::setTunnelHeight(int newSize)
{
    if (newSize < 0) newSize = 0;
    this->tunnelHeight = newSize;
//    this->tunnelHeightSlider->setValue(newSize);
    computeTunnelPreview();
}
