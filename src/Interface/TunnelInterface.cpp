#include "TunnelInterface.h"

#include "Interface/InterfaceUtils.h"
#include "Utils/BSpline.h"
#include "TerrainModification/UnderwaterErosion.h"

TunnelInterface::TunnelInterface(QWidget *parent) : ActionInterface("tunnel", parent), startingShape(KarstHolePredefinedShapes::TUBE), endingShape(KarstHolePredefinedShapes::TUBE)
{

}

void TunnelInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;

//    this->currentTunnelPoints = {Vector3(0, 0, 0), Vector3(10, 0, 0)};
//    createTunnel();
    //    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );
}

void TunnelInterface::display()
{
    for (auto& ctrl : this->controlPoints)
        ctrl->display();
    if (controlPoints.size() > 1) {
        if (this->tunnelPreview.shader != nullptr) {
            this->tunnelPreview.shader->setVector("color", std::vector<float>({0.1f, 0.2f, 0.7f, 0.6f}));
        }
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
        float width = parameters.at("width").get<float>();
        float height = parameters.at("height").get<float>();
        float erosionStrength = parameters.at("erosion_strength").get<float>();
        BSpline path = json_to_bspline(parameters.at("path"));

        UnderwaterErosion erod(this->voxelGrid, 0, erosionStrength, 0);
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

void TunnelInterface::mouseClickInWorldEvent(Vector3 mousePosInWorld, bool mouseInMap, QMouseEvent* event)
{
    if (this->isVisible() && mouseInMap && event->button() == Qt::MouseButton::LeftButton)
        this->addCurvesControlPoint(mousePosInWorld);
}

void TunnelInterface::show()
{
    for (auto& ctrl : this->controlPoints)
        ctrl->show();
    this->tunnelPreview.show();
    CustomInteractiveObject::show();
}

QLayout* TunnelInterface::createGUI()
{
    this->tunnelLayout = new QHBoxLayout;

    addControlPointButton = new QPushButton("Ajouter un point de control");
    tunnelClearControlPointButton = new QPushButton("Tout retirer");
    tunnelWidthSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 30);
    tunnelHeightSlider = new FancySlider(Qt::Orientation::Vertical, 1, 30);
    tunnelStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    tunnelCreateMatter = new QPushButton("Arche");
    tunnelRemoveMatter = new QPushButton("Tunnel");
//    tunnelCreateCrack = new QPushButton("Faille");
    tunnelDisplayButton = new QCheckBox("Afficher");
    startingShapeCombobox = new QComboBox();
    endingShapeCombobox = new QComboBox();

    QIcon tubeIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_tube.png");
    QIcon solubleIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_soluble_bed.png");
    QIcon keyholeIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_keyhole.png");
    QIcon canyonIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_canyon.png");
    QIcon fractureIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_fracture.png");
    QIcon flatCrackIcon = QIcon(":/tunnels/src/assets/tunnels_icons/tunnel_type_fracture_flat.png");

    startingShapeCombobox->addItem(tubeIcon, "Tube", KarstHolePredefinedShapes::TUBE);
    startingShapeCombobox->addItem(solubleIcon, "Soluble bed", KarstHolePredefinedShapes::SOLUBLE_BED);
    startingShapeCombobox->addItem(keyholeIcon, "Keyhole", KarstHolePredefinedShapes::KEYHOLE);
    startingShapeCombobox->addItem(canyonIcon, "Canyon", KarstHolePredefinedShapes::CANYON);
    startingShapeCombobox->addItem(fractureIcon, "Fracture", KarstHolePredefinedShapes::CRACK);
    startingShapeCombobox->addItem(flatCrackIcon, "Fond plat", KarstHolePredefinedShapes::STAR);

    endingShapeCombobox->addItem(tubeIcon, "Tube", KarstHolePredefinedShapes::TUBE);
    endingShapeCombobox->addItem(solubleIcon, "Soluble bed", KarstHolePredefinedShapes::SOLUBLE_BED);
    endingShapeCombobox->addItem(keyholeIcon, "Keyhole", KarstHolePredefinedShapes::KEYHOLE);
    endingShapeCombobox->addItem(canyonIcon, "Canyon", KarstHolePredefinedShapes::CANYON);
    endingShapeCombobox->addItem(fractureIcon, "Fracture", KarstHolePredefinedShapes::CRACK);
    endingShapeCombobox->addItem(flatCrackIcon, "Fond plat", KarstHolePredefinedShapes::STAR);

    tunnelLayout->addWidget(createVerticalGroup({tunnelCreateMatter, tunnelRemoveMatter/*, tunnelCreateCrack*/}));
    tunnelLayout->addWidget(createVerticalGroup({/*addControlPointButton, */tunnelClearControlPointButton}));
    tunnelLayout->addWidget(createSliderGroup("Largeur", tunnelWidthSlider));
    tunnelLayout->addWidget(createSliderGroup("Hauteur", tunnelHeightSlider));
    tunnelLayout->addWidget(createVerticalGroup({new QLabel("Entrée"), startingShapeCombobox}));
    tunnelLayout->addWidget(createVerticalGroup({new QLabel("Sortie"), endingShapeCombobox}));
//    tunnelLayout->addWidget(createSliderGroup("Force", tunnelStrengthSlider));
//    tunnelLayout->addWidget(tunnelDisplayButton);


    QObject::connect(tunnelWidthSlider, &FancySlider::valueChanged, this, &TunnelInterface::setTunnelWidth);
    QObject::connect(tunnelHeightSlider, &FancySlider::valueChanged, this, &TunnelInterface::setTunnelHeight);
    QObject::connect(tunnelStrengthSlider, &FancySlider::floatValueChanged, this, &TunnelInterface::setErosionStrength);
    QObject::connect(tunnelCreateMatter, &QPushButton::pressed, this, [=](){ this->createTunnel(false); } );
    QObject::connect(tunnelRemoveMatter, &QPushButton::pressed, this, [=](){ this->createTunnel(true);  } );
//    QObject::connect(tunnelCreateCrack, &QPushButton::pressed, this, [=](){ this->createCrack(true); } );
//    QObject::connect(addControlPointButton, &QPushButton::pressed, this, [=](){this->setCurvesErosionConstructionMode(true); });
    QObject::connect(tunnelClearControlPointButton, &QPushButton::pressed, this, [=](){this->clearTunnelPoints(); computeTunnelPreview(); });
//    QObject::connect(tunnelDisplayButton, &QCheckBox::toggled, this, &TunnelInterface::setVisibility );
    QObject::connect(startingShapeCombobox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [&](int x) {this->updateStartingShape(); });
    QObject::connect(endingShapeCombobox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [&](int x) {this->updateEndingShape(); });


    this->tunnelWidthSlider->setValue(tunnelWidth);
    this->tunnelHeightSlider->setValue(tunnelHeight);
    this->tunnelStrengthSlider->setfValue(erosionStrength);

    this->updateStartingShape();
    this->updateEndingShape();

    return this->tunnelLayout;
}


void TunnelInterface::addCurvesControlPoint(Vector3 pos, bool justUpdatePath)
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
            this->controlPoints.push_back(std::make_unique<ControlPoint>(pos, 5.f, INACTIVE));
            std::unique_ptr<ControlPoint>& newCtrl = this->controlPoints.back();
            newCtrl->allowAllAxisTranslation(true);
            newCtrl->displayOnTop = true;
            QObject::connect(newCtrl.get(), &ControlPoint::modified,
                             this, [&](){ this->addCurvesControlPoint(Vector3(), true); });
        }
    }
    this->currentTunnelPoints.clear();
//    bool atLeastOnePointIsManipulated = false;
    for (auto& control : this->controlPoints) {
        this->currentTunnelPoints.push_back(control->getPosition());
        if (control->isManipulated() && this->voxelGrid->contains(control->getPosition())) {
//            atLeastOnePointIsManipulated = true;
            Q_EMIT this->needToClipView(
                        control->getFluidMovement(),
                        control->getPosition(),
                        true
                        );
            QObject::connect(control.get(), &ControlPoint::afterModified,
                             this, [&]() -> void { Q_EMIT this->needToClipView(Vector3(), Vector3(), false); });
        }
    }
    this->computeTunnelPreview();

    Q_EMIT updated();
}

void TunnelInterface::updateStartingShape()
{
    startingShape = static_cast<KarstHolePredefinedShapes>(this->startingShapeCombobox->currentData().toInt());
    this->computeTunnelPreview();
}

void TunnelInterface::updateEndingShape()
{
    endingShape = static_cast<KarstHolePredefinedShapes>(this->endingShapeCombobox->currentData().toInt());
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

    UnderwaterErosion erod(this->voxelGrid, 0, erosionStrength, 0);
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
                                              {"path", bspline_to_json(path)}
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
    this->tunnelWidthSlider->setValue(newSize);
    computeTunnelPreview();
}
void TunnelInterface::setTunnelHeight(int newSize)
{
    if (newSize < 0) newSize = 0;
    this->tunnelHeight = newSize;
    this->tunnelHeightSlider->setValue(newSize);
    computeTunnelPreview();
}
