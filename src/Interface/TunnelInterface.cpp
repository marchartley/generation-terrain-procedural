#include "TunnelInterface.h"

#include "Interface/InterfaceUtils.h"
#include "Utils/BSpline.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Karst/KarstHole.h"

TunnelInterface::TunnelInterface(QWidget *parent) : CustomInteractiveObject(parent)
{

}

void TunnelInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    //    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );
}

void TunnelInterface::display()
{
    for (auto& ctrl : this->controlPoints)
        ctrl->display();
    this->tunnelPreview.display(); //GL_LINES, 5.f);
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
    tunnelCreateCrack = new QPushButton("Faille");
    tunnelDisplayButton = new QCheckBox("Afficher");
    tunnelLayout->addWidget(createVerticalGroup({tunnelCreateMatter, tunnelRemoveMatter, tunnelCreateCrack}));
    tunnelLayout->addWidget(createVerticalGroup({/*addControlPointButton, */tunnelClearControlPointButton}));
    tunnelLayout->addWidget(createSliderGroup("Largeur", tunnelWidthSlider));
    tunnelLayout->addWidget(createSliderGroup("Hauteur", tunnelHeightSlider));
//    tunnelLayout->addWidget(createSliderGroup("Force", tunnelStrengthSlider));
//    tunnelLayout->addWidget(tunnelDisplayButton);

    this->tunnelWidthSlider->setValue(tunnelWidth);
    this->tunnelHeightSlider->setValue(tunnelHeight);
    this->tunnelStrengthSlider->setfValue(erosionStrength);

    QObject::connect(tunnelWidthSlider, &FancySlider::valueChanged, this, &TunnelInterface::setTunnelWidth);
    QObject::connect(tunnelHeightSlider, &FancySlider::valueChanged, this, &TunnelInterface::setTunnelHeight);
    QObject::connect(tunnelStrengthSlider, &FancySlider::floatValueChanged, this, &TunnelInterface::setErosionStrength);
    QObject::connect(tunnelCreateMatter, &QPushButton::pressed, this, [=](){ this->createTunnel(false); } );
    QObject::connect(tunnelRemoveMatter, &QPushButton::pressed, this, [=](){ this->createTunnel(true);  } );
    QObject::connect(tunnelCreateCrack, &QPushButton::pressed, this, [=](){ this->createCrack(true); } );
//    QObject::connect(addControlPointButton, &QPushButton::pressed, this, [=](){this->setCurvesErosionConstructionMode(true); });
    QObject::connect(tunnelClearControlPointButton, &QPushButton::pressed, this, [=](){this->clearTunnelPoints(); computeTunnelPreview(); });
//    QObject::connect(tunnelDisplayButton, &QCheckBox::toggled, this, &TunnelInterface::setVisibility );

    return this->tunnelLayout;
}


void TunnelInterface::addCurvesControlPoint(Vector3 pos, bool justUpdatePath)
{
    if (!justUpdatePath)
    {
        bool addTheNewPoint = true;
        for (auto& controls : this->controlPoints) {
            if (controls->manipFrame.isManipulated()) {
                addTheNewPoint = false;
                break;
            }
        }
        if (addTheNewPoint) {
            this->controlPoints.push_back(std::make_unique<ControlPoint>(pos, 5.f, INACTIVE));
            QObject::connect(this->controlPoints.back().get(), &ControlPoint::modified,
                             this, [&](){ this->addCurvesControlPoint(Vector3(), true); });
        }
    }
    this->currentTunnelPoints.clear();
    for (auto& controls : this->controlPoints) {
        this->currentTunnelPoints.push_back(controls->getPosition());
//        controls->onUpdate([=]{ this->addCurvesControlPoint(Vector3(), true); });
    }
    this->computeTunnelPreview();

    Q_EMIT updated();
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
    KarstHole hole(path, this->tunnelWidth, this->tunnelHeight);
    this->tunnelPreview.fromArray(erod.CreateTunnel(hole, !removingMatter, true));
    this->currentTunnelPoints.clear();
    this->controlPoints.clear();
    this->tunnelPreview.update();

    Q_EMIT updated();
}
void TunnelInterface::createCrack(bool removingMatter)
{
    if (this->currentTunnelPoints.size() < 2) return;

    UnderwaterErosion erod(this->voxelGrid, 0, erosionStrength, 0);
    BSpline path(this->currentTunnelPoints);
    KarstHole hole(path, this->tunnelWidth, this->tunnelHeight, CRACK, CRACK);
    this->tunnelPreview.fromArray(erod.CreateTunnel(hole, removingMatter, true));
    this->currentTunnelPoints.clear();
    this->controlPoints.clear();
    this->tunnelPreview.update();

    Q_EMIT updated();
}

void TunnelInterface::computeTunnelPreview() {
    if (this->currentTunnelPoints.size() > 1) {
        BSpline path(this->currentTunnelPoints);
        KarstHole previewHole(path, this->tunnelWidth, this->tunnelHeight);
        std::vector<std::vector<Vector3>> vertices = previewHole.generateMesh();
        std::vector<Vector3> meshVertices;
        for (const auto& triangle : vertices) {
            meshVertices.push_back(triangle[0]);
            meshVertices.push_back(triangle[1]);
            meshVertices.push_back(triangle[2]);
        }
        this->tunnelPreview.fromArray(meshVertices);
        if (this->tunnelPreview.shader != nullptr) {
            this->tunnelPreview.shader->setVector("color", std::vector<float>({0.1f, 0.2f, 0.7f, 0.6f}));
        }
    } else {
        this->tunnelPreview.clear();
    }
    this->tunnelPreview.update();
    Q_EMIT this->updated();
}

void TunnelInterface::wheelEvent(QWheelEvent* event)
{
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
