#include "FaultSlipInterface.h"

#include "Utils/Utils.h"


FaultSlipInterface::FaultSlipInterface(QWidget *parent) : ActionInterface("faultSlip", parent)
{
    this->faultSlip = FaultSlip();
    this->firstSlipControlPoint = std::make_unique<ControlPoint>(Vector3(), 5.f);
    this->slipVector = std::make_unique<InteractiveVector>(Vector3(), Vector3());
    this->createGUI();
}

void FaultSlipInterface::display()
{
    if (voxelGrid != nullptr) {
        this->firstSlipControlPoint->display();
        this->slipVector->display();
        if (this->planeMesh.shader != nullptr)
            this->planeMesh.shader->setVector("color", std::vector<float>({.1f, .5f, 1.f, .5f}));
        this->planeMesh.display();
    }
}

void FaultSlipInterface::remesh()
{
    Vector3 minAABB = Vector3(0, 0, 0);
    Vector3 maxAABB = Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ);

    Vector3 p1 = this->firstSlipControlPoint->getPosition();
    Vector3 p2 = this->slipVector->getStartingVector();
    Vector3 p3 = p1 + this->slipVector->getResultingVector();

    Vector3 largeur  = (p2 - p1).normalize() * maxAABB;
    Vector3 hauteur = (p3 - p1).normalize() * maxAABB;

    this->planeMesh.fromArray({
                                  p1 - largeur - hauteur, p2 + largeur - hauteur, p1 - largeur + hauteur,
                                  p1 - largeur + hauteur, p2 + largeur - hauteur, p2 + largeur + hauteur,

                                  p2 + largeur - hauteur, p1 - largeur - hauteur, p1 - largeur + hauteur,
                                  p2 + largeur - hauteur, p1 - largeur + hauteur, p2 + largeur + hauteur
                              });
    this->planeMesh.update();
}

void FaultSlipInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->firstSlipControlPoint->move(Vector3(voxelGrid->sizeX / 2.f, 0, voxelGrid->sizeZ));
    this->slipVector->setPositions(Vector3(0, voxelGrid->sizeY / 2.f, voxelGrid->sizeZ), Vector3(0, voxelGrid->sizeY / 2.f, voxelGrid->sizeZ / 2.f));
    this->voxelGrid = voxelGrid;
    this->remesh();

    this->setBindings();
}

void FaultSlipInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        Vector3 slippingDirection = json_to_vec3(parameters.at("slipping_direction"));
        float slippingDistance = parameters.at("slipping_distance").get<float>();
        Vector3 firstPointPos = json_to_vec3(parameters.at("first_point_pos"));
        Vector3 secondPointPos = json_to_vec3(parameters.at("second_point_pos"));
        bool positiveSideFalling = parameters.at("positive_side").get<bool>();

        this->faultSlip.slippingDirection = slippingDirection;
        this->faultSlip.slippingDistance = slippingDistance;
        this->faultSlip.firstPointFault = firstPointPos;
        this->faultSlip.secondPointFault = secondPointPos;
        this->faultSlip.positiveSideFalling = positiveSideFalling;
        if (this->voxelGrid != nullptr)
        {
            this->faultSlip.Apply(this->voxelGrid, false);
        }
    }
}
void FaultSlipInterface::updateSlipVector(Vector3 newSlipVector)
{
    this->remesh();
}

void FaultSlipInterface::updatePoints()
{
    this->remesh();
}

void FaultSlipInterface::computeFaultSlip()
{
    Vector3 slippingDirection = this->slipVector->getResultingVector();
    float slippingDistance = 1.f;
    Vector3 firstPointPos = this->firstSlipControlPoint->getPosition();
    Vector3 secondPointPos = this->slipVector->getStartingVector();
    bool positiveSideFalling = this->faultSideApplied->isChecked();
    this->faultSlip.slippingDirection = slippingDirection;
    this->faultSlip.slippingDistance = slippingDistance;
    this->faultSlip.firstPointFault = firstPointPos;
    this->faultSlip.secondPointFault = secondPointPos;
    this->faultSlip.positiveSideFalling = positiveSideFalling;
    if (this->voxelGrid != nullptr)
    {
        this->faultSlip.Apply(this->voxelGrid, true);
        this->addTerrainAction(nlohmann::json({
                                               {"slipping_direction", vec3_to_json(slippingDirection)},
                                               {"slipping_distance", slippingDistance},
                                               {"first_point_pos", vec3_to_json(firstPointPos)},
                                               {"second_point_pos", vec3_to_json(secondPointPos)},
                                               {"positive_side", positiveSideFalling}
                                              }));

        Q_EMIT this->faultSlipApplied();
    }
}

void FaultSlipInterface::setSideAffected(bool isRightSide)
{
    this->faultSlip.positiveSideFalling = isRightSide;
}

void FaultSlipInterface::hide()
{
    this->firstSlipControlPoint->hide();
    this->slipVector->hide();
    this->planeMesh.hide();
    CustomInteractiveObject::hide();
}

void FaultSlipInterface::show()
{
    this->firstSlipControlPoint->show();
    this->slipVector->show();
    this->planeMesh.show();
    CustomInteractiveObject::show();
}


QLayout *FaultSlipInterface::createGUI()
{
    if (this->faultSlipLayout != nullptr) return faultSlipLayout;

    this->faultSlipLayout = new QHBoxLayout;

    this->faultApplyButton = new QPushButton("Chuter");
    this->faultSideApplied = new QCheckBox("Partie de droite chute");
    this->faultDisplayButton = new QCheckBox("Afficher");
    faultSlipLayout->addWidget(faultApplyButton);
    faultSlipLayout->addWidget(faultSideApplied);


    faultSideApplied->setChecked(this->faultSlip.positiveSideFalling);

    return this->faultSlipLayout;
}


void FaultSlipInterface::setBindings()
{
//    if (this->voxelGrid != nullptr)
//    {
        QObject::connect(faultApplyButton, &QPushButton::pressed, this, [=](){ this->computeFaultSlip(); } );
        QObject::connect(faultSideApplied, &QCheckBox::toggled, this, &FaultSlipInterface::setSideAffected);
        QObject::connect(faultDisplayButton, &QCheckBox::toggled, this, &FaultSlipInterface::setVisibility);

        QObject::connect(slipVector.get(), &InteractiveVector::modified, this, &FaultSlipInterface::updateSlipVector);
        QObject::connect(firstSlipControlPoint.get(), &ControlPoint::modified, this, &FaultSlipInterface::updatePoints);
//    }
}
