#include "FaultSlipInterface.h"

#include "Utils/Utils.h"


FaultSlipInterface::FaultSlipInterface(QWidget *parent) : CustomInteractiveObject(parent)
{

    this->faultSlip = FaultSlip();
    this->firstSlipControlPoint = std::make_unique<ControlPoint>(Vector3(), 5.f);
    this->slipVector = std::make_unique<InteractiveVector>(Vector3(), Vector3());
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
    /*
    Vector3 x1 = intersectionPointRayPlane(p1, p1 - p2, minAABB, Vector3(1, 0, 0));
    Vector3 x2 = intersectionPointRayPlane(p1, p1 - p2, maxAABB, Vector3(1, 0, 0));
    Vector3 y1 = intersectionPointRayPlane(p1, p1 - p2, minAABB, Vector3(0, 1, 0));
    Vector3 y2 = intersectionPointRayPlane(p1, p1 - p2, maxAABB, Vector3(0, 1, 0));
    Vector3 z1 = intersectionPointRayPlane(p1, p1 - p2, minAABB, Vector3(0, 0, 1));
    Vector3 z2 = intersectionPointRayPlane(p1, p1 - p2, maxAABB, Vector3(0, 0, 1));*/

    Vector3 largeur  = (p2 - p1).normalize() * maxAABB;
    Vector3 hauteur = (p3 - p1).normalize() * maxAABB;

    this->planeMesh.fromArray({
                                  p1 - largeur - hauteur, p2 + largeur - hauteur, p1 - largeur + hauteur,
                                  p1 - largeur + hauteur, p2 + largeur - hauteur, p2 + largeur + hauteur
                              });
    this->planeMesh.update();
}

void FaultSlipInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->firstSlipControlPoint->move(Vector3(voxelGrid->sizeX / 2.f, 0, voxelGrid->sizeZ));
    this->slipVector->setPositions(Vector3(0, voxelGrid->sizeY / 2.f, voxelGrid->sizeZ), Vector3(0, voxelGrid->sizeY / 2.f, 0));
    this->voxelGrid = voxelGrid;
    this->remesh();

    this->setBindings();
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
    this->faultSlip.slippingDirection = this->slipVector->getResultingVector();
    this->faultSlip.slippingDistance = 1.f;
    this->faultSlip.firstPointFault = this->firstSlipControlPoint->getPosition();
    this->faultSlip.secondPointFault = this->slipVector->getStartingVector();
    this->faultSlip.positiveSideFalling = this->faultSideApplied->isChecked();
    if (this->voxelGrid != nullptr)
    {
        this->faultSlip.Apply(this->voxelGrid, true);
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
    this->faultSlipLayout = new QHBoxLayout;

    this->faultApplyButton = new QPushButton("Chuter");
    this->faultSideApplied = new QCheckBox("Partie de droite chute");
    this->faultDisplayButton = new QCheckBox("Afficher");
    faultSlipLayout->addWidget(faultApplyButton);
    faultSlipLayout->addWidget(faultSideApplied);
    // With new interface, we display/hide directly with buttons
//    faultSlipLayout->addWidget(faultDisplayButton);


    faultSideApplied->setChecked(this->faultSlip.positiveSideFalling);

    this->setBindings();

    return this->faultSlipLayout;
}


void FaultSlipInterface::setBindings()
{
    if (this->voxelGrid != nullptr)
    {
        QObject::connect(faultApplyButton, &QPushButton::pressed, this, [=](){ this->computeFaultSlip(); } );
        QObject::connect(faultSideApplied, &QCheckBox::toggled, this, &FaultSlipInterface::setSideAffected);
        QObject::connect(faultDisplayButton, &QCheckBox::toggled, this, &FaultSlipInterface::setVisibility);

        QObject::connect(slipVector.get(), &InteractiveVector::modified, this, &FaultSlipInterface::updateSlipVector);
        QObject::connect(firstSlipControlPoint.get(), &ControlPoint::modified, this, &FaultSlipInterface::updatePoints);
    }
}
