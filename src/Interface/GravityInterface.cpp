#include "GravityInterface.h"

GravityInterface::GravityInterface(QWidget *parent) : CustomInteractiveObject(parent)
{
    this->createGUI();
}

void GravityInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}

void GravityInterface::display()
{

}


bool GravityInterface::createGlobalGravity()
{
    this->voxelGrid->makeItFall();
    Q_EMIT updated();
    return false;
    /*
    this->startAnimation();
    this->applyLetItFall = !this->applyLetItFall;
//    if (this->applyLetItFall)
        // this->displayMessage( "Gravity is making his job!" );
//    else
        // this->displayMessage( "Gravity stopped caring" );
    update();
    return this->applyLetItFall;*/
}

bool GravityInterface::createSandGravity()
{
    this->voxelGrid->letGravityMakeSandFall(true);
    Q_EMIT updated();
    return false;
    /*
    this->startAnimation();
    this->applyLetSandFall = !this->applyLetSandFall;
//    if (this->applyLetSandFall)
        // this->displayMessage( "Sand is falling!" );
//    else
        // this->displayMessage( "Sand stopped falling" );
    update();
    return this->applyLetSandFall;*/
}

void GravityInterface::hide()
{
    CustomInteractiveObject::hide();
}

void GravityInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* GravityInterface::createGUI()
{
    this->gravityLayout = new QHBoxLayout;

    gravityComputeButton = new QPushButton("Calculer");
//    gravityDisplayButton = new QCheckBox("Afficher");
    gravityLayout->addWidget(gravityComputeButton);
//    gravityLayout->addWidget(gravityDisplayButton);

//    gravityDisplayButton->setChecked(this->visible);

//    QObject::connect(gravityDisplayButton, &QCheckBox::toggled, this, &GravityInterface::setVisibility);
    QObject::connect(gravityComputeButton, &QPushButton::pressed, this, &GravityInterface::createSandGravity);

    return this->gravityLayout;
}
