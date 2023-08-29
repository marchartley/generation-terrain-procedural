#include "SmoothInterface.h"

SmoothInterface::SmoothInterface(QWidget *parent)
    : ActionInterface("smooth", "Smooth voxels", "", "", "", parent) // No menu
{
}

void SmoothInterface::display(const Vector3& camPos)
{

}

void SmoothInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        this->applySmooth();
    }
}


bool SmoothInterface::applySmooth()
{
    this->voxelGrid->smoothVoxels();
    this->heightmap->fromVoxelGrid(*voxelGrid);

    this->addTerrainAction(nlohmann::json({}));
    Q_EMIT updated();
    return false;
}
void SmoothInterface::hide()
{
    CustomInteractiveObject::hide();
}

void SmoothInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* SmoothInterface::createGUI()
{
    this->smoothLayout = new QHBoxLayout;

    smoothComputeButton = new QPushButton("Calculer");
//    smoothDisplayButton = new QCheckBox("Afficher");
    smoothLayout->addWidget(smoothComputeButton);
//    smoothLayout->addWidget(smoothDisplayButton);

//    smoothDisplayButton->setChecked(this->visible);

//    QObject::connect(smoothDisplayButton, &QCheckBox::toggled, this, &SmoothInterface::setVisibility);
    QObject::connect(smoothComputeButton, &QPushButton::pressed, this, &SmoothInterface::applySmooth);

    return this->smoothLayout;
}

