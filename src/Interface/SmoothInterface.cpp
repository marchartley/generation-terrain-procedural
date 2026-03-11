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

InterfaceUI* SmoothInterface::createGUI()
{
    auto UI = new InterfaceUI();

    auto smoothComputeButton = new ButtonElement("Calculer", [=]() { this->applySmooth(); });
    UI->add(smoothComputeButton);

    return UI;
}

