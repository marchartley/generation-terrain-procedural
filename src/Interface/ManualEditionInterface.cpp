#include "ManualEditionInterface.h"

#include "Interface/InterfaceUtils.h"

ManualEditionInterface::ManualEditionInterface()
{

}

void ManualEditionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}

void ManualEditionInterface::show()
{
    CustomInteractiveObject::show();
}

void ManualEditionInterface::hide()
{
    CustomInteractiveObject::hide();
}


QLayout *ManualEditionInterface::createGUI()
{
    this->manualEditLayout = new QHBoxLayout();
    // Manual rock erosion layout
    this->manualEditSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->manualEditStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->addingModeButton = new QRadioButton("Ajouter de la matière");
    this->suppressModeButton = new QRadioButton("Detruire de la matière");
    manualEditLayout->addWidget(createSliderGroup("Taille", manualEditSizeSlider));
    manualEditLayout->addWidget(createSliderGroup("Force", manualEditStrengthSlider));
    manualEditLayout->addWidget(createVerticalGroup({addingModeButton, suppressModeButton}));

    return this->manualEditLayout;
}
