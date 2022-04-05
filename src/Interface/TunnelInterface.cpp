#include "TunnelInterface.h"

#include "Interface/InterfaceUtils.h"

TunnelInterface::TunnelInterface()
{

}

void TunnelInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
//    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );
}

void TunnelInterface::hide()
{
    CustomInteractiveObject::hide();
}

void TunnelInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* TunnelInterface::createGUI()
{
    this->tunnelLayout = new QHBoxLayout;

    addControlPointButton = new QPushButton("Ajouter un point de control");
    tunnelClearControlPointButton = new QPushButton("Tout retirer");
    tunnelSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 30);
    tunnelStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    tunnelCreateMatter = new QPushButton("Creer un pont");
    tunnelRemoveMatter = new QPushButton("Creer un tunnel");
    tunnelCreateCrack = new QPushButton("Creer une faille");
    tunnelDisplayButton = new QCheckBox("Afficher");
    tunnelLayout->addWidget(createVerticalGroup({tunnelCreateMatter, tunnelRemoveMatter, tunnelCreateCrack}));
    tunnelLayout->addWidget(createVerticalGroup({addControlPointButton, tunnelClearControlPointButton}));
    tunnelLayout->addWidget(createSliderGroup("Taille", tunnelSizeSlider));
    tunnelLayout->addWidget(createSliderGroup("Force", tunnelStrengthSlider));
    tunnelLayout->addWidget(tunnelDisplayButton);

    return this->tunnelLayout;
}

