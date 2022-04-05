#include "FlowFieldInterface.h"

FlowFieldInterface::FlowFieldInterface()
{
    this->createGUI();
}

void FlowFieldInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );
}

void FlowFieldInterface::hide()
{
    CustomInteractiveObject::hide();
}

void FlowFieldInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* FlowFieldInterface::createGUI()
{
    this->flowFieldLayout = new QHBoxLayout;

    flowFieldComputeButton = new QPushButton("Calculer");
    flowFieldDisplayButton = new QCheckBox("Afficher");
    flowFieldLayout->addWidget(flowFieldComputeButton);
    flowFieldLayout->addWidget(flowFieldDisplayButton);

    QObject::connect(flowFieldDisplayButton, &QCheckBox::toggled, this, &FlowFieldInterface::setVisibility);
    flowFieldDisplayButton->setChecked(this->visible);

    if (voxelGrid != nullptr)
        QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, [=](){ voxelGrid->computeFlowfield(); } );

    return this->flowFieldLayout;
}
