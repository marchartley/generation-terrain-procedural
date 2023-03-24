#include "FlowFieldInterface.h"

FlowFieldInterface::FlowFieldInterface(QWidget *parent) : ActionInterface("flowfield", parent)
{
//    this->createGUI();
}

/*void FlowFieldInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}*/

void FlowFieldInterface::display()
{
    if (this->flowMesh.shader != nullptr)
        this->flowMesh.shader->setVector("color", std::vector<float>({0.f, 0.f, 1.f, .4f}));
    this->flowMesh.display(GL_LINES, 5.f);
}

void FlowFieldInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        this->voxelGrid->computeFlowfield();
    }
}

void FlowFieldInterface::recomputeFlowfield()
{
    this->voxelGrid->computeFlowfield();
    this->updateFlowfieldDebugMesh();

    this->addTerrainAction(nlohmann::json({
                                              {}
                                          }));
    Q_EMIT this->updated();
}

void FlowFieldInterface::updateFlowfieldDebugMesh()
{
    Matrix3<Vector3> flowNormalized = this->voxelGrid->getFlowfield();
    float max = -1, min = std::numeric_limits<float>::max();
    for (auto& v : flowNormalized){
        max = std::max(max, v.norm2());
        min = std::min(min, v.norm2());
    }
    for (auto& v : flowNormalized)
        v /= max;
    std::vector<Vector3> normals;
    for (int x = this->voxelGrid->fluidSimRescale; x < this->voxelGrid->getSizeX()-1; x+= this->voxelGrid->fluidSimRescale) {
        for (int y = this->voxelGrid->fluidSimRescale; y < this->voxelGrid->getSizeY()-1; y+= this->voxelGrid->fluidSimRescale) {
            for (int z = this->voxelGrid->fluidSimRescale; z < this->voxelGrid->getSizeZ() - 1; z+= this->voxelGrid->fluidSimRescale) {
                normals.push_back(Vector3(x, y, z) + Vector3(.5f, .5f, .5f)); // - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
                normals.push_back(Vector3(x, y, z) + (flowNormalized.at(x, y, z) * (float)voxelGrid->fluidSimRescale) + Vector3(.5f, .5f, .5f)); // - Vector3(this->voxelGrid->sizeX/2.0, this->voxelGrid->sizeY/2.0));
            }
        }
    }
    this->flowMesh.fromArray(normals);
    this->flowMesh.update();

//    Q_EMIT this->updated();
}

void FlowFieldInterface::hide()
{
    this->flowMesh.hide();
    CustomInteractiveObject::hide();
}

void FlowFieldInterface::show()
{
    this->flowMesh.show();
    CustomInteractiveObject::show();
}

QLayout* FlowFieldInterface::createGUI()
{
    this->flowFieldLayout = new QHBoxLayout;

    flowFieldComputeButton = new QPushButton("Calculer");
    flowFieldDisplayButton = new QCheckBox("Afficher");
    flowFieldLayout->addWidget(flowFieldComputeButton);
//    flowFieldLayout->addWidget(flowFieldDisplayButton);

    flowFieldDisplayButton->setChecked(this->visible);

    QObject::connect(flowFieldDisplayButton, &QCheckBox::toggled, this, &FlowFieldInterface::setVisibility);
    QObject::connect(flowFieldComputeButton, &QPushButton::pressed, this, &FlowFieldInterface::recomputeFlowfield );

    return this->flowFieldLayout;
}
