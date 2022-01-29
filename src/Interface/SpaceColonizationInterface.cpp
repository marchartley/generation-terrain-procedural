#include "SpaceColonizationInterface.h"

#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/InterfaceUtils.h"

SpaceColonizationInterface::SpaceColonizationInterface()
{

}

void SpaceColonizationInterface::display()
{
    if (!isHidden)
    {
        for (ControlPoint*& ctrl : this->controlPoints) {
            ctrl->display();
        }
        this->pathsMeshes.shader->setVector("color", std::vector<float>({0/255.f, 255/255.f, 0/255.f, 1.f}));
        this->pathsMeshes.display(GL_LINES, 5.f);
    }
}

void SpaceColonizationInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    Matrix3<int> availableGrid(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ, 0);
    for (int x = 0; x < availableGrid.sizeX; x++) {
        for (int y = 0; y < availableGrid.sizeY; y++) {
            for (int z = 0; z < availableGrid.sizeZ; z++) {
                float voxelVal = voxelGrid->getVoxelValue(x, y, z);
                if (voxelVal > 0) {
                    availableGrid.at(x, y, z) = 1;
                }
            }
        }
    }
    FastPoissonGraph<TreeColonisationAlgo::NODE_TYPE> poissonGraph(availableGrid, 20.f);
    int nb_special_nodes = 10;
    std::vector<Vector3> keyPoints;
    for (int i = 0; i < nb_special_nodes; i++) {
        keyPoints.push_back(poissonGraph.nodes[i * poissonGraph.nodes.size() / (float)nb_special_nodes]->pos);
    }
    this->colonizer = new TreeColonisationAlgo::TreeColonisation(keyPoints, Vector3(0, 0, 0), 10.f);
    this->colonizer->nodeMinDistance = this->voxelGrid->blockSize;
    this->colonizer->nodeMaxDistance = this->voxelGrid->blockSize * this->voxelGrid->chunkSize * 5;

    for (size_t i = 0; i < keyPoints.size(); i++) {
        this->controlPoints.push_back(new ControlPoint(keyPoints[i], 5.f));
        QObject::connect(this->controlPoints[i], &ControlPoint::modified, this, &SpaceColonizationInterface::computeKarst);
    }
}

void SpaceColonizationInterface::initSpaceColonizer()
{
    std::vector<Vector3> newNodes;
    for (ControlPoint*& ctrl : this->controlPoints)
        newNodes.push_back(ctrl->position);
    this->colonizer->reset(newNodes);
}

void SpaceColonizationInterface::computeKarst()
{
    this->initSpaceColonizer();
    this->colonizer->process();
    this->updateKarstPath();
}

void SpaceColonizationInterface::updateKarstPath()
{
    this->karstPaths.clear();
    std::vector<Vector3> pathPositions;
    std::vector<std::vector<Vector3>> allPaths = this->colonizer->simplifyPaths();
    std::cout << "For " << this->colonizer->segments.size() << " segments, " << allPaths.size() << " paths created" << std::endl;
    for (const auto& path : allPaths)
    {
        pathPositions.insert(pathPositions.end(), path.begin(), path.end());
        this->karstPaths.push_back(BSpline(path));
    }
    this->pathsMeshes.fromArray(pathPositions);
}

void SpaceColonizationInterface::createKarst()
{
    this->updateKarstPath();
    UnderwaterErosion erod(this->voxelGrid, 20.f, 2.f, 10);
//    std::cout << this->karstPaths.size() << " tunnels to create." << std::endl;
    erod.CreateMultipleTunnels(this->karstPaths);
}

QHBoxLayout *SpaceColonizationInterface::createGUI()
{
    this->spaceColonizationLayout = new QHBoxLayout;

    QPushButton* spaceColonizerPreviewButton = new QPushButton("Calculer");
    QPushButton* spaceColonizerConfirmButton = new QPushButton("Creer le karst");
    QCheckBox* spaceColonizerDisplay = new QCheckBox("Afficher");
    FancySlider* spaceColonizerSegmentSize = new FancySlider(Qt::Orientation::Horizontal, 1.0, 40.0, 1.0);
    FancySlider* spaceColonizerRandomness = new FancySlider(Qt::Orientation::Horizontal, 0.0, 1.0, 0.1);
    this->spaceColonizationLayout->addWidget(createVerticalGroup({spaceColonizerPreviewButton, spaceColonizerConfirmButton}));
    this->spaceColonizationLayout->addWidget(createVerticalGroup({
                                                           createSliderGroup("Taille des segments", spaceColonizerSegmentSize),
                                                           createSliderGroup("Aleatoire", spaceColonizerRandomness)
                                                       }));
    this->spaceColonizationLayout->addWidget(spaceColonizerDisplay);

    QObject::connect(spaceColonizerPreviewButton, &QPushButton::pressed, this, &SpaceColonizationInterface::computeKarst);
    QObject::connect(spaceColonizerConfirmButton, &QPushButton::pressed, this, &SpaceColonizationInterface::createKarst); // [=](){ this->viewer->createTunnelFromSpaceColonizer(); } );
    QObject::connect(spaceColonizerRandomness, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->randomness = val; } );
    QObject::connect(spaceColonizerSegmentSize, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->segmentLength = val; } );
    QObject::connect(spaceColonizerDisplay, &QCheckBox::toggled, this, [=](bool display){ this->isHidden = !display; } );

    spaceColonizerDisplay->setChecked(!isHidden);

    return this->spaceColonizationLayout;
}
