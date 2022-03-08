#include "SpaceColonizationInterface.h"

#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/InterfaceUtils.h"
#include <QGLViewer/manipulatedCameraFrame.h>

SpaceColonizationInterface::SpaceColonizationInterface()
{

}

SpaceColonizationInterface::~SpaceColonizationInterface()
{
}

void SpaceColonizationInterface::display()
{
    if (!isHidden)
    {
        for (ControlPoint*& ctrl : this->controlPoints) {
            ctrl->display();
        }
        this->startingPoint->mesh.shader->setVector("color", std::vector<float>({100/255.f, 10/255.f, 255/255.f, 1.f}));
        this->startingPoint->display();
        this->pathsMeshes.shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.f}));
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
    int nb_special_nodes = 8;
    std::vector<Vector3> keyPoints;
    for (int i = 0; i < nb_special_nodes; i++) {
        keyPoints.push_back(poissonGraph.nodes[i * poissonGraph.nodes.size() / (float)nb_special_nodes]->pos);
    }
    Vector3 startPos(0, 0, 0);
//    startPos = Vector3(10, 10, 0); // TO REMOVE
//    for (auto& pt : keyPoints)
//        pt.z = 80; // TO REMOVE
    this->colonizer = new TreeColonisationAlgo::TreeColonisation(keyPoints, startPos, 10.f);
    this->colonizer->nodeMinDistance = this->voxelGrid->blockSize;
    this->colonizer->nodeMaxDistance = this->voxelGrid->blockSize * this->voxelGrid->chunkSize * 5;

    for (size_t i = 0; i < keyPoints.size(); i++) {
        this->controlPoints.push_back(new ControlPoint(keyPoints[i], 5.f));
        QObject::connect(this->controlPoints[i], &ControlPoint::modified, this, &SpaceColonizationInterface::computeKarst);
    }
    this->startingPoint = new ControlPoint(startPos, 5.f);
    QObject::connect(this->startingPoint, &ControlPoint::modified, this, &SpaceColonizationInterface::computeKarst);
}

void SpaceColonizationInterface::initSpaceColonizer()
{
    std::vector<Vector3> newNodes;
    for (ControlPoint*& ctrl : this->controlPoints)
        newNodes.push_back(ctrl->position);
    this->colonizer->startPosition = this->startingPoint->position;
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
    if (!this->visitingCamera)
        this->visitingCamera = new qglviewer::Camera();
    this->visitingCamera->setFieldOfView(3.141592 / 3.f);
    this->visitingCamera->setUpVector(Vector3(0, 0, 1));
    this->karstPaths.clear();
    std::vector<Vector3> pathPositions;
    std::vector<std::vector<Vector3>> allPaths = this->colonizer->simplifyPaths();
    std::vector<BSpline> simplifiedKarstParths;

    for (const auto& path : allPaths)
    {
        BSpline simplifiedPath = BSpline(path); //.simplifyByRamerDouglasPeucker(3.0);
        this->karstPaths.push_back(BSpline(path));
        simplifiedKarstParths.push_back(simplifiedPath);
        for (size_t i = 0; i < path.size() - 1; i++) {
            pathPositions.push_back(path[i]);
            pathPositions.push_back(path[i+1]);
        }
    }
    this->pathsMeshes.fromArray(pathPositions);
    if (allPaths.size() > 0) {
        this->cameraConstraint = new PathCameraConstraint(this->visitingCamera, simplifiedKarstParths);
        this->visitingCamera->frame()->setConstraint(this->cameraConstraint);
        this->visitingCamera->setPosition(pathPositions[0]);
        this->visitingCamera->lookAt(pathPositions[1]);
    }

    Q_EMIT this->karstPathUpdated();
}

void SpaceColonizationInterface::createKarst(bool usingSpheres)
{
    if (this->karstPaths.empty())
        this->updateKarstPath();
    UnderwaterErosion erod(this->voxelGrid, 15.f, 3.f, 10);
    erod.CreateMultipleTunnels(this->karstPaths, false, usingSpheres);

    std::string nodes, links;
    std::tie(nodes, links) = this->colonizer->toFile();
    std::ofstream nodeFile;
    nodeFile.open("../karst_nodes.dat", std::ios::out);
    nodeFile << nodes;
    nodeFile.close();
    std::ofstream linkFile;
    linkFile.open("../karst_links.dat", std::ios::out);
    linkFile << links;
    linkFile.close();

}

QHBoxLayout *SpaceColonizationInterface::createGUI()
{
    this->spaceColonizationLayout = new QHBoxLayout;

    QPushButton* spaceColonizerPreviewButton = new QPushButton("Calculer");
    QPushButton* spaceColonizerConfirmButton = new QPushButton("Creer le karst");
    QPushButton* spaceColonizerQuickConfirmButton = new QPushButton("Tunnel rond");
    QCheckBox* spaceColonizerDisplay = new QCheckBox("Afficher");
    QCheckBox* useAsMainCamera = new QCheckBox("Observer l'interieur");
    FancySlider* spaceColonizerSegmentSize = new FancySlider(Qt::Orientation::Horizontal, 1.0, 40.0, 1.0);
    FancySlider* spaceColonizerRandomness = new FancySlider(Qt::Orientation::Horizontal, 0.0, 1.0, 0.1);
    this->spaceColonizationLayout->addWidget(createVerticalGroup({spaceColonizerPreviewButton, spaceColonizerConfirmButton, spaceColonizerQuickConfirmButton}));
    this->spaceColonizationLayout->addWidget(createVerticalGroup({
                                                           createSliderGroup("Taille des segments", spaceColonizerSegmentSize),
                                                           createSliderGroup("Aleatoire", spaceColonizerRandomness)
                                                       }));
    this->spaceColonizationLayout->addWidget(createVerticalGroup({spaceColonizerDisplay, useAsMainCamera}));

    QObject::connect(spaceColonizerPreviewButton, &QPushButton::pressed, this, &SpaceColonizationInterface::computeKarst);
    QObject::connect(spaceColonizerConfirmButton, &QPushButton::pressed, [=](){ this->createKarst(false); } );
    QObject::connect(spaceColonizerQuickConfirmButton, &QPushButton::pressed, this, [=](){ this->createKarst(true); } );
    QObject::connect(spaceColonizerRandomness, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->randomness = val; this->computeKarst(); } );
    QObject::connect(spaceColonizerSegmentSize, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->segmentLength = val; this->computeKarst(); } );
    QObject::connect(spaceColonizerDisplay, &QCheckBox::toggled, this, [=](bool display){ this->isHidden = !display; } );
    QObject::connect(useAsMainCamera, &QCheckBox::toggled, this, [=](bool display){
        Q_EMIT this->useAsMainCamera(this->visitingCamera, display);
    } );

    spaceColonizerDisplay->setChecked(!isHidden);
    useAsMainCamera->setChecked(false);

    return this->spaceColonizationLayout;
}
