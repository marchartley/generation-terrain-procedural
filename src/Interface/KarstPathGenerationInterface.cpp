#include "KarstPathGenerationInterface.h"

#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/InterfaceUtils.h"


KarstPathGenerationInterface::KarstPathGenerationInterface()
{
    this->sourceControlPoint = new ControlPoint(Vector3(0, 0, 0), 5.f);
    this->fractureVector = new InteractiveVector();
    this->waterHeightSlider = new Slider3D();
    this->waterHeightSlider->sliderControlPoint->setGrabberStateColor({
                                                                          {INACTIVE, {0/255.f, 0/255.f, 180/255.f, 1.f}},
                                                                          {ACTIVE, {0/255.f, 0/255.f, 255/255.f, 1.f}},
                                                                      });
    QObject::connect(this->sourceControlPoint, &ControlPoint::modified, this, &KarstPathGenerationInterface::computeKarst);
}

void KarstPathGenerationInterface::display()
{
    this->sourceControlPoint->display();
    this->fractureVector->display();
    this->waterLevelMesh.shader->setVector("color", std::vector<float>({0/255.f, 0/255.f, 255/255.f, 1.f}));
    this->waterLevelMesh.display(GL_LINES, 5.f);
    this->pathsMeshes.shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.f}));
    this->pathsMeshes.display(GL_LINES, 5.f);
    this->waterHeightSlider->display();
}

void KarstPathGenerationInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    this->AABBoxMinPos = Vector3(0, 0, 0);
    this->AABBoxMaxPos = Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ) * voxelGrid->blockSize;
    this->fractureVector = new InteractiveVector(AABBoxMinPos, Vector3(AABBoxMinPos.x, AABBoxMinPos.y, AABBoxMaxPos.z));
    this->waterHeightSlider = new Slider3D(Vector3(AABBoxMinPos.x - 10, AABBoxMinPos.y - 10, AABBoxMinPos.z), Vector3(AABBoxMinPos.x - 10, AABBoxMinPos.y - 10, AABBoxMaxPos.z), 0.f, AABBoxMinPos.z, AABBoxMaxPos.z);
    this->waterHeightSlider->sliderControlPoint->setGrabberStateColor({
                                                                          {INACTIVE, {0/255.f, 0/255.f, 180/255.f, 1.f}},
                                                                          {ACTIVE, {0/255.f, 0/255.f, 255/255.f, 1.f}},
                                                                      });

    QObject::connect(this->fractureVector, &InteractiveVector::modified, this, &KarstPathGenerationInterface::updateFracture);
    QObject::connect(this->waterHeightSlider, &Slider3D::valueChanged, this, &KarstPathGenerationInterface::updateWaterHeight);

    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->position.z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });


    std::cout << "Error during the available matrix?" << std::endl;
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
    std::cout << "No. Error during init karst paths gen?" << std::endl;
    this->karstCreator = new KarstPathsGeneration(availableGrid, Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ), 10.f);
    std::cout << "No." << std::endl;
}

void KarstPathGenerationInterface::updateFracture(Vector3 newFractureDir)
{
    if (karstCreator->fracturesDirections.empty())
        this->karstCreator->addFractureDirection(FractureDirection(Vector3(0, 0, 0), 1.f));
    karstCreator->fracturesDirections.back().direction = newFractureDir;
}

void KarstPathGenerationInterface::updateWaterHeight(float newHeight)
{
    if (karstCreator->waterHeights.empty())
        this->karstCreator->addWaterHeight(WaterHeight(0.f, 1.f));
    karstCreator->waterHeights.back().height = newHeight;
    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->position.z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });
}

void KarstPathGenerationInterface::computeKarst()
{
    std::cout << "Error during the porosity matrix?" << std::endl;
    Matrix3<float> porosityMap(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ, 1.f);
    for (int x = 0; x < porosityMap.sizeX; x++) {
        for (int y = 0; y < porosityMap.sizeY; y++) {
            for (int z = 0; z < porosityMap.sizeZ; z++) {
                float voxelVal = voxelGrid->getVoxelValue(x, y, z);
                if (voxelVal > 0) {
                    porosityMap.at(x, y, z) = voxelVal;
                }
            }
        }
    }
    this->karstCreator->porosityMap = porosityMap;
    int nb_special_nodes = 10;
    std::cout << "No. Error during reset special?" << std::endl;
    this->karstCreator->resetSpecialNodes();
    this->karstCreator->setSpecialNode(0, SOURCE);
    this->karstCreator->graph.nodes[0]->pos = this->sourceControlPoint->position;

    std::cout << "No. Error during new special?" << std::endl;
    for (int i = 0; i < nb_special_nodes - 1; i++)
        this->karstCreator->setSpecialNode((i + 1) * this->karstCreator->graph.nodes.size() / nb_special_nodes, FORCED_POINT);

    std::cout << "No. Error during edges?" << std::endl;
    this->karstCreator->createEdges(10, 50.f);
    std::cout << "No. Error during paths?" << std::endl;
    this->karstCreator->computeAllPathsBetweenSpecialNodes();

    std::cout << "No. Error during update?" << std::endl;
    this->updateKarstPath();
    std::cout << "No." << std::endl;

}

void KarstPathGenerationInterface::updateKarstPath()
{
    if (!this->visitingCamera)
        this->visitingCamera = new qglviewer::Camera();
    this->visitingCamera->setFieldOfView(3.141592 / 3.f);
    this->visitingCamera->setUpVector(Vector3(0, 0, 1));
    this->karstPaths.clear();
    std::vector<Vector3> pathPositions;
    for (const auto& path : this->karstCreator->finalPaths)
    {
        for (size_t i = 0; i < path.size() - 1; i++) {
            pathPositions.push_back(this->karstCreator->getNodePos(path[i    ]));
            pathPositions.push_back(this->karstCreator->getNodePos(path[i + 1]));
            karstPaths.push_back(BSpline(std::vector<Vector3>({this->karstCreator->getNodePos(path[i    ]),
                                                               this->karstCreator->getNodePos(path[i + 1])})));
        }
    }
    this->pathsMeshes.fromArray(pathPositions);
    if (this->karstCreator->finalPaths.size() > 0) {
        this->cameraConstraint = new PathCameraConstraint(this->visitingCamera, karstPaths);
        this->visitingCamera->frame()->setConstraint(this->cameraConstraint);
        this->visitingCamera->setPosition(pathPositions[0]);
        this->visitingCamera->lookAt(pathPositions[1]);
    }
    Q_EMIT this->karstPathUpdated();
}

void KarstPathGenerationInterface::createKarst()
{
    UnderwaterErosion erod(this->voxelGrid, 20.f, 2.f, 10);
    erod.CreateMultipleTunnels(this->karstPaths);
    //    Q_EMIT this->update();
}

void KarstPathGenerationInterface::hide()
{
    this->sourceControlPoint->hide();
    this->fractureVector->hide();
    this->waterHeightSlider->hide();
    this->waterLevelMesh.hide();
    this->pathsMeshes.hide();
    CustomInteractiveObject::hide();
}

void KarstPathGenerationInterface::show()
{
    this->sourceControlPoint->show();
    this->fractureVector->show();
    this->waterHeightSlider->show();
    this->waterLevelMesh.show();
    this->pathsMeshes.show();
    CustomInteractiveObject::show();
}



QHBoxLayout *KarstPathGenerationInterface::createGUI()
{
    this->karstCreationLayout = new QHBoxLayout;

    QPushButton* karstCreationPreviewButton = new QPushButton("Calculer");
    QPushButton* karstCreationConfirmButton = new QPushButton("Creer le karst");
    QCheckBox* karstCreationDisplay = new QCheckBox("Afficher");
    QCheckBox* karstCreationChangeCam = new QCheckBox("Observer l'interieur");
    FancySlider* karstCreationDistanceWeights = new FancySlider(Qt::Orientation::Horizontal, 0, 10, 0.1);
    FancySlider* karstCreationFractureWeights = new FancySlider(Qt::Orientation::Horizontal, 0, 10, 0.1);
    FancySlider* karstCreationWaterWeights = new FancySlider(Qt::Orientation::Horizontal, 0, 10, 0.1);
    FancySlider* karstCreationPorosityWeights = new FancySlider(Qt::Orientation::Horizontal, 0, 10, 0.1);
    FancySlider* karstCreationGamma = new FancySlider(Qt::Orientation::Horizontal, 1, 10, 0.1);
    FancySlider* karstCreationTortuosity = new FancySlider(Qt::Orientation::Horizontal, 0, 10, 0.1);
    karstCreationLayout->addWidget(createVerticalGroup({karstCreationPreviewButton, karstCreationConfirmButton}));
    karstCreationLayout->addWidget(createMultipleSliderGroup({
                                                           {"Distances", karstCreationDistanceWeights},
                                                           {"Fractures", karstCreationFractureWeights},
                                                           {"Porosite", karstCreationPorosityWeights},
                                                           {"Niveau d'eau", karstCreationWaterWeights},
                                                       }));
    karstCreationLayout->addWidget(createVerticalGroup({
                                                           createSliderGroup("Gamma", karstCreationGamma),
                                                           createSliderGroup("Tortuosite (m)", karstCreationTortuosity)
                                                       }));
    karstCreationLayout->addWidget(createVerticalGroup({karstCreationDisplay, karstCreationChangeCam}));

    QObject::connect(karstCreationPreviewButton, &QPushButton::pressed, this, [=](){ this->computeKarst(); } );
    QObject::connect(karstCreationConfirmButton, &QPushButton::pressed, this, [=](){ this->createKarst(); } );
    QObject::connect(karstCreationDistanceWeights, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->distanceWeight = val; } );
    QObject::connect(karstCreationFractureWeights, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->fractureWeight = val; } );
    QObject::connect(karstCreationWaterWeights, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->waterHeightWeight = val; } );
    QObject::connect(karstCreationPorosityWeights, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->porosityWeight = val; } );
    QObject::connect(karstCreationGamma, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->gamma = val; } );
    QObject::connect(karstCreationTortuosity, &FancySlider::floatValueChanged, this, [=](float val){ this->karstCreator->updateTortuosity(val, {0}); this->updateKarstPath(); } );
    QObject::connect(karstCreationDisplay, &QCheckBox::toggled, this, &KarstPathGenerationInterface::setVisibility);
    QObject::connect(karstCreationChangeCam, &QCheckBox::toggled, this, [=](bool display){
        Q_EMIT this->useAsMainCamera(this->visitingCamera, display);
    } );


    if (this->karstCreator) {
        karstCreationDistanceWeights->setfValue(this->karstCreator->distanceWeight);
        karstCreationFractureWeights->setfValue(this->karstCreator->fractureWeight);
        karstCreationWaterWeights->setfValue(this->karstCreator->waterHeightWeight);
        karstCreationPorosityWeights->setfValue(this->karstCreator->porosityWeight);
        karstCreationGamma->setfValue(this->karstCreator->gamma);
    }
    karstCreationTortuosity->setfValue(0.f);
    karstCreationDisplay->setChecked(!this->isHidden());
    karstCreationChangeCam->setChecked(false);

    return this->karstCreationLayout;
}
