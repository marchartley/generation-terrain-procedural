#include "KarstPathGenerationInterface.h"

#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/InterfaceUtils.h"

KarstPathGenerationInterface::KarstPathGenerationInterface(QWidget *parent)
    : ActionInterface("karstpeytavie", "Karts system (Peytavie)", "digging", "Create karsts with graphs", "karst_peytavie_button.png", parent)
{
    this->karstCreator = new KarstPathsGeneration();
    this->sourceControlPoint = std::make_unique<ControlPoint>(Vector3(), 5.f);
    this->fractureVector = std::make_unique<InteractiveVector>();
    this->waterHeightSlider = std::make_unique<Slider3D>(Vector3(), 1.0, 0.0, 0.0, 1.0, Slider3DOrientation::Z);
    this->waterHeightSlider->sliderControlPoint->setGrabberStateColor({
                                                                          {INACTIVE, {0/255.f, 0/255.f, 180/255.f, 1.f}},
                                                                          {ACTIVE, {0/255.f, 0/255.f, 255/255.f, 1.f}},
                                                                      });
    QObject::connect(this->sourceControlPoint.get(), &ControlPoint::pointModified, this, &KarstPathGenerationInterface::computeKarst);
    QObject::connect(this->fractureVector.get(), &InteractiveVector::modified, this, &KarstPathGenerationInterface::updateFracture);
    QObject::connect(this->waterHeightSlider.get(), &Slider3D::valueChanged, this, &KarstPathGenerationInterface::updateWaterHeight);
}

KarstPathGenerationInterface::~KarstPathGenerationInterface()
{
    delete karstCreator;
    sourceControlPoint.reset();
    fractureVector.reset();
    waterHeightSlider.reset();
}

void KarstPathGenerationInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    this->sourceControlPoint->display();
    this->fractureVector->display();
    if (this->waterLevelMesh.shader != nullptr)
        this->waterLevelMesh.shader->setVector("color", std::vector<float>({0/255.f, 0/255.f, 255/255.f, 1.f}));
    this->waterLevelMesh.display(GL_LINES, 5.f);
    if (this->pathsMeshes.shader != nullptr)
        this->pathsMeshes.shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.f}));
    this->pathsMeshes.display(GL_LINES, 5.f);
    this->waterHeightSlider->display();
}

void KarstPathGenerationInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        // Nothing for now
    }
}

void KarstPathGenerationInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
//    this->voxelGrid = voxelGrid;
    this->AABBoxMinPos = Vector3(0, 0, 0);
    this->AABBoxMaxPos = voxelGrid->getDimensions(); // Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ); // * voxelGrid->blockSize;
    this->fractureVector->setPositions(AABBoxMinPos, Vector3(AABBoxMinPos.x, AABBoxMinPos.y, AABBoxMaxPos.z));
    this->waterHeightSlider->setPositions(AABBoxMinPos + Vector3(-10, -10, 0), Vector3(AABBoxMinPos.x, AABBoxMinPos.y, AABBoxMaxPos.z) + Vector3(-10, -10, 0));
    this->waterHeightSlider->sliderControlPoint->setGrabberStateColor({
                                                                          {INACTIVE, {0/255.f, 0/255.f, 180/255.f, 1.f}},
                                                                          {ACTIVE, {0/255.f, 0/255.f, 255/255.f, 1.f}},
                                                                      });

    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->getPosition().z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });


    GridF voxels = voxelGrid->getVoxelValues();
    GridI availableGrid(voxelGrid->getDimensions(), 0.f); // (voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ, 0);
    for (int x = 0; x < availableGrid.sizeX; x++) {
        for (int y = 0; y < availableGrid.sizeY; y++) {
            for (int z = 0; z < availableGrid.sizeZ; z++) {
                float voxelVal = voxels.at(x, y, z);
                if (voxelVal > 0) {
                    availableGrid.at(x, y, z) = 1;
                }
            }
        }
    }
    this->karstCreator = new KarstPathsGeneration(availableGrid, voxelGrid->getDimensions(), 20.f); //Vector3(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ), 20.f);
}

void KarstPathGenerationInterface::updateFracture(const Vector3& newFractureDir)
{
    if (karstCreator->fracturesDirections.empty())
        this->karstCreator->addFractureDirection(FractureDirection(Vector3(0, 0, 0), 1.f));
    karstCreator->fracturesDirections.back().direction = newFractureDir;
    /// For demonstration purpose only
    this->computeKarst();
}

void KarstPathGenerationInterface::updateWaterHeight(float newHeight)
{
    this->AABBoxMinPos = Vector3(0, 0, 0);
    this->AABBoxMaxPos = voxelGrid->getDimensions();
    if (karstCreator->waterHeights.empty())
        this->karstCreator->addWaterHeight(WaterHeight(0.f, 20.f, 2.0));
    karstCreator->waterHeights.back().height = newHeight;
    Vector3 offX = Vector3(AABBoxMaxPos.x - AABBoxMinPos.x, 0, 0);
    Vector3 offY = Vector3(0, AABBoxMaxPos.y - AABBoxMinPos.y, 0);
    Vector3 offZ = Vector3(0, 0, this->waterHeightSlider->sliderControlPoint->getPosition().z);
    this->waterLevelMesh.fromArray({
                                       {AABBoxMinPos + offZ, AABBoxMinPos + offX + offZ,
                                        AABBoxMinPos + offX + offZ, AABBoxMinPos + offX + offY + offZ,
                                        AABBoxMinPos + offX + offY + offZ, AABBoxMinPos + offY + offZ,
                                        AABBoxMinPos + offY + offZ, AABBoxMinPos + offZ}
                                   });
    /// For demonstration purpose only
    this->computeKarst();
}

void KarstPathGenerationInterface::computeKarst()
{
    if (karstCreator->graph.nodes.empty()) return;
    GridF porosityMap(voxelGrid->getDimensions(), 1.f); // voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ, 1.f);
    GridF voxelValues = voxelGrid->getVoxelValues();
    for (int x = 0; x < porosityMap.sizeX; x++) {
        for (int y = 0; y < porosityMap.sizeY; y++) {
            for (int z = 0; z < porosityMap.sizeZ; z++) {
                float voxelVal = voxelValues.at(x, y, z);
                if (voxelVal > 0) {
                    porosityMap.at(x, y, z) = voxelVal;
                }
            }
        }
    }
    this->karstCreator->porosityMap = porosityMap;
    int nb_special_nodes = 10;
    this->karstCreator->resetSpecialNodes();
    this->karstCreator->setSpecialNode(0, SOURCE);
    this->karstCreator->graph.nodes[0]->pos = this->sourceControlPoint->getPosition();

    for (int i = 0; i < nb_special_nodes - 1; i++)
        this->karstCreator->setSpecialNode((i + 1) * this->karstCreator->graph.nodes.size() / nb_special_nodes, FORCED_POINT);

    this->karstCreator->createEdges(10, 50.f);
    this->karstCreator->computeAllPathsBetweenSpecialNodes();

    this->updateKarstPath();

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
    if (this->karstCreator->finalPaths.size() > 0 && pathPositions.size() > 1) {
        this->cameraConstraint = new PathCameraConstraint(this->visitingCamera, karstPaths);
        this->visitingCamera->frame()->setConstraint(this->cameraConstraint);
        this->visitingCamera->setPosition(pathPositions[0]);
        this->visitingCamera->lookAt(pathPositions[1]);
    }
    Q_EMIT this->updated();
}

void KarstPathGenerationInterface::createKarst()
{
    UnderwaterErosion erod(this->voxelGrid.get(), 10.f, 2.f, 10);
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



QLayout *KarstPathGenerationInterface::createGUI()
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
    karstCreationLayout->addWidget(createVerticalGroup({/*karstCreationDisplay, */karstCreationChangeCam}));

    QObject::connect(karstCreationPreviewButton, &QPushButton::pressed, this, [=](){ this->computeKarst(); } );
    QObject::connect(karstCreationConfirmButton, &QPushButton::pressed, this, [=](){ this->createKarst(); } );
    QObject::connect(karstCreationDistanceWeights, &FancySlider::floatValueChanged, this, [=](float val){
        this->karstCreator->distanceWeight = val; this->computeKarst(); } );
    QObject::connect(karstCreationFractureWeights, &FancySlider::floatValueChanged, this, [=](float val){
        this->karstCreator->fractureWeight = val; this->computeKarst(); } );
    QObject::connect(karstCreationWaterWeights, &FancySlider::floatValueChanged, this, [=](float val){
        this->karstCreator->waterHeightWeight = val; this->computeKarst(); } );
    QObject::connect(karstCreationPorosityWeights, &FancySlider::floatValueChanged, this, [=](float val){
        this->karstCreator->porosityWeight = val; this->computeKarst(); } );
    QObject::connect(karstCreationGamma, &FancySlider::floatValueChanged, this, [=](float val){
        this->karstCreator->gamma = val; this->computeKarst(); } );
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
