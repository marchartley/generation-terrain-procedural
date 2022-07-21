#include "SpaceColonizationInterface.h"

#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/InterfaceUtils.h"
#include <QGLViewer/manipulatedCameraFrame.h>

SpaceColonizationInterface::SpaceColonizationInterface(QWidget *parent) : ActionInterface("space_colonization", parent)
{
    this->startingPoint = std::make_unique<ControlPoint>();
//    this->
}

SpaceColonizationInterface::~SpaceColonizationInterface()
{
}

void SpaceColonizationInterface::display()
{
    if (this->isVisible())
    {
        // Hide control points when visiting
        if (visitingCamera->isVisiting) {
            for (auto& ctrl : this->controlPoints) {
                ctrl->hide();
            }
            this->startingPoint->hide();
        } else {
            for (auto& ctrl : this->controlPoints) {
                ctrl->show();
            }
            this->startingPoint->show();
        }


        for (auto& ctrl : this->controlPoints) {
            ctrl->display();
        }
        if (this->startingPoint->mesh.shader != nullptr)
            this->startingPoint->mesh.shader->setVector("color", std::vector<float>({100/255.f, 10/255.f, 255/255.f, 1.f}));
        this->startingPoint->display();
        if (this->pathsMeshes.shader != nullptr)
            this->pathsMeshes.shader->setVector("color", std::vector<float>({255/255.f, 0/255.f, 0/255.f, 1.f}));
        this->pathsMeshes.display(GL_LINES, 5.f);
    } else {
        // Just make sure that they are not interacting when hidden
        for (auto& ctrl : this->controlPoints) {
            ctrl->hide();
        }
    }
}

void SpaceColonizationInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        Vector3 startingPoint = json_to_vec3(parameters.at("starting_point"));
        std::vector<Vector3> control_points;
        for (auto& ctrl_json : parameters.at("control_points"))
            control_points.push_back(json_to_vec3(ctrl_json));
        float width = parameters.at("width").get<float>();
        float randomness = parameters.at("randomness").get<float>();
        float segmentLength = parameters.at("segment_length").get<float>();
        bool usingSpheres = parameters.at("using_spheres").get<bool>();
        float minDistance = parameters.at("min_distance").get<float>();
        float maxDistance = parameters.at("max_distance").get<float>();

        this->colonizer = new TreeColonisationAlgo::TreeColonisation(control_points, startingPoint, segmentLength, randomness);
        colonizer->nodeMinDistance = minDistance;
        colonizer->nodeMaxDistance = maxDistance;
        this->colonizer->process();
        this->updateKarstPath();
        UnderwaterErosion erod(this->voxelGrid, width, 3.f, 10);
        erod.CreateMultipleTunnels(this->karstPaths, false, usingSpheres);

    }
}

void SpaceColonizationInterface::hide()
{
    for (auto& ctrl : this->controlPoints) {
        ctrl->hide();
    }
    this->startingPoint->hide();
    this->pathsMeshes.hide();
    CustomInteractiveObject::hide();
}

void SpaceColonizationInterface::show()
{
    for (auto& ctrl : this->controlPoints) {
        ctrl->show();
    }
    this->startingPoint->show();
    this->pathsMeshes.show();
    CustomInteractiveObject::show();
}

void SpaceColonizationInterface::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Up) {
        if (visitingCamera->isVisiting)
            visitingCamera->moveForward(1.f);
    } else if (event->key() == Qt::Key_Down) {
        if (visitingCamera->isVisiting)
            visitingCamera->moveBackward(1.f);
    }
    ActionInterface::keyPressEvent(event);
}

void SpaceColonizationInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    Matrix3<float> voxels = voxelGrid->getVoxelValues();
    Matrix3<int> availableGrid(voxelGrid->sizeX, voxelGrid->sizeY, voxelGrid->sizeZ, 0);
    for (int x = 0; x < availableGrid.sizeX; x++) {
        for (int y = 0; y < availableGrid.sizeY; y++) {
            for (int z = 0; z < availableGrid.sizeZ; z++) {
                float voxelVal = voxels(x, y, z);
                if (voxelVal > 0) {
                    availableGrid.at(x, y, z) = 1;
                }
            }
        }
    }
    FastPoissonGraph<TreeColonisationAlgo::NODE_TYPE> poissonGraph(availableGrid, 20.f);
    int nb_special_nodes = std::min(8, int(poissonGraph.nodes.size()));
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
        this->controlPoints.push_back(std::make_unique<ControlPoint>(keyPoints[i], 5.f));
        this->controlPoints.back()->allowAllAxisTranslation(true);
        QObject::connect(this->controlPoints.back().get(), &ControlPoint::modified,
                         this, &SpaceColonizationInterface::computeKarst);
    }
    this->startingPoint = std::make_unique<ControlPoint>(startPos, 5.f);
    QObject::connect(this->startingPoint.get(), &ControlPoint::modified, this, &SpaceColonizationInterface::computeKarst);

    this->visitingCamera = new VisitingCamera();
}

void SpaceColonizationInterface::initSpaceColonizer()
{
    std::vector<Vector3> newNodes;
    for (auto& ctrl : this->controlPoints)
        newNodes.push_back(ctrl->getPosition());
    this->colonizer->startPosition = this->startingPoint->getPosition();
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
    this->visitingCamera->setZNearCoefficient(0.001);
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
//        this->cameraConstraint = new PathCameraConstraint(this->visitingCamera, simplifiedKarstParths);
//        this->visitingCamera->frame()->setConstraint(this->cameraConstraint);
        this->visitingCamera->setPosition(pathPositions[0]);
        this->visitingCamera->lookAt(pathPositions[1]);
        this->visitingCamera->paths = this->karstPaths;
    }

    Q_EMIT this->karstPathUpdated();
}

void SpaceColonizationInterface::createKarst(bool usingSpheres)
{
    if (this->karstPaths.empty())
        this->updateKarstPath();
    UnderwaterErosion erod(this->voxelGrid, this->karstWidth, 3.f, 10);
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

    std::vector<nlohmann::json> controlPointsPos;
    for (auto& ctrl : this->controlPoints) controlPointsPos.push_back(vec3_to_json(ctrl->getPosition()));
    this->addTerrainAction(nlohmann::json({
                                              {"starting_point", vec3_to_json(startingPoint->getPosition()) },
                                              {"control_points", controlPointsPos},
                                              {"width", karstWidth},
                                              {"randomness", colonizer->randomness},
                                              {"segment_length", colonizer->segmentLength},
                                              {"using_spheres", usingSpheres},
                                              {"min_distance", colonizer->nodeMinDistance},
                                              {"max_distance", colonizer->nodeMaxDistance}
                                          }));
}

QLayout *SpaceColonizationInterface::createGUI()
{
    this->spaceColonizationLayout = new QHBoxLayout;
    QPushButton* spaceColonizerPreviewButton = new QPushButton("Calculer");
    QPushButton* spaceColonizerConfirmButton = new QPushButton("Creer le karst");
    QPushButton* spaceColonizerQuickConfirmButton = new QPushButton("Tunnel rond");
    QCheckBox* spaceColonizerDisplay = new QCheckBox("Afficher");
    QCheckBox* useAsMainCamera = new QCheckBox("Observer l'interieur");
    FancySlider* spaceColonizerSegmentSize = new FancySlider(Qt::Orientation::Horizontal, 1.0, 40.0, 1.0);
    FancySlider* spaceColonizerRandomness = new FancySlider(Qt::Orientation::Horizontal, 0.0, 1.0, 0.1);
    FancySlider* spaceColonizerTunnelWidth = new FancySlider(Qt::Orientation::Horizontal, 0.0, 30.0);
    this->spaceColonizationLayout->addWidget(createVerticalGroup({spaceColonizerPreviewButton, spaceColonizerConfirmButton, spaceColonizerQuickConfirmButton}));
    this->spaceColonizationLayout->addWidget(createVerticalGroup({
                                                           createSliderGroup("Grossier", spaceColonizerSegmentSize),
                                                           createSliderGroup("Tortuosité", spaceColonizerRandomness),
                                                           createSliderGroup("Largeur", spaceColonizerTunnelWidth)
                                                       }));
    this->spaceColonizationLayout->addWidget(createVerticalGroup({/*spaceColonizerDisplay, */useAsMainCamera}));

    spaceColonizerTunnelWidth->setfValue(this->karstWidth);

    QObject::connect(spaceColonizerPreviewButton, &QPushButton::pressed, this, &SpaceColonizationInterface::computeKarst);
    QObject::connect(spaceColonizerConfirmButton, &QPushButton::pressed, [=](){ this->createKarst(false); } );
    QObject::connect(spaceColonizerQuickConfirmButton, &QPushButton::pressed, this, [=](){ this->createKarst(true); } );
    QObject::connect(spaceColonizerRandomness, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->randomness = val; this->computeKarst(); } );
    QObject::connect(spaceColonizerSegmentSize, &FancySlider::floatValueChanged, this, [=](float val){ this->colonizer->segmentLength = val; this->computeKarst(); } );
    QObject::connect(spaceColonizerTunnelWidth, &FancySlider::floatValueChanged, this, [=](float val){ this->karstWidth = val; } );
    QObject::connect(spaceColonizerDisplay, &QCheckBox::toggled, this, [=](bool display){ this->setVisibility(display); } );
    QObject::connect(useAsMainCamera, &QCheckBox::toggled, this, [=](bool display){
        Q_EMIT this->useAsMainCamera(this->visitingCamera, display);
    } );

    spaceColonizerDisplay->setChecked(this->isVisible());
    useAsMainCamera->setChecked(false);

    return this->spaceColonizationLayout;
}
