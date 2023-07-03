#include "Interface/Interface.h"
#include "Interface/InterfaceUtils.h"

ViewerInterface::ViewerInterface() {
    this->setWindowFlag(Qt::WindowType::WindowMaximizeButtonHint);
    this->setWindowFlag(Qt::WindowType::WindowMinimizeButtonHint);
    this->viewer = new Viewer(this);

    this->actionsOnMap = std::make_shared<std::vector<nlohmann::json>>();
    std::string historyFilename = "MyChanges.json";
    std::shared_ptr<std::fstream> actionsHistoryFile = std::make_shared<std::fstream>(historyFilename);


    this->karstPathGeneration = std::make_shared<KarstPathGenerationInterface>(this);
    this->spaceColonization = std::make_shared<SpaceColonizationInterface>(this);
    this->faultSlip = std::make_shared<FaultSlipInterface>(this);
    this->flowField = std::make_shared<FlowFieldInterface>(this);
    this->tunnelInterface = std::make_shared<TunnelInterface>(this);
    this->manualEditionInterface = std::make_shared<ManualEditionInterface>(this);
    this->gravityInterface = std::make_shared<GravityInterface>(this);
    this->undoRedoInterface = std::make_shared<UndoRedoInterface>(this);
    this->terrainGenerationInterface = std::make_shared<TerrainGenerationInterface>(this);
    this->erosionInterface = std::make_shared<ErosionInterface>(this);
    this->erosionInterface->viewer = this->viewer;
    this->heightmapErosionInterface = std::make_shared<HeightmapErosionInterface>(this);
    this->biomeInterface = std::make_shared<BiomeInterface>(this);
    this->smoothInterface = std::make_shared<SmoothInterface>(this);
    this->patchesInterface = std::make_shared<PrimitivePatchesInterface>(this);
    this->savingInterface = std::make_shared<TerrainSavingInterface>(this);
    this->meshInstanceAmplificationInterface = std::make_shared<MeshInstanceAmplificationInterface>(this);
    this->sphSimulationInterface = std::make_shared<SPHSimulationInterface>(this);
    this->flipSimulationInterface = std::make_shared<FLIPSimulationInterface>(this);
    this->coralIslandGeneratorInterface = std::make_shared<CoralIslandGeneratorInterface>(this);

    this->actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>(
                {
                    { "spaceColonization", spaceColonization},
                    { "karstPathGeneration", karstPathGeneration},
                    { "faultSlip", faultSlip},
                    { "flowField", flowField},
                    { "tunnelInterface", tunnelInterface},
                    { "manualEditionInterface", manualEditionInterface},
                    { "gravityInterface", gravityInterface},
                    { "undoRedoInterface", undoRedoInterface},
                    { "terrainGenerationInterface", terrainGenerationInterface},
                    { "erosionInterface", erosionInterface},
                    { "heightmapErosionInterface", heightmapErosionInterface},
                    { "biomeInterface", biomeInterface},
                    { "smoothInterface", smoothInterface},
                    { "terrainSavingInterface", savingInterface},
                    { "meshInstanceAmplificationInterface", meshInstanceAmplificationInterface},
                    { "primitivePatchInterface", patchesInterface},
                    { "SPHSimulation", sphSimulationInterface },
                    { "FLIPSimulation", flipSimulationInterface },
                    { "coralIslandGeneratorInterface", coralIslandGeneratorInterface}
                });
    viewer->interfaces = this->actionInterfaces;
    for (auto& actionInterface : this->actionInterfaces) {
        actionInterface.second->hide();
        actionInterface.second->affectSavingFile(actionsOnMap, actionsHistoryFile, historyFilename);
    }

//    QObject::connect(terrainGenerationInterface.get(), &TerrainGenerationInterface::updated, erosionInterface.get(), &ErosionInterface::computePredefinedRocksLocations);
//    QObject::connect(terrainGenerationInterface.get(), &TerrainGenerationInterface::updated, sphSimulationInterface.get(), &SPHSimulationInterface::updateParticles);

    QObject::connect(this->viewer, &Viewer::viewerInitialized, this, [&](){
//        this->terrainGenerationInterface->createTerrainFromNoise(3, 3, 2, 1.0, 0.3);

//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/biomes/mayotte.json");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/rock_begin.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/coral_long.png");
        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/map1.png");

//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope_original.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/coral_base.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxel_grids/midCave");

//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/heightmap.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/river.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxel_grids/overhang.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/vase.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/trench.json");

//        this->terrainGenerationInterface->prepareShader();
        this->viewer->voxelGrid = this->terrainGenerationInterface->voxelGrid;
        this->viewer->heightmap = this->terrainGenerationInterface->heightmap;
        this->viewer->layerGrid = this->terrainGenerationInterface->layerGrid;
        this->viewer->implicitTerrain = this->terrainGenerationInterface->implicitTerrain;

        for (auto& actionInterface : this->actionInterfaces) {
            actionInterface.second->affectTerrains(viewer->heightmap, viewer->voxelGrid, viewer->layerGrid, viewer->implicitTerrain);
            actionInterface.second->reloadShaders();

            QObject::connect(actionInterface.second.get(), &ActionInterface::updated, this, [&](){
                this->viewer->update();
                qApp->processEvents();
            });
            QObject::connect(this->viewer, &Viewer::mouseClickOnMap, actionInterface.second.get(), &ActionInterface::mouseClickedOnMapEvent);
            QObject::connect(this->viewer, &Viewer::mouseMovedOnMap, actionInterface.second.get(), &ActionInterface::mouseMovedOnMapEvent);
            QObject::connect(this->viewer, &Viewer::mouseDoubleClickedOnMap, actionInterface.second.get(), &ActionInterface::mouseDoubleClickedOnMapEvent);

            for (auto& otherActionInterface : this->actionInterfaces) {
                QObject::connect(actionInterface.second.get(), &ActionInterface::terrainUpdated, this, [&]() {
                    QObject::blockSignals(true);
                    otherActionInterface.second->afterTerrainUpdated();
                    QObject::blockSignals(false);
                    this->viewer->update();
                });
            }
        }

        if (terrainGenerationInterface->biomeGenerationNeeded) {
            this->biomeInterface->biomeModel = BiomeModel::fromJson(terrainGenerationInterface->biomeGenerationModelData);
            this->biomeInterface->generateBiomes();
        }
        viewer->setSceneCenter(viewer->voxelGrid->getDimensions() / 2.f);

        QObject::connect(this->biomeInterface.get(), &BiomeInterface::terrainViewModified, this->terrainGenerationInterface.get(), &TerrainGenerationInterface::updateDisplayedView);
    });

    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::needToClipView,
                     this->viewer, &Viewer::clipViewTemporarily);
    QObject::connect(this->spaceColonization.get(), &SpaceColonizationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->karstPathGeneration.get(), &KarstPathGenerationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::tunnelCreated, this->biomeInterface.get(), &BiomeInterface::addTunnel);

    QObject::connect(qApp, &QApplication::focusChanged, this, [=](QWidget*, QWidget*) {
        this->setFocus(Qt::OtherFocusReason);
//        viewer->setFocus(Qt::OtherFocusReason);
    });
    Plotter::init(nullptr, this);
//    this->installEventFilter(viewer);
    viewer->installEventFilter(this);
    setupUi();
}

ViewerInterface::~ViewerInterface()
{
    for (auto& action : actionInterfaces)
        action.second->setParent(nullptr);
}

void ViewerInterface::setupUi()
{
    QToolBar* toolbar = new QToolBar("Main tools");
    this->addToolBar(Qt::ToolBarArea::TopToolBarArea, toolbar);

    QMenu* fileMenu = new QMenu("File");
    QMenu* editMenu = new QMenu("Edit");
    QMenu* viewMenu = new QMenu("Affichage");
    QMenu* modelMenu = new QMenu("Model");
    QMenu* recordingMenu = new QMenu("Enregistrements");
    QMenu* physicsMenu = new QMenu("Physiques");
    QMenu* diggingMenu = new QMenu("Creuser");
    QMenuBar* menu = new QMenuBar(this);
    menu->addMenu(fileMenu);
    menu->addMenu(editMenu);
    menu->addMenu(viewMenu);
    menu->addMenu(recordingMenu);
    menu->addMenu(physicsMenu);
    menu->addMenu(diggingMenu);
    menu->addMenu(modelMenu);
    this->setMenuBar(menu);

    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, QMenu*>> interfacesToOpen = {
//         Interface to open                      Button image                    Description                             Menu
//     Called with "createGUI()"
        {terrainGenerationInterface,            "terrain_button.png",           "Terrain generation",                   modelMenu},
        {faultSlip,                             "fault-slip_button.png",        "Fault slips",                          diggingMenu},
        {flowField,                             "flowfield_button.png",         "Currents simulation",                  physicsMenu},
        {tunnelInterface,                       "tunnel_button.png",            "Tunnels creation",                     diggingMenu},
        {manualEditionInterface,                "manual-edit_button.png",       "Manual editing",                       diggingMenu},
        {gravityInterface,                      "gravity_button.png",           "Gravity",                              physicsMenu},
        {heightmapErosionInterface,             "heightmap-erosion.png",        "Erosion on heightmap",                 physicsMenu},
        {biomeInterface,                        "biomes.png",                   "Biotopes system",                      modelMenu},
        {patchesInterface,                      "feature_primitive_button.png", "Feature primitives",                   diggingMenu},
        {savingInterface,                       "save_geometry.png",            "Save the geometry",                    viewMenu},
        {erosionInterface,                      "erosion.png",                  "Erosion on 3D terrain",                physicsMenu},
        {meshInstanceAmplificationInterface,    "amplification_instances.png",  "Amplify the terrain with meshes",      viewMenu},
        {spaceColonization,                     "karst_button.png",             "Create karsts using space colonization",diggingMenu},
        {karstPathGeneration,                   "karst_peytavie_button.png",    "Create karsts with graphs",            diggingMenu},
        {sphSimulationInterface,                "",                             "SPH Simulation",                       physicsMenu},
        {flipSimulationInterface,               "",                             "FLIP Simulation",                      physicsMenu},
        {coralIslandGeneratorInterface,         "",                             "Create coral island",                  diggingMenu}
    };
    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, QMenu* , std::function<void(void)>>> actionsToUse = {
//         Main interface     Button image                        Description                         Menu            Function to call
        {undoRedoInterface,             "undo_button.png",                  "Undo last action",                 editMenu,       [=]() { undoRedoInterface->undo(); }},
        {undoRedoInterface,             "redo_button.png",                  "Redo last action",                 editMenu,       [=]() { undoRedoInterface->redo(); }},
//        {terrainGenerationInterface,    "reload_button.png",                "Reload terrain",                   fileMenu,       [=]() { terrainGenerationInterface->reloadTerrain(this->actionInterfaces); }},
//        {terrainGenerationInterface,    "open_button.png",                  "Open an existing map",             fileMenu,       [=]() { terrainGenerationInterface->openMapUI(); }},
//        {terrainGenerationInterface,    "save_button.png",                  "Save the current map",             fileMenu,       [=]() { terrainGenerationInterface->saveMapUI(); }},
        {nullptr,                       "recording_button.png",             "Start / stop recording the screen",recordingMenu,  [=]() { this->viewer->startStopRecording(); }},
        {nullptr,                       "snapshot.png",                     "Take snapshot of screen",          recordingMenu,  [=]() { this->viewer->screenshot(); }},
        {nullptr,                       "display-marching-cubes_button.png","Use the Marching Cubes algorithm", viewMenu,       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES); } },
        {nullptr,                       "display-voxels_button.png",        "Use raw voxels",                   viewMenu,       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE); } },
        {nullptr,                       "heightmap_button.png",             "Display using heightmap structure",modelMenu,      [=]() { this->viewer->setMapMode(MapMode::GRID_MODE); } },
        {nullptr,                       "layerbased_button.png",            "Display using layered structure",  modelMenu,      [=]() { this->viewer->setMapMode(MapMode::LAYER_MODE); } },
        {nullptr,                       "voxelgrid_button.png",             "Display using voxels structure",   modelMenu,      [=]() { this->viewer->setMapMode(MapMode::VOXEL_MODE);} },
        {nullptr,                       "implicit_button.png",              "Display using implicit structure", modelMenu,      [=]() { this->viewer->setMapMode(MapMode::IMPLICIT_MODE);} },
        {nullptr,                       "no_display.png",                   "Hide map",                         viewMenu,       [=]() { this->viewer->setViewerMode(NO_DISPLAY);} },
        {nullptr,                       "wire_mode.png",                    "Display wireframe",                viewMenu,       [=]() { this->viewer->setViewerMode(WIRE_MODE); } },
        {nullptr,                       "fill_mode.png",                    "Display solid",                    viewMenu,       [=]() { this->viewer->setViewerMode(FILL_MODE); } },
        {smoothInterface,               "smooth_button.png",                "Smooth the terrain",               modelMenu,      [=]() { smoothInterface->applySmooth(); } },
        {terrainGenerationInterface,    "",                                 "Use heightmap as reference",       modelMenu,      [=]() { terrainGenerationInterface->heightmapToAll(); } },
        {terrainGenerationInterface,    "",                                 "Use voxels as reference",          modelMenu,      [=]() { terrainGenerationInterface->voxelsToAll(); } },
        {terrainGenerationInterface,    "",                                 "Use layers as reference",          modelMenu,      [=]() { terrainGenerationInterface->layersToAll(); } },
        {terrainGenerationInterface,    "",                                 "Use implicit as reference",        modelMenu,      [=]() { terrainGenerationInterface->implicitToAll(); } },
    };

    for (auto& informations : interfacesToOpen) {
        auto& interfaceObject = std::get<0>(informations);
        auto& iconFilename = std::get<1>(informations);
        auto& description = std::get<2>(informations);
        auto& menu = std::get<3>(informations);
        QIcon icon(QString::fromStdString(":/icons/src/assets/" + iconFilename));
        QAction* action = new QAction(icon, QString::fromStdString(description));
        menu->addAction(action);
        QObject::connect(action, &QAction::triggered, this, [=]() { this->openInterface(interfaceObject); });
    }
    for (auto& informations : actionsToUse) {
        auto& iconFilename = std::get<1>(informations);
        auto& description = std::get<2>(informations);
        auto& menu = std::get<3>(informations);
        auto& function = std::get<4>(informations);
        QIcon icon(QString::fromStdString(":/icons/src/assets/" + iconFilename));
        QAction* action = new QAction(icon, QString::fromStdString(description));
        menu->addAction(action);
        QObject::connect(action, &QAction::triggered, this, function);
    }

    toolbar->addActions(fileMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(editMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(viewMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(recordingMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(physicsMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(diggingMenu->actions());
    toolbar->addSeparator();
    toolbar->addActions(modelMenu->actions());

    QStatusBar* status = new QStatusBar(this);
    this->setStatusBar(status);

    QDockWidget* displayOptionWidget = new QDockWidget("Affichage", this);
    displayOptionWidget->setFeatures(QDockWidget::DockWidgetFloatable);
    this->addDockWidget(Qt::DockWidgetArea::RightDockWidgetArea, displayOptionWidget);

    this->viewer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    mainLayout = new QGridLayout;
    mainLayout->setColumnStretch(0, 99);
    mainLayout->setColumnStretch(1,  1);

    this->displayModeLayout = new QHBoxLayout;
    this->displayModeBox = new Spoiler("Affichage");

//    this->LoDChooserLayout = new QHBoxLayout;
    this->LoDChooserBox = new Spoiler("Niveau de détail");
    this->mapSliceSliderX = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderY = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderZ = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSmooth = new QCheckBox("Shrink on borders");
    mapSliceSmooth->setChecked(viewer->voxelsSmoothedOnBorders > 1);
    QObject::connect(mapSliceSmooth, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->voxelsSmoothedOnBorders = 5;}
        else { this->viewer->voxelsSmoothedOnBorders = 1; }
        this->viewer->update();
    });
    QCheckBox* sliderXactivation = new QCheckBox("Activer");
    sliderXactivation->setChecked(true);
    QObject::connect(sliderXactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapX = this->mapSliceSliderX->min_value(); this->viewer->maxSliceMapX = this->mapSliceSliderX->max_value(); }
        else        { this->viewer->minSliceMapX = -10.f;                              this->viewer->maxSliceMapX = 10.f; }
        this->viewer->update();
    });
    QCheckBox* sliderYactivation = new QCheckBox("Activer");
    sliderYactivation->setChecked(true);
    QObject::connect(sliderYactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapY = this->mapSliceSliderY->min_value(); this->viewer->maxSliceMapY = this->mapSliceSliderY->max_value(); }
        else        { this->viewer->minSliceMapY = -10.f;                              this->viewer->maxSliceMapY = 10.f; }
        this->viewer->update();
    });
    QCheckBox* sliderZactivation = new QCheckBox("Activer");
    sliderZactivation->setChecked(true);
    QObject::connect(sliderZactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapZ = this->mapSliceSliderZ->min_value(); this->viewer->maxSliceMapZ = this->mapSliceSliderZ->max_value(); }
        else        { this->viewer->minSliceMapZ = -10.f;                              this->viewer->maxSliceMapZ = 10.f; }
        this->viewer->update();
    });
    displayModeLayout->addWidget(createVerticalGroup({createMultipleSliderGroupWithCheckbox({
                                                               {"X", mapSliceSliderX, sliderXactivation},
                                                               {"Y", mapSliceSliderY, sliderYactivation},
                                                               {"Z", mapSliceSliderZ, sliderZactivation}
                                                           }),
                                                      mapSliceSmooth}));

    this->isolevelSelectionSlider = new RangeSlider(Qt::Orientation::Vertical, 0.f, 3.f, 0.1f);
    QCheckBox* isolevelSelectionActivation = new QCheckBox("Activer");
//    LoDChooserLayout->addWidget(createMultipleSliderGroupWithCheckbox({
//                                                                          {"Densité", {isolevelSelectionSlider, isolevelSelectionActivation}}
//                                                                      }));
    isolevelSelectionActivation->setChecked(true);
    QObject::connect(isolevelSelectionActivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->terrainGenerationInterface->minIsoLevel = isolevelSelectionSlider->min_value(); this->terrainGenerationInterface->maxIsoLevel = isolevelSelectionSlider->max_value();}
        else { this->terrainGenerationInterface->minIsoLevel = -1000.f; this->terrainGenerationInterface->maxIsoLevel = 1000.f; }
        this->viewer->update();
    });
    FancySlider* waterLevelSlider = new FancySlider(Qt::Orientation::Vertical, 0.f, 1.f, 0.01f);
    QObject::connect(waterLevelSlider, &FancySlider::floatValueChanged, terrainGenerationInterface.get(), &TerrainGenerationInterface::setWaterLevel);

    mainLayout->addWidget(viewer, 1, 0);

    QPushButton* reloadShadersButton = new QPushButton("Recharger tous les shaders");
    QVBoxLayout* displayOptionLayout = new QVBoxLayout();

    for (int i = 0; i < 4; i++) {
        QPushButton* viewerSetupExperimental = new QPushButton("Setup view " + QString::fromStdString(std::to_string(i)));
        QPushButton* viewerSetupExperimental_save = new QPushButton("Save");
        std::string associatedFilename = "experiments_state" + std::to_string(i) + ".xml";
        QObject::connect(viewerSetupExperimental, &QPushButton::pressed, this, [=]() { this->viewer->setupViewFromFile(associatedFilename); });
        QObject::connect(viewerSetupExperimental_save, &QPushButton::pressed, this, [=]() { this->viewer->saveViewToFile(associatedFilename); });

        displayOptionLayout->addWidget(createHorizontalGroup({viewerSetupExperimental, viewerSetupExperimental_save}));
    }
    displayOptionLayout->addItem(displayModeLayout);
    displayOptionLayout->addWidget(createVerticalGroup({
                                                           createHorizontalGroup({
                                                               createMultipleSliderGroupWithCheckbox({
                                                                   {"Density", isolevelSelectionSlider, isolevelSelectionActivation}
                                                               }),
                                                               createSliderGroup("Water", waterLevelSlider)
                                                           }),
                                                           reloadShadersButton
                                                 }));
    QGroupBox* displayOptionBox = new QGroupBox();
    displayOptionBox->setLayout(displayOptionLayout);
    displayOptionWidget->setWidget(displayOptionBox);


    QWidget* mainFrame = new QWidget(this);
    mainFrame->setLayout(mainLayout);
    this->setCentralWidget(mainFrame);


    frame = new StickyFrame(this->viewer, 0, -1, -1, 1, false);
    frame->hide();

    QObject::connect(reloadShadersButton, &QPushButton::pressed, this, [&]() {
        this->viewer->reloadAllShaders();
    });

    this->setAllValuesToFitViewerDefaults(this->viewer);
    this->setupBindings();
    this->retranslateUi();
} // setupUi

void ViewerInterface::setupBindings()
{
    QObject::connect(mapSliceSliderX, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapX = min; viewer->maxSliceMapX = max; viewer->update(); });
    QObject::connect(mapSliceSliderY, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapY = min; viewer->maxSliceMapY = max; viewer->update(); });
    QObject::connect(mapSliceSliderZ, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapZ = min; viewer->maxSliceMapZ = max; viewer->update(); });

    QObject::connect(isolevelSelectionSlider, &RangeSlider::alt_valueChanged, this, [=](float min, float max){
        this->terrainGenerationInterface->minIsoLevel = min;
        terrainGenerationInterface->maxIsoLevel = max;
        viewer->update();
    });

    QMetaObject::connectSlotsByName(this);
} //setupBindings

void ViewerInterface::retranslateUi()
{
    this->setWindowTitle(QString("Simulation"));
} // retranslateUi

void ViewerInterface::setAllValuesToFitViewerDefaults(Viewer* viewer)
{

}

void ViewerInterface::closeEvent(QCloseEvent *e)
{
    viewer->closeEvent(e);
}

void ViewerInterface::resizeEvent(QResizeEvent *e)
{
    QMainWindow::resizeEvent(e);
}

void ViewerInterface::focusOutEvent(QFocusEvent *e)
{
    QMainWindow::focusOutEvent(e);
}


void ViewerInterface::openInterface(std::shared_ptr<ActionInterface> object)
{
    this->hideAllInteractiveParts();
    if (lastOpenedInterface && lastOpenedInterface == object) {
        lastOpenedInterface = nullptr;
        this->frame->hide();
    } else {
        lastOpenedInterface = object;
        object->show();
        this->frame->setContent(object->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::hideAllInteractiveParts()
{
    for (auto& actionInterface : this->actionInterfaces)
        actionInterface.second->hide();

    this->viewer->update();
    this->update();
}

//void ViewerInterface::openMapUI()
//{
//    QString q_filename = QFileDialog::getOpenFileName(this, QString("Ouvrir une carte"), QString::fromStdString(this->mapSavingFolder));
//    terrainGenerationInterface->createTerrainFromFile(q_filename.toStdString(), this->actionInterfaces);
//    viewer->setSceneCenter(viewer->voxelGrid->getDimensions() / 2.f);
////    this->terrainGenerationInterface->prepareShader(true);
//}

//void ViewerInterface::saveMapUI()
//{
//    QString q_filename = QFileDialog::getSaveFileName(this, QString("Enregistrer la carte"), QString::fromStdString(this->mapSavingFolder));
//    terrainGenerationInterface->saveTerrain(q_filename.toStdString());
//}

