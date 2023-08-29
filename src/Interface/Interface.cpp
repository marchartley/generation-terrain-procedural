#include "Interface/Interface.h"
#include "Interface/InterfaceUtils.h"

ViewerInterface::ViewerInterface() {
    this->setWindowFlag(Qt::WindowType::WindowMaximizeButtonHint);
    this->setWindowFlag(Qt::WindowType::WindowMinimizeButtonHint);
    this->viewer = new Viewer(this);

    this->actionsOnMap = std::make_shared<std::vector<nlohmann::json>>();
    std::string historyFilename = "MyChanges.json";
    std::shared_ptr<std::fstream> actionsHistoryFile = std::make_shared<std::fstream>(historyFilename);


    std::vector<std::shared_ptr<ActionInterface>> interfacesList = {
        std::make_shared<KarstPathGenerationInterface>(this),
        std::make_shared<SpaceColonizationInterface>(this),
        std::make_shared<FaultSlipInterface>(this),
        std::make_shared<FlowFieldInterface>(this),
        std::make_shared<TunnelInterface>(this),
        std::make_shared<ManualEditionInterface>(this),
        std::make_shared<GravityInterface>(this),
        std::make_shared<UndoRedoInterface>(this),
        std::make_shared<TerrainGenerationInterface>(this),
        std::make_shared<ErosionInterface>(this),
        std::make_shared<HeightmapErosionInterface>(this),
        std::make_shared<BiomeInterface>(this),
        std::make_shared<SmoothInterface>(this),
        std::make_shared<PrimitivePatchesInterface>(this),
        std::make_shared<TerrainSavingInterface>(this),
        std::make_shared<MeshInstanceAmplificationInterface>(this),
        std::make_shared<SPHSimulationInterface>(this),
        std::make_shared<FLIPSimulationInterface>(this),
        std::make_shared<WarpFluidSimulationInterface>(this),
        std::make_shared<LBMFluidSimulationInterface>(this),
        std::make_shared<CoralIslandGeneratorInterface>(this),
        std::make_shared<SpheroidalErosionInterface>(this)
    };

    this->actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>();
    for (auto& interf : interfacesList) {
        interf->viewer = viewer;
        this->actionInterfaces[interf->actionType] = interf;
    }

    std::shared_ptr<TerrainGenerationInterface> terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);

    /*
    this->erosionInterface->viewer = this->viewer;
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
    this->warpFluidSimulationInterface = std::make_shared<WarpFluidSimulationInterface>(this);
    this->LbmFluidSimulationInterface = std::make_shared<LBMFluidSimulationInterface>(this);
    this->coralIslandGeneratorInterface = std::make_shared<CoralIslandGeneratorInterface>(this);
    this->spheroidalErosionInterface = std::make_shared<SpheroidalErosionInterface>(this);

    this->actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>(
                {
                    { "spaceColonization", spaceColonization},
                    { "karstPathGeneration", karstPathGeneration},
                    { "faultSlip", faultSlip},
//                    { "flowField", flowField},
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
//                    { "SPHSimulation", sphSimulationInterface },
//                    { "FLIPSimulation", flipSimulationInterface },
//                    { "WarpFluidSimulation", warpFluidSimulationInterface },
//                    { "LBMFluidSimulation", LbmFluidSimulationInterface },
                    { "coralIslandGeneratorInterface", coralIslandGeneratorInterface},
                    { "spheroidalErosionInterface", spheroidalErosionInterface}
                });
    */
    viewer->interfaces = this->actionInterfaces;
    for (auto& actionInterface : this->actionInterfaces) {
        actionInterface.second->hide();
        actionInterface.second->affectSavingFile(actionsOnMap, actionsHistoryFile, historyFilename);
    }

    QObject::connect(this->viewer, &Viewer::viewerInitialized, this, [=](){
//        this->terrainGenerationInterface->createTerrainFromNoise(3, 3, 2, 1.0, 0.3);

//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/biomes/mayotte.json");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/rock_begin.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/goblin_test.jpg");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/trench.json");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/flat (copy).png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/map1.png");

//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope_original.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/coral_base.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/cube.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/CubeTunnel.data");
        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/one_slope_noisy_reinforced.data");


//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/heightmap.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/river.png");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/voxel_grids/overhang.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/vase.data");
//        this->terrainGenerationInterface->createTerrainFromFile("saved_maps/trench.json");
//        terrainGenerationInterface->createTerrainFromNoise(30, 30, 30, true);

//        this->terrainGenerationInterface->prepareShader();
        this->viewer->voxelGrid = terrainGenerationInterface->voxelGrid;
        this->viewer->heightmap = terrainGenerationInterface->heightmap;
        this->viewer->layerGrid = terrainGenerationInterface->layerGrid;
        this->viewer->implicitTerrain = terrainGenerationInterface->implicitTerrain;

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

        auto biomeInterface = std::static_pointer_cast<BiomeInterface>(actionInterfaces["biomes"]);
        if (terrainGenerationInterface->biomeGenerationNeeded) {
            biomeInterface->biomeModel = BiomeModel::fromJson(terrainGenerationInterface->biomeGenerationModelData);
            biomeInterface->generateBiomes();
        }
        viewer->setSceneCenter(viewer->voxelGrid->getDimensions() / 2.f);

        QObject::connect(biomeInterface.get(), &BiomeInterface::terrainViewModified, terrainGenerationInterface.get(), &TerrainGenerationInterface::updateDisplayedView);
    });

    /*
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::needToClipView,
                     this->viewer, &Viewer::clipViewTemporarily);
    QObject::connect(this->spaceColonization.get(), &SpaceColonizationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->karstPathGeneration.get(), &KarstPathGenerationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::tunnelCreated, this->biomeInterface.get(), &BiomeInterface::addTunnel);
*/


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

    QMenuBar* menu = new QMenuBar(this);
    this->setMenuBar(menu);

    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, QMenu*>> interfacesToOpen;
    std::map<std::string, QMenu*> menus;
    for (auto& [id, interface] : actionInterfaces) {
        if (interface->interfaceType != "") continue;
        if (menus.count(interface->interfaceType) == 0) {
            menus[interface->interfaceType] = new QMenu(QString::fromStdString(toCapitalize(interface->interfaceType)));
            menu->addMenu(menus[interface->interfaceType]);
        }
        interfacesToOpen.push_back({interface, interface->mainActionButtonLogo, interface->mainActionDescription, menus[interface->interfaceType]});
    }
    auto undoRedoInterface = std::static_pointer_cast<UndoRedoInterface>(actionInterfaces["undoredo"]);
    auto smoothInterface = std::static_pointer_cast<SmoothInterface>(actionInterfaces["smooth"]);
    auto terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);

    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, std::string, std::function<void(void)>>> actionsToUse = {
//         Main interface     Button image                        Description                         Menu            Function to call
        {undoRedoInterface,             "undo_button.png",                  "Undo last action",                 "edit",       [=]() { undoRedoInterface->undo(); }},
        {undoRedoInterface,             "redo_button.png",                  "Redo last action",                 "edit",       [=]() { undoRedoInterface->redo(); }},
        {terrainGenerationInterface,    "reload_button.png",                "Reload terrain",                   "file",       [=]() { terrainGenerationInterface->reloadTerrain(this->actionInterfaces); }},
//        {terrainGenerationInterface,    "open_button.png",                  "Open an existing map",             "file",       [=]() { terrainGenerationInterface->openMapUI(); }},
//        {terrainGenerationInterface,    "save_button.png",                  "Save the current map",             "file",       [=]() { terrainGenerationInterface->saveMapUI(); }},
        {nullptr,                       "recording_button.png",             "Start / stop recording the screen","recording",  [=]() { this->viewer->startStopRecording(); }},
        {nullptr,                       "snapshot.png",                     "Take snapshot of screen",          "recording",  [=]() { this->viewer->screenshot(); }},
        {nullptr,                       "display-marching-cubes_button.png","Use the Marching Cubes algorithm", "view",       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES); } },
        {nullptr,                       "display-voxels_button.png",        "Use raw voxels",                   "view",       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE); } },
        {nullptr,                       "heightmap_button.png",             "Display using heightmap structure","model",      [=]() { this->viewer->setMapMode(MapMode::GRID_MODE); } },
        {nullptr,                       "layerbased_button.png",            "Display using layered structure",  "model",      [=]() { this->viewer->setMapMode(MapMode::LAYER_MODE); } },
        {nullptr,                       "voxelgrid_button.png",             "Display using voxels structure",   "model",      [=]() { this->viewer->setMapMode(MapMode::VOXEL_MODE);} },
        {nullptr,                       "implicit_button.png",              "Display using implicit structure", "model",      [=]() { this->viewer->setMapMode(MapMode::IMPLICIT_MODE);} },
        {nullptr,                       "no_display.png",                   "Hide map",                         "view",       [=]() { this->viewer->setViewerMode(NO_DISPLAY);} },
        {nullptr,                       "wire_mode.png",                    "Display wireframe",                "view",       [=]() { this->viewer->setViewerMode(WIRE_MODE); } },
        {nullptr,                       "fill_mode.png",                    "Display solid",                    "view",       [=]() { this->viewer->setViewerMode(FILL_MODE); } },
        {smoothInterface,               "smooth_button.png",                "Smooth the terrain",               "model",      [=]() { smoothInterface->applySmooth(); } },
        {terrainGenerationInterface,    "heightmap_to_all_button.png",      "Use heightmap as reference",       "model",      [=]() { terrainGenerationInterface->heightmapToAll(); } },
        {terrainGenerationInterface,    "voxels_to_all_button.png",         "Use voxels as reference",          "model",      [=]() { terrainGenerationInterface->voxelsToAll(); } },
        {terrainGenerationInterface,    "layers_to_all_button.png",         "Use layers as reference",          "model",      [=]() { terrainGenerationInterface->layersToAll(); } },
        {terrainGenerationInterface,    "implicit_to_all_button.png",       "Use implicit as reference",        "model",      [=]() { terrainGenerationInterface->implicitToAll(); } },
        {terrainGenerationInterface,    "reset_button.png",                 "Reload terrain from file",         "model",      [=]() { terrainGenerationInterface->reloadTerrain(this->actionInterfaces); } },
        {terrainGenerationInterface,    "",                                 "Reinforce borders",                "digging",    [=]() { terrainGenerationInterface->reinforceVoxels(); } },
        {nullptr,                       "",                                 "Erosion tests",                    "edit",       [=]() { this->erosionsTests(); } },
        {nullptr,                       "",                                 "Save depo/erod",                   "model",      [=]() { terrainGenerationInterface->saveErosionDepositionTextureMasksOnMultiple(); }}
    };

    for (auto& [action, logo, descr, menuName, function] : actionsToUse) {
        if (menus.count(menuName) == 0)
            menus[menuName] = new QMenu(QString::fromStdString(toCapitalize(menuName)));
    }

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
        auto& menu = menus[std::get<3>(informations)];
        auto& function = std::get<4>(informations);
        QIcon icon(QString::fromStdString(":/icons/src/assets/" + iconFilename));
        QAction* action = new QAction(icon, QString::fromStdString(description));
        menu->addAction(action);
        QObject::connect(action, &QAction::triggered, this, function);
    }

    for (auto& [menuName, menu] : menus) {
        toolbar->addActions(menu->actions());
        toolbar->addSeparator();
    }

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
        if (active) { terrainGenerationInterface->minIsoLevel = isolevelSelectionSlider->min_value(); terrainGenerationInterface->maxIsoLevel = isolevelSelectionSlider->max_value();}
        else { terrainGenerationInterface->minIsoLevel = -1000.f; terrainGenerationInterface->maxIsoLevel = 1000.f; }
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
    auto terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);
    QObject::connect(mapSliceSliderX, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapX = min; viewer->maxSliceMapX = max; viewer->update(); });
    QObject::connect(mapSliceSliderY, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapY = min; viewer->maxSliceMapY = max; viewer->update(); });
    QObject::connect(mapSliceSliderZ, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapZ = min; viewer->maxSliceMapZ = max; viewer->update(); });

    QObject::connect(isolevelSelectionSlider, &RangeSlider::alt_valueChanged, this, [=](float min, float max){
        terrainGenerationInterface->minIsoLevel = min;
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
        this->frame->clearContent();
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


std::map<std::string, float> generateRandomValuesFrom(std::map<std::string, std::vector<float>> variables, std::vector<std::string> unlockedVariables = {}, std::vector<std::string> lockedVariables = {}, float t = -1)
{
//    std::vector<int> indices;
    /*
    if (lockedVariables.empty() && unlockedVariables.size() > 0) {
        for (const auto& var : variables)
            if (!isIn(var.first, unlockedVariables)) lockedVariables.push_back(var.first);
    }
    else if (unlockedVariables.empty() && lockedVariables.size() > 0) {
        for (const auto& var : variables)
            if (!isIn(var.first, lockedVariables)) unlockedVariables.push_back(var.first);
    }
    else if (unlockedVariables.empty() && lockedVariables.empty()) {
        for (const auto& var : variables)
            unlockedVariables.push_back(var.first);
    }
    std::map<std::string, float> results;
    for (const auto& [var, val] : variables) {
        auto [mini, maxi, defaultVal] = val;
        if (isIn(var, unlockedVariables)) {
            if (t == -1) {
                results[var] = random_gen::generate(mini, maxi); // Forget about the step
            } else {
                results[var] = interpolation::inv_linear(t, mini, maxi);
            }
        } else {
            results[var] = defaultVal;
        }
    }
    return results;
    */
    return {};
}

void ViewerInterface::erosionsTests()
{
    auto erosionInterface = std::static_pointer_cast<ErosionInterface>(actionInterfaces["erosion"]);
    auto terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);
    auto savingInterface = std::static_pointer_cast<TerrainSavingInterface>(actionInterfaces["terrainsaving"]);

    std::map<std::string, std::vector<float>> parameters = {
        {"Strength", {0.1, 0.2, 0.5}},
        {"Qtt", {320, 100, 1000}},
        {"B.coef", {0.01, 0.1, 0.2, 1.0}},
        {"B.ness", {1.0, 0.7, 0.3}},
    };

    srand(1);

//    UnderwaterErosion::EROSION_APPLIED terrainType = this->applyOn; //UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP; // 0 = voxels, 1 = heightmap, 2 = implicit, ...
//    if (terrainType == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS)
//        viewer->setMapMode(MapMode::VOXEL_MODE);
//    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP)
//        viewer->setMapMode(MapMode::GRID_MODE);
//    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN)
//        viewer->setMapMode(MapMode::LAYER_MODE);
//    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN)
//        viewer->setMapMode(MapMode::IMPLICIT_MODE);
//    srand(43);

//    TerrainGenerationInterface* terrainInterface = static_cast<TerrainGenerationInterface*>(this->viewer->interfaces["terrainGenerationInterface"].get());
//    ErosionInterface* erosionInterface = static_cast<ErosionInterface*>(this->viewer->interfaces["erosionInterface"].get());

    VoxelGrid initialVoxelGrid = *(terrainGenerationInterface->voxelGrid);
//    Heightmap initialHeightmap = *heightmap;
//    LayerBasedGrid initialLayerGrid = *layerGrid;

//    this->viewer->restoreFromFile("experiments_state2.xml");
//    this->viewer->setStateFileName("experiments_state2.xml");
//    this->viewer->saveStateToFile();

//    Vector3 cameraPosition = Vector3(this->viewer->camera()->position());
//    Vector3 cameraDirection = Vector3(this->viewer->camera()->viewDirection());

    std::map<std::string, std::vector<float>> varyingVariables =
    {
        // Name,            min, max, step
        {"bouncingCoefficient",         {/*0.1, */0.2, 1.0} },
        {"bounciness",                  {1.0, 0.7, 0.3} },
        {"maxCapacityFactor",           {1.0, 0.1, 10.0} },
        {"erosionFactor",               {0.5, 1.0, 5.0} },
        {"depositFactor",               {0.5, 1.0, 5.0} },
        {"matterDensity",               {500.0, 1662.0} },
        {"airForce",                    {0.f, 0.1} },
        {"waterForce",                  {/*0.f, */0.1} },
        {"particleSize",                {8.0, 12.0/*, 4.0*/} },
        {"astrength",                   {0.2, 0.5} },
        {"nbParticles",                 {320} },
        {"waterLevel",                  {0.0, 0.25} },
        {"from",                        {0, 1} },
        {"resistance",                  {0, 1} }
    };

    int nbTests = 1;
    std::vector<std::string> params;
    for (const auto& var : varyingVariables) {
        params.push_back(var.first);
        nbTests *= var.second.size();
    }

    std::vector<int> paramsCurrentIndices(varyingVariables.size());
    paramsCurrentIndices[0] = -1;

    std::string mainFolder = "";
#ifdef linux
    mainFolder = "/data/erosionsTests/";
#else
    mainFolder = "erosionsTests/";
#endif
    makedir(mainFolder);
    // Create the general CSV
    // Check if already exists
    std::string CSVname = mainFolder + "allData.csv";
    bool exists = checkPathExists(CSVname);
    std::fstream mainCSVfile;
    std::vector<std::map<std::string, float>> alreadyTestedParameters;
    if (exists) {
        mainCSVfile.open(CSVname, std::ios_base::in);
        mainCSVfile.seekg(0);
        mainCSVfile.seekp(0);
        std::string header, lineContent;
        std::vector<std::string> headerValues, sLineValues;
        std::getline(mainCSVfile, header); // get header
        headerValues = split(header, ";");
        while (std::getline(mainCSVfile, lineContent)) {
            sLineValues = split(lineContent, ";");
            std::map<std::string, float> testedParameters;
            for (size_t i = 0; i < sLineValues.size(); i++) {
                float val;
                try {
                    val = std::stof(sLineValues[i]);
                    testedParameters[headerValues[i]] = val;
                } catch (std::exception e) {
                    try {
                        val = std::stof(replaceInString(sLineValues[i], ",", "."));
                        testedParameters[headerValues[i]] = val;
                    } catch (std::exception e2) {

                    }
                 }
            }
            alreadyTestedParameters.push_back(testedParameters);
        }
        mainCSVfile.close();
        mainCSVfile.open(CSVname, std::ios_base::out | std::ios_base::app);

    } else {
        mainCSVfile.open(CSVname, std::ios_base::out | std::ios_base::app);
        for (const auto& param : params)
            mainCSVfile << param << ";";
        mainCSVfile << "folder_name" << std::endl;
        mainCSVfile.close();
    }

    std::vector<std::vector<std::string>> testedVariables = {
//        { "particleSize" },
//        { "strength" },
//        { "erosionFactor" },
//        { "criticalShearStress" },
//        { "shearingStressConstantK" },
//        { "nbParticles" },
//        { "maxCapacityFactor" },
//        { "matterDensity" },
//        { "depositFactor" },
//        { "waterLevel" },
//        { "waterForce" },
//        { "matterDensity" },
    };

    savingInterface->saveVoxels = true;
    savingInterface->saveHeightmap = false;
    savingInterface->saveLayers = false;

    int filenameResultIndex = int(alreadyTestedParameters.size());

    int stopAfterStep = -1; // Set to -1 for infinite tries
    for (int iCombination = 0; iCombination < nbTests; iCombination++) {
        mainCSVfile.open(CSVname, std::ios_base::out | std::ios_base::app);
        bool theresARest = true;
        int indexToIncrement = 0;
//        while (theresARest) {
//            theresARest = false;
//            paramsCurrentIndices[indexToIncrement]++;
//            if (paramsCurrentIndices[indexToIncrement] >= int(varyingVariables[params[indexToIncrement]].size())) {
//                theresARest = true;
//                paramsCurrentIndices[indexToIncrement] = 0;
//                indexToIncrement++;
//            }
//            if (indexToIncrement >= int(params.size()))
//                return;
//        }
        while (theresARest) {
            if (paramsCurrentIndices[indexToIncrement] == int(varyingVariables[params[indexToIncrement]].size()) - 1) {
                theresARest = true;
                paramsCurrentIndices[indexToIncrement] = 0;
                indexToIncrement++;
            } else {
                theresARest = false;
                paramsCurrentIndices[indexToIncrement] = int(varyingVariables[params[indexToIncrement]].size()) - 1;
                indexToIncrement++;
            }
        }

        std::map<std::string, float> variables;
        for (int i = 0; i < params.size(); i++) {
            variables[params[i]] = varyingVariables[params[i]][paramsCurrentIndices[i]];
        }
//        std::cout << "Testing " << join(testedVariables[iCombination], " and ") << std::endl;
//        for (int i = 0; i < stopAfterStep || stopAfterStep == -1; i++) {
            // Initialize
//            float t = (testedVariables[iCombination].size() == 1 ? float(i)/float(stopAfterStep - 1) : -1);
//            auto variables = generateRandomValuesFrom(varyingVariables, testedVariables[iCombination], {}, t);
            bool tested = false;
            for (auto& x : alreadyTestedParameters) {
                tested = true;
                for (auto& [paramName, paramVal] : variables) {
                    if (!startsWith(paramName, "cam") && replaceInString(std::to_string(variables[paramName]), ",", ".") != replaceInString(std::to_string(x[paramName]), ",", ".")) {
                        tested = false;
                    }
                }
                if (tested)
                    break;
            }
            if (tested) {
                continue;
            }
            alreadyTestedParameters.push_back(variables);

            *terrainGenerationInterface->voxelGrid = initialVoxelGrid;

            terrainGenerationInterface->setWaterLevel(variables["waterLevel"]);

            std::cout << "Test " << filenameResultIndex << ": \n";
            for (const auto& var : variables)
                std::cout << "\t- \"" << var.first << "\" = " << var.second << "\n";
            std::cout << std::flush;

            /*

        {"bouncingCoefficient",         {0.01, 0.1, 0.2, 1.0} },
        {"bounciness",                  {1.0, 0.7, 0.3} },
        {"maxCapacityFactor",           {1.0, 0.1, 3.0, 10.0} },
        {"erosionFactor",               {0.5, 1.0, 5.0} },
        {"depositFactor",               {0.5, 1.0, 5.0} },
        {"matterDensity",               {500.0, 1662.0} },
        {"airForce",                    {0.f, 0.03, 0.1, 1.0} },
        {"waterForce",                  {0.f, 0.03, 0.1, 1.0} },
        {"particleSize",                {4.0, 8.0, 12.0} },
        {"strength",                    {0.01, 0.1, 0.3, 0.5} },
        {"nbParticles",                 {320} },
        {"waterLevel",                  {0.0, 0.25, 0.5} },
        {"from",                        {0, 1} },
        {"resistance",                  {0, 1} }
             */
            erosionInterface->bouncingCoefficient = variables["bouncingCoefficient"];
            erosionInterface->bounciness = variables["bounciness"];
            erosionInterface->maxCapacityFactor = variables["maxCapacityFactor"];
            erosionInterface->erosionFactor = variables["erosionFactor"];
            erosionInterface->depositFactor = variables["depositFactor"];
            erosionInterface->matterDensity = variables["matterDensity"];
            erosionInterface->airForce = variables["airForce"];
            erosionInterface->waterForce = variables["waterForce"];
            erosionInterface->erosionSize = variables["particleSize"];
            erosionInterface->erosionStrength = variables["astrength"];
            erosionInterface->erosionQtt = variables["nbParticles"];

            ErosionInterface::PARTICLE_INITIAL_LOCATION loc;
            if (variables["from"] == 0) {
                loc = ErosionInterface::PARTICLE_INITIAL_LOCATION::JUST_ABOVE_VOXELS;
            } else {
                loc = ErosionInterface::PARTICLE_INITIAL_LOCATION::RANDOM;
            }

            UnderwaterErosion::DENSITY_TYPE densType;
            if (variables["resistance"] == 0) {
                densType = UnderwaterErosion::DENSITY_TYPE::NATIVE;
            } else {
                densType = UnderwaterErosion::DENSITY_TYPE::RANDOM_DENSITY;
            }

            erosionInterface->densityUsed = densType;
            erosionInterface->numberOfIterations = 50;
            erosionInterface->throwFrom(loc);

//            this->savingInterface->mainFilename = "erosionsTests/param_";
            std::string subFolderName = savingInterface->saveTerrainGeometry("erosionsTests/param_" + std::to_string(filenameResultIndex))[0];
            filenameResultIndex++;

            std::string parametersToCSV = "";
            for (const auto& param : params)
                parametersToCSV += std::to_string(variables[param]) + ";";
            parametersToCSV += subFolderName;
            mainCSVfile << parametersToCSV << std::endl;

//        }
            mainCSVfile.close();
    }
    mainCSVfile.close();

    std::cout << "Unbelivable but true, it's finished!" << std::endl;
}
