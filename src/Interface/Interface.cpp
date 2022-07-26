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
                                                                                  { "biomeInterface", biomeInterface}
                                                                              });
    viewer->interfaces = this->actionInterfaces;
    for (auto& actionInterface : this->actionInterfaces) {
        actionInterface.second->hide();
        actionInterface.second->affectSavingFile(actionsOnMap, actionsHistoryFile, historyFilename);
        QObject::connect(actionInterface.second.get(), &ActionInterface::updated, this, [&]() { this->viewer->update(); });
//        this->viewer->installEventFilter(actionInterface.second.get());
    }

    QObject::connect(this->viewer, &Viewer::viewerInitialized, this, [&](){
//        this->terrainGenerationInterface->createTerrainFromNoise(3, 3, 2, 1.0, 0.3);
#ifdef linux
        this->terrainGenerationInterface->createTerrainFromFile("/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/saved_maps/biomes/mayotte.json");
//        this->terrainGenerationInterface->createTerrainFromFile("/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/saved_maps/heightmaps/one_slope.png");
#else
        this->terrainGenerationInterface->createTerrainFromFile("C:/codes/Qt/generation-terrain-procedural/saved_maps/biomes/mayotte.json");
//        this->terrainGenerationInterface->createTerrainFromFile("C:/codes/Qt/generation-terrain-procedural/saved_maps/heightmaps/one_slope.png");
#endif
        this->terrainGenerationInterface->prepareShader();
        this->viewer->voxelGrid = this->terrainGenerationInterface->voxelGrid;
        this->viewer->grid = this->terrainGenerationInterface->heightmapGrid;

        for (auto& actionInterface : this->actionInterfaces) {
            actionInterface.second->affectVoxelGrid(this->viewer->voxelGrid);
            actionInterface.second->affectHeightmap(this->viewer->grid);
        }

        this->heightmapErosionInterface->heightmap = this->terrainGenerationInterface->heightmapGrid;
//        this->biomeInterface->affectHeightmap(this->terrainGenerationInterface->heightmapGrid);
        if (terrainGenerationInterface->biomeGenerationNeeded) {
            this->biomeInterface->biomeModel = BiomeModel::fromJson(terrainGenerationInterface->biomeGenerationModelData);
            this->biomeInterface->generateBiomes(); // std::async([&]() -> void {this->biomeInterface->generateBiomes(); });
        }
        viewer->setSceneCenter(viewer->voxelGrid->getDimensions() * viewer->voxelGrid->getBlockSize() / 2.f);

        QObject::connect(this->biomeInterface.get(), &BiomeInterface::terrainViewModified, this->terrainGenerationInterface.get(), &TerrainGenerationInterface::updateDisplayedView);
    });

    QObject::connect(this->karstPathGeneration.get(), &KarstPathGenerationInterface::karstPathUpdated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->spaceColonization.get(), &SpaceColonizationInterface::karstPathUpdated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->faultSlip.get(), &FaultSlipInterface::faultSlipApplied,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->viewer, &Viewer::mouseClickOnMap,
                     this->tunnelInterface.get(), &TunnelInterface::mouseClickInWorldEvent);
    QObject::connect(this->viewer, &Viewer::mouseClickOnMap,
                     this->manualEditionInterface.get(), &ManualEditionInterface::mouseClickedOnMapEvent);
    QObject::connect(this->viewer, &Viewer::mouseMovedOnMap,
                     this->manualEditionInterface.get(), &ManualEditionInterface::setPosition);
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::needToClipView,
                     this->viewer, &Viewer::clipViewTemporarily);
    QObject::connect(this->viewer, &Viewer::mouseClickOnMap,
                     this->biomeInterface.get(), &BiomeInterface::mouseClickedOnMapEvent);
    QObject::connect(this->spaceColonization.get(), &SpaceColonizationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->karstPathGeneration.get(), &KarstPathGenerationInterface::useAsMainCamera, this->viewer, &Viewer::swapCamera);
    QObject::connect(this->viewer, &Viewer::mouseDoubleClickedOnMap, this->biomeInterface.get(), &BiomeInterface::mouseDoubleClickOnMapEvent);

//    QObject::connect(qApp, &QApplication::focusChanged, this, [=](QWidget*, QWidget*) {
//        this->setFocus(Qt::OtherFocusReason);
//    });
    setupUi();
}

ViewerInterface::~ViewerInterface()
{
    for (auto& action : actionInterfaces)
        action.second->setParent(nullptr);
//    this->setUpdatesEnabled(false);
//    qDeleteAll(this->findChildren<QWidget*>("", Qt::FindDirectChildrenOnly));
//    this->setUpdatesEnabled(true);
}

void ViewerInterface::setupUi()
{
    // Icons
    QIcon openIcon(":/icons/src/assets/open_button.png");
    QIcon faultSlipIcon(":/icons/src/assets/fault-slip_button.png");
    QIcon flowfieldIcon(":/icons/src/assets/flowfield_button.png");
    QIcon gravityIcon(":/icons/src/assets/gravity_button.png");
    QIcon karstIcon(":/icons/src/assets/karst_button.png");
    QIcon karstPeytavieIcon(":/icons/src/assets/karst_peytavie_button.png");
    QIcon recordingIcon(":/icons/src/assets/recording_button.png");
    QIcon savingIcon(":/icons/src/assets/save_button.png");
    QIcon tunnelIcon(":/icons/src/assets/tunnel_button.png");
    QIcon manualEditIcon(":/icons/src/assets/manual-edit_button.png");
    QIcon undoIcon(":/icons/src/assets/undo_button.png");
    QIcon redoIcon(":/icons/src/assets/redo_button.png");
    QIcon marchingCubesIcon(":/icons/src/assets/display-marching-cubes_button.png");
    QIcon rawVoxelsIcon(":/icons/src/assets/display-voxels_button.png");
    QIcon heightmapIcon(":/icons/src/assets/heightmap_button.png");
    QIcon layerBasedIcon(":/icons/src/assets/layerbased_button.png");
    QIcon voxelGridIcon(":/icons/src/assets/voxelgrid_button.png");
    QIcon noDisplayIcon(":/icons/src/assets/no_display.png");
    QIcon wireModeIcon(":/icons/src/assets/wire_mode.png");
    QIcon fillModeIcon(":/icons/src/assets/fill_mode.png");
    QIcon erosionIcon(":/icons/src/assets/erosion.png");
    QIcon heightmapErosionIcon(":/icons/src/assets/heightmap-erosion.png");
    QIcon biomeIcon(":/icons/src/assets/biomes.png");
    // Actions
    QAction *openAction = new QAction(openIcon, "Ouvrir une map existante");
    openAction->setShortcuts(QKeySequence::Open);

    QAction *faultSlipAction = new QAction(faultSlipIcon, "Glissements de terrains");

    QAction *flowfieldAction = new QAction(flowfieldIcon, "Simulation de courant");

    QAction *gravityAction = new QAction(gravityIcon, "Gravité");

    QAction *karstAction = new QAction(karstIcon, "Création de karsts par colonization d'espace");

    QAction *karstPeytavieAction = new QAction(karstPeytavieIcon, "Création de karsts par graphe");

    QAction *recordingAction = new QAction(recordingIcon, "Enregistrement vidéo");

    QAction *savingAction = new QAction(savingIcon, "Sauvegarder la map");
    savingAction->setShortcuts(QKeySequence::Save);

    QAction *tunnelAction = new QAction(tunnelIcon, "Création de tunnels");

    QAction *manualEditAction = new QAction(manualEditIcon, "Manipulations manuelles");

    QAction *undoAction = new QAction(undoIcon, "Undo");
    undoAction->setShortcuts(QKeySequence::Undo);

    QAction *redoAction = new QAction(redoIcon, "Redo");
    redoAction->setShortcuts(QKeySequence::Redo);

    QAction *marchingCubesAction = new QAction(marchingCubesIcon, "Affichage lisse (Marching cubes)");
    marchingCubesAction->setCheckable(true);

    QAction *rawVoxelsAction = new QAction(rawVoxelsIcon, "Affichage brut (voxels)");
    rawVoxelsAction->setCheckable(true);

    QAction *heightmapAction = new QAction(heightmapIcon, "Grille 2D (carte de hauteur)");
    heightmapAction->setCheckable(true);

    QAction *layerBasedAction = new QAction(layerBasedIcon, "Stacks (Layer based)");
    layerBasedAction->setCheckable(true);

    QAction *voxelGridAction = new QAction(voxelGridIcon, "Grille 3D (voxel grid)");
    voxelGridAction->setCheckable(true);
    voxelGridAction->setChecked(true);

    QAction *wireModeAction = new QAction(wireModeIcon, "Affichage fil de fer");
    wireModeAction->setCheckable(true);

    QAction *fillModeAction = new QAction(fillModeIcon, "Affichage normal");
    fillModeAction->setCheckable(true);

    QAction *noDisplayAction = new QAction(noDisplayIcon, "Cacher la carte");
    noDisplayAction->setCheckable(true);
    fillModeAction->setChecked(true);

    QAction *erosionAction = new QAction(erosionIcon, "Eroder la carte");

    QAction *heightmapErosionAction = new QAction(heightmapErosionIcon, "Eroder la heightmap");

    QAction *biomeAction = new QAction(biomeIcon, "Gérer les biomes");



    QMenu* fileMenu = new QMenu("File");
    fileMenu->addActions({openAction, savingAction});

    QMenu* editMenu = new QMenu("Edit");
    editMenu->addActions({undoAction, redoAction});

    QMenu* viewMenu = new QMenu("Affichage");

    QActionGroup* terrainSmoothGroup = new QActionGroup(viewMenu);
    terrainSmoothGroup->addAction(marchingCubesAction);
    terrainSmoothGroup->addAction(rawVoxelsAction);
    viewMenu->addActions(terrainSmoothGroup->actions());

    QActionGroup* terrainDisplayGroup = new QActionGroup(viewMenu);
    terrainDisplayGroup->addAction(wireModeAction);
    terrainDisplayGroup->addAction(fillModeAction);
    terrainDisplayGroup->addAction(noDisplayAction);
    viewMenu->addActions(terrainDisplayGroup->actions());

    QMenu* modelMenu = new QMenu("Model");

    QActionGroup* terrainTypeGroup = new QActionGroup(modelMenu);
    terrainTypeGroup->addAction(voxelGridAction);
    terrainTypeGroup->addAction(heightmapAction);
    terrainTypeGroup->addAction(layerBasedAction);
    modelMenu->addActions(terrainTypeGroup->actions());
    modelMenu->addActions({biomeAction});

    QMenu* recordingMenu = new QMenu("Enregistrements");
    recordingMenu->addActions({recordingAction});

    QMenu* physicsMenu = new QMenu("Physiques");
    physicsMenu->addActions({gravityAction, flowfieldAction, erosionAction, heightmapErosionAction});

    QMenu* diggingMenu = new QMenu("Creuser");
    diggingMenu->addActions({tunnelAction, karstAction, karstPeytavieAction, manualEditAction, faultSlipAction});

    QMenuBar* menu = new QMenuBar(this);
    menu->addMenu(fileMenu);
    menu->addMenu(editMenu);
    menu->addMenu(viewMenu);
    menu->addMenu(recordingMenu);
    menu->addMenu(physicsMenu);
    menu->addMenu(diggingMenu);
    menu->addMenu(modelMenu);

    this->setMenuBar(menu);

    QToolBar* toolbar = new QToolBar("Main tools");
    toolbar->addActions({openAction, savingAction});
    toolbar->addSeparator();
    toolbar->addActions({undoAction, redoAction});
    toolbar->addSeparator();
    toolbar->addActions({marchingCubesAction, rawVoxelsAction});
    toolbar->addSeparator();
    toolbar->addActions({heightmapAction, voxelGridAction, layerBasedAction});
    toolbar->addSeparator();
    toolbar->addActions({biomeAction});
    toolbar->addSeparator();
    toolbar->addActions({recordingAction});
    toolbar->addSeparator();
    toolbar->addActions({gravityAction, flowfieldAction, erosionAction, heightmapErosionAction});
    toolbar->addSeparator();
    toolbar->addActions({tunnelAction, karstAction, karstPeytavieAction, manualEditAction, faultSlipAction});
    toolbar->addSeparator();
    toolbar->addActions({wireModeAction, fillModeAction, noDisplayAction});

    this->addToolBar(Qt::ToolBarArea::TopToolBarArea, toolbar);

    QStatusBar* status = new QStatusBar(this);
    this->setStatusBar(status);

    QDockWidget* displayOptionWidget = new QDockWidget("Affichage", this);
    displayOptionWidget->setFeatures(QDockWidget::DockWidgetFloatable);
    displayOptionWidget->setLayout(new QVBoxLayout());
    this->addDockWidget(Qt::DockWidgetArea::BottomDockWidgetArea, displayOptionWidget);

    this->viewer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    mainLayout = new QGridLayout;
    mainLayout->setColumnStretch(0, 99);
    mainLayout->setColumnStretch(1,  1);

    controlLayout = new QVBoxLayout;
    this->displayModeLayout = new QHBoxLayout;
    this->displayModeBox = new Spoiler("Affichage");
    this->LoDChooserLayout = new QHBoxLayout;
    this->LoDChooserBox = new Spoiler("Niveau de détail");
    this->mapSliceSliderX = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderY = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderZ = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    QCheckBox* sliderXactivation = new QCheckBox("Activer");
    sliderXactivation->setChecked(true);
    QObject::connect(sliderXactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapX = this->mapSliceSliderX->min_value(); this->viewer->maxSliceMapX = this->mapSliceSliderX->max_value(); }
        else        { this->viewer->minSliceMapX = 0.f;                                this->viewer->maxSliceMapX = 1.f; }
        this->viewer->update();
    });
    QCheckBox* sliderYactivation = new QCheckBox("Activer");
    sliderYactivation->setChecked(true);
    QObject::connect(sliderYactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapY = this->mapSliceSliderY->min_value(); this->viewer->maxSliceMapY = this->mapSliceSliderY->max_value(); }
        else        { this->viewer->minSliceMapY = 0.f;                                this->viewer->maxSliceMapY = 1.f; }
        this->viewer->update();
    });
    QCheckBox* sliderZactivation = new QCheckBox("Activer");
    sliderZactivation->setChecked(true);
    QObject::connect(sliderZactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapZ = this->mapSliceSliderZ->min_value(); this->viewer->maxSliceMapZ = this->mapSliceSliderZ->max_value(); }
        else        { this->viewer->minSliceMapZ = 0.f;                                this->viewer->maxSliceMapZ = 1.f; }
        this->viewer->update();
    });
    displayModeLayout->addWidget(createMultipleSliderGroupWithCheckbox({
                                                               {"X", {mapSliceSliderX, sliderXactivation}},
                                                               {"Y", {mapSliceSliderY, sliderYactivation}},
                                                               {"Z", {mapSliceSliderZ, sliderZactivation}}
                                                           }));

    this->isolevelSelectionSlider = new RangeSlider(Qt::Orientation::Horizontal, 0.f, 3.f, 0.1f);
    QCheckBox* isolevelSelectionActivation = new QCheckBox("Activer");
    LoDChooserLayout->addWidget(createMultipleSliderGroupWithCheckbox({
                                                                          {"Densité", {isolevelSelectionSlider, isolevelSelectionActivation}}
                                                                      }));
    isolevelSelectionActivation->setChecked(true);
    QObject::connect(isolevelSelectionActivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->terrainGenerationInterface->minIsoLevel = isolevelSelectionSlider->min_value(); this->terrainGenerationInterface->maxIsoLevel = isolevelSelectionSlider->max_value();}
        else { this->terrainGenerationInterface->minIsoLevel = -1000.f; this->terrainGenerationInterface->maxIsoLevel = 1000.f; }
        this->viewer->update();
    });

    mainLayout->addWidget(viewer, 1, 0);


    QHBoxLayout* displayOptionLayout = new QHBoxLayout();
    displayOptionLayout->addItem(displayModeLayout);
    displayOptionLayout->addWidget(createVerticalGroup({
                                                           createMultipleSliderGroupWithCheckbox({
                                                               {"Densité", {isolevelSelectionSlider, isolevelSelectionActivation}}
                                                           })
                                                 }));
    QGroupBox* displayOptionBox = new QGroupBox();
    displayOptionBox->setLayout(displayOptionLayout);
    displayOptionWidget->setWidget(displayOptionBox);


    QWidget* mainFrame = new QWidget(this);
    mainFrame->setLayout(mainLayout);
    this->setCentralWidget(mainFrame);


    frame = new StickyFrame(this->viewer, 0, -1, -1, 1, false);
    frame->hide();

    QObject::connect(openAction, &QAction::triggered, this, &ViewerInterface::openMapUI);
    QObject::connect(savingAction, &QAction::triggered, this, &ViewerInterface::saveMapUI);
    QObject::connect(faultSlipAction, &QAction::triggered, this, &ViewerInterface::openFaultSlipInterface);
    QObject::connect(flowfieldAction, &QAction::triggered, this, &ViewerInterface::openFlowfieldInterface);
    QObject::connect(karstAction, &QAction::triggered, this, &ViewerInterface::openKarstInterface);
    QObject::connect(karstPeytavieAction, &QAction::triggered, this, &ViewerInterface::openKarstPeytavieInterface);
    QObject::connect(recordingAction, &QAction::triggered, this->viewer, &Viewer::startStopRecording);
    QObject::connect(gravityAction, &QAction::triggered, this, &ViewerInterface::openGravityInterface);
    QObject::connect(tunnelAction, &QAction::triggered, this, &ViewerInterface::openTunnelInterface);
    QObject::connect(manualEditAction, &QAction::triggered, this, &ViewerInterface::openManualEditionInterface);
    QObject::connect(erosionAction, &QAction::triggered, this, &ViewerInterface::openErosionInterface);
    QObject::connect(heightmapErosionAction, &QAction::triggered, this, &ViewerInterface::openHeightmapErosionInterface);
    QObject::connect(biomeAction, &QAction::triggered, this, &ViewerInterface::openBiomeInterface);
    QObject::connect(undoAction, &QAction::triggered, this->undoRedoInterface.get(), &UndoRedoInterface::undo);
    QObject::connect(redoAction, &QAction::triggered, this->undoRedoInterface.get(), &UndoRedoInterface::redo);
    QObject::connect(marchingCubesAction, &QAction::triggered, this, [&]() {
        this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES);
    });
    QObject::connect(rawVoxelsAction, &QAction::triggered, this, [&]() {
        this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE);
    });
    QObject::connect(heightmapAction, &QAction::triggered, this, [&]() {
        this->viewer->setMapMode(MapMode::GRID_MODE);
    });
    QObject::connect(layerBasedAction, &QAction::triggered, this, [&]() {
        this->viewer->setMapMode(MapMode::LAYER_MODE);
    });
    QObject::connect(voxelGridAction, &QAction::triggered, this, [&]() {
        this->viewer->setMapMode(MapMode::VOXEL_MODE);
    });
    QObject::connect(wireModeAction, &QAction::triggered, this, [&]() {
        this->viewer->setViewerMode(WIRE_MODE);
    });
    QObject::connect(fillModeAction, &QAction::triggered, this, [&]() {
        this->viewer->setViewerMode(FILL_MODE);
    });
    QObject::connect(noDisplayAction, &QAction::triggered, this, [&]() {
        this->viewer->setViewerMode(NO_DISPLAY);
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

void ViewerInterface::openInterface(std::string interfaceName)
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == interfaceName) {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = interfaceName;
        this->faultSlip->show();
        this->frame->setContent(this->actionInterfaces[interfaceName]->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openFaultSlipInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "FaultSlip") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "FaultSlip";
        this->faultSlip->show();
        this->frame->setContent(this->faultSlip->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openFlowfieldInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "FlowField") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "FlowField";
        this->flowField->show();
        this->frame->setContent(this->flowField->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openKarstInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "Karsts") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "Karsts";
        this->spaceColonization->show();
        this->frame->setContent(this->spaceColonization->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}
void ViewerInterface::openKarstPeytavieInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "Karsts-Peytavie") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "Karsts-Peytavie";
        this->karstPathGeneration->show();
        this->frame->setContent(this->karstPathGeneration->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}
void ViewerInterface::openTunnelInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "Tunnels") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "Tunnels";
        this->tunnelInterface->show();
        this->frame->setContent(this->tunnelInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}
void ViewerInterface::openManualEditionInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "ManualEdition") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "ManualEdition";
        this->manualEditionInterface->show();
        this->frame->setContent(this->manualEditionInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}
void ViewerInterface::openGravityInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "Gravity") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "Gravity";
        this->gravityInterface->show();
        this->frame->setContent(this->gravityInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openErosionInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "Erosion") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "Erosion";
        this->erosionInterface->show();
        this->frame->setContent(this->erosionInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openHeightmapErosionInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "HeightmapErosion") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "HeightmapErosion";
        this->heightmapErosionInterface->show();
        this->frame->setContent(this->heightmapErosionInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::openBiomeInterface()
{
    this->hideAllInteractiveParts();
    if (lastPanelOpenedByStickyFrame == "BiomeGeneration") {
        lastPanelOpenedByStickyFrame = "";
        this->frame->hide();
    } else {
        lastPanelOpenedByStickyFrame = "BiomeGeneration";
        this->biomeInterface->show();
        this->frame->setContent(this->biomeInterface->createGUI());
        this->frame->show();
    }
    this->viewer->update();
}

void ViewerInterface::hideAllInteractiveParts()
{
    for (auto& actionInterface : this->actionInterfaces)
        actionInterface.second->hide();
    /*this->karstPathGeneration->hide();
    this->spaceColonization->hide();
    this->faultSlip->hide();
    this->flowField->hide();
    this->tunnelInterface->hide();
    this->manualEditionInterface->hide();
    this->gravityInterface->hide();
    this->erosionInterface->hide();
    this->heightmapErosionInterface->hide();
    this->biomeInterface->hide();*/

    this->viewer->update();
    this->update();
}

void ViewerInterface::openMapUI()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Ouvrir une carte"), QString::fromStdString(this->mapSavingFolder));
    terrainGenerationInterface->createTerrainFromFile(q_filename.toStdString(), this->actionInterfaces);
}

void ViewerInterface::saveMapUI()
{
    QString q_filename = QFileDialog::getSaveFileName(this, QString("Enregistrer la carte"), QString::fromStdString(this->mapSavingFolder));
    terrainGenerationInterface->saveTerrain(q_filename.toStdString());
}

