#include "Interface/Interface.h"
#include "Interface/InterfaceUtils.h"

ViewerInterface::ViewerInterface() {
    this->setWindowFlag(Qt::WindowType::WindowMaximizeButtonHint);
    this->setWindowFlag(Qt::WindowType::WindowMinimizeButtonHint);
    this->viewer = new Viewer(this);
    this->karstPathGeneration = std::make_shared<KarstPathGenerationInterface>(this);
    this->karstPathGeneration->hide();
    this->viewer->karstPathInterface = this->karstPathGeneration;
    this->spaceColonization = std::make_shared<SpaceColonizationInterface>(this);
    this->spaceColonization->hide();
    this->viewer->spaceColonizationInterface = this->spaceColonization;
    this->faultSlip = std::make_shared<FaultSlipInterface>(this);
    this->faultSlip->hide();
    this->viewer->faultSlipInterface = this->faultSlip;
    this->flowField = std::make_shared<FlowFieldInterface>(this);
    this->viewer->flowFieldInterface = flowField;
    this->flowField->hide();
    this->tunnelInterface = std::make_shared<TunnelInterface>(this);
    this->viewer->tunnelInterface = this->tunnelInterface;
    this->tunnelInterface->hide();
    this->manualEditionInterface = std::make_shared<ManualEditionInterface>(this);
    this->viewer->manualEditionInterface = this->manualEditionInterface;
    this->manualEditionInterface->hide();
    this->gravityInterface = std::make_shared<GravityInterface>(this);
    this->gravityInterface->hide();

    QObject::connect(this->viewer, &Viewer::viewerInitialized, this, [&](){
        this->karstPathGeneration->affectVoxelGrid(this->viewer->voxelGrid);
        this->spaceColonization->affectVoxelGrid(this->viewer->voxelGrid);
        this->faultSlip->affectVoxelGrid(this->viewer->voxelGrid);
        this->flowField->affectVoxelGrid(this->viewer->voxelGrid);
        this->tunnelInterface->affectVoxelGrid(this->viewer->voxelGrid);
        this->manualEditionInterface->affectVoxelGrid(this->viewer->voxelGrid);
        this->gravityInterface->affectVoxelGrid(this->viewer->voxelGrid);
    });

    QObject::connect(this->karstPathGeneration.get(), &KarstPathGenerationInterface::karstPathUpdated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->spaceColonization.get(), &SpaceColonizationInterface::karstPathUpdated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->faultSlip.get(), &FaultSlipInterface::faultSlipApplied,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->flowField.get(), &FlowFieldInterface::updated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::updated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->viewer, &Viewer::mouseClickOnMap,
                     this->tunnelInterface.get(), &TunnelInterface::mouseClickInWorldEvent);
    QObject::connect(this->manualEditionInterface.get(), &ManualEditionInterface::updated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->viewer, &Viewer::mouseClickOnMap,
                     this->manualEditionInterface.get(), &ManualEditionInterface::mouseClickedOnMapEvent);
    QObject::connect(this->viewer, &Viewer::mouseMovedOnMap,
                     this->manualEditionInterface.get(), &ManualEditionInterface::setPosition);
    QObject::connect(this->gravityInterface.get(), &GravityInterface::updated,
                     this, [&](){ this->viewer->update(); });
    QObject::connect(this->tunnelInterface.get(), &TunnelInterface::needToClipView,
                     this->viewer, &Viewer::clipViewTemporarily);

    QObject::connect(qApp, &QApplication::focusChanged, this, [=](QWidget*, QWidget*) {
        this->setFocus(Qt::OtherFocusReason);
    });
    setupUi();
}

ViewerInterface::~ViewerInterface()
{
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
    // Actions
    QAction *openAction = new QAction(openIcon, "Ouvrir une map existante");
    openAction->setShortcuts(QKeySequence::Open);
    QAction *faultSlipAction = new QAction(faultSlipIcon, "Glissements de terrains");
//    faultSlipAction->setShortcuts(QKeySequence::Open);
    QAction *flowfieldAction = new QAction(flowfieldIcon, "Simulation de courant");
//    flowfieldAction->setShortcuts(QKeySequence::Open);
    QAction *gravityAction = new QAction(gravityIcon, "Gravité");
//    gravityAction->setShortcuts(QKeySequence::Open);
    QAction *karstAction = new QAction(karstIcon, "Création de karsts");
//    karstAction->setShortcuts(QKeySequence::Open);
    QAction *recordingAction = new QAction(recordingIcon, "Enregistrement vidéo");
//    recordingAction->setShortcuts(QKeySequence::);
    QAction *savingAction = new QAction(savingIcon, "Sauvegarder la map");
    savingAction->setShortcuts(QKeySequence::Save);
    QAction *tunnelAction = new QAction(tunnelIcon, "Création de tunnels");
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *manualEditAction = new QAction(manualEditIcon, "Manipulations manuelles");
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *undoAction = new QAction(undoIcon, "Undo");
    undoAction->setShortcuts(QKeySequence::Undo);
    QAction *redoAction = new QAction(redoIcon, "Redo");
    redoAction->setShortcuts(QKeySequence::Redo);
    QAction *marchingCubesAction = new QAction(marchingCubesIcon, "Affichage lisse (Marching cubes)");
    marchingCubesAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *rawVoxelsAction = new QAction(rawVoxelsIcon, "Affichage brut (voxels)");
    rawVoxelsAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *heightmapAction = new QAction(heightmapIcon, "Grille 2D (carte de hauteur)");
    heightmapAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *layerBasedAction = new QAction(layerBasedIcon, "Stacks (Layer based)");
    layerBasedAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *voxelGridAction = new QAction(voxelGridIcon, "Grille 3D (voxel grid)");
    voxelGridAction->setCheckable(true);
    voxelGridAction->setChecked(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *wireModeAction = new QAction(wireModeIcon, "Affichage fil de fer");
    wireModeAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *fillModeAction = new QAction(fillModeIcon, "Affichage normal");
    fillModeAction->setCheckable(true);
//    tunnelAction->setShortcuts(QKeySequence::Open);
    QAction *noDisplayAction = new QAction(noDisplayIcon, "Cacher la carte");
    noDisplayAction->setCheckable(true);
    //    tunnelAction->setShortcuts(QKeySequence::Open);
    fillModeAction->setChecked(true);


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
    QMenu* recordingMenu = new QMenu("Enregistrements");
    recordingMenu->addActions({recordingAction});
    QMenu* physicsMenu = new QMenu("Physiques");
    physicsMenu->addActions({gravityAction, flowfieldAction});
    QMenu* diggingMenu = new QMenu("Creuser");
    diggingMenu->addActions({tunnelAction, karstAction, manualEditAction, faultSlipAction});

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
    toolbar->addActions({recordingAction});
    toolbar->addSeparator();
    toolbar->addActions({gravityAction, flowfieldAction});
    toolbar->addSeparator();
    toolbar->addActions({tunnelAction, karstAction, manualEditAction, faultSlipAction});
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
//    this->toolbox = new QMenuBar(this);
//    this->toolbox->addActions({openAction, savingAction, manualEditAction});

    mainLayout->setColumnStretch(0, 99);
    mainLayout->setColumnStretch(1,  1);

    controlLayout = new QVBoxLayout;

    // Create the controls layouts
    this->loadSaveLayout = new QHBoxLayout;
    this->loadSaveBox = new Spoiler("Enregistrer / Charger");
    this->randomRocksLayout = new QHBoxLayout;
    this->randomRocksBox = new Spoiler("Erosion");
    this->currentSimulationLayout = new QHBoxLayout;
    this->currentSimulationBox = new Spoiler("Simulation de courants");
    this->manualRocksLayout = new QHBoxLayout;
    this->manualRocksBox = new Spoiler("Modifs manuelles");
    this->curvesErosionLayout = new QHBoxLayout;
    this->curvesErosionBox = new Spoiler("Erosion par courbes aléatoires");
//    this->karstCreationLayout = new QHBoxLayout;
    this->karstCreationBox = new Spoiler("Creation de karsts");
    this->spaceColonizerLayout = new QHBoxLayout;
    this->spaceColonizerBox = new Spoiler("Creation de karsts par colonisation");
    this->faultSlipLayout = new QHBoxLayout;
    this->faultSlipBox = new Spoiler("Fracture");
    this->gravityLayout = new QHBoxLayout;
    this->gravityBox = new Spoiler("Gravité");
    this->recordingLayout = new QHBoxLayout;
    this->recordingBox = new Spoiler("Enregistrements");
    this->displayModeLayout = new QHBoxLayout;
    this->displayModeBox = new Spoiler("Affichage");
    this->algorithmLayout = new QHBoxLayout;
    this->algorithmBox = new Spoiler("Algorithme de mesh");
    this->displayTypeLayout = new QHBoxLayout;
    this->displayTypeBox = new Spoiler("Type de terrain");
    this->LoDChooserLayout = new QHBoxLayout;
    this->LoDChooserBox = new Spoiler("Niveau de détail");


    // random rocks layout
    this->saveButton = new QPushButton("Sauvegarder...");
    this->loadButton = new QPushButton("Charger...");
    loadSaveLayout->addWidget(saveButton);
    loadSaveLayout->addWidget(loadButton);
    loadSaveBox->setContentLayout(*loadSaveLayout);
    controlLayout->addWidget(loadSaveBox);


    // random rocks layout
    this->randomRocksSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->randomRocksStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->randomRocksQuantitySlider = new FancySlider(Qt::Orientation::Horizontal, 1, 1000);
    this->sendRandomRocksFromCam = new QPushButton("Send from cam");
    this->sendRandomRocks = new QPushButton("Send random");
    this->displayRandomRocks = new QCheckBox("Afficher");
    this->displayFailedRandomRocks = new QCheckBox("Echecs");
    randomRocksLayout->addWidget(createSliderGroup("Taille", randomRocksSizeSlider));
    randomRocksLayout->addWidget(createSliderGroup("Force", randomRocksStrengthSlider));
    randomRocksLayout->addWidget(createSliderGroup("Quantite", randomRocksQuantitySlider));
    randomRocksLayout->addWidget(createVerticalGroup({sendRandomRocksFromCam, sendRandomRocks}));
    randomRocksLayout->addWidget(createVerticalGroup({displayRandomRocks, displayFailedRandomRocks}));
    randomRocksBox->setContentLayout(*randomRocksLayout);
    controlLayout->addWidget(randomRocksBox);


    // Manual rock erosion layout
    this->currentSimulationFlowfieldStrengthSlider= new FancySlider(Qt::Orientation::Horizontal, -1.0, 1.0, 0.01);
    this->currentSimulationStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->currentSimulationRandomDirectionSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 0.5, 0.01);
    this->displayFlowfield = new QCheckBox("Afficher le courant");
    this->recomputeFlowfieldButton = new QPushButton("Recalculer le courant");
    currentSimulationLayout->addWidget(createSliderGroup("Facteur d'esquive des collisions", currentSimulationFlowfieldStrengthSlider));
    currentSimulationLayout->addWidget(createSliderGroup("Force du courant", currentSimulationStrengthSlider));
    currentSimulationLayout->addWidget(createSliderGroup("Aleatoire", currentSimulationRandomDirectionSlider));
    currentSimulationLayout->addWidget(createVerticalGroup({displayFlowfield, recomputeFlowfieldButton}));
    currentSimulationBox->setContentLayout(*currentSimulationLayout);
    controlLayout->addWidget(currentSimulationBox);


    // Manual rock erosion layout
    this->manualRocksSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->manualRocksStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->addingMode = new QRadioButton("Ajouter de la matière"); addingMode->setChecked(true);
    this->suppressMode = new QRadioButton("Detruire de la matière");
    manualRocksLayout->addWidget(createSliderGroup("Taille", manualRocksSizeSlider));
    manualRocksLayout->addWidget(createSliderGroup("Force", manualRocksStrengthSlider));
    manualRocksLayout->addWidget(createVerticalGroup({addingMode, suppressMode}));
    manualRocksBox->setContentLayout(*manualRocksLayout);
    controlLayout->addWidget(manualRocksBox);


    // Curved erosion layout
    this->curvesAddControlPointButton = new QPushButton("Ajouter un point de control");
    this->curvesClearControlPointButton = new QPushButton("Tout retirer");
    this->curvesErosionSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 30);
    this->curvesErosionStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->curvesErosionCreateMatter = new QPushButton("Creer un pont");
    this->curvesErosionRemoveMatter = new QPushButton("Creer un tunnel");
    this->curvesErosionCreateCrack = new QPushButton("Creer une faille");
    this->displayCurvesErosion = new QCheckBox("Afficher");
    curvesErosionLayout->addWidget(createVerticalGroup({curvesErosionCreateMatter, curvesErosionRemoveMatter, curvesErosionCreateCrack}));
    curvesErosionLayout->addWidget(createVerticalGroup({curvesAddControlPointButton, curvesClearControlPointButton}));
    curvesErosionLayout->addWidget(createSliderGroup("Taille", curvesErosionSizeSlider));
    curvesErosionLayout->addWidget(createSliderGroup("Force", curvesErosionStrengthSlider));
    curvesErosionLayout->addWidget(displayCurvesErosion);
    curvesErosionBox->setContentLayout(*curvesErosionLayout);
    controlLayout->addWidget(curvesErosionBox);

    // Karst creation layout
    this->karstCreationLayout = this->karstPathGeneration->createGUI();
    karstCreationBox->setContentLayout(*this->karstCreationLayout);//*karstCreationLayout);
    controlLayout->addWidget(karstCreationBox);

    // Space colonization layout
    this->spaceColonizerLayout = this->spaceColonization->createGUI();
    spaceColonizerBox->setContentLayout(*spaceColonizerLayout);
    controlLayout->addWidget(spaceColonizerBox);

    // Fault slip layout
    this->faultSlipLayout = this->faultSlip->createGUI();
    faultSlipBox->setContentLayout(*faultSlipLayout);
    controlLayout->addWidget(faultSlipBox);

    // Gravity controls
    this->gravityGlobalButton = new QPushButton("Chutes de pierres");
//    this->gravityGlobalButton->setCheckable(true);
    this->gravitySandButton = new QPushButton("Chutes de sable");
//    this->gravitySandButton->setCheckable(true);
    gravityLayout->addWidget(gravityGlobalButton);
    gravityLayout->addWidget(gravitySandButton);
    gravityBox->setContentLayout(*gravityLayout);
    controlLayout->addWidget(gravityBox);

    // Recording control
    this->startStopRecording = new QPushButton("Enregistrer / stopper l'enregistrement");
    this->startStopRecording->setCheckable(true);
    recordingLayout->addWidget(startStopRecording);
    recordingBox->setContentLayout(*recordingLayout);
    controlLayout->addWidget(recordingBox);

    // Display mode control
    this->wireModeButton = new QRadioButton("Wire mode");
    this->fillModeButton = new QRadioButton("Fill mode");//  fillModeButton->setChecked(true);
    this->invisibleModeButton = new QRadioButton("No display");
    this->mapSliceSliderX = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderY = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    this->mapSliceSliderZ = new RangeSlider(Qt::Orientation::Horizontal, 0, 1, 0.01f);
    displayModeLayout->addWidget(wireModeButton);
    displayModeLayout->addWidget(fillModeButton);
    displayModeLayout->addWidget(invisibleModeButton);


    /*QGroupBox* optionalXGroup = new QGroupBox();
    QHBoxLayout* optionalXLayout = new QHBoxLayout();*/
    QCheckBox* sliderXactivation = new QCheckBox("Activer");
    sliderXactivation->setChecked(true);
    QObject::connect(sliderXactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapX = this->mapSliceSliderX->min_value(); this->viewer->maxSliceMapX = this->mapSliceSliderX->max_value(); }
        else        { this->viewer->minSliceMapX = 0.f;                                this->viewer->maxSliceMapX = 1.f; }
        this->viewer->update();
    });/*
    optionalXLayout->addWidget(createSliderGroup("X", mapSliceSliderX, true));
    optionalXLayout->addWidget(sliderXactivation);
    optionalXGroup->setLayout(optionalXLayout);


    QGroupBox* optionalYGroup = new QGroupBox();
    QHBoxLayout* optionalYLayout = new QHBoxLayout();*/
    QCheckBox* sliderYactivation = new QCheckBox("Activer");
    sliderYactivation->setChecked(true);
    QObject::connect(sliderYactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapY = this->mapSliceSliderY->min_value(); this->viewer->maxSliceMapY = this->mapSliceSliderY->max_value(); }
        else        { this->viewer->minSliceMapY = 0.f;                                this->viewer->maxSliceMapY = 1.f; }
        this->viewer->update();
    });/*
    optionalYLayout->addWidget(createSliderGroup("Y", mapSliceSliderY, true));
    optionalYLayout->addWidget(sliderYactivation);
    optionalYGroup->setLayout(optionalYLayout);


    QGroupBox* optionalZGroup = new QGroupBox();
    QHBoxLayout* optionalZLayout = new QHBoxLayout();*/
    QCheckBox* sliderZactivation = new QCheckBox("Activer");
    sliderZactivation->setChecked(true);
    QObject::connect(sliderZactivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { this->viewer->minSliceMapZ = this->mapSliceSliderZ->min_value(); this->viewer->maxSliceMapZ = this->mapSliceSliderZ->max_value(); }
        else        { this->viewer->minSliceMapZ = 0.f;                                this->viewer->maxSliceMapZ = 1.f; }
        this->viewer->update();
    });/*
    optionalZLayout->addWidget(createSliderGroup("Z", mapSliceSliderZ, true));
    optionalZLayout->addWidget(sliderZactivation);
    optionalZGroup->setLayout(optionalZLayout);
    */

//    displayModeLayout->addWidget(createVerticalGroup({mapSliceSliderX, optionalYGroup, optionalZGroup}));
                                                         /*
                         createOptionalSlider(mapSliceSliderX, "X", true, [&](bool activate, RangeSlider* slider) {
                             if (activate) {
                                 this->viewer->minSliceMapX = slider->min_value();
                                 this->viewer->maxSliceMapX = slider->max_value();
                             } else {
                                 this->viewer->minSliceMapX = 0.f;
                                 this->viewer->maxSliceMapX = 1.f;
                             }
                             this->viewer->update();
                         }),
                         createOptionalSlider(mapSliceSliderY, "Y", true, [&](bool activate, RangeSlider* slider) {
                             if (activate) {
                                 this->viewer->minSliceMapY = slider->min_value();
                                 this->viewer->maxSliceMapY = slider->max_value();
                             } else {
                                 this->viewer->minSliceMapY = 0.f;
                                 this->viewer->maxSliceMapY = 1.f;
                             }
                             this->viewer->update();
                         }),
                         createOptionalSlider(mapSliceSliderZ, "Z", true, [&](bool activate, RangeSlider* slider) {
                             if (activate) {
                                 this->viewer->minSliceMapZ = slider->min_value();
                                 this->viewer->maxSliceMapZ = slider->max_value();
                             } else {
                                 this->viewer->minSliceMapZ = 0.f;
                                 this->viewer->maxSliceMapZ = 1.f;
                             }
                             this->viewer->update();
                         })
                     }
                     ));*/
/*
    displayModeLayout->addWidget(createMultipleSliderGroup({
                                                               {"X", mapSliceSliderX},
                                                               {"Y", mapSliceSliderY},
                                                               {"Z", mapSliceSliderZ}
                                                           })); */
    displayModeLayout->addWidget(createMultipleSliderGroupWithCheckbox({
                                                               {"X", {mapSliceSliderX, sliderXactivation}},
                                                               {"Y", {mapSliceSliderY, sliderYactivation}},
                                                               {"Z", {mapSliceSliderZ, sliderZactivation}}
                                                           }));
//    displayModeBox->setContentLayout(*displayModeLayout);
//    controlLayout->addWidget(displayModeBox);

    // Display algorithm controls
    this->noAlgorithmButton = new QRadioButton("Voxels");
    this->marchingCubesButton = new QRadioButton("Marching cubes"); marchingCubesButton->setChecked(true);
    algorithmLayout->addWidget(noAlgorithmButton);
    algorithmLayout->addWidget(marchingCubesButton);
    algorithmBox->setContentLayout(*algorithmLayout);
//    controlLayout->addWidget(algorithmBox);

    // Display type controls
    this->gridModeButton = new QRadioButton("Heightmap");
    this->voxelsModeButton = new QRadioButton("3D grid"); voxelsModeButton->setChecked(true);
    this->layerModeButton = new QRadioButton("Layer stacks");
    displayTypeLayout->addWidget(gridModeButton);
    displayTypeLayout->addWidget(voxelsModeButton);
    displayTypeLayout->addWidget(layerModeButton);
    displayTypeBox->setContentLayout(*displayTypeLayout);
//    controlLayout->addWidget(displayTypeBox);

    // Display Level of Detail controls
    this->LoDChooserSlider = new FancySlider(Qt::Orientation::Horizontal, 0, this->viewer->getMaxLoDAvailable(), 1);
    this->LoDChooserConfirmButton = new QPushButton("Regénérer");
    LoDChooserLayout->addWidget(createSliderGroup("Niveau de détail", LoDChooserSlider));
    LoDChooserLayout->addWidget(LoDChooserConfirmButton);
    LoDChooserBox->setContentLayout(*LoDChooserLayout);
    controlLayout->addWidget(LoDChooserBox); // TODO : Connect the slider and button

    for (auto& child : controlLayout->children())
    {
        if (child->isWidgetType())
        {
            controlLayout->setAlignment((QWidget*) child, Qt::AlignmentFlag::AlignTop);
        }
    }
    controlLayout->addStretch();

//    mainLayout->addWidget(toolbox, 0, 0, 1, 2);
    mainLayout->addWidget(viewer, 1, 0);
//    mainLayout->addLayout(controlLayout, 1, 1);


    QHBoxLayout* displayOptionLayout = new QHBoxLayout();
//    displayOptionLayout->addWidget(algorithmBox);
//    displayOptionLayout->addWidget(displayTypeBox);
    displayOptionLayout->addItem(displayModeLayout);
    displayOptionLayout->addWidget(createVerticalGroup({
                                                     createSliderGroup("Niveau de détail", LoDChooserSlider),
                                                     LoDChooserConfirmButton
                                                 }));
    QGroupBox* displayOptionBox = new QGroupBox();
    displayOptionBox->setLayout(displayOptionLayout);
    displayOptionWidget->setWidget(displayOptionBox);


    QWidget* mainFrame = new QWidget(this);
    mainFrame->setLayout(mainLayout);
    this->setCentralWidget(mainFrame);


    frame = new StickyFrame(this->viewer, 0, -1, -1, 1, false);
    frame->hide();

    QObject::connect(openAction, &QAction::triggered, this->viewer, &Viewer::loadMapUI);
    QObject::connect(savingAction, &QAction::triggered, this->viewer, &Viewer::saveMapUI);
    QObject::connect(faultSlipAction, &QAction::triggered, this, &ViewerInterface::openFaultSlipInterface);
    QObject::connect(flowfieldAction, &QAction::triggered, this, &ViewerInterface::openFlowfieldInterface);
    QObject::connect(karstAction, &QAction::triggered, this, &ViewerInterface::openKarstInterface);
    QObject::connect(recordingAction, &QAction::triggered, this->viewer, &Viewer::startStopRecording);
    QObject::connect(gravityAction, &QAction::triggered, this, &ViewerInterface::openGravityInterface);
    QObject::connect(tunnelAction, &QAction::triggered, this, &ViewerInterface::openTunnelInterface);
    QObject::connect(manualEditAction, &QAction::triggered, this, &ViewerInterface::openManualEditionInterface);
    QObject::connect(undoAction, &QAction::triggered, this->viewer, &Viewer::undo);
    QObject::connect(redoAction, &QAction::triggered, this->viewer, &Viewer::redo);
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

    QObject::connect(loadButton, &QCheckBox::pressed, this, [=](){this->viewer->loadMapUI(); } );
    QObject::connect(saveButton, &QCheckBox::pressed, this, [=](){this->viewer->saveMapUI(); } );

    QObject::connect(randomRocksSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setErosionRocksSize(int)));
    QObject::connect(randomRocksStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionRocksStrength(float)));
    QObject::connect(randomRocksQuantitySlider, SIGNAL(valueChanged(int)), viewer, SLOT(setErosionRocksQuantity(int)));
    QObject::connect(sendRandomRocksFromCam, &QPushButton::pressed, this, [=](){this->viewer->erodeMap(true); } );
    QObject::connect(sendRandomRocks, &QCheckBox::pressed, this, [=](){this->viewer->erodeMap(false); } );
    QObject::connect(displayRandomRocks, &QCheckBox::toggled, this, [=](bool display){this->viewer->debugMeshes[ROCK_TRAILS].isDisplayed = display; viewer->update(); } );
    QObject::connect(displayFailedRandomRocks, &QCheckBox::toggled, this, [=](bool display){this->viewer->debugMeshes[FAILED_ROCKS].isDisplayed = display; viewer->update(); } );

    QObject::connect(currentSimulationFlowfieldStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionFlowfieldFactor(float)));
//    QObject::connect(currentSimulationStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setManualErosionRocksStrength(float)));
    QObject::connect(currentSimulationRandomDirectionSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionFlowfieldRandomness(float)));
    QObject::connect(displayFlowfield, &QCheckBox::toggled, this, [=](bool display){this->viewer->debugMeshes[FLOW_TRAILS].isDisplayed = display; viewer->update(); } );
//    QObject::connect(recomputeFlowfieldButton, &QPushButton::pressed, this, [=](){this->viewer->recomputeFlowfield(); } );

    QObject::connect(manualRocksSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setManualErosionRocksSize(int)));
    QObject::connect(manualRocksStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setManualErosionRocksStrength(float)));
    QObject::connect(addingMode, &QRadioButton::clicked, this, [=](){this->viewer->setAddingMatterMode(true); } );
    QObject::connect(suppressMode, &QRadioButton::clicked, this, [=](){this->viewer->setAddingMatterMode(false); } );

//    QObject::connect(curvesErosionSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setCurvesErosionSize(int)));
//    QObject::connect(curvesErosionStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setCurvesErosionStrength(float)));
//    QObject::connect(curvesErosionCreateMatter, &QPushButton::pressed, this, [=](){ this->viewer->createTunnel(false); } );
//    QObject::connect(curvesErosionRemoveMatter, &QPushButton::pressed, this, [=](){ this->viewer->createTunnel(true); } );
//    QObject::connect(curvesErosionCreateCrack, &QPushButton::pressed, this, [=](){ this->viewer->createCrack(true); } );
//    QObject::connect(curvesAddControlPointButton, &QPushButton::pressed, this, [=](){this->viewer->setCurvesErosionConstructionMode(true); });
//    QObject::connect(curvesClearControlPointButton, &QPushButton::pressed, this, [=](){this->viewer->clearTunnelPoints(); });
//    QObject::connect(displayCurvesErosion, &QCheckBox::toggled, this, [=](bool display){ this->viewer->debugMeshes[TUNNEL_PATHS].isDisplayed = display; viewer->update(); } );

    QObject::connect(gravityGlobalButton, &QPushButton::pressed, this, [=](){ this->viewer->createGlobalGravity(); } );
    QObject::connect(gravitySandButton  , &QPushButton::pressed, this, [=](){ this->viewer->createSandGravity(); } );

    QObject::connect(startStopRecording  , &QPushButton::pressed, this, [=](){ this->viewer->startStopRecording(); } );

    QObject::connect(wireModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setViewerMode(ViewerMode::WIRE_MODE); } );
    QObject::connect(fillModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setViewerMode(ViewerMode::FILL_MODE); } );
    QObject::connect(invisibleModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setViewerMode(ViewerMode::NO_DISPLAY); } );
    QObject::connect(mapSliceSliderX, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapX = min; viewer->maxSliceMapX = max; viewer->update(); });
    QObject::connect(mapSliceSliderY, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapY = min; viewer->maxSliceMapY = max; viewer->update(); });
    QObject::connect(mapSliceSliderZ, &RangeSlider::alt_valueChanged, this, [=](float min, float max){this->viewer->minSliceMapZ = min; viewer->maxSliceMapZ = max; viewer->update(); });

    QObject::connect(noAlgorithmButton, &QRadioButton::clicked, this, [=](){this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE); } );
    QObject::connect(marchingCubesButton, &QRadioButton::clicked, this, [=](){this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES); } );

    QObject::connect(gridModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::GRID_MODE); } );
    QObject::connect(voxelsModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::VOXEL_MODE); } );
    QObject::connect(layerModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::LAYER_MODE); } );

    QObject::connect(LoDChooserSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setLoD(int)));
    QObject::connect(LoDChooserConfirmButton, &QPushButton::pressed, viewer, &Viewer::computeLoD);


    QObject::connect(karstPathGeneration.get(), &KarstPathGenerationInterface::karstPathUpdated, this, [&](){ this->viewer->update(); });
    QObject::connect(spaceColonization.get(), &SpaceColonizationInterface::karstPathUpdated, this, [&](){ this->viewer->update(); });

    QMetaObject::connectSlotsByName(this);
} //setupBindings

void ViewerInterface::retranslateUi()
{
    this->setWindowTitle(QString("Simulation"));
} // retranslateUi

void ViewerInterface::setAllValuesToFitViewerDefaults(Viewer* viewer)
{

    randomRocksSizeSlider->setValue(viewer->erosionSize);
    randomRocksStrengthSlider->setfValue(viewer->erosionStrength);
    randomRocksQuantitySlider->setValue(viewer->erosionQtt);
//    displayRandomRocks->setChecked(viewer->debugMeshes[ROCK_TRAILS].isDisplayed);
//    displayFailedRandomRocks->setChecked(viewer->debugMeshes[FAILED_ROCKS].isDisplayed);

    currentSimulationFlowfieldStrengthSlider->setfValue(viewer->erosionFlowfieldFactor);
//    manualRocksStrengthSlider->setfValue(viewer->manualErosionStrength);
    currentSimulationRandomDirectionSlider->setfValue(viewer->erosionFlowfieldRandomness);
//    displayFlowfield->setChecked(viewer->debugMeshes[FLOW_TRAILS].isDisplayed);

    manualRocksSizeSlider->setValue(viewer->manualErosionSize);
    manualRocksStrengthSlider->setfValue(viewer->manualErosionStrength);

    curvesErosionSizeSlider->setValue(viewer->curvesErosionSize);
    curvesErosionStrengthSlider->setfValue(viewer->curvesErosionStrength);
//    displayCurvesErosion->setChecked(viewer->debugMeshes[TUNNEL_PATHS].isDisplayed);

//    karstCreationDistanceWeights->setfValue(viewer->karstPathCreator.distanceWeight);
//    karstCreationFractureWeights->setfValue(viewer->karstPathCreator.fractureWeight);
//    karstCreationWaterWeights->setfValue(viewer->karstPathCreator.waterHeightWeight);
//    karstCreationPorosityWeights->setfValue(viewer->karstPathCreator.porosityWeight);
//    karstCreationGamma->setfValue(viewer->karstPathCreator.gamma);
//    karstCreationTortuosity->setfValue(0.f);
//    karstCreationDisplay->setChecked(viewer->debugMeshes[KARST_PATHS].isDisplayed);

//    spaceColonizerDisplay->setChecked(viewer->debugMeshes[SPACE_COLONI].isDisplayed);

    wireModeButton->setChecked(viewer->viewerMode == WIRE_MODE);
    fillModeButton->setChecked(viewer->viewerMode == FILL_MODE);
    invisibleModeButton->setChecked(viewer->viewerMode == NO_DISPLAY);

    noAlgorithmButton->setChecked(viewer->algorithm == SmoothingAlgorithm::NONE);
    marchingCubesButton->setChecked(viewer->algorithm == SmoothingAlgorithm::MARCHING_CUBES);

    gridModeButton->setChecked(viewer->mapMode == MapMode::GRID_MODE);
    voxelsModeButton->setChecked(viewer->mapMode == MapMode::VOXEL_MODE);
    layerModeButton->setChecked(viewer->mapMode == MapMode::LAYER_MODE);
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
        /*this->karstPathGeneration->show();
        this->spaceColonization->show();
        QVBoxLayout* mix = new QVBoxLayout();
        mix->addLayout(this->karstPathGeneration->createGUI());
        mix->addLayout(this->spaceColonization->createGUI());
        this->frame->setContent(mix);*/
        this->spaceColonization->show();
        this->frame->setContent(this->spaceColonization->createGUI());
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

void ViewerInterface::hideAllInteractiveParts()
{
    this->karstPathGeneration->hide();
    this->spaceColonization->hide();
    this->faultSlip->hide();
    this->flowField->hide();
    this->tunnelInterface->hide();
    this->manualEditionInterface->hide();

    this->viewer->update();
    this->update();
}

