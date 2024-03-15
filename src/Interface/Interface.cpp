#include "Interface/Interface.h"
#include "Interface/InterfaceUtils.h"

#include "Interface/KarstPathGenerationInterface.h"
#include "Interface/SpaceColonizationInterface.h"
#include "Interface/FaultSlipInterface.h"
#include "Interface/StableFluidSimulationInterface.h"
#include "Interface/TunnelInterface.h"
#include "Interface/ManualEditionInterface.h"
#include "Interface/GravityInterface.h"
#include "Interface/UndoRedoInterface.h"
#include "Interface/TerrainGenerationInterface.h"
#include "Interface/ErosionInterface.h"
#include "Interface/HeightmapErosionInterface.h"
#include "Interface/BiomeInterface.h"
#include "Interface/SmoothInterface.h"
#include "Interface/PrimitivePatchesInterface.h"
#include "Interface/TerrainSavingInterface.h"
#include "Interface/MeshInstanceAmplificationInterface.h"
#include "Interface/SPHSimulationInterface.h"
#include "Interface/FLIPSimulationInterface.h"
#include "Interface/WarpFluidSimulationInterface.h"
#include "Interface/LBMFluidSimulationInterface.h"
#include "Interface/CoralIslandGeneratorInterface.h"
#include "Interface/SpheroidalErosionInterface.h"
#include "Interface/EnvObjsInterface.h"
#include "Interface/EnvObjectFluidSimulation.h"
#include "Interface/TerrainComparatorInterface.h"

#include "Interface/CommonInterface.h"

#include "FluidSimulation/OpenFoamParser.h"

ViewerInterface::ViewerInterface() {
    Plotter::init(nullptr, nullptr);
    this->setWindowFlag(Qt::WindowType::WindowMaximizeButtonHint);
    this->setWindowFlag(Qt::WindowType::WindowMinimizeButtonHint);
    this->viewer = new Viewer(this);

    this->actionsOnMap = std::make_shared<std::vector<nlohmann::json>>();
    std::string historyFilename = "MyChanges.json";
    std::shared_ptr<std::fstream> actionsHistoryFile = std::make_shared<std::fstream>(historyFilename);


    std::vector<std::shared_ptr<ActionInterface>> interfacesList = {
        std::make_shared<KarstPathGenerationInterface>(this),
        std::make_shared<SpaceColonizationInterface>(this),
//        std::make_shared<FaultSlipInterface>(this),
        std::make_shared<StableFluidSimulationInterface>(this),
//        std::make_shared<TunnelInterface>(this),
        std::make_shared<ManualEditionInterface>(this),
//        std::make_shared<GravityInterface>(this),
        std::make_shared<UndoRedoInterface>(this),
        std::make_shared<TerrainGenerationInterface>(this),
        std::make_shared<ErosionInterface>(this),
//        std::make_shared<HeightmapErosionInterface>(this),
        std::make_shared<BiomeInterface>(this),
        std::make_shared<SmoothInterface>(this),
        std::make_shared<PrimitivePatchesInterface>(this),
        std::make_shared<TerrainSavingInterface>(this),
//        std::make_shared<MeshInstanceAmplificationInterface>(this),
//        std::make_shared<SPHSimulationInterface>(this),
//        std::make_shared<FLIPSimulationInterface>(this),
//        std::make_shared<WarpFluidSimulationInterface>(this),
//        std::make_shared<LBMFluidSimulationInterface>(this),
//        std::make_shared<CoralIslandGeneratorInterface>(this),
//        std::make_shared<SpheroidalErosionInterface>(this),
        std::make_shared<EnvObjsInterface>(this),
        std::make_shared<EnvObjectFluidSimulation>(this),
        std::make_shared<TerrainComparatorInterface>(this)
    };

    this->actionInterfaces = std::map<std::string, std::shared_ptr<ActionInterface>>();
    for (auto& interf : interfacesList) {
        interf->viewer = viewer;
        this->actionInterfaces[interf->actionType] = interf;
    }

    std::shared_ptr<TerrainGenerationInterface> terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);
//    std::shared_ptr<CoralIslandGeneratorInterface> coralGenerationInterface = std::static_pointer_cast<CoralIslandGeneratorInterface>(actionInterfaces["coralisland"]);
    std::shared_ptr<EnvObjsInterface> envObjectsInterface = std::static_pointer_cast<EnvObjsInterface>(actionInterfaces["envobjects"]);

    viewer->interfaces = this->actionInterfaces;
    for (auto& actionInterface : this->actionInterfaces) {
        actionInterface.second->hide();
        actionInterface.second->affectSavingFile(actionsOnMap, actionsHistoryFile, historyFilename);
    }

    QObject::connect(this->viewer, &Viewer::viewerInitialized, this, [=](){
        envObjectsInterface->setDefinitionFile("saved_maps/primitives.json");

//        terrainGenerationInterface->createTerrainFromNoise(3, 3, 2, 1.0, 0.3);

//        terrainGenerationInterface->createTerrainFromFile("saved_maps/biomes/mayotte.json");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/rock_begin.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/test_openfoam.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/trench.json");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/flat.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/map1.png");

//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmap2.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/slope_with_hole.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/Geometry/_ToClassify/map_2023-08-19__20-31-35-voxels.stl");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/labyrinthe.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/gaussian.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/volcano3_2.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/test.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/Geometry/Pipes/map_2023-08-19__10-23-28-voxels.stl");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/volcano.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope_original.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/coral_base.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/cube.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/CubeTunnel.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/one_slope_noisy_reinforced.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxels/corridor.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/Coral_basis.json");


//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/heightmap.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/heightmaps/new_one_slope.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/river.png");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/voxel_grids/overhang.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/vase.data");
//        terrainGenerationInterface->createTerrainFromFile("saved_maps/trench.json");
        terrainGenerationInterface->createTerrainFromNoise(100, 100, 40, true, 0.f, 0.f);

//        terrainGenerationInterface->prepareShader();
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

                QObject::connect(actionInterface.second.get(), &ActionInterface::waterLevelChanged, this, [&](float newLevel) {
                    QObject::blockSignals(true);
                    otherActionInterface.second->afterWaterLevelChanged();
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

        envObjectsInterface->fromGanUI();
//        float time = timeIt([&]() {
//            std::string simFolder = "OpenFoam/OF_Sim_Marcos/"; //"OpenFOAM/simple/";
//            OpenFoamParser::createSimulationFile(simFolder, viewer->voxelGrid->getVoxelValues().resize(Vector3(10, 10, 5)));
//        });
//        std::cout << "Time for mesh definition: " << showTime(time) << std::endl;

        viewer->setSceneCenter(viewer->voxelGrid->getDimensions() / 2.f);

        QObject::connect(biomeInterface.get(), &BiomeInterface::terrainViewModified, terrainGenerationInterface.get(), &TerrainGenerationInterface::updateDisplayedView);
    });
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

    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, QMenu*>> interfacesToOpen;
    std::map<std::string, QMenu*> menus;
    for (auto& [id, interface] : actionInterfaces) {
        if (interface->interfaceType == "") continue;
        if (menus.count(interface->interfaceType) == 0) {
            menus[interface->interfaceType] = new QMenu(QString::fromStdString(toCapitalize(interface->interfaceType)));
            menu->addMenu(menus[interface->interfaceType]);
        }
        interfacesToOpen.push_back({interface, interface->mainActionButtonLogo, interface->mainActionDescription, menus[interface->interfaceType]});
    }
    auto undoRedoInterface = std::static_pointer_cast<UndoRedoInterface>(actionInterfaces["undoredo"]);
    auto smoothInterface = std::static_pointer_cast<SmoothInterface>(actionInterfaces["smooth"]);
    auto terrainGenerationInterface = std::static_pointer_cast<TerrainGenerationInterface>(actionInterfaces["terraingeneration"]);
    auto erosionInterface = std::static_pointer_cast<ErosionInterface>(actionInterfaces["erosion"]);

    std::vector<std::tuple<std::shared_ptr<ActionInterface>, std::string, std::string, std::string, std::function<void(void)>>> actionsToUse = {
//         Main interface     Button image                        Description                         Menu            Function to call
        {undoRedoInterface,             "undo_button.png",                  "Undo last action",                 "edit",       [=]() { undoRedoInterface->undo(); }},
        {undoRedoInterface,             "redo_button.png",                  "Redo last action",                 "edit",       [=]() { undoRedoInterface->redo(); }},
        {nullptr,                       "recording_button.png",             "Start / stop recording the screen","recording",  [=]() { this->viewer->startStopRecording(); }},
        {nullptr,                       "snapshot.png",                     "Take snapshot of screen",          "recording",  [=]() { this->viewer->screenshot(); }},
        {nullptr,                       "display-marching-cubes_button.png","Use the Marching Cubes algorithm", "view",       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES); } },
        {nullptr,                       "display-voxels_button.png",        "Use raw voxels",                   "view",       [=]() { this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE); } },
        {nullptr,                       "heightmap_button.png",             "Display using heightmap structure","model",      [=]() { this->viewer->setMapMode(MapMode::GRID_MODE); } },
        {nullptr,                       "layerbased_button.png",            "Display using layered structure",  "model",      [=]() { this->viewer->setMapMode(MapMode::LAYER_MODE); } },
        {nullptr,                       "voxelgrid_button.png",             "Display using voxels structure",   "model",      [=]() { this->viewer->setMapMode(MapMode::VOXEL_MODE);} },
        {nullptr,                       "implicit_representation_button.png","Display using implicit structure", "model",      [=]() { this->viewer->setMapMode(MapMode::IMPLICIT_MODE);} },
        {nullptr,                       "no_display.png",                   "Hide map",                         "view",       [=]() { this->viewer->setViewerMode(NO_DISPLAY);} },
        {nullptr,                       "wire_mode.png",                    "Display wireframe",                "view",       [=]() { this->viewer->setViewerMode(WIRE_MODE); } },
        {nullptr,                       "fill_mode.png",                    "Display solid",                    "view",       [=]() { this->viewer->setViewerMode(FILL_MODE); } },
        {smoothInterface,               "smooth_button.png",                "Smooth the terrain",               "model",      [=]() { smoothInterface->applySmooth(); } },
        {terrainGenerationInterface,    "heightmap_to_all_button.png",      "Use heightmap as reference",       "model",      [=]() { terrainGenerationInterface->heightmapToAll(); } },
        {terrainGenerationInterface,    "voxels_to_all_button.png",         "Use voxels as reference",          "model",      [=]() { terrainGenerationInterface->voxelsToAll(); } },
        {terrainGenerationInterface,    "layers_to_all_button.png",         "Use layers as reference",          "model",      [=]() { terrainGenerationInterface->layersToAll(); } },
        {terrainGenerationInterface,    "implicit_to_all_button.png",       "Use implicit as reference",        "model",      [=]() { terrainGenerationInterface->implicitToAll(); } },
        {terrainGenerationInterface,    "reset_button.png",                 "Reload terrain from file",         "model",      [=]() { terrainGenerationInterface->reloadTerrain(this->actionInterfaces); } },
        {terrainGenerationInterface,    "reinforce_borders_button.png",     "Reinforce borders",                "digging",    [=]() { terrainGenerationInterface->reinforceVoxels(); } },
        {nullptr,                       "",                                 "Erosion tests",                    "edit",       [=]() { this->erosionsTests(); } },
        {nullptr,                       "save_erod_depo_button.png",        "Save depo/erod",                   "model",      [=]() { terrainGenerationInterface->saveErosionDepositionTextureMasksOnMultiple(); }},
        {erosionInterface,              "",                                 "Erosion Coast Paris 2019",         "digging",    [=]() { erosionInterface->ErosionParis2019SeaErosion(); }},
        {erosionInterface,              "",                                 "Erosion Perco Paris 2019",         "digging",    [=]() { erosionInterface->ErosionParis2019InvasionPercolation(); }}
    };

    for (auto& [action, logo, descr, menuName, function] : actionsToUse) {
        if (menus.count(menuName) == 0) {
            menus[menuName] = new QMenu(QString::fromStdString(toCapitalize(menuName)));
            menu->addMenu(menus[menuName]);
        }
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

    this->setMenuBar(menu);

    QStatusBar* status = new QStatusBar(this);
    this->setStatusBar(status);

    QDockWidget* displayOptionWidget = new QDockWidget("Affichage", this);
    displayOptionWidget->setFeatures(QDockWidget::DockWidgetFloatable);
    this->addDockWidget(Qt::DockWidgetArea::RightDockWidgetArea, displayOptionWidget);
//    return;

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
    isolevelSelectionActivation->setChecked(true);
    QObject::connect(isolevelSelectionActivation, &QCheckBox::toggled, this, [&](bool active) {
        if (active) { terrainGenerationInterface->minIsoLevel = isolevelSelectionSlider->min_value(); terrainGenerationInterface->maxIsoLevel = isolevelSelectionSlider->max_value();}
        else { terrainGenerationInterface->minIsoLevel = -1000.f; terrainGenerationInterface->maxIsoLevel = 1000.f; }
        this->viewer->update();
    });
    mainLayout->addWidget(viewer, 1, 0);


    QPushButton* reloadShadersButton = new QPushButton("Recharger tous les shaders");

    InterfaceUI* displayOptionUI = new InterfaceUI(new QVBoxLayout, "Options");

    for (int i = 0; i < 4; i++) {
        ButtonElement* viewerSetupExperimental = new ButtonElement("Setup view " + std::to_string(i));
        ButtonElement* viewerSetupExperimental_save = new ButtonElement("Save");
        std::string associatedFilename = "experiments_state" + std::to_string(i) + ".xml";
        viewerSetupExperimental->setOnClick([=]() { this->viewer->setupViewFromFile(associatedFilename); });
        viewerSetupExperimental_save->setOnClick([=]() { this->viewer->saveViewToFile(associatedFilename); });

        displayOptionUI->add(createHorizontalGroupUI({viewerSetupExperimental, viewerSetupExperimental_save}));
    }

    CheckboxElement* displayAsComparisonTerrainButton = new CheckboxElement("Comparison terrain");
    displayAsComparisonTerrainButton->setChecked(terrainGenerationInterface->displayAsComparativeMode);
    displayAsComparisonTerrainButton->setOnChecked([=](bool check) { terrainGenerationInterface->changeDisplayToComparativeMode(check);});

    CheckboxElement* displayWaterDepthButton = new CheckboxElement("Depth");
    displayWaterDepthButton->setChecked(terrainGenerationInterface->displayDepth);
    displayWaterDepthButton->setOnChecked([=](bool check) { terrainGenerationInterface->changeDisplayDepthMode(check);});

    CheckboxElement* displayShadowsButton = new CheckboxElement("Shadows");
    displayShadowsButton->setChecked(terrainGenerationInterface->displayShadows);
    displayShadowsButton->setOnChecked([=](bool check) { terrainGenerationInterface->changeDisplayShadowsMode(check);});

    SliderElement* waterLevelSlider = new SliderElement("Water", 0.f, 1.f, 0.01f, terrainGenerationInterface->waterLevel);
    waterLevelSlider->setOnValueChanged([=](float newValue) { terrainGenerationInterface->setWaterLevel(newValue); });
    SliderElement* ambiantOcclusionSlider = new SliderElement("AO", 0.f, 1.f, 0.01f, terrainGenerationInterface->ambiantOcclusionFactor);
    ambiantOcclusionSlider->setOnValueChanged([=](float newValue) { terrainGenerationInterface->setAmbiantOcclusion(newValue); });
    SliderElement* heightFactorSlider = new SliderElement("Height", 0.f, 4.f, 0.01f, terrainGenerationInterface->heightFactor);
    heightFactorSlider->setOnValueChanged([=](float newValue) { terrainGenerationInterface->setHeightFactor(newValue); });

    displayOptionUI->add(displayModeLayout);
    displayOptionUI->add(createVerticalGroupUI({
                                                createHorizontalGroupUI({
                                                    //createMultipleSliderGroupWithCheckboxUI({
                                                    //    {"Density", isolevelSelectionSlider, isolevelSelectionActivation}
                                                    //}),
                                                    waterLevelSlider,
                                                    ambiantOcclusionSlider,
                                                    heightFactorSlider
                                                }),
                                                //reloadShadersButton,
                                                   createHorizontalGroupUI({
                                                       createVerticalGroupUI({
                                                           displayAsComparisonTerrainButton,
                                                           new CheckboxElement("Animated?", [&](bool checked) { (checked ? viewer->startAnimation() : viewer->stopAnimation()); })
                                                       }),
                                                       createVerticalGroupUI({
                                                           displayShadowsButton,
                                                           displayWaterDepthButton
                                                       })
                                                   })
                                            }));
    displayOptionWidget->setWidget(displayOptionUI->getWidget());


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
    this->update();
}

void ViewerInterface::hideAllInteractiveParts()
{
    for (auto& actionInterface : this->actionInterfaces)
        actionInterface.second->hide();

//    this->viewer->update();
//    this->update();
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

//    EROSION_APPLIED terrainType = this->applyOn; //EROSION_APPLIED::HEIGHTMAP; // 0 = voxels, 1 = heightmap, 2 = implicit, ...
//    if (terrainType == EROSION_APPLIED::DENSITY_VOXELS)
//        viewer->setMapMode(MapMode::VOXEL_MODE);
//    else if (terrainType == EROSION_APPLIED::HEIGHTMAP)
//        viewer->setMapMode(MapMode::GRID_MODE);
//    else if (terrainType == EROSION_APPLIED::LAYER_TERRAIN)
//        viewer->setMapMode(MapMode::LAYER_MODE);
//    else if (terrainType == EROSION_APPLIED::IMPLICIT_TERRAIN)
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

            DENSITY_TYPE densType;
            if (variables["resistance"] == 0) {
                densType = DENSITY_TYPE::NATIVE;
            } else {
                densType = DENSITY_TYPE::RANDOM_DENSITY;
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
