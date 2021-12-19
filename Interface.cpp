#include "Interface.h"

FancySlider::FancySlider(float multiplier, QWidget *parent)
        : QSlider(parent), multiplier(multiplier)
{
    connect(this, SIGNAL(valueChanged(int)),
                this, SLOT(notifyValueChanged(int)));
}
FancySlider::FancySlider(Qt::Orientation orientation, float min, float max, float multiplier, QWidget *parent)
        : QSlider(orientation, parent), multiplier(multiplier)
{
    connect(this, SIGNAL(valueChanged(int)),
                this, SLOT(notifyValueChanged(int)));
    this->setRange(std::round(min / multiplier), std::round(max / multiplier));
}
FancySlider::~FancySlider()
{

}
void FancySlider::setfValue(float val)
{
    int a = int(std::round(val / multiplier));
    this->setValue(a);
}
void FancySlider::setfRange(float min, float max)
{
    this->setRange(int(std::round(min / multiplier)), int(std::round(max / multiplier)));
}
void FancySlider::notifyValueChanged(int value) {
    double doubleValue = value * multiplier;
    Q_EMIT floatValueChanged(doubleValue);
}
void FancySlider::sliderChange(SliderChange change)
{
    QSlider::sliderChange(change);

    if (change == QAbstractSlider::SliderValueChange )
    {
        QStyleOptionSlider opt;
        initStyleOption(&opt);

        QRect sr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
        QPoint bottomRightCorner = sr.bottomLeft();

        QToolTip::showText(mapToGlobal( QPoint( bottomRightCorner.x(), bottomRightCorner.y() ) ), QString::number((double)value() * multiplier), this);
    }
}

QGroupBox* createSliderGroup(std::string label, QSlider* slider)
{
    QLabel* lab = new QLabel(QString::fromStdString(label));
    QVBoxLayout* layout = new QVBoxLayout;
    QGroupBox* group = new QGroupBox;
    layout->addWidget(lab);
    layout->addWidget(slider);
    group->setLayout(layout);

//    QObject::connect(slider, SIGNAL(changedValue(int)))
    return group;
}
QGroupBox* createVerticalGroup(std::vector<QWidget*> widgets)
{
    QVBoxLayout* layout = new QVBoxLayout;
    QGroupBox* group = new QGroupBox;
    for (QWidget*& w : widgets)
        layout->addWidget(w);
    group->setLayout(layout);
    return group;
}

ViewerInterface::ViewerInterface() {
    setupUi(this);
}

void ViewerInterface::setupUi(QDialog *Dialog)
{
    this->setWindowFlag(Qt::WindowType::WindowMaximizeButtonHint);
    this->setWindowFlag(Qt::WindowType::WindowMinimizeButtonHint);
    this->viewer = new Viewer();
    this->viewer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // Add the viewer
    mainLayout = new QGridLayout;
    mainLayout->addWidget(viewer, 0, 0);
    mainLayout->setColumnStretch(0, 8);
    mainLayout->setColumnStretch(1, 2);

    controlLayout = new QVBoxLayout;

    // Create the controls layouts
    this->loadSaveLayout = new QHBoxLayout;
    this->loadSaveBox = new QGroupBox("Enregistrer / Charger");
    this->randomRocksLayout = new QHBoxLayout;
    this->randomRocksBox = new QGroupBox("Erosion");
    this->currentSimulationLayout = new QHBoxLayout;
    this->currentSimulationBox = new QGroupBox("Simulation de courants");
    this->manualRocksLayout = new QHBoxLayout;
    this->manualRocksBox = new QGroupBox("Modifs manuelles");
    this->curvesErosionLayout = new QHBoxLayout;
    this->curvesErosionBox = new QGroupBox("Erosion par courbes aléatoires");
    this->gravityLayout = new QHBoxLayout;
    this->gravityBox = new QGroupBox("Gravité");
    this->recordingLayout = new QHBoxLayout;
    this->recordingBox = new QGroupBox("Enregistrements");
    this->displayModeLayout = new QHBoxLayout;
    this->displayModeBox = new QGroupBox("Affichage");
    this->algorithmLayout = new QHBoxLayout;
    this->algorithmBox = new QGroupBox("Algorithme de mesh");
    this->displayTypeLayout = new QHBoxLayout;
    this->displayTypeBox = new QGroupBox("Type de terrain");
    this->LoDChooserLayout = new QHBoxLayout;
    this->LoDChooserBox = new QGroupBox("Niveau de détail");


    // random rocks layout
    this->saveButton = new QPushButton("Sauvegarder...");
    this->loadButton = new QPushButton("Charger...");
    loadSaveLayout->addWidget(saveButton);
    loadSaveLayout->addWidget(loadButton);
    loadSaveBox->setLayout(loadSaveLayout);
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
    randomRocksBox->setLayout(randomRocksLayout);
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
    currentSimulationBox->setLayout(currentSimulationLayout);
    controlLayout->addWidget(currentSimulationBox);


    // Manual rock erosion layout
    this->manualRocksSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->manualRocksStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->addingMode = new QRadioButton("Ajouter de la matière"); addingMode->setChecked(true);
    this->suppressMode = new QRadioButton("Detruire de la matière");
    manualRocksLayout->addWidget(createSliderGroup("Taille", manualRocksSizeSlider));
    manualRocksLayout->addWidget(createSliderGroup("Force", manualRocksStrengthSlider));
    manualRocksLayout->addWidget(createVerticalGroup({addingMode, suppressMode}));
    manualRocksBox->setLayout(manualRocksLayout);
    controlLayout->addWidget(manualRocksBox);


    // Curved erosion layout
    this->curvesAddControlPointButton = new QPushButton("Ajouter un point de control");
    this->curvesErosionSizeSlider = new FancySlider(Qt::Orientation::Horizontal, 1, 200);
    this->curvesErosionStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.0, 3.0, 0.1);
    this->curvesErosionCreateMatter = new QPushButton("Creer un pont");
    this->curvesErosionRemoveMatter = new QPushButton("Creer un tunnel");
    this->displayCurvesErosion = new QCheckBox("Afficher");
    curvesErosionLayout->addWidget(createVerticalGroup({curvesErosionCreateMatter, curvesErosionRemoveMatter}));
    curvesErosionLayout->addWidget(curvesAddControlPointButton);
    curvesErosionLayout->addWidget(createSliderGroup("Taille", curvesErosionSizeSlider));
    curvesErosionLayout->addWidget(createSliderGroup("Force", curvesErosionStrengthSlider));
    curvesErosionLayout->addWidget(displayCurvesErosion);
    curvesErosionBox->setLayout(curvesErosionLayout);
    controlLayout->addWidget(curvesErosionBox);

    // Gravity controls
    this->gravityGlobalButton = new QPushButton("Chutes de pierres");
//    this->gravityGlobalButton->setCheckable(true);
    this->gravitySandButton = new QPushButton("Chutes de sable");
//    this->gravitySandButton->setCheckable(true);
    gravityLayout->addWidget(gravityGlobalButton);
    gravityLayout->addWidget(gravitySandButton);
    gravityBox->setLayout(gravityLayout);
    controlLayout->addWidget(gravityBox);

    // Recording control
    this->startStopRecording = new QPushButton("Enregistrer / stopper l'enregistrement");
    this->startStopRecording->setCheckable(true);
    recordingLayout->addWidget(startStopRecording);
    recordingBox->setLayout(recordingLayout);
    controlLayout->addWidget(recordingBox);

    // Display mode control
    this->wireModeButton = new QRadioButton("Wire mode");
    this->fillModeButton = new QRadioButton("Fill mode"); fillModeButton->setChecked(true);
    displayModeLayout->addWidget(wireModeButton);
    displayModeLayout->addWidget(fillModeButton);
    displayModeBox->setLayout(displayModeLayout);
    controlLayout->addWidget(displayModeBox);

    // Display algorithm controls
    this->noAlgorithmButton = new QRadioButton("Voxels");
    this->marchingCubesButton = new QRadioButton("Marching cubes"); marchingCubesButton->setChecked(true);
    algorithmLayout->addWidget(noAlgorithmButton);
    algorithmLayout->addWidget(marchingCubesButton);
    algorithmBox->setLayout(algorithmLayout);
    controlLayout->addWidget(algorithmBox);

    // Display type controls
    this->gridModeButton = new QRadioButton("Heightmap");
    this->voxelsModeButton = new QRadioButton("3D grid"); voxelsModeButton->setChecked(true);
    this->layerModeButton = new QRadioButton("Layer stacks");
    displayTypeLayout->addWidget(gridModeButton);
    displayTypeLayout->addWidget(voxelsModeButton);
    displayTypeLayout->addWidget(layerModeButton);
    displayTypeBox->setLayout(displayTypeLayout);
    controlLayout->addWidget(displayTypeBox);

    // Display Level of Detail controls
    this->LoDChooserSlider = new FancySlider(Qt::Orientation::Horizontal, 0, this->viewer->getMaxLoDAvailable(), 1);
    this->LoDChooserConfirmButton = new QPushButton("Regénérer");
    LoDChooserLayout->addWidget(createSliderGroup("Niveau de détail", LoDChooserSlider));
    LoDChooserLayout->addWidget(LoDChooserConfirmButton);
    LoDChooserBox->setLayout(LoDChooserLayout);
    controlLayout->addWidget(LoDChooserBox); // TODO : Connect the slider and button

    mainLayout->addLayout(controlLayout, 0, 1);

    Dialog->setLayout(mainLayout);

    this->setAllValuesToFitViewerDefaults(this->viewer);
    this->setupBindings(this);
    this->retranslateUi(this);
} // setupUi

void ViewerInterface::setupBindings(QDialog* Dialog)
{

    QObject::connect(loadButton, &QCheckBox::pressed, this, [=](){this->viewer->loadMapUI(); } );
    QObject::connect(saveButton, &QCheckBox::pressed, this, [=](){this->viewer->saveMapUI(); } );

    QObject::connect(randomRocksSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setErosionRocksSize(int)));
    QObject::connect(randomRocksStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionRocksStrength(float)));
    QObject::connect(randomRocksQuantitySlider, SIGNAL(valueChanged(int)), viewer, SLOT(setErosionRocksQuantity(int)));
    QObject::connect(sendRandomRocksFromCam, &QPushButton::pressed, this, [=](){this->viewer->erodeMap(true); } );
    QObject::connect(sendRandomRocks, &QCheckBox::pressed, this, [=](){this->viewer->erodeMap(false); } );
    QObject::connect(displayRandomRocks, &QCheckBox::toggled, this, [=](bool display){this->viewer->displayRockTrajectories = display; viewer->update(); } );
    QObject::connect(displayFailedRandomRocks, &QCheckBox::toggled, this, [=](bool display){this->viewer->displayFailedRockTrajectories = display; viewer->update(); } );

    QObject::connect(currentSimulationFlowfieldStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionFlowfieldFactor(float)));
//    QObject::connect(currentSimulationStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setManualErosionRocksStrength(float)));
    QObject::connect(currentSimulationRandomDirectionSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setErosionFlowfieldRandomness(float)));
    QObject::connect(displayFlowfield, &QCheckBox::toggled, this, [=](bool display){this->viewer->displayFlowfield = display; viewer->update(); } );
    QObject::connect(recomputeFlowfieldButton, &QPushButton::pressed, this, [=](){this->viewer->recomputeFlowfield(); } );

    QObject::connect(manualRocksSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setManualErosionRocksSize(int)));
    QObject::connect(manualRocksStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setManualErosionRocksStrength(float)));
    QObject::connect(addingMode, &QRadioButton::clicked, this, [=](){this->viewer->setAddingMatterMode(true); } );
    QObject::connect(suppressMode, &QRadioButton::clicked, this, [=](){this->viewer->setAddingMatterMode(false); } );

    QObject::connect(curvesErosionSizeSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setCurvesErosionSize(int)));
    QObject::connect(curvesErosionStrengthSlider, SIGNAL(floatValueChanged(float)), viewer, SLOT(setCurvesErosionStrength(float)));
    QObject::connect(curvesErosionCreateMatter, &QPushButton::pressed, this, [=](){ this->viewer->createTunnel(false); } );
    QObject::connect(curvesErosionRemoveMatter, &QPushButton::pressed, this, [=](){ this->viewer->createTunnel(true); } );
    QObject::connect(curvesAddControlPointButton, &QPushButton::pressed, this, [=](){this->viewer->setCurvesErosionConstructionMode(true); });
    QObject::connect(displayCurvesErosion, &QCheckBox::toggled, this, [=](bool display){ this->viewer->displayTunnelsPath = display; viewer->update(); } );

    QObject::connect(gravityGlobalButton, &QPushButton::pressed, this, [=](){ this->viewer->createGlobalGravity(); } );
    QObject::connect(gravitySandButton  , &QPushButton::pressed, this, [=](){ this->viewer->createSandGravity(); } );

    QObject::connect(startStopRecording  , &QPushButton::pressed, this, [=](){ this->viewer->startStopRecording(); } );

    QObject::connect(wireModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setViewerMode(ViewerMode::WIRE_MODE); } );
    QObject::connect(fillModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setViewerMode(ViewerMode::FILL_MODE); } );

    QObject::connect(noAlgorithmButton, &QRadioButton::clicked, this, [=](){this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::NONE); } );
    QObject::connect(marchingCubesButton, &QRadioButton::clicked, this, [=](){this->viewer->setSmoothingAlgorithm(SmoothingAlgorithm::MARCHING_CUBES); } );

    QObject::connect(gridModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::GRID_MODE); } );
    QObject::connect(voxelsModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::VOXEL_MODE); } );
    QObject::connect(layerModeButton, &QRadioButton::clicked, this, [=](){this->viewer->setMapMode(MapMode::LAYER_MODE); } );

    QObject::connect(LoDChooserSlider, SIGNAL(valueChanged(int)), viewer, SLOT(setLoD(int)));
    QObject::connect(LoDChooserConfirmButton, &QPushButton::pressed, viewer, &Viewer::computeLoD);

    QMetaObject::connectSlotsByName(Dialog);
} //setupBindings

void ViewerInterface::retranslateUi(QDialog *Dialog)
{
    Dialog->setWindowTitle(QCoreApplication::translate("Dialog", "Dialog", nullptr));
} // retranslateUi

void ViewerInterface::setAllValuesToFitViewerDefaults(Viewer* viewer)
{

    randomRocksSizeSlider->setValue(viewer->erosionSize);
    randomRocksStrengthSlider->setfValue(viewer->erosionStrength);
    randomRocksQuantitySlider->setValue(viewer->erosionQtt);
    displayRandomRocks->setChecked(viewer->displayRockTrajectories);
    displayFailedRandomRocks->setChecked(viewer->displayFailedRockTrajectories);

    currentSimulationFlowfieldStrengthSlider->setfValue(viewer->erosionFlowfieldFactor);
//    manualRocksStrengthSlider->setfValue(viewer->manualErosionStrength);
    currentSimulationRandomDirectionSlider->setfValue(viewer->erosionFlowfieldRandomness);
    displayFlowfield->setChecked(viewer->displayFlowfield);

    manualRocksSizeSlider->setValue(viewer->manualErosionSize);
    manualRocksStrengthSlider->setfValue(viewer->manualErosionStrength);

    curvesErosionSizeSlider->setValue(viewer->curvesErosionSize);
    curvesErosionStrengthSlider->setfValue(viewer->curvesErosionStrength);
    displayCurvesErosion->setChecked(viewer->displayTunnelsPath);

    wireModeButton->setChecked(viewer->viewerMode == WIRE_MODE);
    fillModeButton->setChecked(viewer->viewerMode == FILL_MODE);

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
