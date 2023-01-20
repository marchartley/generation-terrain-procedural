#include "HeightmapErosionInterface.h"

#include "Interface/InterfaceUtils.h"

HeightmapErosionInterface::HeightmapErosionInterface(QWidget *parent)
    : ActionInterface("heightmap-erosion", parent)
{
    hydraulicMesh = Mesh({}, true, GL_LINES);
    hydraulicMesh.useIndices = false;

    windDirectionSelector = std::make_unique<InteractiveVector>();
}

void HeightmapErosionInterface::affectTerrains(std::shared_ptr<Grid> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid);
    this->windDirectionSelector->setPositions(Vector3(0, 0, voxelGrid->getSizeZ()), (windDirection.normalized() * voxelGrid->getSizeX() / 10.f) * Vector3(1, 1, 0) + Vector3(0, 0, voxelGrid->getSizeZ()));
    QObject::connect(windDirectionSelector.get(), &InteractiveVector::modified, [&](Vector3 newVal) { this->windDirection = newVal.normalized() * 2.f; } );
}

/*void HeightmapErosionInterface::affectHeightmap(std::shared_ptr<Grid> heightmap)
{
    this->heightmap = heightmap;
}*/

void HeightmapErosionInterface::display()
{
    if (hydraulicMesh.shader != nullptr)
        hydraulicMesh.shader->setVector("color", std::vector<float>({1.f, .0f, .0f, 1.f}));
    hydraulicMesh.display(GL_LINES);

    windDirectionSelector->display();
}

void HeightmapErosionInterface::replay(nlohmann::json action)
{
//    if (this->isConcerned(action)) {
//        auto& parameters = action.at("parameters");
//        Vector3 pos = json_to_vec3(parameters.at("position"));
//        Vector3 dir = json_to_vec3(parameters.at("direction"));
//        float size = parameters.at("size").get<float>();
//        int qtt = parameters.at("quantity").get<int>();
//        float strength = parameters.at("strength").get<float>();
//        float randomness = parameters.at("randomness").get<float>();
//        UnderwaterErosion erod(this->voxelGrid, size, strength, qtt);
//        erod.Apply(pos, dir, 0.f, 0.f, randomness, false);
//    }
}

void HeightmapErosionInterface::hydraulicErosion()
{
    if (this->heightmap != nullptr) {
        std::vector<std::vector<Vector3>> traces = heightmap->hydraulicErosion(
                    hydraulicNumIterations,
                    hydraulicErosionRadius,
                    hydraulicMaxDropletLifetime,
                    hydraulicErodeSpeed,
                    hydraulicDepositSpeed,
                    hydraulicEvaporateSpeed,
                    hydraulicGravity,
                    hydraulicInertia,
                    hydraulicSedimentCapacityFactor,
                    hydraulicApplyDeposit);
        std::vector<Vector3> segments;
        for (const auto& trace : traces) {
            if (trace.size() < 2) continue;
            for (size_t i = 0; i < trace.size() - 1; i++) {
                segments.push_back(trace[i]);
                segments.push_back(trace[i + 1]);
            }
        }
        hydraulicMesh.fromArray(segments);
        Q_EMIT updated();
    }
}

void HeightmapErosionInterface::thermalErosion()
{
    if (this->heightmap != nullptr) {
        for (int i = 0; i < 20; i++) {
            heightmap->thermalErosion(thermalErosionFactor, thermalMinSlope);
        }
        Q_EMIT updated();
    }
}

void HeightmapErosionInterface::windErosion()
{
    if (this->heightmap != nullptr) {
        std::vector<std::vector<Vector3>> traces = heightmap->windErosion(
                    windNumberOfParticles,
                    windDirection,
                    windBedrocksProportionInGround,
                    windSuspension,
                    windAbrasion,
                    windRoughness,
                    windSettling,
                    windScale,
                    windDt);
        std::vector<Vector3> segments;
        for (const auto& trace : traces) {
            if (trace.size() < 2) continue;
            for (size_t i = 0; i < trace.size() - 1; i++) {
                segments.push_back(trace[i]);
                segments.push_back(trace[i + 1]);
            }
        }
        hydraulicMesh.fromArray(segments);
        Q_EMIT updated();
    }
}

QGroupBox *HeightmapErosionInterface::createHydraulicErosionGUI()
{
    QVBoxLayout* layout = new QVBoxLayout;

    FancySlider* numIterationsSlider = new FancySlider(Qt::Horizontal, 0, 5000);
    FancySlider* erosionRadiusSlider = new FancySlider(Qt::Horizontal, 0, 20);
    FancySlider* dropletLifetimeSlider = new FancySlider(Qt::Horizontal, 0, 100);
    FancySlider* erodingSpeedSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* depositSpeedSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* evaporationSpeedSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* gravitySlider = new FancySlider(Qt::Horizontal, 0, 10, .1f);
    FancySlider* inertiaSlider = new FancySlider(Qt::Horizontal, 0, 1, .001f);
    FancySlider* sedimentCapacityFactorSlider = new FancySlider(Qt::Horizontal, 0, 5, .1f);
    QCheckBox* applyDepositCheckbox = new QCheckBox("With deposit");
    QPushButton* hydraulicErosionButton = new QPushButton("Erosion hydraulique");

    layout->addWidget(createMultipleSliderGroup({
                                  { "Iterations", numIterationsSlider },
                                  { "Rayon", erosionRadiusSlider },
                                  { "Duree", dropletLifetimeSlider},
                                  { "Erosion", erodingSpeedSlider },
                                  { "Depot", depositSpeedSlider },
                                  { "Evaporation", evaporationSpeedSlider },
                                  { "Gravite", gravitySlider },
                                  { "Inertie", inertiaSlider },
                                  { "Capacite", sedimentCapacityFactorSlider }
                              }));
    layout->addWidget(applyDepositCheckbox);
    layout->addWidget(hydraulicErosionButton);

    numIterationsSlider->setfValue(hydraulicNumIterations);
    erosionRadiusSlider->setfValue(hydraulicErosionRadius);
    dropletLifetimeSlider->setfValue(hydraulicMaxDropletLifetime);
    erodingSpeedSlider->setfValue(hydraulicErodeSpeed);
    depositSpeedSlider->setfValue(hydraulicDepositSpeed);
    evaporationSpeedSlider->setfValue(hydraulicEvaporateSpeed);
    gravitySlider->setfValue(hydraulicGravity);
    inertiaSlider->setfValue(hydraulicInertia);
    sedimentCapacityFactorSlider->setfValue(hydraulicSedimentCapacityFactor);
    applyDepositCheckbox->setChecked(hydraulicApplyDeposit);

    QObject::connect(numIterationsSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicNumIterations = val; });
    QObject::connect(erosionRadiusSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicErosionRadius = val; });
    QObject::connect(dropletLifetimeSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicMaxDropletLifetime = val; });
    QObject::connect(erodingSpeedSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicErodeSpeed = val; });
    QObject::connect(depositSpeedSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicDepositSpeed = val; });
    QObject::connect(evaporationSpeedSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicEvaporateSpeed = val; });
    QObject::connect(gravitySlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicGravity = val; });
    QObject::connect(inertiaSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicInertia = val; });
    QObject::connect(sedimentCapacityFactorSlider, &FancySlider::floatValueChanged, [&](float val) { hydraulicSedimentCapacityFactor = val; });
    QObject::connect(applyDepositCheckbox, &QCheckBox::toggled, [&](bool val) { hydraulicApplyDeposit = val; });
    QObject::connect(hydraulicErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::hydraulicErosion);

    QGroupBox* box = new QGroupBox("Erosion Hydraulique");
    box->setLayout(layout);
    return box;
}

QGroupBox *HeightmapErosionInterface::createThermicErosionGUI()
{
    QVBoxLayout* layout = new QVBoxLayout;

    FancySlider* erosionFactorSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* minSlopeSlider = new FancySlider(Qt::Horizontal, 0, 1, .001f);
    QPushButton* thermalErosionButton = new QPushButton("Erosion thermique");

    layout->addWidget(createMultipleSliderGroup({
                                  { "Erosion", erosionFactorSlider },
                                  { "Pente min", minSlopeSlider }
                              }));
    layout->addWidget(thermalErosionButton);

    erosionFactorSlider->setfValue(thermalErosionFactor);
    minSlopeSlider->setfValue(thermalMinSlope);

    QObject::connect(erosionFactorSlider, &FancySlider::floatValueChanged, [&](float val) { thermalErosionFactor = val; });
    QObject::connect(minSlopeSlider, &FancySlider::floatValueChanged, [&](float val) { thermalMinSlope = val; });
    QObject::connect(thermalErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::thermalErosion);

    QGroupBox* box = new QGroupBox("Erosion Thermique");
    box->setLayout(layout);
    return box;
}

QGroupBox *HeightmapErosionInterface::createWindErosionGUI()
{
    QVBoxLayout* layout = new QVBoxLayout;

    FancySlider* numParticlesSlider = new FancySlider(Qt::Horizontal, 0, 5000);
    FancySlider* bedrockSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* suspensionSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* abrasionSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* roughnessSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* settlingSlider = new FancySlider(Qt::Horizontal, 0, 1, .01f);
    FancySlider* scaleSlider = new FancySlider(Qt::Horizontal, 0, 100, .1f);
    FancySlider* dtSlider = new FancySlider(Qt::Horizontal, 0, 1, .001f);
    QPushButton* windErosionButton = new QPushButton("Erosion de vent");

    layout->addWidget(createMultipleSliderGroup({
                                  { "Particules", numParticlesSlider },
                                  { "Ratio roche/sable", bedrockSlider },
                                  { "Suspension", suspensionSlider},
                                  { "Abrasion", abrasionSlider },
                                  { "Rugosite", roughnessSlider },
                                  { "Decantation", settlingSlider },
                                  { "Echelle", scaleSlider },
                                  { "delta-temps", dtSlider },
                              }));
    layout->addWidget(windErosionButton);

    numParticlesSlider->setfValue(windNumberOfParticles);
    bedrockSlider->setfValue(windBedrocksProportionInGround);
    suspensionSlider->setfValue(windSuspension);
    abrasionSlider->setfValue(windAbrasion);
    roughnessSlider->setfValue(windRoughness);
    settlingSlider->setfValue(windSettling);
    scaleSlider->setfValue(windScale);
    dtSlider->setfValue(windDt);

    QObject::connect(numParticlesSlider, &FancySlider::floatValueChanged, [&](float val) { windNumberOfParticles = val; });
    QObject::connect(bedrockSlider, &FancySlider::floatValueChanged, [&](float val) { windBedrocksProportionInGround = val; });
    QObject::connect(suspensionSlider, &FancySlider::floatValueChanged, [&](float val) { windSuspension = val; });
    QObject::connect(abrasionSlider, &FancySlider::floatValueChanged, [&](float val) { windAbrasion = val; });
    QObject::connect(roughnessSlider, &FancySlider::floatValueChanged, [&](float val) { windRoughness = val; });
    QObject::connect(settlingSlider, &FancySlider::floatValueChanged, [&](float val) { windSettling = val; });
    QObject::connect(scaleSlider, &FancySlider::floatValueChanged, [&](float val) { windScale = val; });
    QObject::connect(dtSlider, &FancySlider::floatValueChanged, [&](float val) { windDt = val; });
    QObject::connect(windErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::windErosion);

    QGroupBox* box = new QGroupBox("Erosion Aérienne");
    box->setLayout(layout);
    return box;
}

QLayout *HeightmapErosionInterface::createGUI()
{
//    if (this->erosionLayout != nullptr) return erosionLayout;
    this->erosionLayout = new QHBoxLayout();

//    QPushButton* hydraulicErosionButton = new QPushButton("Erosion hydraulique");
//    QPushButton* thermalErosionButton = new QPushButton("Erosion thermique");
//    QPushButton* windErosionButton = new QPushButton("Erosion de vent");
    erosionLayout->addWidget(createHydraulicErosionGUI());
    erosionLayout->addWidget(createThermicErosionGUI());
    erosionLayout->addWidget(createWindErosionGUI());
    /*createVerticalGroup({
                                                     hydraulicErosionButton,
                                                     thermalErosionButton,
                                                     windErosionButton
                                                 }));*/

//    QObject::connect(thermalErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::thermalErosion);
//    QObject::connect(windErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::windErosion);

    return erosionLayout;
}

void HeightmapErosionInterface::show()
{
    this->hydraulicMesh.show();

    this->windDirectionSelector->show();
    ActionInterface::show();
}

void HeightmapErosionInterface::hide()
{
    this->hydraulicMesh.hide();
    this->windDirectionSelector->hide();
    ActionInterface::hide();
}
