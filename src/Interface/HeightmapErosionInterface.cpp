#include "HeightmapErosionInterface.h"

#include "GUIElements/InterfaceUtils.h"

HeightmapErosionInterface::HeightmapErosionInterface(QWidget *parent)
    : ActionInterface("heightmaperosion", "Erosion on heightmaps", "physics", "Erosion on heightmap", "heightmap-erosion.png", parent)
{
    hydraulicMesh = Mesh({}, true, GL_LINES);
    hydraulicMesh.useIndices = false;

    windDirectionSelector = std::make_shared<InteractiveVector>();
}

void HeightmapErosionInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    this->windDirectionSelector->setPositions(Vector3(0, 0, voxelGrid->getSizeZ()), (windDirection.normalized() * voxelGrid->getSizeX() / 10.f) * Vector3(1, 1, 0) + Vector3(0, 0, voxelGrid->getSizeZ()));
    QObject::connect(windDirectionSelector.get(), &InteractiveVector::modified, [&](const Vector3& newVal) { this->windDirection = newVal.normalized() * 2.f; } );
}

/*void HeightmapErosionInterface::affectHeightmap(std::shared_ptr<Grid> heightmap)
{
    this->heightmap = heightmap;
}*/

void HeightmapErosionInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

    if (hydraulicMesh.shader != nullptr)
        hydraulicMesh.shader->setVector("color", std::vector<float>({1.f, .0f, .0f, 1.f}));
    hydraulicMesh.display(GL_LINES);

    windDirectionSelector->display();
}

void HeightmapErosionInterface::replay(nlohmann::json action)
{
//    if (this->isConcerned(action)) {
//        auto& parameters = action.at("parameters");
//        Vector3 pos = parameters.at("position");
//        Vector3 dir = parameters.at("direction");
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

InterfaceUI* HeightmapErosionInterface::createHydraulicErosionGUI()
{
    auto UI = new InterfaceUI();

    auto numIterationsSlider = new SliderElement("Iterations", 0, 5000, hydraulicNumIterations);
    auto erosionRadiusSlider = new SliderElement("Rayon", 0, 20, hydraulicErosionRadius);
    auto dropletLifetimeSlider = new SliderElement("Duree", 0, 100, hydraulicMaxDropletLifetime);
    auto erodingSpeedSlider = new SliderElement("Erosion", 0, 1, .01f, hydraulicErodeSpeed);
    auto depositSpeedSlider = new SliderElement("Depot", 0, 1, .01f, hydraulicDepositSpeed);
    auto evaporationSpeedSlider = new SliderElement("Evaporation", 0, 1, .01f, hydraulicEvaporateSpeed);
    auto gravitySlider = new SliderElement("Gravite", 0, 10, .1f, hydraulicGravity);
    auto inertiaSlider = new SliderElement("Inertie", 0, 1, .001f, hydraulicInertia);
    auto sedimentCapacityFactorSlider = new SliderElement("Capacite", 0, 5, .1f, hydraulicSedimentCapacityFactor);
    auto applyDepositCheckbox = new CheckboxElement("With deposit", hydraulicApplyDeposit);
    auto hydraulicErosionButton = new ButtonElement("Erosion hydraulique", [=]() { this->hydraulicErosion(); });


    UI->add(std::vector<UIElement*>{
        numIterationsSlider,
        erosionRadiusSlider,
        dropletLifetimeSlider,
        erodingSpeedSlider,
        depositSpeedSlider,
        evaporationSpeedSlider,
        gravitySlider,
        inertiaSlider,
        sedimentCapacityFactorSlider,
        applyDepositCheckbox,
        hydraulicErosionButton
    });


    return UI;
}

InterfaceUI* HeightmapErosionInterface::createThermicErosionGUI()
{
    auto UI = new InterfaceUI();

    auto erosionFactorSlider = new SliderElement("Erosion", 0, 1, .01f, thermalErosionFactor);
    auto minSlopeSlider = new SliderElement("Pente Min", 0, 1, .001f, thermalMinSlope);
    auto thermalErosionButton = new ButtonElement("Erosion thermique", [=]() { this->thermalErosion(); });


    UI->add(std::vector<UIElement*>{
        erosionFactorSlider,
        minSlopeSlider,
        thermalErosionButton
    });


    return UI;
}

InterfaceUI* HeightmapErosionInterface::createWindErosionGUI()
{
    auto UI = new InterfaceUI();

    auto numParticlesSlider = new SliderElement("Particules", 0, 5000, windNumberOfParticles);
    auto bedrockSlider = new SliderElement("Ratio roche/sable", 0, 1, .01f, windBedrocksProportionInGround);
    auto suspensionSlider = new SliderElement("Suspension", 0, 1, .01f, windSuspension);
    auto abrasionSlider = new SliderElement("Abrasion", 0, 1, .01f, windAbrasion);
    auto roughnessSlider = new SliderElement("Rugosite", 0, 1, .01f, windRoughness);
    auto settlingSlider = new SliderElement("Decantation", 0, 1, .01f, windSettling);
    auto scaleSlider = new SliderElement("Echelle", 0, 100, .1f, windScale);
    auto dtSlider = new SliderElement("delta-temps", 0, 1, .001f, windDt);
    auto windErosionButton = new ButtonElement("Erosion de vent", [=]() { this->windErosion(); });


    UI->add(std::vector<UIElement*>{
            numParticlesSlider,
            bedrockSlider,
            suspensionSlider,
            abrasionSlider,
            roughnessSlider,
            settlingSlider,
            scaleSlider,
            dtSlider,
            windErosionButton
    });

    return UI;
}

InterfaceUI* HeightmapErosionInterface::createGUI()
{
//    if (this->erosionLayout != nullptr) return erosionLayout;
    auto UI = new InterfaceUI();

//    QPushButton* hydraulicErosionButton = new QPushButton("Erosion hydraulique");
//    QPushButton* thermalErosionButton = new QPushButton("Erosion thermique");
//    QPushButton* windErosionButton = new QPushButton("Erosion de vent");
    auto hydrauUI = createHydraulicErosionGUI();
    auto thermalUI = createThermicErosionGUI();
    auto windUI = createWindErosionGUI();
    // UI->add(std::vector<UIElement*>{
        // std::move(hydrauUI),
        // std::move(thermalUI),
        // std::move(windUI)
    // });
    UI->add(std::vector<UIElement*>{
         hydrauUI,
         thermalUI,
         windUI
    });

//    QObject::connect(thermalErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::thermalErosion);
//    QObject::connect(windErosionButton, &QPushButton::pressed, this, &HeightmapErosionInterface::windErosion);

    return UI;
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
