#include "ErosionInterface.h"

#include "Interface/InterfaceUtils.h"

#include <chrono>

ErosionInterface::ErosionInterface(QWidget *parent)
    : ActionInterface("rock-throwing", parent)
{

}

void ErosionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
    this->erosion = std::make_shared<UnderwaterErosion>(voxelGrid, erosionSize, erosionStrength, erosionQtt);

    const char* vNoShader = ":/src/Shaders/no_shader.vert";
    const char* fNoShader = ":/src/Shaders/no_shader.frag";

    this->rocksPathFailure = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->rocksPathFailure.shader->setVector("color", std::vector<float>({.7f, .2f, .1f, .5f}));
    this->rocksPathSuccess = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->rocksPathSuccess.shader->setVector("color", std::vector<float>({.1f, .7f, .2f, .5f}));
}

void ErosionInterface::display()
{
    this->rocksPathSuccess.display(GL_LINES, 3.f);
    this->rocksPathFailure.display(GL_LINES, 3.f);
}

void ErosionInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        Vector3 pos = json_to_vec3(parameters.at("position")) + Vector3::random(0.f, 20.f);
        Vector3 dir = json_to_vec3(parameters.at("direction")) + Vector3::random();
        float size = parameters.at("size").get<float>() + random_gen::generate(0.f, 3.f);
        int qtt = parameters.at("quantity").get<int>() + random_gen::generate(0.f, 100.f);
        float strength = parameters.at("strength").get<float>() + random_gen::generate(0.f, 1.f);
        float randomness = parameters.at("randomness").get<float>() + random_gen::generate(0.f, .1f);
        UnderwaterErosion erod(this->voxelGrid, size, strength, qtt);
        erod.Apply(pos, dir, randomness);
    }
}

void ErosionInterface::throwFromSky()
{
    UnderwaterErosion erod = UnderwaterErosion(voxelGrid, this->erosionSize, this->erosionStrength, this->erosionQtt);

    std::vector<std::vector<Vector3>> lastRocksLaunched, lastFailedRocksLaunched;
    std::tie(lastRocksLaunched, lastFailedRocksLaunched) = erod.Apply(Vector3(false), Vector3(false), this->rockRandomness, true,
                                                                      gravity,
                                                                      bouncingCoefficient,
                                                                      bounciness,
                                                                      minSpeed,
                                                                      maxSpeed,
                                                                      maxCapacityFactor,
                                                                      erosionFactor,
                                                                      depositFactor,
                                                                      matterDensity + .1f,
                                                                      materialImpact,
                                                                      airFlowfieldRotation,
                                                                      waterFlowfieldRotation,
                                                                      airForce,
                                                                      waterForce
                                                                      );

    std::vector<Vector3> asOneVector;
    for(std::vector<Vector3>& points : lastRocksLaunched) {
        BSpline path = BSpline(points).simplifyByRamerDouglasPeucker(0.1f);
        for (size_t i = 0; i < path.points.size() - 1; i++) {
            asOneVector.push_back(path.points[i]);
            asOneVector.push_back(path.points[i + 1]);
        }
    }
    this->rocksPathSuccess.fromArray(asOneVector);
    asOneVector.clear();
    for(std::vector<Vector3>& points : lastFailedRocksLaunched) {
        BSpline path = BSpline(points).simplifyByRamerDouglasPeucker(0.1f);
        for (size_t i = 0; i < path.points.size() - 1; i++) {
            asOneVector.push_back(path.points[i]);
            asOneVector.push_back(path.points[i + 1]);
        }
    }
    this->rocksPathFailure.fromArray(asOneVector);

    Q_EMIT this->updated();
}
void ErosionInterface::throwFromCam()
{
    Vector3 pos;
    Vector3 dir;
    pos = this->viewer->camera()->position();
    dir = this->viewer->camera()->viewDirection();
    this->throwFrom(pos, dir);
}
void ErosionInterface::throwFromSide()
{
    Vector3 pos;
    Vector3 dir;
    pos = this->viewer->camera()->position();
    dir = this->viewer->camera()->viewDirection();

    pos = Vector3(this->voxelGrid->getSizeX() * -1.f, this->voxelGrid->getSizeY() * .5f, this->voxelGrid->getSizeZ() * .5f);
    dir = Vector3(1.f, 0.f, 0.f);

    this->throwFrom(pos, dir);
}
void ErosionInterface::throwFrom(Vector3 pos, Vector3 dir)
{
    UnderwaterErosion erod = UnderwaterErosion(voxelGrid, this->erosionSize, this->erosionStrength, this->erosionQtt);

    std::vector<std::vector<Vector3>> lastRocksLaunched, lastFailedRocksLaunched;
    std::tie(lastRocksLaunched, lastFailedRocksLaunched) = erod.Apply(pos, dir, this->rockRandomness, false,
                                                                      gravity,
                                                                      bouncingCoefficient,
                                                                      bounciness,
                                                                      minSpeed,
                                                                      maxSpeed,
                                                                      maxCapacityFactor,
                                                                      erosionFactor,
                                                                      depositFactor,
                                                                      matterDensity + 0.1f,
                                                                      materialImpact,
                                                                      airFlowfieldRotation,
                                                                      waterFlowfieldRotation,
                                                                      airForce,
                                                                      waterForce
                                                                      );
    std::vector<Vector3> asOneVector;
    for(std::vector<Vector3>& points : lastRocksLaunched) {
        BSpline path = BSpline(points).simplifyByRamerDouglasPeucker(0.1f);
        for (size_t i = 0; i < path.points.size() - 1; i++) {
            asOneVector.push_back(path.points[i]);
            asOneVector.push_back(path.points[i + 1]);
        }
    }
    this->rocksPathSuccess.fromArray(asOneVector);
    asOneVector.clear();
    for(std::vector<Vector3>& points : lastFailedRocksLaunched) {
        BSpline path = BSpline(points).simplifyByRamerDouglasPeucker(0.1f);
        for (size_t i = 0; i < path.points.size() - 1; i++) {
            asOneVector.push_back(path.points[i]);
            asOneVector.push_back(path.points[i + 1]);
        }
    }
    this->rocksPathFailure.fromArray(asOneVector);

    this->addTerrainAction(nlohmann::json({
                                              {"position", vec3_to_json(pos) },
                                              {"direction", vec3_to_json(dir) },
                                              {"size", erosionSize},
                                              {"strength", erosionStrength},
                                              {"randomness", rockRandomness},
                                              {"quantity", erosionQtt}
                                          }));

    Q_EMIT this->updated();
}

QLayout *ErosionInterface::createGUI()
{
    this->erosionLayout = new QHBoxLayout();

    FancySlider* rockSizeSlider = new FancySlider(Qt::Horizontal, 0.f, 100.f);
    FancySlider* rockStrengthSlider = new FancySlider(Qt::Horizontal, 0.f, 3.f, .1f);
    FancySlider* rockQttSlider = new FancySlider(Qt::Horizontal, 0.f, 1000.f);
    FancySlider* rockRandomnessSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* gravitySlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* bouncingCoefficientSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* bouncinessSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* minSpeedSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* maxSpeedSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* maxCapacityFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 10.f, .01f);
    FancySlider* erosionFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 5.f, .01f);
    FancySlider* depositFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 5.f, .01f);
    FancySlider* matterDensitySlider = new FancySlider(Qt::Horizontal, 0.f, 2000.f, .25f);
    FancySlider* materialImpactSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);

    FancySlider* airFlowfieldRotationSlider = new FancySlider(Qt::Horizontal, 0.f, 360.f, 45.f);
    FancySlider* waterFlowfieldRotationSlider = new FancySlider(Qt::Horizontal, 0.f, 360.f, 45.f);
    FancySlider* airForceSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .1f);
    FancySlider* waterForceSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .1f);

    QPushButton* confirmButton = new QPushButton("Envoyer");
    QPushButton* confirmFromCamButton = new QPushButton("Camera");
    QPushButton* confirmFromSkyButton = new QPushButton("Pluie");

    erosionLayout->addWidget(createMultipleSliderGroup({
                                                           {"Taille", rockSizeSlider},
                                                           {"Strength", rockStrengthSlider},
                                                           {"Quantity", rockQttSlider},
                                                           {"gravitySlider", gravitySlider},
                                                           {"bouncingCoefficientSlider", bouncingCoefficientSlider},
                                                           {"bouncinessSlider", bouncinessSlider},
                                                           {"minSpeedSlider", minSpeedSlider},
                                                           {"maxSpeedSlider", maxSpeedSlider},
                                                           {"maxCapacityFactorSlider", maxCapacityFactorSlider},
                                                           {"erosionFactorSlider", erosionFactorSlider},
                                                           {"depositFactorSlider", depositFactorSlider},
                                                           {"matterDensitySlider", matterDensitySlider},
                                                           {"materialImpactSlider", materialImpactSlider},
                                                           {"airRotation", airFlowfieldRotationSlider},
                                                           {"waterRotation", waterFlowfieldRotationSlider},
                                                           {"airForce", airForceSlider},
                                                           {"waterForce", waterForceSlider}
                                                       }));
    erosionLayout->addWidget(confirmButton);
    erosionLayout->addWidget(confirmFromCamButton);
    erosionLayout->addWidget(confirmFromSkyButton);

    QObject::connect(rockSizeSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionSize = newVal; });
    QObject::connect(rockStrengthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionStrength = newVal; });
    QObject::connect(rockQttSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionQtt = newVal; });
    QObject::connect(rockRandomnessSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->rockRandomness = newVal; });

    QObject::connect(gravitySlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->gravity = newVal; });
    QObject::connect(bouncingCoefficientSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->bouncingCoefficient = newVal; });
    QObject::connect(bouncinessSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->bounciness = newVal; });
    QObject::connect(minSpeedSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->minSpeed = newVal; });
    QObject::connect(maxSpeedSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->maxSpeed = newVal; });
    QObject::connect(maxCapacityFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->maxCapacityFactor = newVal; });
    QObject::connect(erosionFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionFactor = newVal; });
    QObject::connect(depositFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->depositFactor = newVal; });
    QObject::connect(matterDensitySlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->matterDensity = newVal; });
    QObject::connect(materialImpactSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->materialImpact = newVal; });
    QObject::connect(airFlowfieldRotationSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->airFlowfieldRotation = newVal; });
    QObject::connect(waterFlowfieldRotationSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->waterFlowfieldRotation = newVal; });
    QObject::connect(airForceSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->airForce = newVal; });
    QObject::connect(waterForceSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->waterForce = newVal; });
    QObject::connect(confirmButton, &QPushButton::pressed, this, &ErosionInterface::throwFromSide);
    QObject::connect(confirmFromCamButton, &QPushButton::pressed, this, &ErosionInterface::throwFromCam);
    QObject::connect(confirmFromSkyButton, &QPushButton::pressed, this, &ErosionInterface::throwFromSky);

    rockSizeSlider->setfValue(this->erosionSize);
    rockStrengthSlider->setfValue(this->erosionStrength);
    rockQttSlider->setfValue(this->erosionQtt);
    rockRandomnessSlider->setfValue(this->rockRandomness);
    gravitySlider->setfValue(this->gravity);
    bouncingCoefficientSlider->setfValue(this->bouncingCoefficient);
    bouncinessSlider->setfValue(this->bounciness);
    minSpeedSlider->setfValue(this->minSpeed);
    maxSpeedSlider->setfValue(this->maxSpeed);
    maxCapacityFactorSlider->setfValue(this->maxCapacityFactor);
    erosionFactorSlider->setfValue(this->erosionFactor);
    depositFactorSlider->setfValue(this->depositFactor);
    matterDensitySlider->setfValue(this->matterDensity);
    materialImpactSlider->setfValue(this->materialImpact);
    airFlowfieldRotationSlider->setfValue(this->airFlowfieldRotation);
    waterFlowfieldRotationSlider->setfValue(this->waterFlowfieldRotation);
    airForceSlider->setfValue(this->airForce);
    waterForceSlider->setfValue(this->waterForce);

    return erosionLayout;
}

void ErosionInterface::show()
{
    this->rocksPathSuccess.show();
    this->rocksPathFailure.show();
    ActionInterface::show();
}

void ErosionInterface::hide()
{
    this->rocksPathSuccess.hide();
    this->rocksPathFailure.hide();
    ActionInterface::hide();
}
