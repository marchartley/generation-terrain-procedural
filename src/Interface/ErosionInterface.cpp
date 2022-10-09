#include "ErosionInterface.h"

#include "Interface/InterfaceUtils.h"

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
        Vector3 pos = json_to_vec3(parameters.at("position"));
        Vector3 dir = json_to_vec3(parameters.at("direction"));
        float size = parameters.at("size").get<float>();
        int qtt = parameters.at("quantity").get<int>();
        float strength = parameters.at("strength").get<float>();
        float randomness = parameters.at("randomness").get<float>();
        UnderwaterErosion erod(this->voxelGrid, size, strength, qtt);
        erod.Apply(pos, dir, 0.f, 0.f, randomness, false);
    }
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
    std::tie(lastRocksLaunched, lastFailedRocksLaunched) = erod.Apply(pos, dir, 0.f, 0.f, this->rockRandomness, true);
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
    QPushButton* confirmButton = new QPushButton("Envoyer");
    QPushButton* confirmFromCamButton = new QPushButton("Camera");

    erosionLayout->addWidget(createMultipleSliderGroup({
                                                           {"Taille", rockSizeSlider},
                                                           {"Strength", rockStrengthSlider},
                                                           {"Quantity", rockQttSlider},
                                                           {"Randomness", rockRandomnessSlider}
                                                       }));
    erosionLayout->addWidget(confirmButton);
    erosionLayout->addWidget(confirmFromCamButton);

    QObject::connect(rockSizeSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionSize = newVal; });
    QObject::connect(rockStrengthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionStrength = newVal; });
    QObject::connect(rockQttSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionQtt = newVal; });
    QObject::connect(rockRandomnessSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->rockRandomness = newVal; });
    QObject::connect(confirmButton, &QPushButton::pressed, this, &ErosionInterface::throwFromSide);
    QObject::connect(confirmFromCamButton, &QPushButton::pressed, this, &ErosionInterface::throwFromCam);

    rockSizeSlider->setfValue(this->erosionSize);
    rockStrengthSlider->setfValue(this->erosionStrength);
    rockQttSlider->setfValue(this->erosionQtt);
    rockRandomnessSlider->setfValue(this->rockRandomness);

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
