#include "CoralIslandGeneratorInterface.h"
#include "Interface/FancySlider.h"
#include "Interface/InterfaceUtils.h"
#include "Interface/RangeSlider.h"

CoralIslandGeneratorInterface::CoralIslandGeneratorInterface(QWidget *parent)
    : ActionInterface("coralisland", "Coral island generator", "digging", "Create coral island", "", parent)
{

}

void CoralIslandGeneratorInterface::display(const Vector3& camPos)
{
    return ActionInterface::display(camPos);
}

void CoralIslandGeneratorInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);

    this->startingHeightmap = *this->heightmap;
}

void CoralIslandGeneratorInterface::replay(nlohmann::json action)
{
    // return ActionInterface::replay(action);
}

void CoralIslandGeneratorInterface::mouseMoveEvent(QMouseEvent *event)
{
    return ActionInterface::mouseMoveEvent(event);
}

void CoralIslandGeneratorInterface::keyPressEvent(QKeyEvent *event)
{
    return ActionInterface::keyPressEvent(event);
}

void CoralIslandGeneratorInterface::keyReleaseEvent(QKeyEvent *event)
{
    return ActionInterface::keyReleaseEvent(event);
}

void CoralIslandGeneratorInterface::wheelEvent(QWheelEvent *event)
{
    return ActionInterface::wheelEvent(event);
}

void CoralIslandGeneratorInterface::mousePressEvent(QMouseEvent *event)
{
    return ActionInterface::mousePressEvent(event);
}

QLayout *CoralIslandGeneratorInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout;
    QPushButton* applyButton = new QPushButton("Apply");
    FancySlider* subsidenceSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
    RangeSlider* coralLevelsSlider = new RangeSlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
    FancySlider* verticalScaleSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
    FancySlider* horizontalScaleSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
    FancySlider* alphaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);

    layout->addWidget(createMultipleSliderGroup({
                                                    {"Subsidence", subsidenceSlider},
                                                    {"Coral", coralLevelsSlider},
                                                    {"Vertical", verticalScaleSlider},
                                                    {"Horizontal", horizontalScaleSlider},
                                                    {"Alpha", alphaSlider}
                                                }));

    layout->addWidget(applyButton);

    subsidenceSlider->setfValue(subsidence);
    coralLevelsSlider->setMinValue(coralLevelMin);
    coralLevelsSlider->setMaxValue(coralLevelMax);
    verticalScaleSlider->setfValue(vScale);
    horizontalScaleSlider->setfValue(hScale);
    alphaSlider->setfValue(alpha);

    QObject::connect(applyButton, &QPushButton::pressed, this, &CoralIslandGeneratorInterface::validateTerrainChange);
    QObject::connect(subsidenceSlider, &FancySlider::floatValueChanged, this, &CoralIslandGeneratorInterface::setSubsidence);
    QObject::connect(coralLevelsSlider, &RangeSlider::minValueChanged, this, &CoralIslandGeneratorInterface::setCoralLevelMin);
    QObject::connect(coralLevelsSlider, &RangeSlider::maxValueChanged, this, &CoralIslandGeneratorInterface::setCoralLevelMax);
    QObject::connect(verticalScaleSlider, &FancySlider::floatValueChanged, this, &CoralIslandGeneratorInterface::setVScale);
    QObject::connect(horizontalScaleSlider, &FancySlider::floatValueChanged, this, &CoralIslandGeneratorInterface::setHScale);
    QObject::connect(alphaSlider, &FancySlider::floatValueChanged, this, &CoralIslandGeneratorInterface::setAlpha);

    return layout;
}

void CoralIslandGeneratorInterface::hide()
{
    ActionInterface::hide();
//    if (this->heightmap) {
//        *this->heightmap = this->startingHeightmap;
//    }
}

void CoralIslandGeneratorInterface::show()
{
    if (this->heightmap) {
        this->startingHeightmap = *this->heightmap;
    }

    return ActionInterface::show();
}

void CoralIslandGeneratorInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    if (this->isVisible()) {
        this->voxelGrid->setVoxelValues(coralBoulderGen.volume);
        coralBoulderGen.step();
        this->voxelGrid->setVoxelValues(coralBoulderGen.volume);
    }
    return ActionInterface::mouseClickedOnMapEvent(mousePosInMap, mouseInMap, event, model);
}

void CoralIslandGeneratorInterface::setSubsidence(float newVal)
{
    this->subsidence = newVal;
    updateCoral();
}

void CoralIslandGeneratorInterface::setCoralLevelMin(float newCoral)
{
    this->coralLevelMin = newCoral;
    updateCoral();
}

void CoralIslandGeneratorInterface::setCoralLevelMax(float newCoral)
{
    this->coralLevelMax = newCoral;
    updateCoral();
}

void CoralIslandGeneratorInterface::setVScale(float newVal)
{
    vScale = newVal;
    updateCoral();
}

void CoralIslandGeneratorInterface::setHScale(float newVal)
{
    hScale = newVal;
    updateCoral();
}

void CoralIslandGeneratorInterface::setAlpha(float newVal)
{
    alpha = newVal;
    updateCoral();
}

void CoralIslandGeneratorInterface::validateTerrainChange()
{
    if (this->heightmap) {
        this->startingHeightmap.heights = this->heightmap->heights;
    }
}

void CoralIslandGeneratorInterface::updateCoral()
{
    if (this->heightmap) {
        float waterLevel = voxelGrid->storedWaterLevel;
        this->heightmap->heights = CoralIslandGenerator::generate(this->startingHeightmap.heights, 1.f - subsidence, waterLevel, coralLevelMin/* * waterLevel*/, coralLevelMax/* * waterLevel*/, vScale, hScale, alpha);
        Q_EMIT updated();
    }
}

void CoralIslandGeneratorInterface::afterTerrainUpdated()
{
    this->startingHeightmap = *this->heightmap;
}
