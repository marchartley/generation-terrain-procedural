#include "SpheroidalErosionInterface.h"

SpheroidalErosionInterface::SpheroidalErosionInterface(QWidget *parent)
    : ActionInterface("SpheroidalErosionInterface", parent)
{

}

void SpheroidalErosionInterface::display(Vector3 camPos)
{
    return ActionInterface::display(camPos);
}

void SpheroidalErosionInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
}

void SpheroidalErosionInterface::replay(nlohmann::json action)
{
    // return ActionInterface::replay(action);
}

void SpheroidalErosionInterface::mouseMoveEvent(QMouseEvent *event)
{
    return ActionInterface::mouseMoveEvent(event);
}

void SpheroidalErosionInterface::keyPressEvent(QKeyEvent *event)
{
    return ActionInterface::keyPressEvent(event);
}

void SpheroidalErosionInterface::keyReleaseEvent(QKeyEvent *event)
{
    return ActionInterface::keyReleaseEvent(event);
}

void SpheroidalErosionInterface::wheelEvent(QWheelEvent *event)
{
    return ActionInterface::wheelEvent(event);
}

void SpheroidalErosionInterface::mousePressEvent(QMouseEvent *event)
{
    return ActionInterface::mousePressEvent(event);
}

QLayout *SpheroidalErosionInterface::createGUI()
{
    QLayout* layout = new QVBoxLayout;
    QPushButton* applyButton = new QPushButton("Apply");
//    FancySlider* subsidenceSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
//    RangeSlider* coralLevelsSlider = new RangeSlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
//    FancySlider* verticalScaleSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
//    FancySlider* horizontalScaleSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);
//    FancySlider* alphaSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);

//    layout->addWidget(createMultipleSliderGroup({
//                                                    {"Subsidence", subsidenceSlider},
//                                                    {"Coral", coralLevelsSlider},
//                                                    {"Vertical", verticalScaleSlider},
//                                                    {"Horizontal", horizontalScaleSlider},
//                                                    {"Alpha", alphaSlider}
//                                                }));

    layout->addWidget(applyButton);

//    subsidenceSlider->setfValue(subsidence);
//    coralLevelsSlider->setMinValue(coralLevelMin);
//    coralLevelsSlider->setMaxValue(coralLevelMax);
//    verticalScaleSlider->setfValue(vScale);
//    horizontalScaleSlider->setfValue(hScale);
//    alphaSlider->setfValue(alpha);

    QObject::connect(applyButton, &QPushButton::pressed, this, &SpheroidalErosionInterface::applyWeatheringErosion);
//    QObject::connect(subsidenceSlider, &FancySlider::floatValueChanged, this, &SpheroidalErosionInterface::setSubsidence);
//    QObject::connect(coralLevelsSlider, &RangeSlider::minValueChanged, this, &SpheroidalErosionInterface::setCoralLevelMin);
//    QObject::connect(coralLevelsSlider, &RangeSlider::maxValueChanged, this, &SpheroidalErosionInterface::setCoralLevelMax);
//    QObject::connect(verticalScaleSlider, &FancySlider::floatValueChanged, this, &SpheroidalErosionInterface::setVScale);
//    QObject::connect(horizontalScaleSlider, &FancySlider::floatValueChanged, this, &SpheroidalErosionInterface::setHScale);
//    QObject::connect(alphaSlider, &FancySlider::floatValueChanged, this, &SpheroidalErosionInterface::setAlpha);

    return layout;
}

void SpheroidalErosionInterface::hide()
{
    ActionInterface::hide();
}

void SpheroidalErosionInterface::show()
{
    return ActionInterface::show();
}

void SpheroidalErosionInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
{
    return ActionInterface::mouseClickedOnMapEvent(mousePosInMap, mouseInMap, event, model);
}

void SpheroidalErosionInterface::applyWeatheringErosion()
{
    simu.voxelGrid = voxelGrid;

    simu.applyErosion();

    Q_EMIT updated();
}

void SpheroidalErosionInterface::afterTerrainUpdated()
{

}
