#include "CoralIslandGeneratorInterface.h"
#include "Interface/FancySlider.h"
#include "Interface/InterfaceUtils.h"
#include "Interface/RangeSlider.h"

#include "DataStructure/Image.h"

CoralIslandGeneratorInterface::CoralIslandGeneratorInterface(QWidget *parent)
    : ActionInterface("coralisland", "Coral island generator", "digging", "Create coral island", "coral_generation_button.png", parent)
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

    QPushButton* fromGanButton = new QPushButton("From GAN");

    layout->addWidget(createMultipleSliderGroup({
                                                    {"Subsidence", subsidenceSlider},
                                                    {"Coral", coralLevelsSlider},
//                                                    {"Vertical", verticalScaleSlider},
//                                                    {"Horizontal", horizontalScaleSlider},
                                                    {"Alpha", alphaSlider}
                                                }));

    layout->addWidget(applyButton);
    layout->addWidget(fromGanButton);


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

    QObject::connect(fromGanButton, &QPushButton::pressed, this, &CoralIslandGeneratorInterface::fromGanUI);

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
//    if (this->isVisible()) {
//        this->voxelGrid->setVoxelValues(coralBoulderGen.volume);
//        coralBoulderGen.step();
//        this->voxelGrid->setVoxelValues(coralBoulderGen.volume);
//    }
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
        float waterLevel = voxelGrid->properties->waterLevel;
        this->heightmap->heights = CoralIslandGenerator::generate(this->startingHeightmap.heights, 1.f - subsidence, waterLevel, coralLevelMin * waterLevel, coralLevelMax * waterLevel, vScale, hScale, alpha);
        Q_EMIT updated();
    }
}

void CoralIslandGeneratorInterface::fromGanUI()
{
    std::string path = "Python_tests/test_island_heightmapfeatures/";
    QString q_filename= QString::fromStdString(path + "69.png");  //QFileDialog::getOpenFileName(this, "Open feature map", QString::fromStdString(path), "*", nullptr);
    if (!q_filename.isEmpty()) {
        std::string file = q_filename.toStdString();
        GridV3 img = Image::readFromFile(file).colorImage;

        auto envObjects = CoralIslandGenerator::envObjsFromFeatureMap(img, voxelGrid->getDimensions());
        implicitTerrain->deleteAllChildren();
        for (auto& obj : envObjects)
            implicitTerrain->addChild(obj->createImplicitPatch());
        implicitTerrain->_cached = false;
        voxelGrid->fromImplicit(implicitTerrain.get());
    }
}

void CoralIslandGeneratorInterface::afterTerrainUpdated()
{
    this->startingHeightmap = *this->heightmap;
}

void CoralIslandGeneratorInterface::afterWaterLevelChanged()
{
    if (this->isVisible()) {
        this->updateCoral();
    }
}
