#include "CoralIslandGeneratorInterface.h"
#include "GUIElements/FancySlider.h"
#include "GUIElements/InterfaceUtils.h"
#include "GUIElements/RangeSlider.h"

#include "DataStructure/Image.h"

#include "Interface/EnvObjsInterface.h"

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

void CoralIslandGeneratorInterface::replay([[maybe_unused]] nlohmann::json action)
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

InterfaceUI* CoralIslandGeneratorInterface::createGUI()
{
    auto UI = new InterfaceUI();
    auto applyButton = new ButtonElement("Apply", [=]() { this->validateTerrainChange(); });
    auto subsidenceSlider = new SliderElement("Subsidence", 0.f, 1.f, .01f, subsidence);
    auto coralLevelsSlider = new RangeSliderElement("Coral", 0.f, 1.f, .01f, coralLevelMin, coralLevelMax);
    auto verticalScaleSlider = new SliderElement("Vertical", 0.f, 1.f, .01f, vScale);
    auto horizontalScaleSlider = new SliderElement("Horizontal", 0.f, 1.f, .01f, hScale);
    auto alphaSlider = new SliderElement("Alpha", 0.f, 1.f, .01f, alpha);
    auto fromGanButton = new ButtonElement("From GAN", [=]() { this->fromGanUI(); });

    subsidenceSlider->setOnValueChanged([=](float) { this->updateCoral(); });
    coralLevelsSlider->setOnValueChanged([=](float) { this->updateCoral(); });
    verticalScaleSlider->setOnValueChanged([=](float) { this->updateCoral(); });
    horizontalScaleSlider->setOnValueChanged([=](float) { this->updateCoral(); });
    alphaSlider->setOnValueChanged([=](float) { this->updateCoral(); });

    UI->add(std::vector<UIElement*>{
        applyButton,
        subsidenceSlider,
        coralLevelsSlider,
        verticalScaleSlider,
        horizontalScaleSlider,
        alphaSlider,
        fromGanButton
    });

    return UI;
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

        GridF heights(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 0.f);
        auto scene = dynamic_cast<EnvObjsInterface*>(this->findOtherInterface("envobjects").get())->scene;
        auto envObjects = CoralIslandGenerator::envObjsFromFeatureMap(img, voxelGrid->getDimensions(), scene);
        implicitTerrain->deleteAllChildren();
        for (auto& obj : envObjects)
            implicitTerrain->addChild(obj->createImplicitPatch(heights));
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
