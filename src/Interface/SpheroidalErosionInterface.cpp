#include "SpheroidalErosionInterface.h"

SpheroidalErosionInterface::SpheroidalErosionInterface(QWidget *parent)
    : ActionInterface("spheroidalerosion", "Spheroidal Erosion (Beardall 2010)", "physics", "Spheroidal erosion (Beardall 2010)", "spheroidal_erosion_button.png", parent)
{

}

void SpheroidalErosionInterface::display(const Vector3& camPos)
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

InterfaceUI* SpheroidalErosionInterface::createGUI()
{
    auto UI = new InterfaceUI();
    auto applyButton = new ButtonElement("Apply", [=]() { this->applyWeatheringErosion(); });

    UI->add(applyButton);

    return UI;
}

void SpheroidalErosionInterface::hide()
{
    ActionInterface::hide();
}

void SpheroidalErosionInterface::show()
{
    return ActionInterface::show();
}

void SpheroidalErosionInterface::mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent *event, TerrainModel *model)
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
