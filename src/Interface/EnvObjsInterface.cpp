#include "EnvObjsInterface.h"

EnvObjsInterface::EnvObjsInterface(QWidget *parent)
    : ActionInterface("envobjects", "Environmental Objects", "model", "Management of environmental objects generation", "", parent)
{

}

void EnvObjsInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{

}

void EnvObjsInterface::display(const Vector3 &camPos)
{

}

void EnvObjsInterface::replay(nlohmann::json action)
{

}

QLayout *EnvObjsInterface::createGUI()
{

}

void EnvObjsInterface::show()
{

}

void EnvObjsInterface::hide()
{

}

void EnvObjsInterface::afterTerrainUpdated()
{

}
