#include "ActionInterface.h"

ActionInterface::ActionInterface(std::string actionTypeName,
                                 std::string interfaceName,
                                 std::string interfaceType,
                                 std::string mainActionDescription,
                                 std::string mainActionButtonLogo,  QWidget *parent)
    : CustomInteractiveObject(parent),
      actionType(actionTypeName),
      interfaceName(interfaceName),
      interfaceType(interfaceType),
      mainActionDescription(mainActionDescription),
      mainActionButtonLogo(mainActionButtonLogo)
{
    this->actionType = simplify(actionType);
    this->interfaceType = simplify(interfaceType);
}

void ActionInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid) {
    this->voxelGrid = voxelGrid;
}

void ActionInterface::affectHeightmap(std::shared_ptr<Heightmap> heightmap) {
    this->heightmap = heightmap;
}

void ActionInterface::affectLayerGrid(std::shared_ptr<LayerBasedGrid> layerGrid) {
    this->layerGrid = layerGrid;
}

void ActionInterface::affectImplicitTerrain(std::shared_ptr<ImplicitNaryOperator> implicitPatch) {
    this->implicitTerrain = implicitPatch;
}

void ActionInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch) {
    this->affectHeightmap(heightmap);
    this->affectVoxelGrid(voxelGrid);
    this->affectLayerGrid(layerGrid);
    this->affectImplicitTerrain(implicitPatch);
}

QLayout *ActionInterface::createGUI() {
    return nullptr;
}

void ActionInterface::display(const Vector3 &camPos) {

}

void ActionInterface::reloadShaders() {

}

void ActionInterface::addTerrainAction(nlohmann::json actionParameters) {
    jsonActionsHistory->push_back(nlohmann::json({
                                                     {"type", this->actionType},
                                                     {"parameters", actionParameters}
                                                 }));
    //        saveAllActions();
}

void ActionInterface::affectSavingFile(std::shared_ptr<std::vector<nlohmann::json> > history, std::shared_ptr<std::fstream> file, std::string filename) {
    this->jsonActionsHistory = history;
    this->savingFile = file;
    this->savingFilename = filename;
}

void ActionInterface::saveAllActions(std::string filename) {
    savingFile->close();
    savingFile->open((filename.empty() ? savingFilename : filename), std::fstream::in | std::fstream::out | std::fstream::trunc);
    *savingFile << nlohmann::json({{"actions", *jsonActionsHistory}}).dump(1, '\t');
    savingFile->flush();
}

bool ActionInterface::isConcerned(nlohmann::json &action) { return action.contains("type") && action.at("type").get<std::string>() == this->actionType; }

void ActionInterface::log(const std::string &message) {
    if (this->viewer) {
        std::cout << "[INFO " << this->interfaceName << "]  " << message << std::endl;
//        this->viewer->displayMessage(QString::fromStdString("[INFO]  " + message));
    }
}

void ActionInterface::error(const std::string &message) {
    if (this->viewer) {
        std::cerr << "[ERROR " << this->interfaceName << "]  " << message << std::endl;
//        this->viewer->displayMessage(QString::fromStdString("[ERROR] " + message));
    }
}

std::shared_ptr<ActionInterface> ActionInterface::findOtherInterface(std::string name) const
{
    if (!viewer) {
        std::cerr << "The viewer from this interface (" << this->interfaceName << ") is not defined, cannot fetch interface '" << name << "'..." << std::endl;
        return nullptr;
    }
    if (viewer->interfaces.count(name) == 0) {
        std::cerr << "The viewer has no interface '" << name << "'... Here are all available interfaces :\n";
        for (auto& nameAndInterface : viewer->interfaces)
            std::cerr << "- " << nameAndInterface.first << "\n";
        std::cerr << std::endl;
        return nullptr;
    }
    return viewer->interfaces[name];
}

void ActionInterface::afterTerrainUpdated() {

}

void ActionInterface::afterWaterLevelChanged() {

}

void ActionInterface::mouseClickedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model) {

}

void ActionInterface::mouseDoubleClickedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model) {

}

void ActionInterface::mouseReleasedOnMapEvent(const Vector3 &mouseWorldPosition, bool mouseInMap, QMouseEvent *event, TerrainModel *model) {

}

void ActionInterface::mouseMovedOnMapEvent(const Vector3 &mouseWorldPosition, TerrainModel *model) {

}

//ActionInterface::~ActionInterface()
//{
//}

