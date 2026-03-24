#include "GravityInterface.h"
#include "GUIElements/InterfaceUtils.h"

// #include "serialization/Serializer.h"

GravityInterface::GravityInterface(QWidget *parent)
    : ActionInterface("gravity", "Gravity", "physics", "Gravity", "gravity_button.png", parent)
{
//    this->createGUI();
}

/*void GravityInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}*/

void GravityInterface::display([[maybe_unused]] const Vector3& camPos)
{

}

void GravityInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        bool applyGlobalGravity = parameters.at("global_gravity").get<bool>();
        bool applySandGravity = parameters.at("sand_gravity").get<bool>();
        float erosionStrength = parameters.at("erosion_strength").get<float>();
        Vector3 currentDirection = parameters.at("current_direction");

        if (applyGlobalGravity) {
            this->voxelGrid->makeItFall();
        }
        if (applySandGravity) {
            this->voxelGrid->letGravityMakeSandFall(false);
        }
    }
}


bool GravityInterface::createGlobalGravity()
{
    this->voxelGrid->makeItFall();

    this->addTerrainAction(nlohmann::json({
                                               {"global_gravity", true},
                                               {"sand_gravity", false},
                                               {"erosion_strength", 0.f},
                                               {"current_direction", Vector3(0, 0, 0)}
                                            }));
    Q_EMIT updated();
    return false;
    /*
    this->startAnimation();
    this->applyLetItFall = !this->applyLetItFall;
//    if (this->applyLetItFall)
        // this->displayMessage( "Gravity is making his job!" );
//    else
        // this->displayMessage( "Gravity stopped caring" );
    update();
    return this->applyLetItFall;*/
}

bool GravityInterface::createSandGravity()
{
    this->voxelGrid->letGravityMakeSandFall(true);
    Q_EMIT updated();


    this->addTerrainAction(nlohmann::json({
                                           {"global_gravity", false},
                                           {"sand_gravity", true},
                                           {"erosion_strength", 0.f},
                                           {"current_direction", Vector3(0, 0, 0)}
                                          }));
    return false;
    /*
    this->startAnimation();
    this->applyLetSandFall = !this->applyLetSandFall;
//    if (this->applyLetSandFall)
        // this->displayMessage( "Sand is falling!" );
//    else
        // this->displayMessage( "Sand stopped falling" );
    update();
    return this->applyLetSandFall;*/
}

void GravityInterface::hide()
{
    CustomInteractiveObject::hide();
}

void GravityInterface::show()
{
    CustomInteractiveObject::show();
}

InterfaceUI* GravityInterface::createGUI()
{
    auto UI = new InterfaceUI();

    auto gravityComputeButton = new ButtonElement("Calculer", [=]() { this->createSandGravity(); });
    auto arrangingLayersButton = new ButtonElement("Rearrange layers");
    auto gravityLayersButton = new ButtonElement("Apply gravity on layers");

    UI->add({
        gravityComputeButton,
        gravityLayersButton,
        arrangingLayersButton
    });


    gravityLayersButton->setOnPressed([&]() {
        this->layerGrid->thermalErosion();
        this->voxelGrid->fromLayerBased(*this->layerGrid, this->voxelGrid->getSizeZ());
        this->heightmap->fromLayerGrid(*this->layerGrid);
        Q_EMIT this->terrainUpdated();
    });
    arrangingLayersButton->setOnPressed([&]() {
        this->layerGrid->reorderLayers();
        this->voxelGrid->fromLayerBased(*this->layerGrid, this->voxelGrid->getSizeZ());
        this->heightmap->fromLayerGrid(*this->layerGrid);
        Q_EMIT this->terrainUpdated();
    });

    return UI;
}
