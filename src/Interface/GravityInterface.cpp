#include "GravityInterface.h"
#include "Interface/InterfaceUtils.h"

GravityInterface::GravityInterface(QWidget *parent) : ActionInterface("gravity", parent)
{
//    this->createGUI();
}

/*void GravityInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}*/

void GravityInterface::display(Vector3 camPos)
{

}

void GravityInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        bool applyGlobalGravity = parameters.at("global_gravity").get<bool>();
        bool applySandGravity = parameters.at("sand_gravity").get<bool>();
        float erosionStrength = parameters.at("erosion_strength").get<float>();
        Vector3 currentDirection = json_to_vec3(parameters.at("current_direction"));

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
                                               {"current_direction", vec3_to_json(Vector3(0, 0, 0))}
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
                                           {"current_direction", vec3_to_json(Vector3(0, 0, 0))}
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

QLayout* GravityInterface::createGUI()
{
    this->gravityLayout = new QHBoxLayout;

    gravityComputeButton = new QPushButton("Calculer");
    QPushButton* arrangingLayersButton = new QPushButton("Rearrange layers");
    QPushButton* gravityLayersButton = new QPushButton("Apply gravity on layers");
//    gravityDisplayButton = new QCheckBox("Afficher");
    gravityLayout->addWidget(createVerticalGroup({
                                                     gravityComputeButton,
                                                     gravityLayersButton,
                                                     arrangingLayersButton
                                                 }));
//    gravityLayout->addWidget(gravityDisplayButton);

//    gravityDisplayButton->setChecked(this->visible);

//    QObject::connect(gravityDisplayButton, &QCheckBox::toggled, this, &GravityInterface::setVisibility);
    QObject::connect(gravityComputeButton, &QPushButton::pressed, this, &GravityInterface::createSandGravity);
    QObject::connect(gravityLayersButton, &QPushButton::pressed, this, [&]() {
        this->layerGrid->thermalErosion();
        this->voxelGrid->fromLayerBased(*this->layerGrid, this->voxelGrid->getSizeZ());
//        voxelGrid->fromCachedData();
        this->heightmap->fromLayerGrid(*this->layerGrid);
        Q_EMIT this->terrainUpdated();
    });
    QObject::connect(arrangingLayersButton, &QPushButton::pressed, this, [&]() {
        this->layerGrid->reorderLayers();
        this->voxelGrid->fromLayerBased(*this->layerGrid, this->voxelGrid->getSizeZ());
//        this->voxelGrid->fromCachedData();
        this->heightmap->fromLayerGrid(*this->layerGrid);
        Q_EMIT this->terrainUpdated();
    });

    return this->gravityLayout;
}
