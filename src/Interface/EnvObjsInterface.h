#ifndef ENVOBJSINTERFACE_H
#define ENVOBJSINTERFACE_H

#include <QWidget>
#include "Interface/ActionInterface.h"

#include "EnvObject/EnvObject.h"

class EnvObjsInterface : public ActionInterface
{
    Q_OBJECT
public:
    EnvObjsInterface(QWidget *parent = nullptr);

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

    std::tuple<GridF, GridV3> extractErosionDataOnTerrain();

    void createEnvObjectsFromImplicitTerrain();

public Q_SLOTS:
    void show();
    void hide();
    virtual void afterTerrainUpdated();

    void instantiateObject();

    void recomputeErosionValues();

    void updateEnvironmentFromEnvObjects();

public:
    Mesh velocitiesMesh;
    Mesh highErosionsMesh;
    Mesh highDepositionMesh;

    bool displayVelocities = true;
    bool displayHighErosions = true;
    bool displaySediments = false;
    bool displayHighCurrents = false;


    GridF erosionGrid;
    GridV3 velocitiesGrid;
};


#endif // ENVOBJSINTERFACE_H
