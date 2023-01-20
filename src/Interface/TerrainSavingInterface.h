#ifndef TERRAINSAVINGINTERFACE_H
#define TERRAINSAVINGINTERFACE_H

#include "ActionInterface.h"

class TerrainSavingInterface : public ActionInterface
{
    Q_OBJECT
public:
    TerrainSavingInterface(QWidget *parent = nullptr);

    //    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display(); // Nothing to display
    void replay(nlohmann::json action); // Nothing to replay


    QLayout* createGUI();

public Q_SLOTS:
    void saveTerrainGeometry(std::string filename = "");

public:
    std::string mainFilename = "saved_maps/Geometry/map";
    bool saveHeightmap = true;
    bool saveVoxels = true;
    bool saveLayers = true;
};

#endif // TERRAINSAVINGINTERFACE_H
