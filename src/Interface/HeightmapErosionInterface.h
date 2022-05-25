#ifndef HEIGHTMAPEROSIONINTERFACE_H
#define HEIGHTMAPEROSIONINTERFACE_H


class HeightmapErosionInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/Viewer.h"

#include "Interface/FancySlider.h"

class HeightmapErosionInterface : public ActionInterface
{
    Q_OBJECT
public:
    HeightmapErosionInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();

    void replay(nlohmann::json action);

    QLayout* createGUI();

Q_SIGNALS:
    void updated();

public Q_SLOTS:
    void show();
    void hide();

    void hydraulicErosion();
    void windErosion();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Grid> heightmap;

protected:
    Mesh hydraulicMesh;
    QHBoxLayout* erosionLayout = nullptr;
};

#endif // HEIGHTMAPEROSIONINTERFACE_H
