#ifndef SPHEROIDALEROSIONINTERFACE_H
#define SPHEROIDALEROSIONINTERFACE_H

#include "ActionInterface.h"
#include "TerrainModification/SpheroidalWeathering.h"

class SpheroidalErosionInterface : public ActionInterface
{
    Q_OBJECT
public:
    SpheroidalErosionInterface(QWidget *parent = nullptr);

    void display(Vector3 camPos = Vector3(false));

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void replay(nlohmann::json action);

    void mouseMoveEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mousePressEvent(QMouseEvent* event);

    QLayout* createGUI();

public Q_SLOTS:
    void hide();
    void show();

    void mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event, TerrainModel *model);

    void applyWeatheringErosion();

    void afterTerrainUpdated();

public:
    SpheroidalWeathering simu;
};

#endif // SPHEROIDALEROSIONINTERFACE_H
