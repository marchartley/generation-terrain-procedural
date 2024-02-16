#ifndef TerrainComparatorInterface_H
#define TerrainComparatorInterface_H

#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "Interface/ActionInterface.h"

class TerrainComparatorInterface : public ActionInterface
{
    Q_OBJECT
public:
    TerrainComparatorInterface(QWidget *parent = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void replay(nlohmann::json action);

    void mouseMoveEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mousePressEvent(QMouseEvent* event);

    QLayout* createGUI();

    GridV3 extractDifferencesAsImage();

public Q_SLOTS:
    void hide();
    void show();

    void mouseClickedOnMapEvent(const Vector3& mousePosInMap, bool mouseInMap, QMouseEvent* event, TerrainModel *model);

    void afterTerrainUpdated();
    void afterWaterLevelChanged();

    void interpolate(float t);

    void updateStuff();

public:
    VoxelGrid voxelsFromHeightmap;

    GridF displayedVoxels;
    Mesh terrainMesh;

    bool displayUnion = true;
    bool displayIntersection = false;
    bool displaySubstractionAB = false;
    bool displaySubstractionBA = false;

};

#endif // TerrainComparatorInterface_H
