#ifndef CORALISLANDGENERATORINTERFACE_H
#define CORALISLANDGENERATORINTERFACE_H

#include "TerrainModification/CoralGrowth.h"
#include "TerrainModification/CoralIslandGenerator.h"
#include <QWidget>
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "Interface/ActionInterface.h"

class CoralIslandGeneratorInterface : public ActionInterface
{
    Q_OBJECT
public:
    CoralIslandGeneratorInterface(QWidget *parent = nullptr);

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

    void setSubsidence(float newVal);
    void setCoralLevelMin(float newVal);
    void setCoralLevelMax(float newVal);
    void setVScale(float newVal);
    void setHScale(float newVal);
    void setAlpha(float newVal);

    void validateTerrainChange();

    void updateCoral();

    void afterTerrainUpdated();

public:
    float subsidence = 0.f;
    float waterLevel = 1.f;
    float coralLevelMin = 0.7f;
    float coralLevelMax = 0.9f;

    float vScale = 1.f;
    float hScale = .01f;
    float alpha = .1f;


    CoralIslandGenerator coralGen;

    Heightmap startingHeightmap;

    CoralGrowth coralBoulderGen;

};

#endif // CORALISLANDGENERATORINTERFACE_H
