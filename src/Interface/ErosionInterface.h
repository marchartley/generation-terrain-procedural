#ifndef EROSIONINTERFACE_H
#define EROSIONINTERFACE_H


class ErosionInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/Viewer.h"

#include "Interface/FancySlider.h"

class ErosionInterface : public ActionInterface
{
    Q_OBJECT
public:
    ErosionInterface(QWidget *parent = nullptr);

    void affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid);

    void display();

    void replay(nlohmann::json action);

    QLayout* createGUI();

public Q_SLOTS:
    void show();
    void hide();

    void throwFromCam();
    void throwFromSide();
    void throwFrom(Vector3 pos, Vector3 dir);

public:
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<UnderwaterErosion> erosion;
    Viewer* viewer;

protected:
    Mesh rocksPathSuccess;
    Mesh rocksPathFailure;

    float erosionSize = 2.f;
    float erosionStrength = .5f;
    int erosionQtt = 1000;
    float rockRandomness = .1f;

    QHBoxLayout* erosionLayout = nullptr;
};

#endif // EROSIONINTERFACE_H
