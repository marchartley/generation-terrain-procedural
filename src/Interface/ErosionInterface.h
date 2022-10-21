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
    void throwFromSky();
    void throwFromSide();
    void throwFrom(Vector3 pos, Vector3 dir);

public:
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<UnderwaterErosion> erosion;
    Viewer* viewer;

protected:
    Mesh rocksPathSuccess;
    Mesh rocksPathFailure;

    float erosionSize = 8.f;
    float erosionStrength = .1f;
    int erosionQtt = 500;
    float rockRandomness = .1f;

    float gravity = 1.f;
    float bouncingCoefficient = .5f;
    float bounciness = .5f;
    float minSpeed = .1f;
    float maxSpeed = 1.f;
    float maxCapacityFactor = 1.f;
    float erosionFactor = 1.f;
    float depositFactor = 1.f;
    float matterDensity = 1000.f;
    float materialImpact = 0.f;

    float airFlowfieldRotation = 0.f;
    float waterFlowfieldRotation = 180.f;
    float airForce = 1.f;
    float waterForce = 1.f;

    QHBoxLayout* erosionLayout = nullptr;
};

#endif // EROSIONINTERFACE_H
