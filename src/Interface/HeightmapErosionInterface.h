#ifndef HEIGHTMAPEROSIONINTERFACE_H
#define HEIGHTMAPEROSIONINTERFACE_H


class HeightmapErosionInterface;
#include <QWidget>
#include "Interface/ActionInterface.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainModification/UnderwaterErosion.h"
#include "Interface/Viewer.h"

#include "Interface/FancySlider.h"
#include "Interface/InteractiveVector.h"

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
    void thermalErosion();
    void windErosion();

public:
    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<Grid> heightmap;

protected:
    Mesh hydraulicMesh;
    QHBoxLayout* erosionLayout = nullptr;

    QGroupBox* createHydraulicErosionGUI();
    QGroupBox* createThermicErosionGUI();
    QGroupBox* createWindErosionGUI();

    // Hydraulic Erosion parameters
    int hydraulicNumIterations = 1000;
    int hydraulicErosionRadius = 10;
    int hydraulicMaxDropletLifetime = 30;
    float hydraulicErodeSpeed = .3f;
    float hydraulicDepositSpeed = .3f;
    float hydraulicEvaporateSpeed = .01f;
    float hydraulicGravity = 4;
    float hydraulicInertia = .05f;
    float hydraulicSedimentCapacityFactor = 1;
    bool hydraulicApplyDeposit = true;

    // Thermal Erosion parameters
    float thermalErosionFactor = .1f;
    float thermalMinSlope = .01f;

    // Wind Erosion parameters
    int windNumberOfParticles = 100;
    Vector3 windDirection = Vector3(2.f, 0.f, 0.f);
    float windBedrocksProportionInGround = .0f;
    float windSuspension = .002f;
    float windAbrasion = .01f;
    float windRoughness = .005f;
    float windSettling = .05f;
    float windScale = 40.f;
    float windDt = .1;

    std::unique_ptr<InteractiveVector> windDirectionSelector;
};

#endif // HEIGHTMAPEROSIONINTERFACE_H
