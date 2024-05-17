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

    void affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch = nullptr);

    void display(const Vector3& camPos = Vector3(false));

    void replay(nlohmann::json action);

    QLayout* createGUI();

    enum PARTICLE_INITIAL_LOCATION {SKY, RIVER, RANDOM, RIVER2, UNDERWATER, CENTER_TOP, FROM_X, FROM_BIG_X, EVERYWHERE, JUST_ABOVE_VOXELS, VOLCANO, VOLCANO2, VOLCANO3};

public Q_SLOTS:
//    void keyPressEvent(QKeyEvent *e);

    void show();
    void hide();

    void thermalErosionProcess();
    void hydraulicErosionProcess();
    void gobelinsErosionProcess();
    void desertErosionProcess();
    void coastalErosionProcess();
    void fluidSimErosionProcess();
    void volcanoErosionProcess();
    void implicitErosionProcess();
    void heightmapErosionProcess();
    void voxelsErosionProcess();

    void throwFromCam();
    void throwFromSky();
    void throwFromSide();
    void throwFrom(PARTICLE_INITIAL_LOCATION location);

    void testManyManyErosionParameters();

    virtual void afterTerrainUpdated();

    void browseWaterFlowFromFile();
    void browseAirFlowFromFile();
    void browseDensityFieldFromFile();

    void computePredefinedRocksLocations();
    void recomputeAboveVoxelRocksPositions(TerrainModel *terrain, const SpacePartitioning &boundaries);
    void recomputeUnderwaterRocksPositions(TerrainModel *terrain, const SpacePartitioning &boundaries);
    void recomputeRainingPositions(TerrainModel *terrain, const SpacePartitioning &boundaries);

    std::tuple<float, float, float> computeTerrainBoundaries(TerrainModel *terrain, BVHTree *boundariesTree);

    void ErosionParis2019SeaErosion();
    void ErosionParis2019InvasionPercolation();

public:
//    std::shared_ptr<VoxelGrid> voxelGrid;
    std::shared_ptr<UnderwaterErosion> erosion;
    Viewer* viewer;

//protected:
public:
    std::function<Vector3(Vector3)> computeFlowfieldFunction();
    Mesh rocksPathSuccess;
    Mesh rocksPathFailure;
    Mesh boundariesMesh;

    float erosionSize = 8.f;
    float erosionStrength = 0.5f; //.5; // .35f;
    int erosionQtt = 1000;
    float rockRandomness = .1f;

    int maxParticles = 1000;

    float gravity = .981f;
    float bouncingCoefficient = 1.f; // 0.15f; //0.15f;
    float bounciness = 1.f;
    float minSpeed = .1f;
    float maxSpeed = 5.f;
    float maxCapacityFactor = 1.f;
    float erosionFactor = 1.f;
    float depositFactor = 1.f;
    float matterDensity = 500.f; // -1.f; // 500.f;
    float materialImpact = 1.f;

    float airFlowfieldRotation = 0.f; // 270.f;
    float waterFlowfieldRotation = 90.f;
    float airForce = 0.f;
    float waterForce = 0.f;

    float dt = 1.f;

    float shearingStressConstantK = 1.f;
    float shearingRatePower = .5f;
    float erosionPowerValue = 1.f;
    float criticalShearStress = .8f;

    bool continuousRotation = false;
    bool wrapParticles = false;

    int numberOfIterations = 1; // 200;

    int particleMaxCollisions = -1; // 10; // -1 for infinity

    float initialCapacity = .0f;

    EROSION_APPLIED applyOn = EROSION_APPLIED::DENSITY_VOXELS;

    bool displayTrajectories = true;
    bool displayBoundaries = false;

    std::string waterFlowImagePath = "";
    std::string airFlowImagePath = "";
    std::string densityFieldImagePath = "";

    FluidSimType selectedSimulationType = STABLE;

    FLOWFIELD_TYPE flowfieldUsed = FLOWFIELD_TYPE::BASIC;
    DENSITY_TYPE densityUsed = DENSITY_TYPE::NATIVE;

    std::map<PARTICLE_INITIAL_LOCATION, std::vector<std::vector<std::pair<Vector3, Vector3>>>> initialPositionsAndDirections;

    QHBoxLayout* erosionLayout = nullptr;

    UnderwaterErosion erosionProcess;
    bool currentlyModifyingTerrain = false;

    std::vector<std::pair<Vector3, float>> randomObstacles;

    std::vector<int> countsTrajectories;
};

#endif // EROSIONINTERFACE_H
