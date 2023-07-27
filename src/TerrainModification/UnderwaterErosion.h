#ifndef UNDERWATEREROSION_H
#define UNDERWATEREROSION_H

#include "TerrainGen/VoxelGrid.h"

#include <tuple>
#include "Utils/BSpline.h"
#include "Karst/KarstHoleProfile.h"
#include "Karst/KarstHole.h"

struct ErosionParticle {
    Vector3 pos;
    Vector3 force;
    Vector3 velocity;
    Vector3 dir; // Normalized velocity

    float density;
    float radius;
    float volume;
    float mass;
    float capacity;
    float maxCapacity;
};

struct EnvironmentProperty {
    Vector3 gravity;
    Vector3 flowfield;
    float density;
    float viscosity;
};

struct MaterialProperty {
    Vector3 normal;
    float density;
    float criticalShearStress;
    float resistance;
};

class UnderwaterErosion
{
public:
    UnderwaterErosion();
//    UnderwaterErosion(std::shared_ptr<VoxelGrid> grid, int maxRockSize, float maxRockStrength, int rockAmount);
    UnderwaterErosion(VoxelGrid* grid, int maxRockSize, float maxRockStrength, int rockAmount);


    enum FLOWFIELD_TYPE { FLUID_SIMULATION, FLOWFIELD_IMAGE, BASIC };
    enum DENSITY_TYPE { NATIVE, DENSITY_IMAGE, UNUSED, RANDOM_DENSITY};
    enum EROSION_APPLIED { DENSITY_VOXELS, HEIGHTMAP, IMPLICIT_TERRAIN, LAYER_TERRAIN, BINARY_VOXELS};

//    std::vector<std::vector<Vector3>> Apply(int avoidMatter = -1);
    // ApplyOn : 0 = density-voxels, 1 = heightmap, 2 = implicit, 3 = layers, 4 = binary-voxels
    std::tuple<std::vector<BSpline>, int, int> Apply(EROSION_APPLIED applyOn,
                                                    float& particleSimulationTime, float& terrainModifTime,
                                                    Vector3 startingPoint = Vector3(false),
                                                    Vector3 originalDirection = Vector3(false),
                                                    float randomnessFactor = 0.05,
                                                    bool fallFromSky = false,
                                                    float gravity = 1.f,
                                                    float bouncingCoefficient = .5f,
                                                    float bounciness = .5f,
                                                    float minSpeed = .1f,
                                                    float maxSpeed = 1.f,
                                                    float maxCapacityFactor = 1.f,
                                                    float erosion = 1.f,
                                                    float deposit = 1.f,
                                                    float matterDensity = 1000.f,
                                                    float materialImpact = 0.f,
                                                    float airFlowfieldRotation = 0.f,
                                                    float waterFlowfieldRotation = 0.f,
                                                    float airForce = 0.f, // 1.f,
                                                    float waterForce = 0.f, //1.f,
                                                    float dt = .1f,
                                                    float shearingStressConstantK = 1.f,
                                                    float shearingRatePower = .5f,
                                                    float erosionPowerValue = 1.f, // From Wojtan
                                                    float criticalShearValue = 10.f,
                                                     std::vector<std::pair<Vector3, Vector3>> posAndDirs = {},
//                                                     std::function<Vector3(Vector3)> flowfieldFunction = nullptr
                                                     FLOWFIELD_TYPE flowType = BASIC,
                                                     Matrix3<Vector3> waterFlow = Matrix3<Vector3>(),
                                                     Matrix3<Vector3> airFlow = Matrix3<Vector3>(),
                                                     DENSITY_TYPE densityUsed = DENSITY_TYPE::UNUSED,
                                                     Matrix3<float> densityMap = Matrix3<float>(),
                                                     float initialCapacity = 0.f,
                                                     FluidSimType fluidSimType = FluidSimType::LBM);
/*
    std::tuple<std::vector<BSpline>, int, int> ApplyOnAnyTerrain(TerrainModel* terrain,
                                                                 float& particleSimulationTime,
                                                                 float& terrainModifTime,
                                                                 Vector3 startingPoint = Vector3(false),
                                                                 Vector3 originalDirection = Vector3(false),
                                                                 float randomnessFactor = 0.05,
                                                                 bool fallFromSky = false,
                                                                 float gravity = 1.f,
                                                                 float bouncingCoefficient = .5f,
                                                                 float bounciness = .5f,
                                                                 float minSpeed = .1f,
                                                                 float maxSpeed = 1.f,
                                                                 float maxCapacityFactor = 1.f,
                                                                 float erosion = 1.f,
                                                                 float deposit = 1.f,
                                                                 float matterDensity = 1000.f,
                                                                 float materialImpact = 0.f,
                                                                 float airFlowfieldRotation = 0.f,
                                                                 float waterFlowfieldRotation = 0.f,
                                                                 float airForce = 0.f, // 1.f,
                                                                 float waterForce = 0.f, //1.f,
                                                                 float dt = .1f,
                                                                 float shearingStressConstantK = 1.f,
                                                                 float shearingRatePower = .5f,
                                                                 float erosionPowerValue = 1.f, // From Wojtan
                                                                 float criticalShearValue = 10.f,
                                                                 std::vector<std::pair<Vector3, Vector3>> posAndDirs = {},
             //                                                     std::function<Vector3(Vector3)> flowfieldFunction = nullptr
                                                                  FLOWFIELD_TYPE flowType = BASIC,
                                                                  Matrix3<Vector3> waterFlow = Matrix3<Vector3>(),
                                                                  Matrix3<Vector3> airFlow = Matrix3<Vector3>(),
                                                                 DENSITY_TYPE densityUsed = DENSITY_TYPE::UNUSED,
                                                                 Matrix3<float> densityMap = Matrix3<float>(),
                                                                 float initialCapacity = 0.f);
*/
    std::vector<Vector3> CreateTunnel(int numberPoints = 2, bool addingMatter = false, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<Vector3> CreateTunnel(BSpline path, bool addingMatter = false, bool usingSpheres = true, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<std::vector<Vector3>> CreateMultipleTunnels(std::vector<BSpline> paths, bool addingMatter = false, bool usingSpheres = true, bool applyChanges = true, KarstHolePredefinedShapes startingShape = SOLUBLE_BED, KarstHolePredefinedShapes endingShape = KEYHOLE);
    std::vector<Vector3> CreateCrack(Vector3 start, Vector3 end, bool applyChanges = true);

    // With a predefined tunnel :
    std::vector<Vector3> CreateTunnel(KarstHole& tunnel, bool addingMatter = false, bool applyChanges = true);


//    std::shared_ptr<VoxelGrid> grid;
    VoxelGrid* voxelGrid;
    Heightmap* heightmap;
    LayerBasedGrid* layerBasedGrid;
    ImplicitPatch* implicitTerrain;

    int maxRockSize;
    int rockAmount;
    float maxRockStrength;

};

#endif // UNDERWATEREROSION_H
