#ifndef PARTICLEEROSION_H
#define PARTICLEEROSION_H

class ParticleErosion;

#include "DataStructure/Vector3.h"
#include "DataStructure/Matrix3.h"
#include "TerrainGen/Heightmap.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/ImplicitPatch.h"
#include "EnvObject/EnvObject.h"

struct ParticleProperties {
    float density;
    float radius;
    float volume;
    float mass;
    float maxCapacity;
};

class ErosionParticle {
public:
    Vector3 pos;
    Vector3 force;
    Vector3 velocity;
    Vector3 dir;

//    float density;
//    float radius;
//    float volume;
//    float mass;
    float capacity;
//    float maxCapacity;
    std::shared_ptr<ParticleProperties> properties;

    bool rolling = false;
    bool justStartedRolling = false;


    void update();
    Vector3 predictNextPos(float dt);
    std::vector<Vector3> bounce(float dt, SpacePartitioning* boundaries);

    void addForce(const Vector3& f);
    void addVelocity(const Vector3& v);

    Vector3 eulerIntegration(float dt);
    Vector3 verletIntegration(float dt);
};


struct ErosionPoint {
    Vector3 position;
    float erosionValue;
//    Vector3 velocity;
};

struct ParticleState {
    ParticleState() {}
    ParticleState(ErosionParticle& p) : particleState(p) {}
    ErosionParticle particleState;
    bool terrainModificationApplied = false;
};

struct ParticleHistory {
    void add(const ParticleState& state) { history.push_back(state); }

    std::vector<ParticleState> history;

    std::vector<Vector3> getPositions();
    std::vector<float> getCapacities();
    std::vector<Vector3> getErosionPositions();
    std::vector<float> getErosionValues();

    std::vector<ErosionPoint> getErosionPoints();
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


enum FLOWFIELD_TYPE { FLUID_SIMULATION, FLOWFIELD_IMAGE, BASIC, FLOWFIELD_ENVOBJECTS };
enum DENSITY_TYPE { NATIVE, DENSITY_IMAGE, UNUSED, RANDOM_DENSITY, LAYERED_DENSITY};
enum EROSION_APPLIED { DENSITY_VOXELS, HEIGHTMAP, IMPLICIT_TERRAIN, LAYER_TERRAIN, BINARY_VOXELS};


#include "TerrainModification/UnderwaterErosion.h"


class ParticleErosion {
public:
    ParticleErosion();

    std::vector<ParticleHistory> process();
    void initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity);
    ParticleHistory trackParticlePositions(ErosionParticle& particle);
    float computeErosionValue(const ErosionParticle &particle, float vRel);
    float computeDepositionValue(const ErosionParticle& particle);

    std::vector<ErosionPoint> computeErosionValuesFromHistory(const ParticleHistory& history, SpacePartitioning& boundaries);

    void modifyHeightmap(const std::vector<std::vector<ErosionPoint>> &allErosions);
    void modifyVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions);
    void modifyImplicitTerrain(const std::vector<std::vector<ErosionPoint>>& allErosions);
    void modifyLayerBased(const std::vector<std::vector<ErosionPoint>>& allErosions);


    EROSION_APPLIED applyOn;
    TerrainModel* terrain;
    SpacePartitioning* boundariesTree;
    float particleSimulationTime;
    float terrainModifTime;
    Vector3 startingPoint;
    Vector3 originalDirection;
    float randomnessFactor;
    bool fallFromSky;
    float gravity;
    float bouncingCoefficient;
    float bounciness;
    float minSpeed;
    float maxSpeed;
    float maxCapacityFactor;
    float erosion;
    float deposit;
    float matterDensity;
    float materialImpact;
    float airFlowfieldRotation;
    float waterFlowfieldRotation;
    float airForce;
    float waterForce;
    float dt;
    float shearingStressConstantK;
    float shearingRatePower;
    float erosionPowerValue;
    float criticalShearValue;
    std::vector<std::pair<Vector3, Vector3> > posAndDirs;
    FLOWFIELD_TYPE flowType;
    GridV3 waterFlow;
    GridV3 airFlow;
    DENSITY_TYPE densityUsed;
    GridF densityMap;
    float initialCapacity;
    FluidSimType fluidSimType;
    bool wrapPositions;
    bool applyTheErosion;
    int maxCollisions;
    float strength;
    float size;
    int quantity;

    GridF environmentalDensities;
    GridV3 flowfieldValues;

    Vector3 gravityDefault;
    Vector3 terrainSize;
    float particleSize;
    float depositFactor;
    float erosionFactor;

    VoxelGrid* voxelGrid;
    Heightmap* heightmap;
    ImplicitPatch* implicitTerrain;
    LayerBasedGrid* layerBasedGrid;
};



#endif // PARTICLEEROSION_H
