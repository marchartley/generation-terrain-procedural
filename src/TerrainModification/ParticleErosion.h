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

struct ParticleContact {
    Vector3 position;
    Vector3 normal;
    Vector3 velocity;
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
    float capacity = 0.f;
//    float maxCapacity;
    std::shared_ptr<ParticleProperties> properties;

    bool rolling = false;
    bool justStartedRolling = false;


    void update();
    Vector3 predictNextPos(float dt);
    // force stores acceleration; contacts retain the incoming velocity.
    std::vector<ParticleContact> bounce(float dt, SpacePartitioning* boundaries,
                                       float restitution = 1.f, float damping = 1.f, int maxBounces = 8);

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
    // Positive means detached (terrain loses material); negative means deposited.
    ErosionPoint exchange;
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


struct ErosionBalance {
    double initialSediment = 0.;
    double finalSediment = 0.; // remaining load, including particles leaving the domain
    double requestedTerrainChange = 0.;
    double measuredTerrainChange = 0.;
    double unappliedTerrainChange = 0.;
    double sedimentResidual = 0.;
    double proposedTerrainChange = 0.;
    double rejectedDetachment = 0.;
    double rejectedDeposition = 0.;
    // Scalar density sums and occupancy counts are not geometric volume.
    enum Measure { HEIGHT_VOLUME, LAYER_VOLUME, DENSITY_SUM, OCCUPIED_VOXELS, BINARY_EFFECTIVE_VOLUME, UNMEASURED_FIELD };
    Measure measure = UNMEASURED_FIELD;
    bool terrainMeasured = false;
};

class ParticleErosion {
public:
    ParticleErosion();
    ErosionBalance balance;

    std::vector<ParticleHistory> process();
    // Replay proposed contacts in order against native terrain capacity.
    void reconcileTerrainExchanges(std::vector<ParticleHistory>& tracks);
    void initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity);
    ParticleHistory trackParticlePositions(ErosionParticle& particle);
    float computeErosionValue(const ErosionParticle &particle, float vRel);
    float computeDepositionValue(const ErosionParticle& particle);
    Vector3 computeSettlingVelocity(const ErosionParticle& particle) const;

    std::vector<ErosionPoint> computeErosionValuesFromHistory(const ParticleHistory& history, SpacePartitioning& boundaries);

    void modifyHeightmap(const std::vector<std::vector<ErosionPoint>> &allErosions);
    void modifyBinaryVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions);
    void modifyVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions);
    void modifyImplicitTerrain(const std::vector<std::vector<ErosionPoint>>& allErosions);
    void modifyLayerBased(const std::vector<std::vector<ErosionPoint>>& allErosions);


    EROSION_APPLIED applyOn;
    TerrainModel* terrain = nullptr;
    SpacePartitioning* boundariesTree = nullptr;
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

    VoxelGrid* voxelGrid = nullptr;
    Heightmap* heightmap = nullptr;
    ImplicitPatch* implicitTerrain = nullptr;
    LayerBasedGrid* layerBasedGrid = nullptr;

    // Sediment grains and carrier parcels are distinct from the terrain brush.
    float sedimentRadius = 1e-4f; // meters, used only for settling
    float sedimentDensity = 2650.f; // kg/m^3
    float airViscosity = 1.8e-5f; // Pa s
    float waterViscosity = 1e-3f; // Pa s
    float flowResponseTime = 1.f; // seconds; bounded relaxation toward medium velocity
    // Keep sub-hit precision between batches; +/- threshold toggles occupancy.
    GridF binaryHits;
    GridF* binaryHitAccumulator = nullptr;
    float binaryHitScale = 1000.f;
    float binaryHitThreshold = 1000.f;
    std::shared_ptr<EnvironmentalScene> scene;
};



#endif // PARTICLEEROSION_H
