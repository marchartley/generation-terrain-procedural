#include "DataStructure/BVH.h"
#include "Utils/Globals.h"
#include "UnderwaterErosion.h"
#include "TerrainModification/RockErosion.h"
#include "Utils/BSpline.h"
#include "Karst/KarstHole.h"
#include "Utils/Utils.h"
#include "Graph/Matrix3Graph.h"
#include "EnvObject/EnvObject.h"
#include "Utils/Curve1D.h"

UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(VoxelGrid *grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : voxelGrid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

Vector3 getGravityForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, const Vector3& flowfieldValues, Vector3& gravityDirection) {
    // Gravity + Buoyancy
    float gravityForce = gravity * particleMass; // = gravityfieldValues.at(pos + dir);
    float boyancyForce = -environmentDensity * particleVolume; // B = - rho_fluid * V * g
    float gravityCoefficient = gravityForce + boyancyForce;
//    float gravityCoefficient = gravity * particleVolume * (particleDensity - environmentDensity); // Same as F + B
//                std::cout << environmentDensity << " -> " << particleVolume << " * (" << environmentDensity << " - " << matterDensity << ") => " << gravityCoefficient << "\n";
    // float gravityCoefficient = std::max(1.f - (environmentDensity / matterDensity), -1.f); // Keep it between -1 and 1
    Vector3 acceleration = /*flowfieldValues.at(position) +*/ gravityDirection * (gravityForce * gravityCoefficient);
    return acceleration;
}

Vector3 getSettlingForce(float particleDensity, float environmentDensity, float gravity, float capacity, float maxCapacity, Vector3& gravityDirection) {
    // Settling
    float constant = 2.f/9.f;
    float r = 1.f;
    float capacityHinderingExponent = 5.f;
//                acceleration += std::clamp(constant * r * r * (matterDensity - environmentDensity) * gravity * (std::pow(capacity/maxCapacity, capacityHinderingExponent)), -1.f, 1.f);
    Vector3 acceleration = gravityDirection * (constant * r * r * (particleDensity - environmentDensity) * gravity * (1.f - std::pow(capacity/maxCapacity, capacityHinderingExponent)));
    return acceleration;
}
Vector3 getSettlingVelocity(float particleDensity, float particleSize, float environmentDensity, float gravity, float capacity, float maxCapacity, Vector3& gravityDirection) {
    // Settling
    float constant = 2.f/9.f;
    float r = particleSize;
    float capacityHinderingExponent = 5.f;
//                acceleration += std::clamp(constant * r * r * (matterDensity - environmentDensity) * gravity * (std::pow(capacity/maxCapacity, capacityHinderingExponent)), -1.f, 1.f);
    Vector3 acceleration = gravityDirection * (constant * r * r * (particleDensity - environmentDensity) * gravity * (1.f - std::pow(capacity/maxCapacity, capacityHinderingExponent)));
    return acceleration;
}

Vector3 getParticleVelocity() {
    return Vector3();
}

//Vector3 getExternalForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, GridV3& flowfieldValues, Vector3& gravityDirection, float capacity, float maxCapacity) {
Vector3 getExternalForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, const Vector3& flowfieldValues, Vector3& gravityDirection, float capacity, float maxCapacity) {
    Vector3 gravityForce = getGravityForce(position, gravity, particleMass, particleVolume, particleDensity, environmentDensity, flowfieldValues, gravityDirection);
    Vector3 settlingForce = Vector3(); //getSettlingForce(particleDensity, environmentDensity, gravity, capacity, maxCapacity, gravityDirection);
    return gravityForce + settlingForce;
}
//void leapfrogVerlet(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, GridV3& flowfield, bool applyFlow) {
void leapfrogVerlet(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, const Vector3& flowfield, bool applyFlow) {
    // Verlet : x_n+1 = 2 * x_n - x_n-1 + acc * dt^2
    Vector3 oldPos = (position - velocity * dt); // Yeah... it's more an implicit Euler integration...
//    Vector3 usedVelocity = position - oldPos; // (Should be "position - oldPos;")
//    Vector3 v_halfNext = getParticleVelocity(position) + acceleration * dt;
//    position += v_halfNext * dt;
    Vector3 displacement = (position - oldPos) + 0.5f * acceleration * dt * dt;
    if (applyFlow)
        displacement += flowfield * dt;
    //    displacement += flowfield.at(position) * dt;
//    std::cout << displacement << " " << flowfield.at(position) << std::endl;
    Vector3 nextPosition = position + displacement.maxMagnitude(1.f); //(2.f * position - oldPos) + acceleration * dt * dt;
    velocity = nextPosition - position;
    position = nextPosition;
    acceleration *= 0;
}

Vector3 computeParticleForces(ErosionParticle& particle, EnvironmentProperty& env, MaterialProperty& mat)
{
    Vector3 Fg = particle.properties->mass * particle.properties->volume * env.gravity; // Gravity
    Vector3 Fb = -particle.properties->volume * env.density * env.gravity; // Buoyancy
    Vector3 Fd = -6.f * M_PI * env.viscosity * particle.properties->radius * particle.velocity; // Drag

    return Fg + Fb + Fd;
}
void updateParticlePosition(ErosionParticle& particle, EnvironmentProperty& env, MaterialProperty& mat) {
    // Verlet : x_n+1 = 2 * x_n - x_n-1 + acc * dt^2
    float dt = 0.001f;
//    float dt2 = dt * .5f;

    // Euler :
    particle.velocity += (computeParticleForces(particle, env, mat)) * dt;
    particle.pos += particle.velocity;
    particle.dir = particle.velocity.normalized();
    particle.force *= 0.f;
}

std::pair<float, float> computeErosionDeposition(ErosionParticle& particle, EnvironmentProperty& env, MaterialProperty& mat, float erosionFactor, float depositFactor) {
    // Erosion :
    // Constants
    float shearingStressConstantK = 1.f;
    float shearingRatePower = .5f;
    float erosionPowerValue = 1.f;
    float l = particle.properties->radius * 0.002f;

    // Computations
    float relativeVelocity = (1.f - std::abs(particle.dir.dot(mat.normal))); // velocity relative to the surface
    float shearRate = relativeVelocity / l;
    float shear = shearingStressConstantK * std::pow(shearRate, shearingRatePower);
    float amountToErode = erosionFactor * std::pow(std::max(0.f, shear - mat.criticalShearStress), erosionPowerValue) * (1.f - mat.resistance);
    amountToErode = std::min(amountToErode, particle.properties->maxCapacity - particle.capacity);

    // Deposition :
    // Constants
    float hinderedPower = 1.f;
    float densitiesRatio = particle.properties->density / mat.density;

    // Computations
    float hinderedFunction = 1.f - std::pow(particle.capacity / particle.properties->maxCapacity, hinderedPower);
    Vector3 settlingVelocity = (2.f/9.f) * particle.properties->radius * particle.properties->radius * ((particle.properties->density - env.density) / env.viscosity) * env.gravity * hinderedFunction;
    float amountToDeposit = std::min(particle.capacity, depositFactor * densitiesRatio);

    // Here, maybe we could change the settling depending on the normal, but meh
    return {amountToErode, amountToDeposit};
}

void initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity)
{
    particle.pos = position;
    particle.dir = velocity;
    particle.properties = std::make_shared<ParticleProperties>();
    particle.properties->density = density;
    particle.capacity = initialCapacity;
    particle.properties->maxCapacity = maxCapacity;
    particle.force = Vector3();
    particle.properties->volume = 1.f;
    particle.properties->radius = radius;
    particle.properties->volume = .74 * M_PI * particle.properties->radius * particle.properties->radius;
    particle.properties->mass = particle.properties->density * particle.properties->volume;
}

std::tuple<std::vector<BSpline>, int, int, std::vector<std::vector<std::pair<float, Vector3>>>>
UnderwaterErosion::Apply(EROSION_APPLIED applyOn,
                         TerrainModel* terrain,
                         SpacePartitioning& boundariesTree,
                         float &particleSimulationTime, float &terrainModifTime,
                         Vector3 startingPoint,
                         Vector3 originalDirection,
                         float randomnessFactor,
                         bool fallFromSky,
                         float gravity,
                         float bouncingCoefficient,
                         float bounciness,
                         float minSpeed,
                         float maxSpeed,
                         float maxCapacityFactor,
                         float erosion,
                         float deposit,
                         float matterDensity,
                         float materialImpact,
                         float airFlowfieldRotation,
                         float waterFlowfieldRotation,
                         float airForce,
                         float waterForce,
                         float dt,
                         float shearingStressConstantK,
                         float shearingRatePower,
                         float erosionPowerValue,
                         float criticalShearValue,
                         std::vector<std::pair<Vector3, Vector3> > posAndDirs,
//                         std::function<Vector3(Vector3)> flowfieldFunction)
                         FLOWFIELD_TYPE flowType,
                         GridV3 waterFlow,
                         GridV3 airFlow,
                         DENSITY_TYPE densityUsed,
                         GridF densityMap,
                         float initialCapacity,
                         FluidSimType fluidSimType,
                         bool wrapPositions,
                         bool applyTheErosion,
                         int maxCollisions)
{
    /*
    ParticleErosion erosionProcess;

    erosionProcess.applyOn = applyOn;
    erosionProcess.terrain = terrain;
    erosionProcess.boundariesTree = &boundariesTree;
    erosionProcess.particleSimulationTime = particleSimulationTime;
    erosionProcess.terrainModifTime = terrainModifTime;
    erosionProcess.startingPoint = startingPoint;
    erosionProcess.originalDirection = originalDirection;
    erosionProcess.randomnessFactor = randomnessFactor;
    erosionProcess.fallFromSky = fallFromSky;
    erosionProcess.gravity = gravity;
    erosionProcess.bouncingCoefficient = bouncingCoefficient;
    erosionProcess.bounciness = bounciness;
    erosionProcess.minSpeed = minSpeed;
    erosionProcess.maxSpeed = maxSpeed;
    erosionProcess.maxCapacityFactor = maxCapacityFactor;
    erosionProcess.erosion = erosion;
    erosionProcess.deposit = deposit;
    erosionProcess.matterDensity = matterDensity;
    erosionProcess.materialImpact = materialImpact;
    erosionProcess.airFlowfieldRotation = airFlowfieldRotation;
    erosionProcess.waterFlowfieldRotation = waterFlowfieldRotation;
    erosionProcess.airForce = airForce;
    erosionProcess.waterForce = waterForce;
    erosionProcess.dt = dt;
    erosionProcess.shearingStressConstantK = shearingStressConstantK;
    erosionProcess.shearingRatePower = shearingRatePower;
    erosionProcess.erosionPowerValue = erosionPowerValue;
    erosionProcess.criticalShearValue = criticalShearValue;
    erosionProcess.posAndDirs = posAndDirs;
    erosionProcess.flowType = flowType;
    erosionProcess.waterFlow = waterFlow;
    erosionProcess.airFlow = airFlow;
    erosionProcess.densityUsed = densityUsed;
    erosionProcess.densityMap = densityMap;
    erosionProcess.initialCapacity = initialCapacity;
    erosionProcess.fluidSimType = fluidSimType;
    erosionProcess.wrapPositions = wrapPositions;
    erosionProcess.applyTheErosion = applyTheErosion;
    erosionProcess.maxCollisions = maxCollisions;
    erosionProcess.strength = this->maxRockStrength;
    erosionProcess.size = this->maxRockSize;
    erosionProcess.quantity = this->rockAmount;

    erosionProcess.heightmap = heightmap;
    erosionProcess.voxelGrid = voxelGrid;
    erosionProcess.implicitTerrain = implicitTerrain;
    erosionProcess.layerBasedGrid = layerBasedGrid;

    auto allReturns = erosionProcess.process();
    particleSimulationTime = erosionProcess.particleSimulationTime;
    terrainModifTime = erosionProcess.terrainModifTime;

    int quantity = allReturns.size();
    std::vector<BSpline> tunnels(quantity);
    std::vector<std::vector<std::pair<float, Vector3>>> allErosions(quantity);
    int positions = 0;
    int erosions = 0;
    for (size_t i = 0; i < quantity; i++) {
        tunnels[i] = BSpline(allReturns[i].getPositions());
        std::vector<Vector3> erosionPositions = allReturns[i].getErosionPositions();
        std::vector<float> erosionValues = allReturns[i].getErosionValues();
        positions += tunnels[i].points.size();
        erosions += erosionValues.size();
        allErosions[i].resize(erosionValues.size());
        for (size_t j = 0; j < erosionValues.size(); j++) {
            allErosions[i][j] = {erosionValues[j], erosionPositions[j]};
        }
    }
    return {tunnels, positions, erosions, allErosions};
    */


    matterDensity = std::max(matterDensity, 1e-1f);
    VoxelGrid* asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    Heightmap* asHeightmap = dynamic_cast<Heightmap*>(terrain);
    ImplicitPatch* asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    LayerBasedGrid* asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    Vector3 terrainSize = voxelGrid->getDimensions(); //terrain->getDimensions();
    float starting_distance = pow(terrainSize.maxComp()/2.f, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap

    if (densityMap.size() > 0) {
        densityMap.raiseErrorOnBadCoord = false;
        densityMap.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
    }
    // Hydraulic erosion parameters
    float strengthValue = this->maxRockStrength;
    float maxCapacity = 1.f * strengthValue * maxCapacityFactor;
    float erosionFactor = .01f * strengthValue * erosion;
    float depositFactor = .01f * strengthValue * deposit;

    float particleSize = this->maxRockSize;
    float radius = particleSize * .5f;

    GridF modifications((asHeightmap ? Vector3(terrainSize.x, terrainSize.y, 1) : terrainSize));

    Vector3 gravityDefault = Vector3(0, 0, -gravity);
//    auto environmentalDensities = [&](const Vector3& pos) {
//        if (pos.z < 0) return 1000.f;
//        return 1.f;
//    };
//    auto flowfieldValues = [&](const Vector3& pos) {
//        return Vector3();
//    };
//    GridV3 gravityfieldValues = GridV3(terrainSize, Vector3(0, 0, -gravity));
//    gravityfieldValues.raiseErrorOnBadCoord = false;
//    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);

    GridF environmentalDensities = voxelGrid->getEnvironmentalDensities(); // Here also
    GridV3 flowfieldValues = GridV3(terrainSize);
    if (flowType == FLOWFIELD_TYPE::BASIC) {
        Vector3 airDir = Vector3(0, airForce, 0).rotate(0, 0, ((airFlowfieldRotation + random_gen::generate(-45, 45)) / 180) * PI);
        Vector3 waterDir = Vector3(0, waterForce, 0).rotate(0, 0, ((waterFlowfieldRotation + random_gen::generate(-45, 45)) / 180) * PI);
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            Vector3 pos = flowfieldValues.getCoordAsVector3(i);
//            if (pos.x > flowfieldValues.sizeX / 2) continue;
            if (environmentalDensities(pos) < 100) { // In the air
                flowfieldValues(pos) = airDir;
            } else { // In water
                flowfieldValues(pos) = waterDir;
            }
        }
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
        airFlow = airFlow.resize(terrainSize.x, terrainSize.y, 1) * airForce * 10.f;
        waterFlow = waterFlow.resize(terrainSize.x, terrainSize.y, 1) * waterForce * 10.f;
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            if (environmentalDensities.at(flowfieldValues.getCoordAsVector3(i)) < 100) { // In the air
                flowfieldValues.at(i) = airFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            } else { // In water
                flowfieldValues.at(i) = waterFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            }
        }
    } else if (flowType == FLOWFIELD_TYPE::FLUID_SIMULATION) {
        flowfieldValues = GridV3(voxelGrid->getFlowfield(fluidSimType).resize(terrainSize) * airForce).meanSmooth();
        Vector3 aspectRatio = terrainSize.normalized();
        for (auto& v : flowfieldValues)
            v *= aspectRatio;
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS) {
        flowfieldValues = EnvObject::flowfield.resize(terrainSize).meanSmooth() * airForce;
        for (auto& v : flowfieldValues)
            v = v.xy();
    }
    flowfieldValues.raiseErrorOnBadCoord = false;
    flowfieldValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    /*flowfieldValues.iterateParallel([&](const Vector3& pos) {
        flowfieldValues(pos) = Vector3(1.f, std::cos(pos.x * .08f), 0.f).normalized() * 1.f * airForce;
    });*/

    std::vector<BSpline> tunnels(rockAmount);
    std::vector<int> nbPos(rockAmount), nbErosions(rockAmount);
    std::vector<std::vector<std::pair<float, Vector3>>> allErosions(rockAmount);
    std::vector<std::vector<std::pair<float, Vector3>>> allErosionsIfErosionWasOn(rockAmount);

    // Cache computations of the erosion matrix
    RockErosion::createPrecomputedAttackMask(particleSize);
    RockErosion::createPrecomputedAttackMask(particleSize * 2);
    RockErosion::createPrecomputedAttackMask2D(particleSize);
    RockErosion::createPrecomputedAttackMask2D(particleSize * 2);

    std::vector<ErosionParticle> particles(this->rockAmount);

    float totalCollisionTime = 0.f;

    particleSimulationTime = timeIt([&]() {
        #pragma omp parallel for
        for (int i = 0; i < this->rockAmount; i++)
        {
            bool rolling = false;
            std::vector<std::pair<float, Vector3>> erosionValuesAndPositions;

            BSpline tunnel;
            float capacity = maxCapacity * initialCapacity;

            float flowfieldInfluence = 50.0;
            int maxSteps = 500 / dt; // An estimation of how many step we need
            int steps = maxSteps;
            int maxBounces = (maxCollisions < 0 ? 10000 : maxCollisions);
            auto [initialPos, initialDir] = posAndDirs[i];

            ErosionParticle& particle = particles[i];
            initializeParticle(particle, initialPos, initialDir, .5f, matterDensity, capacity, maxCapacity);

            if (terrain->checkIsInGround(particle.pos))
                continue;

    //        ErosionParticle backupParticle;
            bool continueSimulation = true;
            bool hasBeenAtLeastOnceInside = false;
            Vector3 lastCollisionPoint(false);
            int lastBouncingTime = 10000;
            size_t intersectingTriangleIndex = -1;
            std::vector<Vector3> lastCollisions;
            std::tuple<int, int> firstLastCollisionIndices = {-1, -1};
            while (continueSimulation) {
                Vector3 nextPos = particle.pos + (particle.dir * dt);
                float environmentDensity = environmentalDensities.at(particle.pos);
                float particleVolume = .1f;
                float gravityCoefficient = particleVolume * gravity * (matterDensity - environmentDensity) * 0.01f * (matterDensity < 0 ? 0.f : 1.f);
    //            float gravityCoefficient = 0.f;
                Vector3 justFlow = flowfieldValues(particle.pos); // (flowfieldValues.at(pos) + flowfieldValues.at(nextPos)) * .5f;
                Vector3 justGravity = gravityDefault * gravityCoefficient;
    //            Vector3 justGravity = (gravityfieldValues.at(particle.pos) * gravityCoefficient); // (gravityfieldValues.at(nextPos) * gravityCoefficient); // .maxMagnitude(2.f);
    //            std::cout << "Z=" << particle.pos.z << " -> flow " << justFlow << " -- gravity: " << justGravity << "(" << gravityDefault << " * " << gravityCoefficient << ") density = " << environmentDensity << " (" << environmentalDensities << ")" << std::endl;
                Vector3 flowfield = justFlow + justGravity;

                    if (matterDensity < 0)
                        particle.dir = (particle.dir + flowfield * flowfieldInfluence) * .5f;
                    //            if (!rolling)
                    else
                        particle.dir += flowfield * flowfieldInfluence * dt;

                particle.dir.maxMagnitude(maxSpeed);
                if (Vector3::isInBox(nextPos.xy(), -terrainSize.xy(), terrainSize.xy()))
                    hasBeenAtLeastOnceInside = true;
                auto collisionPointAndTriangleIndex = boundariesTree.getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, (rolling ? intersectingTriangleIndex : -1));
                Vector3 collisionPoint;
                std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
                if (rolling || collisionPoint.isValid()) {
    //                steps = std::min(steps, 500);
                    if (std::get<0>(firstLastCollisionIndices) < 0) std::get<0>(firstLastCollisionIndices) = steps;
                    std::get<1>(firstLastCollisionIndices) = steps;

                    maxBounces --;

                    if (steps < 0 || maxBounces <= 0)
                        continueSimulation = false;

                    bool justStartedRolling = false;

                    Vector3 normal;
                    float vRel;

                    if (!rolling && lastBouncingTime < 5 && collisionPoint.isValid()) {
                        //steps = -1000;
                        rolling = true;
                        justStartedRolling = true;
    //                    particle.pos.z += .1f;
                    } else if (rolling && collisionPoint.isValid()) {
                        rolling = false;
                        lastCollisions.push_back(collisionPoint);
                        if (lastCollisions.size() > 5) {
                            lastCollisions.erase(lastCollisions.begin());
                            Vector3 mean;
                            for (const auto& p : lastCollisions)
                                mean += p;
                            mean /= float(lastCollisions.size());
                            float error = 0.f;
                            for (const auto& p : lastCollisions)
                                error += (p - mean).norm2();
    //                        if (error / float(lastCollisions.size()) < 1e-0)
    //                            continueSimulation = false;
                        }

                    }

                    if (rolling && !justStartedRolling) {
                        vRel = particle.dir.norm();
                    } else {
                        normal = boundariesTree.triangles[intersectingTriangleIndex].normal;
                        float dotNormal = std::abs(particle.dir.normalized().dot(normal));
                        vRel = particle.dir.norm() * (1.f - dotNormal); // velocity relative to the surface

                        if (!continueSimulation && normal.z < 0)
                            continueSimulation = true;
                    }

                    if (justStartedRolling) {
                        particle.dir = (particle.dir - (normal * particle.dir.dot(normal))).setMag(particle.dir.norm());
                        collisionPointAndTriangleIndex = boundariesTree.getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, intersectingTriangleIndex);
                        std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
                        if (collisionPoint.isValid())
                            normal = boundariesTree.triangles[intersectingTriangleIndex].normal;
                    }

                    lastBouncingTime = 0;

                    // Using SPH formula
                    // Erosion part
                    float l = particleSize * 0.001f;
                    float theta = vRel / l;
                    float shear = shearingStressConstantK * std::pow(theta, shearingRatePower);
                    float amountToErode = erosionFactor * std::pow(std::max(0.f, shear), erosionPowerValue); // std::pow(std::max(0.f, shear - criticalShearStress), erosionPowerValue)  * (1.f - resistanceValue * materialImpact);
    //                if (rolling)
    //                    amountToErode = 0.f;
                    if (maxCollisions != -1) {
                        amountToErode *= 1000000.f;
                    }
    //                std::cout << theta << " " << shear << " " << amountToErode << " " << particle.maxCapacity << " " << particle.capacity << std::endl;
                    amountToErode = std::min(amountToErode, particle.properties->maxCapacity - particle.capacity);
                    if (nextPos.isValid() && nextPos.z < 3.f) {
//                        std::cout << "Ignore\n";
                        amountToErode = 0;
                    }


    //                    // Deposition part
                    float densitiesRatio = 1.f; //1000.f / matterDensity; // The ground should have a density, I'll assume the value 1000 for now
    //                float u = (2.f / 9.f) * particle.radius * particle.radius * (matterDensity - environmentDensity) * gravity * (1.f - (particle.capacity / particle.maxCapacity));
                    float amountToDeposit = std::min(particle.capacity, depositFactor * (densitiesRatio * capacity * (particle.capacity / particle.properties->maxCapacity)));

                    float amountToErodeIfErosionWasOn = std::pow(std::max(0.f, shear), erosionPowerValue);
                    float amountToDepositIfErosionWasOn = 0.1f;

                    if (amountToErode - amountToDeposit != 0) {
                        particle.capacity += (amountToErode - amountToDeposit);
                        if (nextPos.z >= 0.f) {
                            erosionValuesAndPositions.push_back({amountToErode - amountToDeposit, nextPos});
                        }
                    }
                    if (collisionPoint.isValid()) {
                        totalCollisionTime += timeIt([&]() {
                            // Continue the rock tracing
                            int maxCollisions = 2;
                            do {
                                maxCollisions --;
                                if (maxCollisions <= 0 || !continueSimulation) {
                                    continueSimulation = false;
                                    break;
                                }
                                Vector3 bounce = particle.dir.reflexion(normal);
                                bounce -= normal * bounce.dot(normal) * (1.f - bounciness);
                                particle.dir = bounce * bouncingCoefficient;
    //                            particle.dir.minMagnitude(.5f);

                                particle.pos = collisionPoint;
                                lastCollisionPoint = collisionPoint;
                                std::tie(collisionPoint, intersectingTriangleIndex) = boundariesTree.getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, intersectingTriangleIndex);
                            } while (collisionPoint.isValid());
                        });
                    } else if (rolling){
                        particle.dir *= bouncingCoefficient;
    //                    particle.dir = (particle.dir - (normal * particle.dir.dot(normal))).setMag(particle.dir.norm());
                    }
                } else {

                }
                lastBouncingTime ++;
                steps --;

                tunnel.points.push_back(particle.pos);

                particle.pos = particle.pos + (particle.dir * dt);
                if (wrapPositions)
                    particle.pos = Vector3::wrap(particle.pos, Vector3(0, 0, -100), Vector3(terrain->getSizeX(), terrain->getSizeY(), 1000));
                particle.dir *= 0.99f;
                if (steps < 0 || (hasBeenAtLeastOnceInside && nextPos.z < -20) || particle.pos.z < -20 || particle.dir.norm2() < 1e-4 || !continueSimulation || maxBounces <= 0) {
                    if (depositFactor > 0.f && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                        while (!terrain->checkIsInGround(particle.pos) && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                            particle.pos.z -= .5f;
                        }
                        Vector3 depositPosition = particle.pos;
                        erosionValuesAndPositions.push_back({-particle.capacity, depositPosition});
                    }
                    continueSimulation = false;
                }
    //            if (steps < 0)
    //                continueSimulation = false;
            }
            tunnels[i] = tunnel;
            nbPos[i] = tunnel.points.size();
            nbErosions[i] = erosionValuesAndPositions.size();

            allErosions[i] = erosionValuesAndPositions;
        }
    });

    int positions = 0, erosions = 0;
    for (int i = 0; i < rockAmount; i++) {
        positions += nbPos[i];
        erosions += nbErosions[i];
    }

    if (!applyTheErosion || strengthValue == 0) {
        return {tunnels, positions, erosions, allErosions};
    }

    terrainModifTime = timeIt([&]() {
        std::vector<ImplicitNaryOperator*> allNary(allErosions.size());
        if (asImplicit) {
            for (size_t i = 0; i < allNary.size(); i++) {
                allNary[i] = new ImplicitNaryOperator;
            }
        }

    //    float summary = 0.f;
        if (asLayers) {
            for (size_t i = 0; i < allErosions.size(); i++) {
                const auto& erosionValuesAndPositions = allErosions[i];
                Vector3 previousPos(false);
                bool isTheOnlyDepositRock = true;
                for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
                    bool lastRock = iRock == erosionValuesAndPositions.size() - 1;
                    auto [val, pos] = erosionValuesAndPositions[iRock];
                    if (std::abs(val) < 1e-4) continue;

                    if (lastRock)
                        continue;
                    previousPos = pos;
                    val *= 1.f;
                    float max = std::abs(val) * std::sqrt(particleSize*particleSize)/particleSize;
                    float radius = max;
                    float radiusXY = max;
                    for (int x = -radiusXY; x < radiusXY; x++) {
                        for (int y = -radiusXY; y < radiusXY; y++) {
                            if (x*x + y*y < radiusXY * radiusXY) {
                                float halfHeight = .5f * std::sqrt(radius*radius - (x*x + y*y))/radius;
                                float startZ = std::max(0.f, pos.z - halfHeight);
                                float endZ = pos.z + halfHeight;
                                asLayers->transformLayer(pos.x + x, pos.y + y, startZ, endZ, (val > 0 ? TerrainTypes::AIR : TerrainTypes::SAND));
                            }
                        }
                    }
                }
            }
        } else {
            bool parallel = false;
            if (parallel) {
                std::vector<GridF> submodifications(rockAmount, GridF(modifications.getDimensions()));
                #pragma omp parallel for
                for (size_t i = 0; i < allErosions.size(); i++) {
                    const auto& erosionValuesAndPositions = allErosions[i];
                    Vector3 previousPos(false);
                    for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
                        auto [val, pos] = erosionValuesAndPositions[iRock];
                        int size = particleSize;
                        val *= 100.f;
                        if (iRock == erosionValuesAndPositions.size() - 1) {
                        }

                        if (asVoxels) {
                            RockErosion(size, val).computeErosionMatrix(submodifications[i], pos - Vector3(.5f, .5f, .5f));
                        } else if (asHeightmap) {
                            RockErosion(size, val).computeErosionMatrix2D(submodifications[i], pos);
                        } else if (asImplicit) {
                            if (iRock == erosionValuesAndPositions.size() - 1) {
                                val *= .1f;
                            }
        //                    if (val <= 0.f) continue; // Fits the paper Paris et al. Terrain Amplification with Implicit 3D Features (2019)
                            float dimensions = size;
                            if (dimensions < 1.f) continue;
                            if (abs(val) < .5f) continue;
                            previousPos = pos;
                            auto sphere = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(dimensions, dimensions, dimensions), std::abs(val)*.1f));
                            sphere->material = (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT);
                            sphere->dimensions = Vector3(dimensions, dimensions, dimensions);
                            sphere->supportDimensions = sphere->dimensions;
                            sphere->position = pos - sphere->dimensions * .5f;
                            allNary[i]->composables.push_back(sphere);
                        }
                    }
                }
                for (const auto& sub : submodifications)
                    modifications += sub;
            } else { // Not parallel
                int countRocks = 0;
                GridF rocks(terrainSize);
                for (size_t i = 0; i < allErosions.size(); i++) {
                    const auto& erosionValuesAndPositions = allErosions[i];
                    Vector3 previousPos(false);
                    for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
                        auto [val, pos] = erosionValuesAndPositions[iRock];
//                        int size = particleSize;
                        val *= 100.f;
                        if (iRock == erosionValuesAndPositions.size() - 1) {
                        }

                        if (asVoxels) {
                            RockErosion(particleSize, val).computeErosionMatrix(modifications, pos - Vector3(.5f, .5f, .5f));
                        } else if (asHeightmap) {
                            RockErosion(particleSize, val).computeErosionMatrix2D(modifications, pos);
                        } else if (asImplicit) {
        //                    if (val <= 0.f) continue; // Fits the paper Paris et al. Terrain Amplification with Implicit 3D Features (2019)
                            if (particleSize < 1.f) continue;
//                            if (abs(val) < .5f) continue;
//                            if (previousPos.isValid() && (pos - previousPos).norm2() < 2*2) continue;
                            countRocks++;
                            rocks(pos) += val;
                            previousPos = pos;/*
                            auto sphere = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(particleSize, particleSize, particleSize), std::abs(val)*.01f));
                            sphere->material = (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT);
                            sphere->dimensions = Vector3(particleSize, particleSize, particleSize);
                            sphere->supportDimensions = sphere->dimensions;
                            sphere->position = pos - sphere->dimensions * .5f;
                            allNary[i]->composables.push_back(sphere);*/
                        }
                    }
                }
                std::cout << "Only " << countRocks << " counted. Min = " << rocks.min() << " Max = " << rocks.max() << std::endl;
                rocks.iterate([&](const Vector3& pos) {
                    float val = rocks(pos);
                    if (std::abs(val) < 10) return;

                    auto sphere = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(particleSize, particleSize, particleSize), std::abs(val)*.01f));
                    sphere->material = (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT);
                    sphere->dimensions = Vector3(particleSize, particleSize, particleSize);
                    sphere->supportDimensions = sphere->dimensions;
                    sphere->position = pos - sphere->dimensions * .5f;
                    allNary[0]->composables.push_back(sphere);
                });

            }
        }


        if (densityMap.size() > 0) {
            for (size_t i = 0; i < modifications.size(); i++) {
                modifications[i] = (modifications[i] > 0 ? modifications[i] : modifications[i] * ((1.f - materialImpact) + (1.f - densityMap[i]) * materialImpact));
            }
        }

    //    std::cout << "Total erosion : " << modifications.sum() << std::endl;
    //    std::cout << "Collisions cost " << showTime(totalCollisionTime) << std::endl;
    //    std::cout << "Other cost " << showTime(totalOtherTime) << std::endl;
        if (asVoxels) {
            for (int x = 0; x < modifications.sizeX; x++)
                for (int y = 0; y < modifications.sizeY; y++)
                    for (int z = 0; z < 2; z++)
                        modifications.at(x, y, z) = 0;
            asVoxels->applyModification(modifications * 0.5f);
            asVoxels->limitVoxelValues(2.f, false);
            asVoxels->saveState();
        } else if (asHeightmap) {
            asHeightmap->heights += modifications * 0.5f;
            asHeightmap->heights.iterateParallel([&](size_t i) {
                asHeightmap->heights[i] = std::max(asHeightmap->heights[i], 0.f);
            });
        } else if (asImplicit) {
            ImplicitNaryOperator* totalErosion = new ImplicitNaryOperator;
            for (auto& op : allNary)
                if (!op->composables.empty())
                    totalErosion->composables.push_back(op);
            if (!totalErosion->composables.empty())
            {
                dynamic_cast<ImplicitNaryOperator*>(implicitTerrain)->addChild(totalErosion);
                implicitTerrain->_cached = false;
            }
        } else if (asLayers) {
            // ...
        }
    });
    return {tunnels, positions, erosions, allErosions};

}

std::tuple<GridF, GridF, GridF> UnderwaterErosion::flatteningErodedTerrain(const GridF &initialTerrain, const GridF &currentTerrain)
{
    Vector3 finalDimensions(initialTerrain.sizeX, initialTerrain.sizeY, 1);
    GridF flatMap(finalDimensions);
    GridF erosionMap(finalDimensions);
    GridF depositionMap(finalDimensions);

    GridF flatInitial(finalDimensions);
    GridF flatCurrent(finalDimensions);

    initialTerrain.iterateParallel([&](const Vector3& pos) {
        if (initialTerrain(pos) > 0) flatInitial(pos.xy()) += 1.f;
        if (currentTerrain(pos) > 0) flatCurrent(pos.xy()) += 1.f;
    });

    flatMap.iterateParallel([&](size_t i) {
        depositionMap[i] = std::max(flatCurrent[i] - flatInitial[i], 0.f);
        erosionMap[i] = std::max(flatInitial[i] - flatCurrent[i], 0.f);
        flatMap[i] = flatCurrent[i] - depositionMap[i];
    });

    return {flatMap, erosionMap, depositionMap};
}











std::vector<Vector3> UnderwaterErosion::CreateTunnel(int numberPoints, bool addingMatter, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    BSpline curve = BSpline(numberPoints); // Random curve
    for (Vector3& coord : curve)
        coord = ((coord + Vector3(1.0, 1.0, 1.0)) / 2.0) * voxelGrid->getDimensions(); //Vector3(grid->sizeX, grid->sizeY, grid->sizeZ);
    return CreateTunnel(curve, addingMatter, true, applyChanges, startingShape, endingShape);
}
std::vector<Vector3> UnderwaterErosion::CreateTunnel(BSpline path, bool addingMatter, bool usingSpheres, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    GridF erosionMatrix(voxelGrid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    bool modificationDoesSomething = true;
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 2.0),
                                                     Vector3(1.0, 1.0)
                                                 }));

    std::vector<Vector3> coords;
    if (usingSpheres) {
        float nb_points_on_path = path.length() / (this->maxRockSize/5.f);
        RockErosion rock(this->maxRockSize, this->maxRockStrength);
        for (const auto& pos : path.getPath(nb_points_on_path)) {
            erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos);
        }
        if (erosionMatrix.abs().max() < 1.e-6) {
            modificationDoesSomething = false;
        } else {
            erosionMatrix = erosionMatrix.abs();
            erosionMatrix.toDistanceMap();
            erosionMatrix.normalize();
            for (float& m : erosionMatrix) {
                m = interpolation::linear(m, 0.f, 1.0) * this->maxRockStrength * (addingMatter ? 1.f : -1.f);
        //        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
            }
        }
    } else {
        KarstHole hole(path, this->maxRockSize, this->maxRockSize, startingShape, endingShape);
        GridF holeMatrix;
        Vector3 anchor;
        std::vector<std::vector<Vector3>> triangles = hole.generateMesh();
        std::tie(holeMatrix, anchor) = hole.generateMask(triangles);
        if (holeMatrix.abs().max() == 0) {
            modificationDoesSomething = false;
        } else {
            holeMatrix = holeMatrix.abs().toDistanceMap();
            holeMatrix.normalize();
            for (float& m : holeMatrix) {
//                m = (m > 0 ? 1.0 : 0.0) * -this->maxRockStrength;
                m = interpolation::linear(m, 0.f, 1.0) * -this->maxRockStrength * (addingMatter ? -1.f : 1.f);
        //        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
            }
            RockErosion rock;
            erosionMatrix = rock.computeErosionMatrix(erosionMatrix, holeMatrix, path.getPoint(0), addingMatter, anchor, true);

            for (const auto& triangle : triangles) {
                coords.push_back(triangle[0]);
                coords.push_back(triangle[1]);
                coords.push_back(triangle[1]);
                coords.push_back(triangle[2]);
                coords.push_back(triangle[2]);
                coords.push_back(triangle[0]);
            }
        }
    }
    if (modificationDoesSomething) {
        voxelGrid->applyModification(erosionMatrix);
    }
//    if (applyChanges)
//        grid->remeshAll();
    return coords;
}

std::vector<std::vector<Vector3> > UnderwaterErosion::CreateMultipleTunnels(std::vector<BSpline> paths, bool addingMatter, bool usingSpheres, bool applyChanges,
                                                                            KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    GridF erosionMatrix(voxelGrid->getDimensions()); // this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 0.5),
                                                     Vector3(1.0, 1.0)
                                                 }));

//    float resolution = 0.01;
    std::vector<std::vector<Vector3>> allCoords;
    for (BSpline& path : paths) {
        bool modificationDoesSomething = true;
        if (usingSpheres) {
            float nb_points_on_path = path.length() / (float)(this->maxRockSize / 3.f);
            std::vector<Vector3> coords;
            RockErosion rock(this->maxRockSize, this->maxRockStrength);
            for (const auto& pos : path.getPath(nb_points_on_path)) {
                erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos);
            }
        } else {
            KarstHole hole(path, this->maxRockSize, this->maxRockSize, startingShape, endingShape);
            GridF holeMatrix;
            Vector3 anchor;
            std::vector<std::vector<Vector3>> triangles = hole.generateMesh();
            std::tie(holeMatrix, anchor) = hole.generateMask(triangles);
            if (holeMatrix.abs().max() == 0) {
                modificationDoesSomething = false;
            } else {
                holeMatrix = holeMatrix.abs().toDistanceMap();
                holeMatrix.normalize();
                for (float& m : holeMatrix) {
                    m = interpolation::linear(m, 0.f, 1.0) * this->maxRockStrength * (addingMatter ? 1.f : -1.f);
            //        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
                }
//                holeMatrix *= -this->maxRockStrength;
                RockErosion rock;
                erosionMatrix = rock.computeErosionMatrix(erosionMatrix, holeMatrix, path.getPoint(0), addingMatter, anchor);

                std::vector<Vector3> coords;
                for (const auto& triangle : triangles) {
                    coords.push_back(triangle[0]);
                    coords.push_back(triangle[1]);
                    coords.push_back(triangle[1]);
                    coords.push_back(triangle[2]);
                    coords.push_back(triangle[2]);
                    coords.push_back(triangle[0]);
                }
                allCoords.push_back(coords);
            }
        }
    }
    voxelGrid->applyModification(erosionMatrix * 20.f);
//    if (applyChanges)
//        grid->remeshAll();
    return allCoords;
}

std::vector<Vector3> UnderwaterErosion::CreateCrack(const Vector3& start, const Vector3& end, bool applyChanges)
{
    float rx = 3.f, ry = 3.f, rz = 3.f;
    Vector3 ratio(rx, ry, rz);
    GridI resizedMap = this->voxelGrid->getVoxelValues().resize(voxelGrid->getDimensions() / ratio).binarize();
//    GridI resizedMap = this->grid->getVoxelValues().resize(this->grid->sizeX / rx, this->grid->sizeY / ry, this->grid->sizeZ / rz).binarize();
    Matrix3Graph graph = Matrix3Graph(resizedMap).computeSurface().randomizeEdges(.5f);
    Vector3 clampedStart = start / ratio;
    clampedStart.x = std::clamp(clampedStart.x, 0.f, resizedMap.sizeX - 1.f);
    clampedStart.y = std::clamp(clampedStart.y, 0.f, resizedMap.sizeY - 1.f);
    clampedStart.z = std::clamp(clampedStart.z, 0.f, resizedMap.sizeZ - 1.f);
    Vector3 clampedEnd = end / ratio;
    clampedEnd.x = std::clamp(clampedEnd.x, 0.f, resizedMap.sizeX - 1.f);
    clampedEnd.y = std::clamp(clampedEnd.y, 0.f, resizedMap.sizeY - 1.f);
    clampedEnd.z = std::clamp(clampedEnd.z, 0.f, resizedMap.sizeZ - 1.f);
    std::vector<Vector3> path = graph.shortestPath(clampedStart, clampedEnd);

    ratio = ((start / clampedStart.rounded()) + (end / clampedEnd.rounded())) / 2.f; // To be continued...
    for (Vector3& point : path)
        point *= ratio;
    std::vector<Vector3> meshPath;
    for (size_t i = 0; i < path.size() - 1; i++) {
        meshPath.push_back(path[i]);
        meshPath.push_back(path[i + 1]);
    }
    return this->CreateTunnel(BSpline(path).simplifyByRamerDouglasPeucker(0.5f), false, false, applyChanges, KarstHolePredefinedShapes::CRACK, KarstHolePredefinedShapes::CRACK);
    //    return meshPath;
}

std::vector<Vector3> UnderwaterErosion::CreateTunnel(KarstHole &tunnel, bool addingMatter, bool applyChanges)
{
    GridF erosionMatrix(voxelGrid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    bool modificationDoesSomething = true;
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 2.0),
                                                     Vector3(1.0, 1.0)
                                                 }));

    std::vector<Vector3> coords;
    GridF holeMatrix;
    Vector3 anchor;
    std::vector<std::vector<Vector3>> triangles = tunnel.generateMesh();
    std::tie(holeMatrix, anchor) = tunnel.generateMask(triangles);
    if (holeMatrix.abs().max() == 0) {
        modificationDoesSomething = false;
    } else {
        holeMatrix = holeMatrix.abs().toDistanceMap();
        holeMatrix.normalize();
        for (float& m : holeMatrix) {
//                m = (m > 0 ? 1.0 : 0.0) * -this->maxRockStrength;
            m = interpolation::linear(m, 0.f, 1.0) * -this->maxRockStrength * (addingMatter ? -1.f : 1.f);
    //        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
        }
        holeMatrix = holeMatrix.meanSmooth();
        RockErosion rock;
        erosionMatrix = rock.computeErosionMatrix(erosionMatrix, holeMatrix, tunnel.path.getPoint(0), addingMatter, anchor, true);

        for (const auto& triangle : triangles) {
            coords.push_back(triangle[0]);
            coords.push_back(triangle[1]);
            coords.push_back(triangle[1]);
            coords.push_back(triangle[2]);
            coords.push_back(triangle[2]);
            coords.push_back(triangle[0]);
        }
    }
    if (modificationDoesSomething) {
        voxelGrid->applyModification(erosionMatrix);
    }
//    if (applyChanges)
//        grid->remeshAll();
    return coords;
}

void UnderwaterErosion::ParisSeaErosion()
{
    Curve1D resistanceCurve = Curve1D(BSpline({Vector3(0, 0, 0), Vector3(1, 0, 0)}).getPath(20));
    for (auto& p : resistanceCurve)
        p.y = random_gen::generate(0.f, 1.f);

    float waterLevel = voxelGrid->properties->waterLevel;
    float previousWater = clamp(waterLevel - 0.1f, 0.f, 1.f);
    float nextWater = clamp(waterLevel + 0.1f, 0.f, 1.f);
    Curve1D seaErosionCurve = Curve1D({Vector3(0, 0, 0), Vector3(previousWater, 0, 0), Vector3(waterLevel, 1, 0), Vector3(nextWater, 0, 0), Vector3(1, 0, 0)});
//    std::cout << "Water level: " << waterLevel << std::endl;
//    std::cout << seaErosionCurve.display1DPlot(200, 20) << std::endl;
//    std::cout << resistanceCurve.display1DPlot(200, 20) << std::endl;

    GridF surface = voxelGrid->getVoxelValues().binarize().isosurface();
    std::vector<Vector3> possiblePositions;
    possiblePositions.reserve(surface.sum());
    surface.iterate([&](Vector3 pos) {
        if (surface(pos) > 0.5) possiblePositions.push_back(pos);
    });
    std::vector<std::pair<Vector3, float>> erosionPositionsAndEnergy;
    size_t nbSamples = 400;

    float beta = 10.f;
    float mini = std::numeric_limits<float>::max(), maxi = std::numeric_limits<float>::lowest();
    float maxZ;
    for (size_t i = 0; i < possiblePositions.size(); i++) {
        auto& p = possiblePositions[i];
        float z = p.z / voxelGrid->getSizeZ();
        float resistance = resistanceCurve.get(z);
        float erosion = seaErosionCurve.get(z);
        float energy = erosion * (1.f - resistance) * (1.f - beta) + erosion * beta;


        mini = std::min(mini, energy);
        if (maxi < energy) maxZ = z;
        maxi = std::max(maxi, energy);
        if (energy > beta * .5f) {
            erosionPositionsAndEnergy.push_back({p, energy});
        }
    }
//    std::cout << "\nMin: " << mini << " , max: " << maxi << " at z = " << maxZ << " , " << erosionPositionsAndEnergy.size() << " voxels affected out of " << possiblePositions.size() << " voxels on surface." << std::endl;

    std::shuffle(erosionPositionsAndEnergy.begin(), erosionPositionsAndEnergy.end(), random_gen::random_generator);
    erosionPositionsAndEnergy.resize(std::min(erosionPositionsAndEnergy.size(), nbSamples));

    std::cout << "Erosion time: " << showTime(timeIt([&]() {
        GridF modifs(voxelGrid->getDimensions());
        int radius = 3;
        for (auto& [p, energy] : erosionPositionsAndEnergy) {
            RockErosion erosion(radius * 2 + 1, energy);
            erosion.computeErosionMatrix(modifs, p);
        }
        voxelGrid->applyModification(modifs);
    })) << std::endl;
}

void UnderwaterErosion::ParisInvasionPercolation()
{

}

