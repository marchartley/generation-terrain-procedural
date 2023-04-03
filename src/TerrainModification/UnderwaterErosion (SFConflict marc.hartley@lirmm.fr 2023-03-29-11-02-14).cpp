#include "Utils/Globals.h"
#include "UnderwaterErosion.h"
#include "TerrainModification/RockErosion.h"
#include "Utils/BSpline.h"
#include "Karst/KarstHole.h"
#include "Utils/Utils.h"
#include "Graph/Matrix3Graph.h"


UnderwaterErosion::UnderwaterErosion()
{

}
UnderwaterErosion::UnderwaterErosion(VoxelGrid *grid, int maxRockSize, float maxRockStrength, int rockAmount)
    : voxelGrid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
{

}

void retroChangeFlowfield(std::vector<Vector3>& coords, std::vector<Vector3>& dirs, std::shared_ptr<VoxelGrid> grid)
{
    if (coords.size() == 0)
        return;
    Vector3& impactZone = coords[coords.size() - 1];
    for (size_t i = 0; i < coords.size(); i++)
    {
        Vector3& coord = coords[i];
        Vector3& dir = dirs[i];
        float alpha_effect = 0.5 * ((coord - impactZone).norm2() < 20.0 ? -2.0 : 1.0); // Inverse and double if it's a choc (last coord)
        grid->affectFlowfieldAround(coord, dir * alpha_effect, 3);
//        grid->affectFlowfieldAround(coord, alpha_effect, 10);
    }
}

//Vector3 getGravityForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, Matrix3<Vector3>& flowfieldValues, Vector3& gravityDirection) {
Vector3 getGravityForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, Vector3 flowfieldValues, Vector3& gravityDirection) {
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

//Vector3 getExternalForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, Matrix3<Vector3>& flowfieldValues, Vector3& gravityDirection, float capacity, float maxCapacity) {
Vector3 getExternalForce(Vector3& position, float gravity, float particleMass, float particleVolume, float particleDensity, float environmentDensity, Vector3 flowfieldValues, Vector3& gravityDirection, float capacity, float maxCapacity) {
    Vector3 gravityForce = getGravityForce(position, gravity, particleMass, particleVolume, particleDensity, environmentDensity, flowfieldValues, gravityDirection);
    Vector3 settlingForce = Vector3(); //getSettlingForce(particleDensity, environmentDensity, gravity, capacity, maxCapacity, gravityDirection);
    return gravityForce + settlingForce;
}
//void leapfrogVerlet(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, Matrix3<Vector3>& flowfield, bool applyFlow) {
void leapfrogVerlet(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, Vector3 flowfield, bool applyFlow) {
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

//void updateParticlePosition(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, TerrainModel* terrain, Matrix3<Vector3>& flowfield, bool applyFlow) {
void updateParticlePosition(Vector3& position, Vector3& velocity, Vector3& acceleration, float dt, TerrainModel* terrain, Vector3 flowfield, bool applyFlow) {
    // Verlet : x_n+1 = 2 * x_n - x_n-1 + acc * dt^2
    Vector3 prevVel = velocity;
    leapfrogVerlet(position, velocity, acceleration, dt, flowfield, applyFlow);
    // Euler : x_n+1 = x_n + vel * dt
//    velocity += (acceleration * dt).maxMagnitude(1.f);
//    position += (velocity * dt).maxMagnitude(1.f);
//    acceleration *= 0;
//    Vector3 previous = position;
//    while (terrain->checkIsInGround(position)) {
//        if (prevVel.norm2() == 0) {
//            break;
//        }
//        position -= prevVel * 0.01f; // Move particle until it's not in the terrain anymore
//    }
//    std::string s = previous.toString() + " -> " + position.toString();
//    std::cout << s << std::endl;
//    position = position;
}

void setInitialPositionAndVelocityAndAccelerationOfParticle(Vector3& position, Vector3& velocity, Vector3& acceleration, const Vector3& startingPoint, float starting_distance, const Vector3& originalDirection, float randomnessFactor, bool fallFromSky, const Vector3 terrainDimensions)
{
    if (!startingPoint.isValid()) {
        position = Vector3(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(0.0, 1.0));
        position.normalize();
        position *= starting_distance;
    } else {
        position = startingPoint;
    }

//        velocity = Vector3::random();
    if (originalDirection.isValid())
        acceleration = (originalDirection + Vector3::random(randomnessFactor)).normalize() * 1000;

    if (fallFromSky) {
        position = Vector3(random_gen::generate(terrainDimensions.x), random_gen::generate(terrainDimensions.y), terrainDimensions.z + 10);
        acceleration = Vector3(0, 0, -1.f);
        acceleration += Vector3::random(randomnessFactor);
    }
    velocity = Vector3();

//    position.x = 30.f;
//    position.y = 30.f;
}

std::tuple<std::vector<BSpline>, int, int>
UnderwaterErosion::Apply(EROSION_APPLIED applyOn, float &particleSimulationTime, float &terrainModifTime, Vector3 startingPoint, Vector3 originalDirection, float randomnessFactor, bool fallFromSky,
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
                         Matrix3<Vector3> waterFlow,
                         Matrix3<Vector3> airFlow,
                         DENSITY_TYPE densityUsed,
                         Matrix3<float> densityMap,
                         float initialCapacity)
{
    // ApplyOn : 0 = density-voxels, 1 = heightmap, 2 = implicit, 3 = layers, 4 = binary-voxels
    if (applyOn == EROSION_APPLIED::DENSITY_VOXELS) {
        return this->ApplyOnAnyTerrain(voxelGrid, particleSimulationTime, terrainModifTime, startingPoint, originalDirection, randomnessFactor, fallFromSky, gravity, bouncingCoefficient, bounciness, minSpeed, maxSpeed, maxCapacityFactor, erosion, deposit, matterDensity, materialImpact, airFlowfieldRotation, waterFlowfieldRotation, airForce, waterForce, dt, shearingStressConstantK, shearingRatePower, erosionPowerValue, criticalShearValue, posAndDirs, flowType, waterFlow, airFlow, densityUsed, densityMap, initialCapacity);
    } else if (applyOn == EROSION_APPLIED::HEIGHTMAP) {
        return this->ApplyOnAnyTerrain(heightmap, particleSimulationTime, terrainModifTime, startingPoint, originalDirection, randomnessFactor, fallFromSky, gravity, bouncingCoefficient, bounciness, minSpeed, maxSpeed, maxCapacityFactor, erosion, deposit, matterDensity, materialImpact, airFlowfieldRotation, waterFlowfieldRotation, airForce, waterForce, dt, shearingStressConstantK, shearingRatePower, erosionPowerValue, criticalShearValue, posAndDirs, flowType, waterFlow, airFlow, densityUsed, densityMap, initialCapacity);
    } else if (applyOn == EROSION_APPLIED::IMPLICIT_TERRAIN) {
        return this->ApplyOnAnyTerrain(implicitTerrain, particleSimulationTime, terrainModifTime, startingPoint, originalDirection, randomnessFactor, fallFromSky, gravity, bouncingCoefficient, bounciness, minSpeed, maxSpeed, maxCapacityFactor, erosion, deposit, matterDensity, materialImpact, airFlowfieldRotation, waterFlowfieldRotation, airForce, waterForce, dt, shearingStressConstantK, shearingRatePower, erosionPowerValue, criticalShearValue, posAndDirs, flowType, waterFlow, airFlow, densityUsed, densityMap, initialCapacity);
    } else if (applyOn == EROSION_APPLIED::LAYER_TERRAIN) {
        return this->ApplyOnAnyTerrain(layerBasedGrid, particleSimulationTime, terrainModifTime, startingPoint, originalDirection, randomnessFactor, fallFromSky, gravity, bouncingCoefficient, bounciness, minSpeed, maxSpeed, maxCapacityFactor, erosion, deposit, matterDensity, materialImpact, airFlowfieldRotation, waterFlowfieldRotation, airForce, waterForce, dt, shearingStressConstantK, shearingRatePower, erosionPowerValue, criticalShearValue, posAndDirs, flowType, waterFlow, airFlow, densityUsed, densityMap, initialCapacity);
    }
    return {{}, -1, -1};
}

std::tuple<std::vector<BSpline>, int, int>
UnderwaterErosion::ApplyOnAnyTerrain(TerrainModel *terrain, float &particleSimulationTime, float &terrainModifTime, Vector3 startingPoint, Vector3 originalDirection, float randomnessFactor, bool fallFromSky,
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
//                                     std::function<Vector3(Vector3)> flowfieldFunction)
                                     FLOWFIELD_TYPE flowType,
                                     Matrix3<Vector3> waterFlow,
                                     Matrix3<Vector3> airFlow,
                                     DENSITY_TYPE densityUsed,
                                     Matrix3<float> densityMap,
                                     float initialCapacity)
{
    VoxelGrid* asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    Heightmap* asHeightmap = dynamic_cast<Heightmap*>(terrain);
    ImplicitPatch* asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    LayerBasedGrid* asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    Vector3 terrainSize = voxelGrid->getDimensions(); //terrain->getDimensions();
    float starting_distance = pow(terrainSize.maxComp()/2.f, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap

    // Hydraulic erosion parameters
    float strengthValue = this->maxRockStrength;
    float maxCapacity = 1.f * strengthValue * maxCapacityFactor;
    float erosionFactor = .01f * strengthValue * erosion;
    float depositFactor = .01f * strengthValue * deposit;

    float particleSize = this->maxRockSize;
    float radius = particleSize * .5f;

    Matrix3<float> initialMapValues;
    if (asVoxels) {
        initialMapValues = asVoxels->getVoxelValues();
    } else if (asHeightmap) {
        initialMapValues = asHeightmap->getHeights();
    } else if (asImplicit) {
        asImplicit->_cached = false;
        initialMapValues = asImplicit->getVoxelized(terrainSize);
    } else if (asLayers) {
        initialMapValues = asLayers->voxelize(terrainSize.z, 3.f);
    }

    initialMapValues.raiseErrorOnBadCoord = false;
    initialMapValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    // Use normals from initial map, it shouldn't change too much...
    Matrix3<Vector3> normals = initialMapValues.gradient(); //initialMapValues.binarize().toDistanceMap().gradient();
    normals.raiseErrorOnBadCoord = false;

    if (asHeightmap) {
        normals = -normals + Vector3(0, 0, 1); // Add the Z-component for heightmaps
    }

    for (auto& n : normals)
        n.normalize();

    Matrix3<float> modifications(initialMapValues.getDimensions());

    Matrix3<Vector3> gravityfieldValues = Matrix3<Vector3>(voxelGrid->getDimensions(), Vector3(0, 0, -gravity));
    gravityfieldValues.raiseErrorOnBadCoord = false;
    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);

    Matrix3<float> environmentalDensities = voxelGrid->getEnvironmentalDensities(); // Here also
    Matrix3<Vector3> flowfieldValues = Matrix3<Vector3>(environmentalDensities.getDimensions());
    if (flowType == FLOWFIELD_TYPE::BASIC) {
        Vector3 airDir = Vector3(0, airForce, 0).rotate(0, 0, (airFlowfieldRotation / 180) * PI);
        Vector3 waterDir = Vector3(0, waterForce, 0).rotate(0, 0, (waterFlowfieldRotation / 180) * PI);
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            if (voxelGrid->environmentalDensities.at(i) < 100) { // In the air
                flowfieldValues.at(i) = airDir;
            } else { // In water
                flowfieldValues.at(i) = waterDir;
            }
        }
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
        airFlow = airFlow.resize(terrainSize.x, terrainSize.y, 1) * airForce * 10.f;
        waterFlow = waterFlow.resize(terrainSize.x, terrainSize.y, 1) * waterForce * 10.f;
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            if (voxelGrid->environmentalDensities.at(i) < 100) { // In the air
                flowfieldValues.at(i) = airFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            } else { // In water
                flowfieldValues.at(i) = waterFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            }
        }
    } else if (flowType == FLOWFIELD_TYPE::FLUID_SIMULATION) {
//        flowfieldValues = voxelGrid->getFlowfield().resize(terrainSize) * airForce;
        for (int x = 0; x < flowfieldValues.sizeX; x++) {
            for (int y = 0; y < flowfieldValues.sizeY; y++) {
                for (int z = 0; z < flowfieldValues.sizeZ; z++) {
                    flowfieldValues.at(x, y, z) = Vector3(0, std::cos(3.f * float(y) / terrainSize.y), 0) * 3.f * airForce;
//                    std::cout << flowfieldValues.at(x, y, z) << "\n";
                }
            }
        }
        for (auto& f : flowfieldValues) {
            f = (f + Vector3::random(.1f)).normalize();
        }
    }
    flowfieldValues.raiseErrorOnBadCoord = false;
    flowfieldValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    std::vector<BSpline> tunnels(rockAmount);
    std::vector<Matrix3<float>> submodifications(rockAmount, Matrix3<float>(modifications.getDimensions()));
    std::vector<int> nbPos(rockAmount), nbErosions(rockAmount);
    std::vector<std::vector<std::pair<float, Vector3>>> allErosions(rockAmount);

    auto startTime = std::chrono::system_clock::now();

//    std::cout << densityMap << std::endl;

    #pragma omp parallel for
    for (int i = 0; i < this->rockAmount; i++)
    {
        std::vector<std::pair<float, Vector3>> erosionValuesAndPositions;

        BSpline tunnel;
//        Matrix3<float> currentMapValues = initialMapValues;
//        currentMapValues.raiseErrorOnBadCoord = false;
//        currentMapValues.defaultValueOnBadCoord = -1000;

//        float capacity = 0.f;
//        float capacity = maxCapacity * .5f; //0.f;
        float capacity = maxCapacity * initialCapacity; //0.f;

        float flowfieldInfluence = 1.0;
        int steps = 10 * starting_distance; // An estimation of how many step we need
        auto [pos, dir] = posAndDirs[i];

//        std::vector<Vector3> coords;

        bool touched = false;
        bool hasBeenAtLeastOnceInside = false;
        bool firstHit = true;
        while (!touched) {
            Vector3 nextPos = pos + (dir * dt).maxMagnitude(1.f);
            if (Vector3::isInBox(nextPos.xy(), Vector3(), terrainSize.xy())) {
                float environmentDensity = environmentalDensities.at(nextPos);
                // TODO !!!!
                // Gravity + Buoyancy
                float particleMass = 1.f, particleVolume = .1f;
                float gravityForce = gravity * particleMass;
                float boyancyForce = -environmentDensity * particleVolume; // B = - rho_fluid * V * g
                float gravityCoefficient = particleVolume * gravity * (matterDensity - environmentDensity); // (gravityForce + boyancyForce);
//                std::cout << matterDensity << " " << environmentDensity << " " << gravityCoefficient << "\n";
//                std::cout  << gravityCoefficient << std::endl;
            //    float gravityCoefficient = gravity * particleVolume * (particleDensity - environmentDensity); // Same as F + B
            //                std::cout << environmentDensity << " -> " << particleVolume << " * (" << environmentDensity << " - " << matterDensity << ") => " << gravityCoefficient << "\n";
                // float gravityCoefficient = std::max(1.f - (environmentDensity / matterDensity), -1.f); // Keep it between -1 and 1
//                Vector3 acceleration = /*flowfieldValues.at(position) +*/ gravityDirection * (gravityForce * gravityCoefficient);
//                float gravityCoefficient = std::max(1.f - (environmentDensity / matterDensity), -1.f); // Keep it between -1 and 1
                Vector3 justFlow = (flowfieldValues.at(pos) + flowfieldValues.at(nextPos)) * .5f;
                Vector3 justGravity = (gravityfieldValues.at(nextPos) * gravityCoefficient).maxMagnitude(2.f);
                Vector3 flowfield = justFlow + justGravity;
                dir += flowfield * flowfieldInfluence * dt;
                hasBeenAtLeastOnceInside = true;
                bool justHit = false;
                if (terrain->checkIsInGround(nextPos)) { //currentMapValues.at(nextPos) > 0.0) { // Hit a wall
                    Vector3 normal;
                    if (asHeightmap)
                        normal = normals.at(nextPos.xy());
                    else
                        normal = normals.at(nextPos);

                    float resistance = 0.f;
                    if (densityUsed == DENSITY_TYPE::DENSITY_IMAGE)
                        resistance = /*1.f - */densityMap.at(nextPos.xy()) * materialImpact;

//                    std::cout << resistance << "\n";

//                    float speedRate = (dir.norm() / maxSpeed);
//                    float coef = 1.f; // std::max(std::abs(normal.x), std::abs(normal.y));
                    // Using SPH formula
                    // Erosion part
                    float l = particleSize * 0.001f;
                    float vRel = dir.norm() * (1.f - std::abs(dir.normalized().dot(normal))); // velocity relative to the surface
                    float theta = vRel / l;
                    float shear = shearingStressConstantK * std::pow(theta, shearingRatePower);
                    float amountToErode = erosionFactor * std::pow(std::max(0.f, shear)/* - criticalShearValue)*/, erosionPowerValue) * (1.f - resistance); // std::pow(std::max(0.f, shear - criticalShearStress), erosionPowerValue)  * (1.f - resistanceValue * materialImpact);
                    amountToErode = std::min(amountToErode, maxCapacity - capacity);

//                    // Deposition part
                    float densitiesRatio = 1.f; //1000.f / matterDensity; // The ground should have a density, I'll assume the value 1000 for now
                    float u = (2.f / 9.f) * radius * radius * (matterDensity - environmentDensity) * gravity * (1.f - (capacity / maxCapacity));
//                    float amountToDeposit = std::min(capacity, depositFactor * u);
                    float amountToDeposit = std::min(capacity, depositFactor * (densitiesRatio * capacity * (capacity / maxCapacity))); //maxCapacity));

                    if (!firstHit || true) {
                        if (amountToErode - amountToDeposit != 0) {
                            capacity += (amountToErode - amountToDeposit);
//                            RockErosion(particleSize, amountToErode - amountToDeposit).computeErosionMatrix(submodifications[i], nextPos);
                            erosionValuesAndPositions.push_back({amountToErode - amountToDeposit, nextPos});
                        }
                    }

                    // Continue the rock tracing
                    if (!justHit) {
                        Vector3 bounce = dir.reflexion(normal);
                        bounce -= normal * bounce.dot(normal) * (1.f - bounciness);
                        dir = bounce * bouncingCoefficient;
                    }
                    justHit = true;
                }
            }
            else {
                dir += Vector3(0, 0, -1) * dt;
                dir += Vector3::random(randomnessFactor);
            }
            steps --;
            tunnel.points.push_back(pos);
//            do {
                pos += dir.clamped(0.f, 1.f);
//            } while (terrain->checkIsInGround(pos));

            if ((hasBeenAtLeastOnceInside && !Vector3::isInBox(nextPos.xy(), Vector3(), terrainSize.xy())) || steps < 0 /*|| dir.norm() < 1e-3*/) {
                if ((steps < 0 || dir.norm2() < 1e-3) && depositFactor > 0.f) {
//                    RockErosion(particleSize, -capacity).computeErosionMatrix(submodifications[i], pos - Vector3(0, 0, particleSize * .2f), false);
                    erosionValuesAndPositions.push_back({-capacity, pos - Vector3(0, 0, particleSize - .2f)});
                }
                break;
            }
        }
        tunnels[i] = tunnel;
//        modifications += submodifications[i];
        nbPos[i] = tunnel.points.size();
        nbErosions[i] = erosionValuesAndPositions.size();

        allErosions[i] = erosionValuesAndPositions;
    }
    auto endParticleTime = std::chrono::system_clock::now();

    std::vector<ImplicitNaryOperator*> allNary(allErosions.size());
    if (asImplicit) {
        for (size_t i = 0; i < allNary.size(); i++) {
            allNary[i] = new ImplicitNaryOperator;
        }
    }

    #pragma omp parallel for
    for (size_t i = 0; i < allErosions.size(); i++) {
        const auto& erosionValuesAndPositions = allErosions[i];
        for (auto& [val, pos] : erosionValuesAndPositions) {
            if (asVoxels) {
                RockErosion(particleSize, val).computeErosionMatrix(submodifications[i], pos);
            } else if (asHeightmap) {
                RockErosion(particleSize, val).computeErosionMatrix2D(submodifications[i], pos);
            } else if (asImplicit) {
                float dimensions = particleSize * val * 5.f;
                auto sphere = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(dimensions, dimensions, dimensions), val));
                sphere->material = (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT);
                sphere->dimensions = Vector3(dimensions, dimensions, dimensions);
                sphere->position = pos - sphere->dimensions * .5f;
                allNary[i]->composables.push_back(sphere);
            } else if (asLayers) {
                for (int x = -particleSize; x < particleSize; x++) {
                    for (int y = -particleSize; y < particleSize; y++) {
//                        for (int z = -particleSize; z < particleSize; z++) {
                        if (x*x + y*y < particleSize * particleSize) {
                            float halfHeight = val * std::sqrt(particleSize*particleSize - (x*x + y*y));
                            asLayers->transformLayer(pos.x + x, pos.y + y, pos.z - halfHeight, pos.z + halfHeight, (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT));
                        }
//                        }
                    }
                }
//                asLayers->cleanLayers();
                asLayers->_historyIndex++;
            }
        }
    }

    for (const auto& sub : submodifications)
        modifications += sub;

//    std::cout << "Total erosion : " << modifications.sum() << std::endl;
    if (asVoxels) {
        for (int x = 0; x < modifications.sizeX; x++)
            for (int y = 0; y < modifications.sizeY; y++)
                for (int z = 0; z < 8; z++)
                    modifications.at(x, y, z) = 0;
        asVoxels->applyModification(modifications);
    } else if (asHeightmap) {
        asHeightmap->heights += modifications;
    } else if (asImplicit) {
        ImplicitNaryOperator* totalErosion = new ImplicitNaryOperator;
        for (auto& op : allNary)
            totalErosion->composables.push_back(op);
        totalErosion->composables.push_back(implicitTerrain->copy());
//        delete implicitTerrain;
        implicitTerrain->~ImplicitPatch();
        implicitTerrain = new (implicitTerrain) ImplicitNaryOperator;
        dynamic_cast<ImplicitNaryOperator*>(implicitTerrain)->composables = totalErosion->composables;
        implicitTerrain->_cached = false;
    } else if (asLayers) {
        // ...
    }

    int positions = 0, erosions = 0;
    for (int i = 0; i < rockAmount; i++) {
        positions += nbPos[i];
        erosions += nbErosions[i];
    }
    auto endModificationTime = std::chrono::system_clock::now();

    particleSimulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(endParticleTime - startTime).count();
    terrainModifTime = std::chrono::duration_cast<std::chrono::milliseconds>(endModificationTime - endParticleTime).count();
    return {tunnels, positions, erosions};







    /*
    VoxelGrid* asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    Heightmap* asHeightmap = dynamic_cast<Heightmap*>(terrain);
    ImplicitPatch* asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    LayerBasedGrid* asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    Matrix3<float> environmentalDensities = voxelGrid->getEnvironmentalDensities(); // Here also
    Matrix3<Vector3> flowfieldValues;
    Vector3 airDir = Vector3(0, airForce, 0).rotate(0, 0, ((airFlowfieldRotation + random_gen::generate(-10, 10)) / 180) * PI);
    Vector3 waterDir = Vector3(0, waterForce, 0).rotate(0, 0, ((waterFlowfieldRotation + random_gen::generate(-10, 10)) / 180) * PI);
    if (flowType == FLOWFIELD_TYPE::FLUID_SIMULATION) {
        flowfieldValues = voxelGrid->getFlowfield() * 30.f;
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
        flowfieldValues = Matrix3<Vector3>(environmentalDensities.getDimensions());
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            flowfieldValues.at(i) = (environmentalDensities.at(i) > 100 ? waterFlow.at(flowfieldValues.getCoordAsVector3(i).xy()) * waterForce : airFlow.at(flowfieldValues.getCoordAsVector3(i).xy()) * airForce);
        }
    } else {
        flowfieldValues = Matrix3<Vector3>(environmentalDensities.getDimensions());
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            flowfieldValues.at(i) = (environmentalDensities.at(i) > 100 ? waterDir * waterForce : airDir * airForce);
        }
    }

//    if (flowfieldFunction == nullptr)
//        flowfieldFunction = [&](Vector3 pos) { return (environmentalDensities.at(pos) < 100 ? airDir : waterDir); };

//    for (size_t i = 0; i < flowfieldValues.size(); i++) {
//        flowfieldValues.at(i) = flowfieldFunction(flowfieldValues.getCoordAsVector3(i));
//    }
    Matrix3<float> terrainDensities = Matrix3<float>(environmentalDensities.getDimensions(), 1.f);
    if (densityUsed == DENSITY_TYPE::NATIVE) {
        if (asVoxels)
            terrainDensities = voxelGrid->getVoxelValues();
        else if (asHeightmap || asImplicit)
            terrainDensities = Matrix3<float>(environmentalDensities.getDimensions(), 0.f);
    } else if (densityUsed == DENSITY_TYPE::DENSITY_IMAGE) {
//        std::cout << "Using image as densities:\n" << densityMap.displayValues() << std::endl;
        for (size_t i = 0; i < terrainDensities.size(); i++)
            terrainDensities.at(i) = densityMap.at(terrainDensities.getCoordAsVector3(i).xy());
    }


    auto startTime = std::chrono::system_clock::now();
    int numberOfParticles = this->rockAmount;
    float particleSize = this->maxRockSize;
    randomnessFactor = .1f;
    float starting_distance = pow(terrain->getDimensions().maxComp()/2.0, 2); // pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ))/2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap

    // Hydraulic erosion parameters
    float particleVolume = (4.f / 3.f) * PI * (particleSize * .5f) * 0.001; // Volume = 4/3 * pi * R
    float particleMass = matterDensity * particleVolume * 0.1f; // rho = m/V -> m = rho * V
//    float strengthVolumeDivisor = particleSize*particleSize*particleSize;
    float strength = this->maxRockStrength; // / particleVolume; // strengthVolumeDivisor;
    float maxCapacity = particleSize * maxCapacityFactor;
    float erosionFactor = .01f * strength * erosion; // / (10.f * particleSize);
    float depositFactor = .01f * strength * deposit; // / (10.f * particleSize);

    Matrix3<float> initialMapValues;
    Matrix3<std::vector<std::pair<TerrainTypes, float>>> initialLayers;
    if (asVoxels) {
        initialMapValues = asVoxels->getVoxelValues();
        for (int x = 0; x < initialMapValues.sizeX; x++)
            for (int y = 0; y < initialMapValues.sizeY; y++)
                initialMapValues.at(x, y, 0) = 1.f;
    } else if (asHeightmap) {
        initialMapValues = asHeightmap->getHeights();
    } else if (asImplicit) {
        initialMapValues = asImplicit->getVoxelized();
    } else if (asLayers) {
        initialMapValues = asLayers->voxelize();
        initialLayers = asLayers->getLayers();
        initialLayers.raiseErrorOnBadCoord = false;
        initialLayers.returned_value_on_outside = REPEAT_VALUE;
    }

    initialMapValues.raiseErrorOnBadCoord = false;
    initialMapValues.returned_value_on_outside = REPEAT_VALUE;
    // Use normals from initial map, it shouldn't change too much...
    Matrix3<Vector3> normals = initialMapValues.gradient(); // initialMapValues.gradient(); //initialMapValues.binarize().toDistanceMap().gradient();
    normals.raiseErrorOnBadCoord = false;

    if (asVoxels) {
    } else if (asHeightmap) {
        for (auto& n : normals) {
            n = (-n + Vector3(0, 0, 1)); // Add the Z component for heightmaps
        }
    } else if (asImplicit) {

    } else if (asLayers) {

    }

    for (auto& n : normals) {
        n.normalize();
    }

    Matrix3<float> modifications;

    if (asVoxels) {
        modifications = Matrix3<float>(terrain->getSizeX(), terrain->getSizeY(), terrain->getSizeZ());
    } else if (asHeightmap) {
        modifications = Matrix3<float>(terrain->getSizeX(), terrain->getSizeY());
    } else if (asImplicit) {
        modifications = Matrix3<float>(terrain->getSizeX(), terrain->getSizeY(), terrain->getSizeZ());
    } else if (asLayers) {
        modifications = Matrix3<float>(terrain->getSizeX(), terrain->getSizeY(), terrain->getSizeZ()); //?
    }
    modifications.raiseErrorOnBadCoord = false;

    Matrix3<Vector3> gravityfieldValues = Matrix3<Vector3>(1, 1, 1, Vector3(0, 0, -gravity)); // Will be applied everywhere
    gravityfieldValues.raiseErrorOnBadCoord = false;
    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);


    environmentalDensities.raiseErrorOnBadCoord = false;
    environmentalDensities.defaultValueOnBadCoord = 1.f;

    std::vector<std::vector<BSpline>> tunnels(numberOfParticles);

    std::vector<std::vector<RockErosion>> erosions(numberOfParticles);
    std::vector<std::vector<float>> erosionsStrength(numberOfParticles);
    std::vector<std::vector<Vector3>> erosionsPositions(numberOfParticles);

    Matrix3<float> currentMapValues = initialMapValues;
    currentMapValues.raiseErrorOnBadCoord = false;
    currentMapValues.defaultValueOnBadCoord = -1000;

    // Given in SPH
    float criticalShearStress = shearingStressConstantK * std::pow(criticalShearValue, shearingRatePower);

    Vector3 gravityForce = Vector3(0, 0, -1.f);
    int steps = 100 * 10 * starting_distance; // / dt; // An estimation of how many step we need

    Vector3 terrainDims = Vector3(terrain->getSizeX(), terrain->getSizeY(), terrain->getSizeZ());

    flowfieldValues.raiseErrorOnBadCoord = false;
    environmentalDensities.raiseErrorOnBadCoord = false;
    terrainDensities.raiseErrorOnBadCoord = false;

    flowfieldValues.returned_value_on_outside = REPEAT_VALUE;

    std::vector<Matrix3<float>> subtotalModifications(numberOfParticles, Matrix3<float>(modifications.getDimensions()));
    #pragma omp parallel for
    for (int i = 0; i < numberOfParticles; i++)
    {
        BSpline tunnel;
//        currentMapValues = initialMapValues + modifications;

        float capacity = 0.f;

        Vector3 acceleration;
        Vector3 velocity;
        Vector3 position;

        position = posAndDirs[i].first;
        velocity = posAndDirs[i].second;
//        setInitialPositionAndVelocityAndAccelerationOfParticle(position, velocity, acceleration, startingPoint, starting_distance, originalDirection, randomnessFactor, fallFromSky, terrainDims);
        Vector3 previousPosition = position;
        tunnel.points.push_back(position);

        bool hasBeenAtLeastOnceInside = false;
        bool atLeastOneBounce = false;

        std::vector<Vector3> lastPositions(10, Vector3());

        for (int step = 0; step < steps; step++) {
            bool earlyStopSimulation = false;
//            std::cout << position << std::endl;
            previousPosition = position;


//            int staticPos = 0;
//            for (const auto& p : lastPositions) {
//                staticPos += ((p-position).norm2() < 1e-3 ? 1 : 0);
//            }
//            lastPositions.erase(lastPositions.begin());
//            lastPositions.push_back(position);

//            Vector3 flowfield = flowfieldValues.at(position); // flowfieldFunction(position);
//            updateParticlePosition(position, velocity, acceleration, dt, terrain, flowfield, true); // atLeastOneBounce);

            if (terrain->contains(position)) {
                float environmentDensity = environmentalDensities.at(position);
                Vector3 flow = flowfieldValues.at(position); // flowfieldFunction(position);
                acceleration = getExternalForce(position,
                                                gravity,
                                                particleMass,
                                                particleVolume,
                                                matterDensity,
                                                environmentDensity,
                                                flow, //flowfieldValues,
                                                gravityForce,
                                                capacity,
                                                maxCapacity);

                hasBeenAtLeastOnceInside = true;

                if (terrain->checkIsInGround(position)) { // Hit a wall
                    atLeastOneBounce = true;
                    float resistanceValue = terrainDensities.at(position);
                    Vector3 averageVelocity(position - lastPositions.back());
                    for (size_t i = 0; i < lastPositions.size() - 1; i++) {
                        averageVelocity += (lastPositions[i+1] - lastPositions[i]);
                    }
                    averageVelocity /= (float)(lastPositions.size());

                    Vector3 normal;
                    if (asVoxels) {
                        normal = (normals.at(position - velocity) + normals.at(position)).normalized();
                    } else if (asHeightmap) {
                        normal = (normals.at((position - velocity).xy()) + normals.at(position.xy())).normalized();
                    } else if (asImplicit) {
                        normal = asImplicit->getNormal(position);
                    } else if (asLayers) {
                        normal = (normals.at(position - velocity) + normals.at(position)).normalized();
                    }
//                    float speedRate = (velocity.norm2() / (maxSpeed * maxSpeed));
//                    float coef = 1.f; // std::max(std::abs(normal.x), std::abs(normal.y));
                    // Using SPH formula
                    // Erosion part
                    float l = particleSize * 0.001f;
                    float vRel = averageVelocity.norm() * (1.f - std::abs(averageVelocity.normalized().dot(normal))); // velocity relative to the surface
//                    float vRel = velocity.norm2() * (1.f - std::abs(velocity.normalized().dot(normal))); // velocity relative to the surface
                    float theta = vRel / l;
                    float shear = shearingStressConstantK * std::pow(theta, shearingRatePower);
                    float amountToErode = erosionFactor * std::pow(std::max(0.f, shear - criticalShearStress), erosionPowerValue)  * (1.f - resistanceValue * materialImpact);
                    amountToErode = std::min(amountToErode, maxCapacity - capacity);
                    if (std::abs(normal.z) == 1.f)
                        amountToErode = 0;

                    // Deposition part
                    float densitiesRatio = terrainDensities.at(position) * 1000.f / matterDensity; // The ground should have a density, I'll assume the value 1000 for now
                    float amountToDeposit = std::min(capacity, depositFactor * (densitiesRatio * capacity * (capacity / maxCapacity))); //maxCapacity));

                    if (std::abs(amountToErode - amountToDeposit) > 0) {
                        capacity += (amountToErode - amountToDeposit);
                        if (amountToErode > 0) {
                            RockErosion erosionRock(particleSize, amountToErode);
                            erosions[i].push_back(erosionRock);
                            erosionsStrength[i].push_back(amountToErode);
                            erosionsPositions[i].push_back(position + velocity);
                        }
                        if (amountToDeposit > 0) {
                            RockErosion depositRock(particleSize, -amountToDeposit);
                            erosions[i].push_back(depositRock);
                            erosionsStrength[i].push_back(-amountToDeposit);
                            erosionsPositions[i].push_back(position - velocity - Vector3(0, 0, particleSize*.5f));
                        }
                    }

                    Vector3 bounce = velocity.normalized().reflexion(normal + Vector3::random(1.f)) * velocity.norm();
                    bounce -= normal * bounce.dot(normal) * (1.f - bounciness);
                    Vector3 prevVel = velocity;
                    velocity = bounce * bouncingCoefficient / dt;
                    acceleration *= 0;
                    while (terrain->checkIsInGround(position)) {
                        if (prevVel.norm2() == 0) {
                            earlyStopSimulation = true;
                            break;
                        }
                        position -= prevVel; // Move particle until it's not in the terrain anymore
                    }
                    previousPosition = position - velocity;
                }
            } else {
                acceleration = gravityForce * gravity;
            }

            int staticPos = 0;
            for (const auto& p : lastPositions) {
                staticPos += ((p-position).norm2() < 1e-3 ? 1 : 0);
            }
            lastPositions.erase(lastPositions.begin());
            lastPositions.push_back(position);

            Vector3 flowfield = flowfieldValues.at(position); // flowfieldFunction(position);
            updateParticlePosition(position, velocity, acceleration, dt, terrain, flowfield, true); // atLeastOneBounce);
            tunnel.points.push_back(position);

//            int _i = 0;
//            Vector3 meanFirsts, meanLasts;
//            for (const auto& p : lastPositions) {
//                if (_i < std::floor(lastPositions.size() / 2)) meanFirsts += p;
//                if (_i > std::ceil(lastPositions.size() / 2)) meanLasts += p;
//                _i++;
//            }
//            if ((meanFirsts - meanLasts).norm2() < 1e-8) earlyStopSimulation = true;

            bool speedTooCloseToNull = earlyStopSimulation || (staticPos > 3); // false; // (consecutiveStepsAtLowSpeed >= 8 && atLeastOneHit); // && previousVelocity.norm() < 1e-8);
            bool goingOutside = (hasBeenAtLeastOnceInside && !terrain->contains(position));
            if (goingOutside || speedTooCloseToNull) { // Is getting out of the map
                if (speedTooCloseToNull && depositFactor > 0.f) {
                    float lastDeposit = -capacity*.25f; // 0;
                    RockErosion rock(particleSize*2, lastDeposit);
                    erosions[i].push_back(rock);
                    erosionsStrength[i].push_back(lastDeposit);
                    erosionsPositions[i].push_back(position - Vector3(0, 0, particleSize*.5f));
                    capacity = 0;
                } else {

                }
                break;
            }
        }
        tunnels[i].push_back(tunnel);
        if (capacity > 0 && depositFactor > 0.f && terrain->contains(position)) {
            RockErosion rock(particleSize, -capacity);
            erosions[i].push_back(rock);
            erosionsStrength[i].push_back(-capacity);
            erosionsPositions[i].push_back(position - Vector3(0, 0, particleSize*.5f));
        }


        for (size_t iRock = 0; iRock < erosionsPositions[i].size(); iRock++) {
            RockErosion& rock = erosions[i][iRock];
            float& strength = erosionsStrength[i][iRock];
            Vector3& pos = erosionsPositions[i][iRock];
            if (asVoxels) {
//                rock.maxStrength *= .5f;
                rock.computeErosionMatrix(subtotalModifications[i], pos);
            } else if (asHeightmap) {
                rock.computeErosionMatrix2D(subtotalModifications[i], pos.xy()); // Bring it in 2D to have an effect on the heightmap
            }
        }

    }

    particleSimulationTime = std::chrono::duration<float>(std::chrono::system_clock::now() - startTime).count();
    int totalNumberOfErosionToCompute = 0;
    float totalErosion = 0.f;
    float totalDeposition = 0.f;

    if (asImplicit) {
        ImplicitNaryOperator* erosionPatch = new ImplicitNaryOperator;
        for (size_t iPath = 0; iPath < erosionsPositions.size(); iPath++) {
                for (size_t iRock = 0; iRock < erosionsPositions[iPath].size(); iRock++) {
                    RockErosion& rock = erosions[iPath][iRock];
                    float& strength = erosionsStrength[iPath][iRock];
                    totalErosion -= std::min(0.f, strength);
                    totalDeposition += std::max(0.f, strength);
                    Vector3& pos = erosionsPositions[iPath][iRock];
                    if (asVoxels) {
                        rock.computeErosionMatrix(modifications, pos);
                    } else if (asHeightmap) {
                        rock.computeErosionMatrix(modifications, pos.xy()); // Bring it in 2D to have an effect on the heightmap
                    } else if (asImplicit) {
                        if (std::abs(strength) > 0.01f && iRock < 30) {
                            float sphereSize = particleSize * std::abs(std::max(1.f, strength * 10.f));
                            ImplicitPrimitive* hole = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::PredefinedShapes::Sphere, Vector3(sphereSize, sphereSize, sphereSize), 0.f));
                            if (erosionsStrength[iPath][iRock] > 0.f) {
                                hole->material = TerrainTypes::AIR;
                            } else {
                                hole->material = TerrainTypes::DIRT;
                            }
                            hole->dimensions = Vector3(sphereSize, sphereSize, sphereSize);
                            hole->position = pos - hole->dimensions * .5f;
                            erosionPatch->composables.push_back(hole);
                        }
                    }
                    totalNumberOfErosionToCompute++;
                }
    //        }
        }

        ImplicitOperator* newTerrain = new ImplicitOperator;
        newTerrain->blendingFactor = 2.f;
        newTerrain->composeFunction = ImplicitPatch::CompositionFunction::BLEND;
        newTerrain->positionalB = ImplicitPatch::PositionalLabel::FIXED_POS;
        newTerrain->composableA = asImplicit->copy();
        newTerrain->composableB = erosionPatch;
        asImplicit = newTerrain;
//        *implicitTerrain = *newTerrain;
    } else if (asLayers) {
        for (size_t iPath = 0; iPath < erosionsPositions.size(); iPath++) {
            for (size_t iRock = 0; iRock < erosionsPositions[iPath].size(); iRock++) {
                RockErosion& rock = erosions[iPath][iRock];
                float& strength = erosionsStrength[iPath][iRock];
                Vector3& position = erosionsPositions[iPath][iRock];

                for (int x = std::floor(position.x - particleSize); x < std::ceil(position.x + particleSize); x++) {
                    for (int y = std::floor(position.y - particleSize); y < std::ceil(position.y + particleSize); y++) {
                        auto& column = asLayers->getLayers().at(x, y);
                        float dist = (Vector3(x, y) - position.xy()).norm();
                        float erosionHeight = std::sqrt(particleSize*particleSize - dist*dist);
                        float minHeight = position.z - erosionHeight;
                        float maxHeight = position.z + erosionHeight;
                        float currentHeight = 0.f;
                        for (auto& [mat, height] : column) {
                            float z0 = std::max(minHeight, currentHeight);
                            float z1 = std::min(maxHeight, currentHeight + height);
//                            float amountToRemove = z1 - z0;
                            asLayers->transformLayer(x, y, z0, z1, (strength > 0 ? AIR : SAND));
                        }
                    }
                }
            }
        }
        asLayers->cleanLayers();
    } else {
        for (auto& modif : subtotalModifications)
            modifications += modif;
    }
    if (asVoxels) {
        for (int x = 0; x < modifications.sizeX; x++)
            for (int y = 0; y < modifications.sizeY; y++)
                for (int z = 0; z < 3; z++)
                    modifications.at(x, y, z) = 0.f;
        this->voxelGrid->applyModification(modifications);
        this->voxelGrid->limitVoxelValues(1.f);
        this->voxelGrid->saveState(); // Just to gain time
    } else if (asHeightmap) {
        asHeightmap->heights += modifications.meanSmooth(3, 3, 1);
    } else if (asImplicit) {
        // Nothing to do
        terrain = asImplicit;
        this->implicitTerrain = asImplicit;
    }

    auto endTime = std::chrono::system_clock::now();
//    std::cout << totalErosion << " " << totalDeposition << std::endl;
    auto flatTunnels = flattenArray(tunnels);
    int totalNumberOfEvaluationPoints = 0;
    for (const auto& path : flatTunnels)
        totalNumberOfEvaluationPoints += path.points.size();
    terrainModifTime = std::chrono::duration<float>(endTime - startTime).count() - particleSimulationTime;

    return {flatTunnels, totalNumberOfEvaluationPoints, totalNumberOfErosionToCompute};*/
}



















std::vector<Vector3> UnderwaterErosion::CreateTunnel(int numberPoints, bool addingMatter, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    BSpline curve = BSpline(numberPoints); // Random curve
    for (Vector3& coord : curve.points)
        coord = ((coord + Vector3(1.0, 1.0, 1.0)) / 2.0) * voxelGrid->getDimensions(); //Vector3(grid->sizeX, grid->sizeY, grid->sizeZ);
    return CreateTunnel(curve, addingMatter, true, applyChanges, startingShape, endingShape);
}
std::vector<Vector3> UnderwaterErosion::CreateTunnel(BSpline path, bool addingMatter, bool usingSpheres, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    Matrix3<float> erosionMatrix(voxelGrid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
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
        Matrix3<float> holeMatrix;
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
    Matrix3<float> erosionMatrix(voxelGrid->getDimensions()); // this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
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
            Matrix3<float> holeMatrix;
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
    voxelGrid->applyModification(erosionMatrix);
//    if (applyChanges)
//        grid->remeshAll();
    return allCoords;
}

std::vector<Vector3> UnderwaterErosion::CreateCrack(Vector3 start, Vector3 end, bool applyChanges)
{
    float rx = 3.f, ry = 3.f, rz = 3.f;
    Vector3 ratio(rx, ry, rz);
    Matrix3<int> resizedMap = this->voxelGrid->getVoxelValues().resize(voxelGrid->getDimensions() / ratio).binarize();
//    Matrix3<int> resizedMap = this->grid->getVoxelValues().resize(this->grid->sizeX / rx, this->grid->sizeY / ry, this->grid->sizeZ / rz).binarize();
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
    Matrix3<float> erosionMatrix(voxelGrid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
    bool modificationDoesSomething = true;
    BSpline width = BSpline(std::vector<Vector3>({
                                                     Vector3(0.0, 1.0),
                                                     Vector3(0.5, 2.0),
                                                     Vector3(1.0, 1.0)
                                                 }));

    std::vector<Vector3> coords;
    Matrix3<float> holeMatrix;
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

