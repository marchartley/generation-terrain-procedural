#include "ParticleErosion.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>

#include "DataStructure/BVH.h"
#include "TerrainModification/RockErosion.h"
#include "EnvObject/EnvironmentalScene.h"

Vector3 ErosionParticle::predictNextPos(float dt)
{
    return pos + velocity * dt + force * (.5f * dt * dt);
}

std::vector<ParticleContact> ErosionParticle::bounce(float dt, SpacePartitioning* boundaries,
                                                   float restitution, float damping, int maxBounces)
{
    std::vector<ParticleContact> contacts;
    Vector3 driftVelocity = velocity + force * (.5f * dt);
    float remaining = dt;
    size_t lastTriangle = size_t(-1);
    restitution = std::clamp(restitution, 0.f, 1.f);
    damping = std::clamp(damping, 0.f, 1.f);
    while (remaining > 0.f) {
        const Vector3 destination = pos + driftVelocity * remaining;
        auto [hit, triangle] = boundaries->getIntersectionAndTriangleIndex(pos, destination, lastTriangle);
        if (!hit.isValid()) { pos = destination; break; }
        Vector3 normal = boundaries->triangles[triangle].normal.normalized();
        if (normal.dot(driftVelocity) > 0.f) normal = -normal;
        const float distance = (destination - pos).norm();
        const float fraction = distance > 0.f ? std::clamp((hit - pos).norm() / distance, 0.f, 1.f) : 0.f;
        remaining *= 1.f - fraction;
        contacts.push_back({hit, normal, driftVelocity});
        driftVelocity = (driftVelocity - normal * ((1.f + restitution) * driftVelocity.dot(normal))) * damping;
        // Separate the ray origin from the surface, without moving the recorded contact.
        pos = hit + normal * 1e-4f;
        lastTriangle = triangle;
        if (int(contacts.size()) >= maxBounces) break;
    }
    velocity = driftVelocity + force * (.5f * dt);
    dir = velocity;
    return contacts;
}

void ErosionParticle::addForce(const Vector3& f) { force += f; }
void ErosionParticle::addVelocity(const Vector3& v) { velocity += v; }

Vector3 ErosionParticle::eulerIntegration(float dt)
{
    velocity += force * dt;
    pos += velocity * dt;
    return pos;
}

Vector3 ErosionParticle::verletIntegration(float dt)
{
    // Velocity Verlet for acceleration frozen during this step.
    pos = predictNextPos(dt);
    velocity += force * dt;
    return pos;
}

ParticleErosion::ParticleErosion()
    : applyOn(DENSITY_VOXELS), particleSimulationTime(0), terrainModifTime(0),
      randomnessFactor(0), fallFromSky(false), gravity(1), bouncingCoefficient(1), bounciness(.5f),
      minSpeed(.1f), maxSpeed(1), maxCapacityFactor(1), erosion(1), deposit(1), matterDensity(1000),
      materialImpact(0), airFlowfieldRotation(0), waterFlowfieldRotation(0), airForce(0), waterForce(0),
      dt(.1f), shearingStressConstantK(1), shearingRatePower(.5f), erosionPowerValue(1), criticalShearValue(0),
      flowType(BASIC), densityUsed(UNUSED), initialCapacity(0), fluidSimType(FluidSimType::LBM),
      wrapPositions(false), applyTheErosion(true), maxCollisions(-1), strength(1), size(1), quantity(0),
      particleSize(1), depositFactor(0), erosionFactor(0)
{}

void ParticleErosion::initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity,
                                         float radius, float density, float initialLoad, float maxCapacity)
{
    particle = ErosionParticle();
    particle.pos = position;
    particle.dir = particle.velocity = velocity;
    particle.force = Vector3();
    particle.capacity = std::clamp(initialLoad, 0.f, maxCapacity);
    particle.properties = std::make_shared<ParticleProperties>();
    particle.properties->density = density;
    particle.properties->maxCapacity = maxCapacity;
    particle.properties->radius = radius;
    particle.properties->volume = (4.f / 3.f) * M_PI * radius * radius * radius;
    particle.properties->mass = density * particle.properties->volume;
}

Vector3 ParticleErosion::computeSettlingVelocity(const ErosionParticle& particle) const
{
    const float density = environmentalDensities.interpolate(particle.pos);
    const float viscosity = density < 100.f ? airViscosity : waterViscosity;
    const float concentration = particle.properties->maxCapacity > 0.f
        ? std::clamp(particle.capacity / particle.properties->maxCapacity, 0.f, 1.f) : 0.f;
    const float hindered = std::pow(1.f - concentration, 5.f);
    // Stokes slip of sediment grains, not the radius/density of the carrier parcel.
    return gravityDefault * ((2.f / 9.f) * sedimentRadius * sedimentRadius *
                            (sedimentDensity - density) / viscosity * hindered);
}

float ParticleErosion::computeErosionValue(const ErosionParticle& particle, float vRel)
{
    const float length = std::max(particle.properties->radius * .002f, 1e-6f);
    const float shear = shearingStressConstantK * std::pow(std::max(0.f, vRel) / length, shearingRatePower);
    const float excess = std::max(0.f, shear - criticalShearValue);
    if (excess == 0.f) return 0.f;
    const float hemisphere = .5f * particle.properties->volume;
    const float detached = erosionFactor * std::pow(excess, erosionPowerValue) * hemisphere * dt;
    return std::clamp(detached, 0.f, std::max(0.f, particle.properties->maxCapacity - particle.capacity));
}

float ParticleErosion::computeDepositionValue(const ErosionParticle& particle)
{
    if (particle.capacity <= 0.f || particle.properties->volume <= 0.f) return 0.f;
    const float area = M_PI * particle.properties->radius * particle.properties->radius;
    // Convert settling flux to a volume exchanged over this contact step.
    const float amount = depositFactor * computeSettlingVelocity(particle).norm() * area * dt *
                         (particle.capacity / particle.properties->volume);
    return std::clamp(amount, 0.f, particle.capacity);
}

ParticleHistory ParticleErosion::trackParticlePositions(ErosionParticle& particle)
{
    ParticleHistory history;
    if (terrain->checkIsInGround(particle.pos)) return history;
    const int maxSteps = int(std::min(500. / dt, double(std::numeric_limits<int>::max() - 1)));
    history.history.reserve(std::min(maxSteps, 1024));
    history.add(particle);
    int collisionsLeft = maxCollisions < 0 ? std::numeric_limits<int>::max() : maxCollisions;
    const float relaxation = std::exp(-.5f * dt / flowResponseTime);
    for (int step = 0; step < maxSteps && collisionsLeft > 0; ++step) {
        const float density = environmentalDensities.interpolate(particle.pos);
        particle.force = gravityDefault * (1.f - density / particle.properties->density);
        const float load = particle.properties->maxCapacity > 0.f
            ? particle.capacity / particle.properties->maxCapacity : 0.f;
        const Vector3 target = flowfieldValues.interpolate(particle.pos) + computeSettlingVelocity(particle) * load;
        // Exact drag relaxation in two half steps around the Verlet drift/kick.
        particle.velocity = target + (particle.velocity - target) * relaxation;
        auto contacts = particle.bounce(dt, boundariesTree, bounciness, bouncingCoefficient,
                                       std::min(8, collisionsLeft));
        particle.velocity = target + (particle.velocity - target) * relaxation;
        for (const auto& contact : contacts) {
            --collisionsLeft;
            const Vector3 endPosition = particle.pos;
            particle.pos = contact.position;
            const Vector3 tangent = contact.velocity - contact.normal * contact.velocity.dot(contact.normal);
            float detached = computeErosionValue(particle, tangent.norm());
            if (!densityMap.empty()) {
                const float resistance = std::clamp(densityMap.interpolate(contact.position), 0.f, 1.f);
                detached *= 1.f - std::clamp(materialImpact, 0.f, 1.f) * resistance;
            }
            const Vector3 settling = computeSettlingVelocity(particle);
            const float orientation = settling.norm2() > 0.f
                ? std::max(0.f, -settling.normalized().dot(contact.normal)) : 0.f;
            const float deposited = computeDepositionValue(particle) * orientation;
            const float exchange = detached - deposited;
            particle.capacity = std::clamp(particle.capacity + exchange, 0.f, particle.properties->maxCapacity);
            ParticleState state(particle);
            state.exchange = {contact.position, exchange};
            state.terrainModificationApplied = exchange != 0.f;
            history.add(state);
            particle.pos = endPosition;
        }
        if (wrapPositions)
            particle.pos = Vector3::wrap(particle.pos, Vector3(0, 0, -100), Vector3(terrainSize.x(), terrainSize.y(), 1000));
        history.add(particle);
        if (!Vector3::isInBox(particle.pos.xy(), Vector3(), terrainSize.xy()) || particle.pos.z() < -20.f ||
            (!contacts.empty() && particle.velocity.norm2() < 1e-8f)) break;
    }
    // Unload remaining sediment only onto upward-facing ground below an in-domain parcel.
    // A parcel leaving the domain exports its remaining load instead.
    if (depositFactor > 0.f && particle.capacity > 0.f &&
        Vector3::isInBox(particle.pos.xy(), Vector3(), terrainSize.xy())) {
        auto [hit, triangle] = boundariesTree->getIntersectionAndTriangleIndex(
            particle.pos, Vector3(particle.pos.x(), particle.pos.y(), -1.f));
        if (hit.isValid() && boundariesTree->triangles[triangle].normal.z() > 0.f) {
            particle.pos = hit;
            const float amount = particle.capacity;
            particle.capacity = 0.f;
            ParticleState state(particle);
            state.exchange = {hit, -amount};
            state.terrainModificationApplied = true;
            history.add(state);
        }
    }
    return history;
}

std::vector<ErosionPoint> ParticleErosion::computeErosionValuesFromHistory(const ParticleHistory& history,
                                                                        SpacePartitioning& /*boundaries*/)
{
    std::vector<ErosionPoint> points;
    for (const auto& state : history.history)
        if (state.terrainModificationApplied) points.push_back(state.exchange);
    return points;
}

void ParticleErosion::modifyHeightmap(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    GridF modifications(heightmap->heights.getDimensions());
    std::vector<GridF> submodifications(allErosions.size(), modifications);

    #pragma omp parallel for
    for (size_t i = 0; i < allErosions.size(); i++) {
        const auto& erosionValuesAndPositions = allErosions[i];
        Vector3 previousPos(false);
        for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
            Vector3 pos = erosionValuesAndPositions[iRock].position;
            float val = erosionValuesAndPositions[iRock].erosionValue;
            int size = particleSize;
            val *= 100.f;
            RockErosion(size, val).computeErosionMatrix2D(submodifications[i], pos);
        }
    }

    for (const auto& sub : submodifications)
        modifications += sub;

    if (densityMap.size() > 0) {
        for (size_t i = 0; i < modifications.size(); i++) {
            modifications[i] = (modifications[i] > 0 ? modifications[i] : modifications[i] * ((1.f - materialImpact) + (1.f - densityMap[i]) * materialImpact));
        }
    }

    heightmap->heights += modifications * 0.5f;
}

void ParticleErosion::modifyVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    GridF modifications(voxelGrid->getDimensions());
    std::vector<GridF> submodifications(allErosions.size(), modifications);

    #pragma omp parallel for
    for (size_t i = 0; i < allErosions.size(); i++) {
        const auto& erosionValuesAndPositions = allErosions[i];
        Vector3 previousPos(false);
        for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
            Vector3 pos = erosionValuesAndPositions[iRock].position;
            float val = erosionValuesAndPositions[iRock].erosionValue;
            int size = particleSize;
            val *= 100.f;
            RockErosion(size, val).computeErosionMatrix(submodifications[i], pos - Vector3(.5f, .5f, .5f));
        }
    }
    for (const auto& sub : submodifications)
        modifications += sub;

    if (densityMap.size() > 0) {
        for (size_t i = 0; i < modifications.size(); i++) {
            modifications[i] = (modifications[i] > 0 ? modifications[i] : modifications[i] * ((1.f - materialImpact) + (1.f - densityMap[i]) * materialImpact));
        }
    }

    for (int x = 0; x < modifications.sizeX; x++)
        for (int y = 0; y < modifications.sizeY; y++)
            for (int z = 0; z < 2; z++)
                modifications.at(x, y, z) = 0;
    voxelGrid->applyModification(modifications * 0.5f);
    voxelGrid->limitVoxelValues(2.f, false);
    voxelGrid->saveState();
}

void ParticleErosion::modifyImplicitTerrain(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    auto* root = dynamic_cast<ImplicitNaryOperator*>(implicitTerrain);
    if (!root || typeid(*root) != typeid(ImplicitNaryOperator))
        throw std::invalid_argument("Implicit erosion requires an ordinary 3D root container");
    if (!std::isfinite(particleSize)) throw std::invalid_argument("Erosion particle size must be finite");
    const float dimensions = std::floor(particleSize);
    if (dimensions < 1.f || allErosions.empty()) return;
    std::vector<std::unique_ptr<ImplicitNaryOperator>> groups;
    groups.reserve(allErosions.size());
    for (size_t i = 0; i < allErosions.size(); ++i) groups.push_back(std::make_unique<ImplicitNaryOperator>());
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    #pragma omp parallel for
    for (size_t i = 0; i < allErosions.size(); ++i) {
        if (failed.load(std::memory_order_relaxed)) continue;
        try {
            const auto& points = allErosions[i];
            for (size_t j = 0; j < points.size(); ++j) {
                const auto& point = points[j];
                if (!point.position.isValid() || !std::isfinite(point.position.x()) || !std::isfinite(point.position.y()) ||
                    !std::isfinite(point.position.z()) || !std::isfinite(point.erosionValue))
                    throw std::invalid_argument("Implicit erosion points must be finite");
                const double value = double(point.erosionValue) * 100. * (j + 1 == points.size() ? .1 : 1.);
                const double amplitude = std::abs(value) * .1;
                if (amplitude > std::numeric_limits<float>::max()) throw std::invalid_argument("Implicit erosion amplitude overflow");
                if (amplitude == 0.) continue;
                std::unique_ptr<ImplicitPrimitive> sphere(ImplicitPatch::createPredefinedShape(
                    ImplicitPatch::Sphere, Vector3(dimensions, dimensions, dimensions), float(amplitude)));
                sphere->material = value > 0. ? AIR : DIRT;
                sphere->setSupportDimensions(sphere->dimensions);
                sphere->position = point.position - sphere->dimensions * .5f;
                groups[i]->addChild(sphere.get());
                sphere.release();
            }
        } catch (...) {
            if (!failed.exchange(true)) error = std::current_exception();
        }
    }
    if (error) std::rethrow_exception(error);
    auto total = std::make_unique<ImplicitNaryOperator>();
    for (auto& group : groups) {
        if (!group->composables.empty()) {
            total->addChild(group.get());
            group.release();
        }
    }
    if (!total->composables.empty()) {
        root->addChild(total.get());
        total.release();
    }
}

void ParticleErosion::modifyLayerBased(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    for (size_t i = 0; i < allErosions.size(); i++) {
        const auto& erosionValuesAndPositions = allErosions[i];
        Vector3 previousPos(false);
        bool isTheOnlyDepositRock = true;
        for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
            bool lastRock = iRock == erosionValuesAndPositions.size() - 1;
            Vector3 pos = erosionValuesAndPositions[iRock].position;
            float val = erosionValuesAndPositions[iRock].erosionValue;
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
                        float startZ = std::max(0.f, pos.z() - halfHeight);
                        float endZ = pos.z() + halfHeight;
                        layerBasedGrid->transformLayer(pos.x() + x, pos.y() + y, startZ, endZ, (val > 0 ? TerrainTypes::AIR : TerrainTypes::SAND));
                    }
                }
            }
        }
    }
}














std::vector<ParticleHistory> ParticleErosion::process()
{
    if (!terrain || !boundariesTree || quantity < 0 || posAndDirs.size() < size_t(quantity))
        throw std::invalid_argument("Erosion requires a terrain, boundaries and one emission per particle");
    for (float value : {dt, size, matterDensity, flowResponseTime, sedimentRadius, sedimentDensity,
                        airViscosity, waterViscosity})
        if (!std::isfinite(value) || value <= 0.f)
            throw std::invalid_argument("Erosion time, sizes, densities and viscosities must be finite and positive");
    for (float value : {strength, maxCapacityFactor, erosion, deposit, gravity, shearingStressConstantK,
                        shearingRatePower, erosionPowerValue, criticalShearValue, initialCapacity})
        if (!std::isfinite(value) || value < 0.f)
            throw std::invalid_argument("Erosion strengths and shear parameters must be finite and nonnegative");
    terrainModifTime = 0.f;
    VoxelGrid* asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    Heightmap* asHeightmap = dynamic_cast<Heightmap*>(terrain);
    ImplicitPatch* asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    LayerBasedGrid* asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    terrainSize = voxelGrid->getDimensions(); //terrain->getDimensions();
    float starting_distance = pow(terrainSize.maxComp()/2.f, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap

    if (densityMap.size() > 0) {
        densityMap.raiseErrorOnBadCoord = false;
        densityMap.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
    }
    // Hydraulic erosion parameters
    float strengthValue = strength;
    float maxCapacity = 1.f * strengthValue * maxCapacityFactor;
    erosionFactor = .01f * strengthValue * erosion;
    depositFactor = .01f * strengthValue * deposit;
    float capacity = maxCapacity * initialCapacity;

    particleSize = size;

//    GridF modifications((asHeightmap ? Vector3(terrainSize.x, terrainSize.y, 1) : terrainSize));

    gravityDefault = Vector3(0, 0, -gravity);
//    auto environmentalDensities = [&](const Vector3& pos) {
//        if (pos.z() < 0) return 1000.f;
//        return 1.f;
//    };
//    auto flowfieldValues = [&](const Vector3& pos) {
//        return Vector3();
//    };
//    GridV3 gravityfieldValues = GridV3(terrainSize, Vector3(0, 0, -gravity));
//    gravityfieldValues.raiseErrorOnBadCoord = false;
//    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);

    environmentalDensities = voxelGrid->getEnvironmentalDensities();
    environmentalDensities.raiseErrorOnBadCoord = false;
    environmentalDensities.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    flowfieldValues = GridV3(terrainSize);
    if (flowType == FLOWFIELD_TYPE::BASIC) {
        Vector3 airDir = Vector3(0, airForce, 0).rotate(0, 0, (airFlowfieldRotation / 180) * PI);
        Vector3 waterDir = Vector3(0, waterForce, 0).rotate(0, 0, (waterFlowfieldRotation / 180) * PI);
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            Vector3 pos = flowfieldValues.getCoordAsVector3(i);
//            if (pos.x() > flowfieldValues.sizeX / 2) continue;
            if (environmentalDensities(pos) < 100) { // In the air
                flowfieldValues(pos) = airDir;
            } else { // In water
                flowfieldValues(pos) = waterDir;
            }
        }
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
        /*
        airFlow = airFlow.resize(terrainSize.x, terrainSize.y, 1) * airForce * 10.f;
        waterFlow = waterFlow.resize(terrainSize.x, terrainSize.y, 1) * waterForce * 10.f;
        for (size_t i = 0; i < flowfieldValues.size(); i++) {
            if (environmentalDensities.at(flowfieldValues.getCoordAsVector3(i)) < 100) { // In the air
                flowfieldValues.at(i) = airFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            } else { // In water
                flowfieldValues.at(i) = waterFlow.at(flowfieldValues.getCoordAsVector3(i).xy());
            }
        }*/
    } else if (flowType == FLOWFIELD_TYPE::FLUID_SIMULATION) {
        flowfieldValues = GridV3(voxelGrid->getFlowfield(fluidSimType).resize(terrainSize) * airForce).meanSmooth();
        Vector3 aspectRatio = terrainSize.normalized();
        for (auto& v : flowfieldValues)
            v *= aspectRatio;
    } else if (flowType == FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS) {
        flowfieldValues = this->scene->flowfield.resize(terrainSize).meanSmooth();
        for (auto& v : flowfieldValues)
            v = v.xy();
    }
    flowfieldValues.raiseErrorOnBadCoord = false;
    flowfieldValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

//    ImageViewer::get().addImage(flowfieldValues.sliceXY(flowfieldValues.sizeZ / 2));
//    ImageViewer::get().show();

//    std::vector<BSpline> tunnels(quantity);
//    std::vector<GridF> submodifications(quantity, GridF(modifications.getDimensions()));
//    std::vector<int> nbPos(quantity), nbErosions(quantity);
//    std::vector<std::vector<std::pair<float, Vector3>>> allErosions(quantity);
//    std::vector<std::vector<std::pair<float, Vector3>>> allErosionsIfErosionWasOn(quantity);

    // Cache computations of the erosion matrix
    RockErosion::createPrecomputedAttackMask(particleSize);
    RockErosion::createPrecomputedAttackMask(particleSize * 2);
    RockErosion::createPrecomputedAttackMask2D(particleSize);
    RockErosion::createPrecomputedAttackMask2D(particleSize * 2);

    std::vector<ParticleHistory> tracks(quantity);
    std::vector<std::vector<ErosionPoint>> allErosions(quantity);
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    particleSimulationTime = timeIt([&]() {
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < quantity; ++i) {
            if (failed.load(std::memory_order_relaxed)) continue;
            try {
                auto [position, velocity] = posAndDirs[i];
                ErosionParticle particle;
                initializeParticle(particle, position, velocity, particleSize, matterDensity, capacity, maxCapacity);
                tracks[i] = trackParticlePositions(particle);
                allErosions[i] = computeErosionValuesFromHistory(tracks[i], *boundariesTree);
            } catch (...) {
                if (!failed.exchange(true)) error = std::current_exception();
            }
        }
    });
    if (error) std::rethrow_exception(error);
    if (!applyTheErosion || strengthValue == 0.f) return tracks;
    terrainModifTime = timeIt([&]() {
        if (asImplicit) modifyImplicitTerrain(allErosions);
        else if (asVoxels) modifyVoxelGrid(allErosions);
        else if (asHeightmap) modifyHeightmap(allErosions);
        else if (asLayers) modifyLayerBased(allErosions);
    });
    return tracks;
}

std::vector<Vector3> ParticleHistory::getPositions()
{
    std::vector<Vector3> positions;
    positions.reserve(history.size());
    for (const auto& state : history) positions.push_back(state.particleState.pos);
    return positions;
}

std::vector<Vector3> ParticleHistory::getErosionPositions()
{
    std::vector<Vector3> positions;
    for (const auto& state : history)
        if (state.terrainModificationApplied) positions.push_back(state.exchange.position);
    return positions;
}

std::vector<float> ParticleHistory::getErosionValues()
{
    std::vector<float> values;
    for (const auto& state : history)
        if (state.terrainModificationApplied) values.push_back(state.exchange.erosionValue);
    return values;
}

std::vector<ErosionPoint> ParticleHistory::getErosionPoints()
{
    std::vector<ErosionPoint> points;
    for (const auto& state : history)
        if (state.terrainModificationApplied) points.push_back(state.exchange);
    return points;
}

std::vector<float> ParticleHistory::getCapacities()
{
    std::vector<float> capacities;
    capacities.reserve(history.size());
    for (const auto& state : history) capacities.push_back(state.particleState.capacity);
    return capacities;
}
