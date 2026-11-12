#include "ParticleErosion.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>

#include "DataStructure/BVH.h"
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
    for (float value : {radius,density,maxCapacity,initialLoad})
        if (!std::isfinite(value) || value<0.f)
            throw std::invalid_argument("Invalid erosion particle properties");
    if (radius<=0.f || density<=0.f || !position.isValid() || !velocity.isValid())
        throw std::invalid_argument("Erosion particles require a positive radius/density and valid launch state");
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
    if (!std::isfinite(particle.properties->mass) || !std::isfinite(particle.properties->volume))
        throw std::invalid_argument("Erosion particle mass or volume overflow");
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
        const float endDensity=environmentalDensities.interpolate(particle.pos);
        const Vector3 endAcceleration=gravityDefault*(1.f-endDensity/particle.properties->density);
        particle.velocity+=(endAcceleration-particle.force)*(.5f*dt);
        const Vector3 endTarget=flowfieldValues.interpolate(particle.pos)+computeSettlingVelocity(particle)*load;
        particle.velocity=endTarget+(particle.velocity-endTarget)*relaxation;
        particle.force=endAcceleration;
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

namespace {
struct BrushEvent {
    ErosionPoint point;
    int x0, x1, y0, y1, z0, z1;
    double normalization = 0.;
};

bool finitePosition(const Vector3& p)
{
    return p.isValid() && std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
}

// Equations 9/13: normalized hemisphere height or linear radial voxel kernel.
// Each XY tile owns its output columns, so overlapping brushes need no locks or
// full-terrain buffers per particle/thread. Event order is stable across threads.
GridF rasterizeErosion(const std::vector<std::vector<ErosionPoint>>& tracks,
                       const Vector3& dimensions, float radius, bool heightField)
{
    if (!std::isfinite(radius) || radius <= 0.f)
        throw std::invalid_argument("Erosion brush radius must be finite and positive");
    GridF result(dimensions);
    if (result.empty()) return result;
    constexpr int tileWidth = 8;
    const int nx = int(result.sizeX), ny = int(result.sizeY), nz = int(result.sizeZ);
    const int tilesX = (nx + tileWidth - 1) / tileWidth;
    const int tilesY = (ny + tileWidth - 1) / tileWidth;
    std::vector<BrushEvent> events;
    for (const auto& track : tracks) for (const auto& point : track) {
        if (!finitePosition(point.position) || !std::isfinite(point.erosionValue))
            throw std::invalid_argument("Erosion points must be finite");
        if (point.erosionValue == 0.f) continue;
        BrushEvent event;
        event.point = point;
        // Clamp in floating point before integer conversion (including remote events).
        auto lower = [radius](float c, int n) { return int(std::clamp(std::ceil(double(c)-radius), 0., double(n))); };
        auto upper = [radius](float c, int n) { return int(std::clamp(std::floor(double(c)+radius)+1., 0., double(n))); };
        event.x0 = lower(point.position.x(),nx); event.x1 = upper(point.position.x(),nx);
        event.y0 = lower(point.position.y(),ny); event.y1 = upper(point.position.y(),ny);
        event.z0 = heightField ? 0 : lower(point.position.z(),nz);
        event.z1 = heightField ? 1 : upper(point.position.z(),nz);
        if (event.x0 >= event.x1 || event.y0 >= event.y1 || event.z0 >= event.z1) continue;
        events.push_back(event);
    }
    auto weight = [radius,heightField](const BrushEvent& event, int x, int y, int z) {
        const double dx = double(x)-event.point.position.x(), dy = double(y)-event.point.position.y();
        const double dz = heightField ? 0. : double(z)-event.point.position.z();
        const double squared = dx*dx+dy*dy+dz*dz;
        if (squared >= double(radius)*radius) return 0.;
        return heightField ? std::sqrt(double(radius)*radius-squared) : 1.-std::sqrt(squared)/radius;
    };
    #pragma omp parallel for schedule(dynamic, 8)
    for (size_t i = 0; i < events.size(); ++i) {
        auto& e = events[i];
        for (int z=e.z0; z<e.z1; ++z) for (int y=e.y0; y<e.y1; ++y) for (int x=e.x0; x<e.x1; ++x)
            e.normalization += weight(e,x,y,z);
    }
    std::vector<std::vector<size_t>> tiles(size_t(tilesX)*tilesY);
    for (size_t i=0; i<events.size(); ++i) {
        const auto& e=events[i];
        if (e.normalization <= 0.) continue;
        for (int ty=e.y0/tileWidth; ty<=(e.y1-1)/tileWidth; ++ty)
            for (int tx=e.x0/tileWidth; tx<=(e.x1-1)/tileWidth; ++tx)
                tiles[size_t(ty)*tilesX+tx].push_back(i);
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t tile=0; tile<tiles.size(); ++tile) {
        const int tx=int(tile%tilesX)*tileWidth, ty=int(tile/tilesX)*tileWidth;
        for (size_t index : tiles[tile]) {
            const auto& e=events[index];
            const double scale=-double(e.point.erosionValue)/e.normalization;
            for (int z=e.z0; z<e.z1; ++z)
                for (int y=std::max(ty,e.y0); y<std::min(ty+tileWidth,e.y1); ++y)
                    for (int x=std::max(tx,e.x0); x<std::min(tx+tileWidth,e.x1); ++x)
                        result(x,y,z) += float(scale*weight(e,x,y,z));
        }
    }
    return result;
}
}

void ParticleErosion::modifyHeightmap(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    heightmap->heights += rasterizeErosion(allErosions, heightmap->heights.getDimensions(), particleSize, true);
}

void ParticleErosion::modifyVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    GridF changes = rasterizeErosion(allErosions, voxelGrid->getDimensions(), particleSize, false);
    auto values = voxelGrid->getVoxelValues();
    // Carry upper density excess upward within each independent column (section 4.4).
    // Excess leaving the top is exported. Lower clamping is a density-grid limitation.
    #pragma omp parallel for collapse(2)
    for (int y=0; y<int(values.sizeY); ++y) for (int x=0; x<int(values.sizeX); ++x) {
        float excess=0.f;
        for (int z=0; z<int(values.sizeZ); ++z) {
            if (changes(x,y,z)==0.f && excess==0.f) continue;
            const float previous=values(x,y,z);
            const float value=previous+changes(x,y,z)+excess;
            // Existing converted terrains may contain signed fields beyond [-1,1].
            // Preserve their range rather than clamping unrelated preexisting matter.
            const float upper=std::max(1.f,previous), lower=std::min(-1.f,previous);
            excess=std::max(0.f,value-upper);
            changes(x,y,z)=std::clamp(value,lower,upper)-previous;
        }
    }
    voxelGrid->applyModification(std::move(changes));
    voxelGrid->saveState();
}

void ParticleErosion::modifyBinaryVoxelGrid(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    if (!std::isfinite(binaryHitScale) || !std::isfinite(binaryHitThreshold) ||
        binaryHitScale<=0.f || binaryHitThreshold<=0.f)
        throw std::invalid_argument("Binary erosion hit scale and threshold must be finite and positive");
    GridF changes=rasterizeErosion(allErosions,voxelGrid->getDimensions(),particleSize,false);
    auto values=voxelGrid->getVoxelValues();
    GridF& hits=binaryHitAccumulator ? *binaryHitAccumulator : binaryHits;
    if (hits.getDimensions()!=values.getDimensions()) hits=GridF(values.getDimensions());
    #pragma omp parallel for
    for (size_t i=0;i<values.size();++i) {
        const bool solid=values[i]>0.f;
        // Retain fractional hits so repeated small brushes eventually affect occupancy.
        // Do not charge already empty cells for further detachment (or solid cells for filling).
        const float accumulated=hits[i]-changes[i]*binaryHitScale;
        hits[i]=solid ? std::max(0.f,accumulated) : std::min(0.f,accumulated);
        float occupancy=solid ? 1.f : -1.f;
        if (solid && hits[i]>=binaryHitThreshold) { occupancy=-1.f; hits[i]=0.f; }
        else if (!solid && hits[i]<=-binaryHitThreshold) { occupancy=1.f; hits[i]=0.f; }
        changes[i]=occupancy-values[i];
    }
    voxelGrid->applyModification(std::move(changes));
    voxelGrid->saveState();
}

void ParticleErosion::modifyImplicitTerrain(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    auto* root = dynamic_cast<ImplicitNaryOperator*>(implicitTerrain);
    if (!root || typeid(*root) != typeid(ImplicitNaryOperator))
        throw std::invalid_argument("Implicit erosion requires an ordinary 3D root container");
    if (!std::isfinite(particleSize)) throw std::invalid_argument("Erosion particle size must be finite");
    const float dimensions = 2.f * particleSize;
    if (dimensions <= 0.f) throw std::invalid_argument("Erosion particle radius must be positive");
    if (allErosions.empty()) return;
    std::vector<std::unique_ptr<ImplicitNaryOperator>> groups;
    groups.reserve(allErosions.size());
    for (size_t i = 0; i < allErosions.size(); ++i) {
        groups.push_back(std::make_unique<ImplicitNaryOperator>());
        groups.back()->additive = true;
    }
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
                const double value = point.erosionValue;
                // Corrected compact linear kernel: integral over its support is |Q|.
                const double amplitude = 3. * std::abs(value) / (M_PI * std::pow(double(particleSize), 3));
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
    total->additive = true;
    for (auto& group : groups) {
        if (!group->composables.empty()) {
            total->addChild(group.get());
            group.release();
        }
    }
    if (!total->composables.empty()) {
        if (!root->additive && !root->composables.empty()) {
            // Allocate everything before publishing. Retain the original tree's
            // composition as a single base, then add signed erosion field offsets.
            auto base = std::make_unique<ImplicitNaryOperator>();
            std::vector<ImplicitPatch*> replacement{base.get(), total.get()};
            base->composables.swap(root->composables);
            for (auto* child : base->composables) if (child) child->parent = base.get();
            base->parent = root;
            total->parent = root;
            root->composables.swap(replacement);
            root->additive = true;
            root->firstChildIsBase = true;
            base.release(); total.release();
            for (ImplicitPatch* node=root; node; node=node->parent) node->_cached=false;
        } else {
            root->addChild(total.get());
            root->additive = true;
            total.release();
        }
    }
}

void ParticleErosion::modifyLayerBased(const std::vector<std::vector<ErosionPoint>>& allErosions)
{
    if (!std::isfinite(particleSize) || particleSize<=0.f)
        throw std::invalid_argument("Erosion brush radius must be finite and positive");
    auto& columns = layerBasedGrid->getLayers();
    // Prepare column-local operations using equation 10's hemisphere footprint.
    // Per-event rasterization is only 2D; no volumetric conversion is needed.
    struct Interval { float center, amount; };
    std::vector<std::vector<Interval>> operations(columns.size());
    for (const auto& track : allErosions) for (const auto& point : track) {
        if (!finitePosition(point.position) || !std::isfinite(point.erosionValue))
            throw std::invalid_argument("Erosion points must be finite");
        if (point.erosionValue == 0.f) continue;
        // Sparse footprint, normalized on the in-domain column samples.
        const int x0=int(std::clamp(std::ceil(double(point.position.x())-particleSize),0.,double(columns.sizeX)));
        const int x1=int(std::clamp(std::floor(double(point.position.x())+particleSize)+1.,0.,double(columns.sizeX)));
        const int y0=int(std::clamp(std::ceil(double(point.position.y())-particleSize),0.,double(columns.sizeY)));
        const int y1=int(std::clamp(std::floor(double(point.position.y())+particleSize)+1.,0.,double(columns.sizeY)));
        double normalization=0.;
        auto height = [&](int x,int y) {
            const double dx=x-double(point.position.x()),dy=y-double(point.position.y());
            return std::sqrt(std::max(0.,double(particleSize)*particleSize-dx*dx-dy*dy));
        };
        for (int y=y0;y<y1;++y) for (int x=x0;x<x1;++x) normalization+=height(x,y);
        if (normalization==0.) continue;
        for (int y=y0;y<y1;++y) for (int x=x0;x<x1;++x) {
            const float amount=float(point.erosionValue*height(x,y)/normalization);
            if (amount!=0.f) operations[columns.getIndex(x,y,0)].push_back({point.position.z(),amount});
        }
    }
    // Different workers own different material stacks. Preserve the order of
    // overlapping operations within a stack and keep all final deposit events.
    #pragma omp parallel for schedule(dynamic, 16)
    for (size_t i=0;i<columns.size();++i) {
        auto& column=columns[i];
        for (const auto& op : operations[i]) {
            const bool eroding=op.amount>0.f;
            float remaining=std::abs(op.amount);
            struct Span { float lo,hi; TerrainTypes material; };
            std::vector<Span> spans;
            float top=0.f;
            for (const auto& layer : column) {
                spans.push_back({top,top+layer.second,layer.first}); top+=layer.second;
            }
            // Explicit air above the stack permits deposition and cavity filling.
            spans.push_back({top,std::max(top,op.center)+remaining+particleSize,AIR});
            std::vector<size_t> candidates;
            for (size_t j=0;j<spans.size();++j) {
                const bool air=LayerBasedGrid::invisibleLayers.count(spans[j].material)!=0;
                if (air!=eroding) candidates.push_back(j);
            }
            auto distance = [&](size_t j) {
                return std::max({spans[j].lo-op.center,op.center-spans[j].hi,0.f});
            };
            std::stable_sort(candidates.begin(),candidates.end(),[&](size_t a,size_t b) { return distance(a)<distance(b); });
            std::vector<Span> replacements;
            for (size_t j : candidates) {
                const auto span=spans[j];
                const float amount=std::min(remaining,span.hi-span.lo);
                if (amount<=0.f) continue;
                const float start=std::clamp(op.center-amount*.5f,span.lo,span.hi-amount);
                replacements.push_back({start,start+amount,eroding?AIR:SAND});
                remaining-=amount;
                if (remaining<=1e-7f) break;
            }
            std::vector<std::pair<TerrainTypes,float>> updated;
            for (const auto& span : spans) {
                const auto replacement=std::find_if(replacements.begin(),replacements.end(),[&](const Span& r) {
                    return r.lo>=span.lo && r.hi<=span.hi;
                });
                if (replacement==replacements.end()) updated.push_back({span.material,span.hi-span.lo});
                else {
                    if (replacement->lo>span.lo) updated.push_back({span.material,replacement->lo-span.lo});
                    updated.push_back({replacement->material,replacement->hi-replacement->lo});
                    if (replacement->hi<span.hi) updated.push_back({span.material,span.hi-replacement->hi});
                }
            }
            column=std::move(updated);
        }
        const auto [x,y,z]=columns.getCoord(i);
        layerBasedGrid->cleanLayer(x,y);
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
    balance = ErosionBalance();
    terrainModifTime = 0.f;
    VoxelGrid* asVoxels = dynamic_cast<VoxelGrid*>(terrain);
    Heightmap* asHeightmap = dynamic_cast<Heightmap*>(terrain);
    ImplicitPatch* asImplicit = dynamic_cast<ImplicitPatch*>(terrain);
    LayerBasedGrid* asLayers = dynamic_cast<LayerBasedGrid*>(terrain);

    if (!asVoxels && !asHeightmap && !asImplicit && !asLayers)
        throw std::invalid_argument("Unsupported particle erosion terrain representation");
    if (applyOn==BINARY_VOXELS && !asVoxels)
        throw std::invalid_argument("Binary erosion requires a voxel terrain");
    // Bind the adapter to the actual terrain instead of an optional companion.
    if (asVoxels) voxelGrid=asVoxels;
    if (asHeightmap) heightmap=asHeightmap;
    if (asImplicit) implicitTerrain=asImplicit;
    if (asLayers) layerBasedGrid=asLayers;
    terrainSize=terrain->getDimensions();
    if (!finitePosition(terrainSize) || terrainSize.minComp()<=0.f)
        throw std::invalid_argument("Erosion terrain dimensions must be positive and finite");
    for (const auto& emission : posAndDirs)
        if (!finitePosition(emission.first) || !finitePosition(emission.second))
            throw std::invalid_argument("Erosion emissions must have finite positions and velocities");
    for (float value : {bounciness,bouncingCoefficient,airForce,waterForce,airFlowfieldRotation,
                        waterFlowfieldRotation,materialImpact})
        if (!std::isfinite(value)) throw std::invalid_argument("Erosion flow and collision parameters must be finite");

    if (!densityMap.empty()) {
        densityMap.raiseErrorOnBadCoord=false;
        densityMap.returned_value_on_outside=RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    }
    const float strengthValue=strength;
    const float maxCapacity=strength*maxCapacityFactor;
    erosionFactor=.01f*strength*erosion;
    depositFactor=.01f*strength*deposit;
    if (!std::isfinite(maxCapacity) || !std::isfinite(erosionFactor) || !std::isfinite(depositFactor))
        throw std::invalid_argument("Erosion strength products overflow");
    const float capacity=maxCapacity*std::clamp(initialCapacity,0.f,1.f);
    particleSize=size;
    gravityDefault=Vector3(0,0,-gravity);
    environmentalDensities=terrain->getEnvironmentalDensities();
    environmentalDensities.raiseErrorOnBadCoord=false;
    environmentalDensities.returned_value_on_outside=RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    for (float value : environmentalDensities)
        if (!std::isfinite(value) || value<0.f) throw std::invalid_argument("Invalid medium density");
    // The medium can extend above a heightmap or beyond the terrain's solid extent.
    const Vector3 fieldSize=environmentalDensities.getDimensions();
    flowfieldValues=GridV3(fieldSize);
    if (flowType==BASIC) {
        const Vector3 airDir=Vector3(0,airForce,0).rotate(0,0,airFlowfieldRotation*PI/180.f);
        const Vector3 waterDir=Vector3(0,waterForce,0).rotate(0,0,waterFlowfieldRotation*PI/180.f);
        #pragma omp parallel for
        for (size_t i=0;i<flowfieldValues.size();++i)
            flowfieldValues[i]=environmentalDensities[i]<100.f ? airDir : waterDir;
    } else if (flowType==FLOWFIELD_IMAGE || flowType==FLOWFIELD_ENVOBJECTS) {
        if (flowType==FLOWFIELD_ENVOBJECTS && scene) {
            flowfieldValues=scene->flowfield.resize(fieldSize);
            #pragma omp parallel for
            for (size_t i=0;i<flowfieldValues.size();++i)
                flowfieldValues[i]*=environmentalDensities[i]<100.f ? airForce : waterForce;
        } else {
            // The interface passes environmental-object flow through these fields too.
            // Preserve 3D fields; image fields extend their sole Z slice vertically.
            GridV3 air=airFlow.empty() ? GridV3() : airFlow.resize(fieldSize.x(),fieldSize.y(),airFlow.sizeZ>1 ? fieldSize.z() : 1);
            GridV3 water=waterFlow.empty() ? GridV3() : waterFlow.resize(fieldSize.x(),fieldSize.y(),waterFlow.sizeZ>1 ? fieldSize.z() : 1);
            #pragma omp parallel for
            for (size_t i=0;i<flowfieldValues.size();++i) {
                const Vector3 p=flowfieldValues.getCoordAsVector3(i);
                const bool inAir=environmentalDensities[i]<100.f;
                const GridV3& source=inAir ? air : water;
                if (!source.empty()) flowfieldValues[i]=source.interpolate(p)*(inAir ? airForce : waterForce);
            }
        }
    } else if (flowType==FLUID_SIMULATION) {
        flowfieldValues=terrain->getFlowfield(fluidSimType).resize(fieldSize)*airForce;
    } else throw std::invalid_argument("Unknown erosion velocity field source");
    for (const auto& velocity : flowfieldValues)
        if (!finitePosition(velocity)) throw std::invalid_argument("Invalid medium velocity");
    flowfieldValues.raiseErrorOnBadCoord=false;
    flowfieldValues.returned_value_on_outside=RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

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
    double exchanges=0.;
    for (const auto& track : tracks) {
        if (track.history.empty()) continue;
        balance.initialSediment+=track.history.front().particleState.capacity;
        balance.finalSediment+=track.history.back().particleState.capacity;
        for (const auto& state : track.history)
            if (state.terrainModificationApplied) exchanges+=state.exchange.erosionValue;
    }
    balance.requestedTerrainChange=-exchanges;
    balance.sedimentResidual=balance.initialSediment+exchanges-balance.finalSediment;
    if (!applyTheErosion || strengthValue == 0.f) return tracks;
    // Measure in the native representation's units; never label a field sum as volume.
    auto measure = [&]() {
        double sum=0.;
        if (asHeightmap) {
            balance.measure=ErosionBalance::HEIGHT_VOLUME;
            for (float value : heightmap->heights) sum+=value;
        } else if (asLayers) {
            balance.measure=ErosionBalance::LAYER_VOLUME;
            for (const auto& column : layerBasedGrid->getLayers())
                for (const auto& layer : column)
                    if (!LayerBasedGrid::invisibleLayers.count(layer.first)) sum+=layer.second;
        } else if (asVoxels) {
            balance.measure=applyOn==BINARY_VOXELS ? ErosionBalance::OCCUPIED_VOXELS : ErosionBalance::DENSITY_SUM;
            for (float value : voxelGrid->getVoxelValues()) sum+=applyOn==BINARY_VOXELS ? double(value>0.f) : value;
        }
        return sum;
    };
    const double before=measure();
    terrainModifTime = timeIt([&]() {
        if (asImplicit) modifyImplicitTerrain(allErosions);
        else if (asVoxels && applyOn==BINARY_VOXELS) modifyBinaryVoxelGrid(allErosions);
        else if (asVoxels) modifyVoxelGrid(allErosions);
        else if (asHeightmap) modifyHeightmap(allErosions);
        else if (asLayers) modifyLayerBased(allErosions);
    });
    balance.measuredTerrainChange=measure()-before;
    balance.terrainMeasured=!asImplicit;
    // Occupancy is quantized; compare to continuous sediment only for volume/density adapters.
    if (balance.terrainMeasured && balance.measure!=ErosionBalance::OCCUPIED_VOXELS)
        balance.unappliedTerrainChange=balance.requestedTerrainChange-balance.measuredTerrainChange;
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
