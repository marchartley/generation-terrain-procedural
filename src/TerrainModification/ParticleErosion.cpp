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
    float restitution, float damping, int maxBounces, const Vector3& mediumVelocity, float responseTime)
{
    std::vector<ParticleContact> contacts;
    float remaining=dt;
    restitution=std::clamp(restitution,0.f,1.f);damping=std::clamp(damping,0.f,1.f);
    while (remaining>0.f) {
        // Exact solution for acceleration and a linear drag target frozen over
        // this segment. This avoids timestep-dependent terminal velocities.
        auto motion=[&](float time) {
            if (responseTime<=0.f) return std::pair<Vector3,Vector3>{velocity*time+force*(.5f*time*time),velocity+force*time};
            const Vector3 terminal=mediumVelocity+force*responseTime;
            const float decay=std::exp(-time/responseTime);
            const float integrated=-responseTime*std::expm1(-time/responseTime);
            return std::pair<Vector3,Vector3>{terminal*time+(velocity-terminal)*integrated,
                                            terminal+(velocity-terminal)*decay};
        };
        auto [displacement,endVelocity]=motion(remaining);
        const Vector3 destination=pos+displacement;
        auto [hit,triangle]=boundaries->getIntersectionAndTriangleIndex(pos,destination);
        if (!hit.isValid()) {pos=destination;velocity=endVelocity;break;}
        Vector3 normal=boundaries->triangles[triangle].normal.normalized();
        if (normal.dot(displacement)>0.f) normal=-normal;
        float lo=0.f,hi=remaining;
        // Refine the hit time against the contacted plane along the curved path.
        for (int i=0;i<20;++i) {
            const float middle=(lo+hi)*.5f;
            if (normal.dot(pos+motion(middle).first-hit)>0.f) lo=middle;else hi=middle;
        }
        const float hitTime=(lo+hi)*.5f;
        auto [toHit,incoming]=motion(hitTime);
        hit=pos+toHit;
        contacts.push_back({hit,normal,incoming});
        remaining-=hitTime;
        velocity=(incoming-normal*((1.f+restitution)*incoming.dot(normal)))*damping;
        pos=hit+normal*1e-4f;
        if (int(contacts.size())>=maxBounces) break;
        // Persistent support balances inward acceleration/drag. Tangential motion
        // continues; do not tunnel through the same triangle during the remainder.
        if (std::abs(velocity.dot(normal))<1e-3f && force.dot(normal)<=0.f && mediumVelocity.dot(normal)<=0.f) {
            Vector3 acceleration=force-normal*force.dot(normal);
            Vector3 tangentTarget=mediumVelocity-normal*mediumVelocity.dot(normal);
            velocity-=normal*velocity.dot(normal);
            if (responseTime>0.f) {
                const Vector3 terminal=tangentTarget+acceleration*responseTime;
                const float decay=std::exp(-remaining/responseTime);
                pos+=terminal*remaining+(velocity-terminal)*(-responseTime*std::expm1(-remaining/responseTime));
                velocity=terminal+(velocity-terminal)*decay;
            } else {pos+=velocity*remaining+acceleration*(.5f*remaining*remaining);velocity+=acceleration*remaining;}
            break;
        }
    }
    dir=velocity;
    return contacts;
}

namespace {
float clearSettlingSpeed(float radius,float grainDensity,float fluidDensity,float viscosity,float gravity)
{
    const double stokes=(2./9.)*radius*radius*std::abs(double(grainDensity)-fluidDensity)*gravity/viscosity;
    const double reynolds=2.*radius*fluidDensity*stokes/viscosity;
    if (reynolds<=.1) return float(stokes);
    // Schiller–Naumann below Re=1000, constant drag coefficient above it.
    double lo=0.,hi=stokes;
    for (int i=0;i<32;++i) {
        const double speed=(lo+hi)*.5, re=2.*radius*fluidDensity*speed/viscosity;
        const double correction=re<1000.?1.+.15*std::pow(re,.687):.44*re/24.;
        if (speed*correction<stokes) lo=speed;else hi=speed;
    }
    return float((lo+hi)*.5);
}
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
    const bool cached=settlingCacheValid && settlingCacheKey==std::array<float,5>{sedimentRadius,sedimentDensity,gravity,airViscosity,waterViscosity};
    const float speed=cached && density==1.f?settlingAir:
        (cached && density==1000.f?settlingWater:
         clearSettlingSpeed(sedimentRadius,sedimentDensity,density,viscosity,gravity));
    const float direction=sedimentDensity>=density?1.f:-1.f;
    return Vector3(0,0,-direction*speed*hindered);
}

float ParticleErosion::computeErosionValue(const ErosionParticle& particle, float vRel, float contactTime)
{
    const float length = std::max(particle.properties->radius * .002f, 1e-6f);
    const float shear = shearingStressConstantK * std::pow(std::max(0.f, vRel) / length, shearingRatePower);
    const float excess = std::max(0.f, shear - criticalShearValue);
    if (excess == 0.f) return 0.f;
    const float hemisphere = .5f * particle.properties->volume;
    const float detached = erosionFactor * std::pow(excess, erosionPowerValue) * hemisphere * (contactTime<0.f?dt:contactTime);
    return std::clamp(detached, 0.f, std::max(0.f, particle.properties->maxCapacity - particle.capacity));
}

float ParticleErosion::computeDepositionValue(const ErosionParticle& particle, float contactTime)
{
    if (particle.capacity <= 0.f || particle.properties->volume <= 0.f) return 0.f;
    const float area = M_PI * particle.properties->radius * particle.properties->radius;
    // Convert settling flux to a volume exchanged over this contact step.
    const float amount = depositFactor * computeSettlingVelocity(particle).norm() * area * (contactTime<0.f?dt:contactTime) *
                         (particle.capacity / particle.properties->volume);
    return std::clamp(amount, 0.f, particle.capacity);
}

ParticleHistory ParticleErosion::trackParticlePositions(ErosionParticle& particle)
{
    ParticleHistory history;
    if (terrain->checkIsInGround(particle.pos)) return history;
    history.history.reserve(1024);
    history.add(particle);
    history.termination=ParticleHistory::TIME_LIMIT;
    int collisionsLeft=maxCollisions<0?std::numeric_limits<int>::max():maxCollisions;
    double elapsed=0.;
    int step=0;
    for (; elapsed<maxSimulationTime && collisionsLeft>0 && step<maxTrackingSteps; ++step) {
        const float density=environmentalDensities.interpolate(particle.pos);
        particle.force=gravityDefault*(1.f-density/particle.properties->density);
        const float load=particle.properties->maxCapacity>0.f?particle.capacity/particle.properties->maxCapacity:0.f;
        const Vector3 target=flowfieldValues.interpolate(particle.pos)+computeSettlingVelocity(particle)*load;
        const float speed=std::max(particle.velocity.norm(),target.norm());
        const float stepDt=float(std::min({double(dt),double(maxSimulationTime)-elapsed,.25,
            double(particle.properties->radius)*.25/std::max(double(speed),1e-6),
            std::sqrt(double(particle.properties->radius)*.5/std::max(double(particle.force.norm()),1e-6))}));
        auto contacts=particle.bounce(stepDt,boundariesTree,bounciness,bouncingCoefficient,
                                    std::min(8,collisionsLeft),target,flowResponseTime);
        elapsed+=stepDt;
        const float contactTime=stepDt/std::max(size_t(1),contacts.size());
        for (const auto& contact : contacts) {
            --collisionsLeft;
            const Vector3 endPosition = particle.pos;
            particle.pos = contact.position;
            const Vector3 tangent = contact.velocity - contact.normal * contact.velocity.dot(contact.normal);
            float detached = computeErosionValue(particle, tangent.norm(),contactTime);
            if (!densityMap.empty()) {
                const float resistance = std::clamp(densityMap.interpolate(contact.position), 0.f, 1.f);
                detached *= 1.f - std::clamp(materialImpact, 0.f, 1.f) * resistance;
            }
            const Vector3 settling = computeSettlingVelocity(particle);
            const float orientation = settling.norm2() > 0.f
                ? std::max(0.f, -settling.normalized().dot(contact.normal)) : 0.f;
            const float deposited = computeDepositionValue(particle,contactTime) * orientation;
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
        if (!Vector3::isInBox(particle.pos.xy(),Vector3(),terrainSize.xy()) || particle.pos.z()<-20.f) {
            history.termination=ParticleHistory::LEFT_DOMAIN;break;
        }
        const bool supported=!contacts.empty() && contacts.back().normal.z()>.5f &&
            particle.force.dot(contacts.back().normal)<=0.f && target.dot(contacts.back().normal)<=0.f;
        if (supported && particle.velocity.norm()<=restSpeed && target.xy().norm()<=restSpeed) {
            // Mechanical rest deposits at the actual support, never at a remote
            // point below a trajectory truncated by a time/collision limit.
            if (depositFactor>0.f && particle.capacity>0.f) {
                ParticleState settled(particle);
                settled.exchange={contacts.back().position,-particle.capacity};
                settled.terrainModificationApplied=true;
                particle.capacity=settled.particleState.capacity=0.f;
                history.add(settled);
            }
            history.termination=ParticleHistory::RESTING;break;
        }
        if (contacts.empty() && particle.velocity.norm()<1e-8f && particle.force.norm()<1e-8f && target.norm()<1e-8f) {
            history.termination=ParticleHistory::STATIONARY;break;
        }
    }
    if (history.termination==ParticleHistory::TIME_LIMIT) {
        if (collisionsLeft<=0) history.termination=ParticleHistory::COLLISION_LIMIT;
        else if (step>=maxTrackingSteps) history.termination=ParticleHistory::STEP_LIMIT;
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
bool finitePosition(const Vector3& p)
{
    return p.isValid() && std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
}
}

namespace {
// Alter only eligible intervals near the contact. Return the volume actually
// exchanged; unavailable bedrock cannot become a particle's sediment.
double exchangeColumn(std::vector<std::pair<TerrainTypes,float>>& column,
                      float center, double request, float radius)
{
    struct Span { float lo,hi; TerrainTypes material; };
    std::vector<Span> spans;
    float top=0.f;
    for (const auto& layer : column) {
        spans.push_back({top,top+layer.second,layer.first}); top+=layer.second;
    }
    const float upper=std::max(top,center+radius);
    if (upper>top) spans.push_back({top,upper,AIR});
    const bool eroding=request>0.;
    std::vector<size_t> candidates;
    for (size_t i=0;i<spans.size();++i) {
        const bool air=LayerBasedGrid::invisibleLayers.count(spans[i].material)!=0;
        if (air!=eroding && spans[i].hi>std::max(0.f,center-radius) && spans[i].lo<center+radius)
            candidates.push_back(i);
    }
    auto distance=[&](size_t i) { return std::max({spans[i].lo-center,center-spans[i].hi,0.f}); };
    std::stable_sort(candidates.begin(),candidates.end(),[&](size_t a,size_t b){ return distance(a)<distance(b); });
    std::vector<Span> replacements;
    double remaining=std::abs(request), accepted=0.;
    for (size_t i : candidates) {
        const float lo=std::max({0.f,spans[i].lo,center-radius});
        const float hi=std::min(spans[i].hi,center+radius);
        const float amount=float(std::min(remaining,double(hi-lo)));
        if (amount<=0.f) continue;
        const float start=std::clamp(center-amount*.5f,lo,hi-amount);
        float end=start+amount;
        if (double(end-start)>remaining) end=std::nextafter(end,start);
        if (end<=start) continue;
        replacements.push_back({start,end,eroding?AIR:SAND});
        accepted+=double(end)-start;
        remaining-=double(end)-start;
        if (remaining<=0.) break;
    }
    if (accepted==0.) return 0.;
    std::vector<std::pair<TerrainTypes,float>> updated;
    auto append=[&](TerrainTypes material,float amount) {
        if (amount<=0.f) return;
        if (!updated.empty() && updated.back().first==material) updated.back().second+=amount;
        else updated.push_back({material,amount});
    };
    for (const auto& span : spans) {
        auto replacement=std::find_if(replacements.begin(),replacements.end(),[&](const Span& r){return r.lo>=span.lo && r.hi<=span.hi;});
        if (replacement==replacements.end()) append(span.material,span.hi-span.lo);
        else {
            append(span.material,replacement->lo-span.lo);
            append(replacement->material,replacement->hi-replacement->lo);
            append(span.material,span.hi-replacement->hi);
        }
    }
    while (!updated.empty() && LayerBasedGrid::invisibleLayers.count(updated.back().first)) updated.pop_back();
    column=std::move(updated);
    return eroding?accepted:-accepted;
}
}

void ParticleErosion::reconcileTerrainExchanges(std::vector<ParticleHistory>& tracks)
{
    const bool layers=layerBasedGrid && terrain==layerBasedGrid;
    const bool heights=heightmap && terrain==heightmap;
    const bool voxels=voxelGrid && terrain==voxelGrid;
    const bool binary=voxels && applyOn==BINARY_VOXELS;
    std::vector<std::vector<size_t>> indices(tracks.size());
    std::vector<size_t> next(tracks.size(),0);
    std::vector<double> loads(tracks.size(),0.), initial(tracks.size(),0.);
    for (size_t i=0;i<tracks.size();++i) {
        if (tracks[i].history.empty()) continue;
        initial[i]=loads[i]=tracks[i].history.front().particleState.capacity;
        for (size_t j=0;j<tracks[i].history.size();++j) if (tracks[i].history[j].terrainModificationApplied) {
            indices[i].push_back(j);
            balance.proposedTerrainChange-=tracks[i].history[j].exchange.erosionValue;
        }
    }
    GridF values, original;
    auto columns=layers ? layerBasedGrid->getLayers() : Matrix3<std::vector<std::pair<TerrainTypes,float>>>();
    GridF* hitStorage=binaryHitAccumulator?binaryHitAccumulator:&binaryHits;
    GridF workingHits;
    const double binaryUnit=double(binaryHitThreshold)/binaryHitScale;
    if (binary && (!std::isfinite(binaryUnit) || binaryUnit<=0. || binaryHitScale<=0.f || binaryHitThreshold<=0.f)) throw std::invalid_argument("Invalid binary hit units");
    if (heights) values=heightmap->heights;
    if (voxels) {
        original=voxelGrid->getVoxelValues(); values=original;
        if (binary) {
            workingHits=hitStorage->getDimensions()==values.getDimensions()?*hitStorage:GridF(values.getDimensions());
            for (size_t i=0;i<values.size();++i)
                values[i]=float(std::clamp((original[i]>0.f?binaryUnit:0.)-workingHits[i]/binaryHitScale,0.,binaryUnit));
        }
    }
    const Vector3 dims=layers?columns.getDimensions():values.getDimensions();
    constexpr int width=8;
    const int nx=int(dims.x()), ny=int(dims.y());
    const int txCount=(nx+width-1)/width, tyCount=(ny+width-1)/width;
    struct Event {
        size_t track,state;
        double request,normalization=0.;
        Vector3 position;
        int x0,x1,y0,y1,z0,z1;
        std::vector<double> accepted;
    };
    const bool surface=heights||layers;
    auto weight=[&](const Event& e,int x,int y,int z) {
        const double dx=x-double(e.position.x()),dy=y-double(e.position.y());
        const double dz=surface?0.:z-double(e.position.z());
        const double d2=dx*dx+dy*dy+dz*dz, r=particleSize;
        return d2>=r*r?0.:(surface?std::sqrt(r*r-d2):1.-std::sqrt(d2)/r);
    };
    while (true) {
        std::vector<Event> events;
        for (size_t i=0;i<tracks.size();++i) if (next[i]<indices[i].size()) {
            const size_t j=indices[i][next[i]];
            const auto& state=tracks[i].history[j];
            const double request=std::clamp(double(state.exchange.erosionValue),-loads[i],
                std::max(0.,double(state.particleState.properties->maxCapacity)-loads[i]));
            Event e;
            e.track=i;e.state=j;e.request=request;e.position=state.exchange.position;
            auto lower=[&](float c,int n){return int(std::clamp(std::ceil(double(c)-particleSize),0.,double(n)));};
            auto upper=[&](float c,int n){return int(std::clamp(std::floor(double(c)+particleSize)+1.,0.,double(n)));};
            e.x0=lower(e.position.x(),nx);e.x1=upper(e.position.x(),nx);
            e.y0=lower(e.position.y(),ny);e.y1=upper(e.position.y(),ny);
            e.z0=surface?0:lower(e.position.z(),int(dims.z()));
            e.z1=surface?1:upper(e.position.z(),int(dims.z()));
            events.push_back(std::move(e));
        }
        if (events.empty()) break;
        // Field offsets have no finite matter reservoir. Preserve their native
        // field units while still preventing deposits from an empty carrier.
        if (!heights && !layers && !voxels) {
            for (auto& e : events) e.accepted.push_back(e.request);
        } else {
            #pragma omp parallel for schedule(dynamic,8)
            for (size_t i=0;i<events.size();++i) {
                auto& e=events[i];
                for (int z=e.z0;z<e.z1;++z) for (int y=e.y0;y<e.y1;++y) for (int x=e.x0;x<e.x1;++x)
                    e.normalization+=weight(e,x,y,z);
            }
            struct Work { size_t event,slot; };
            std::vector<std::vector<Work>> tiles(size_t(txCount)*tyCount);
            for (size_t i=0;i<events.size();++i) {
                auto& e=events[i];
                if (e.normalization<=0. || e.request==0.) continue;
                for (int ty=e.y0/width;ty<=(e.y1-1)/width;++ty) for (int tx=e.x0/width;tx<=(e.x1-1)/width;++tx) {
                    tiles[size_t(ty)*txCount+tx].push_back({i,e.accepted.size()});
                    e.accepted.push_back(0.);
                }
            }
            std::atomic<bool> failed{false};std::exception_ptr error;
            #pragma omp parallel for schedule(dynamic,1)
            for (size_t tile=0;tile<tiles.size();++tile) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    const int tx=int(tile%txCount)*width,ty=int(tile/txCount)*width;
                    for (const auto& work : tiles[tile]) {
                        auto& e=events[work.event];double accepted=0.;
                        for (int z=e.z0;z<e.z1;++z) for (int y=std::max(ty,e.y0);y<std::min(ty+width,e.y1);++y)
                            for (int x=std::max(tx,e.x0);x<std::min(tx+width,e.x1);++x) {
                                const double request=e.request*weight(e,x,y,z)/e.normalization;
                                if (request==0.) continue;
                                if (layers) accepted+=exchangeColumn(columns(x,y),e.position.z(),request,particleSize);
                                else {
                                    double remaining=std::abs(request);
                                    // A saturated density sample can deposit above itself.
                                    const int end=(!surface && !binary && request<0.)?int(values.sizeZ):z+1;
                                    for (int zz=z;zz<end && remaining>0.;++zz) {
                                        float& cell=values(x,y,zz);
                                        const double lo=heights||binary?0.:std::min(-1.f,original(x,y,zz));
                                        const double hi=heights?std::numeric_limits<float>::max():
                                            (binary?binaryUnit:std::max(1.f,original(x,y,zz)));
                                        const double amount=std::min(remaining,request>0.?std::max(0.,double(cell)-lo):std::max(0.,hi-cell));
                                        const float previous=cell;
                                        cell=float(request>0.?double(cell)-amount:double(cell)+amount);
                                        if (std::abs(double(cell)-previous)>amount) cell=std::nextafter(cell,previous);
                                        const double actual=std::abs(double(cell)-previous);
                                        remaining-=actual;
                                        accepted+=request>0.?actual:-actual;
                                    }
                                }
                            }
                        e.accepted[work.slot]=accepted;
                    }
                } catch (...) { if (!failed.exchange(true)) error=std::current_exception(); }
            }
            if (error) std::rethrow_exception(error);
        }
        for (const auto& e : events) {
            double accepted=0.;for (double part : e.accepted) accepted+=part;
            auto& state=tracks[e.track].history[e.state];
            const double proposed=state.exchange.erosionValue;
            if (proposed>0.) balance.rejectedDetachment+=std::max(0.,proposed-accepted);
            else balance.rejectedDeposition+=std::max(0.,accepted-proposed);
            loads[e.track]+=accepted;
            state.exchange.erosionValue=float(accepted);
            state.terrainModificationApplied=accepted!=0.;
            ++next[e.track];
        }
    }
    for (size_t i=0;i<tracks.size();++i) {
        double load=initial[i];
        for (auto& state : tracks[i].history) {
            if (state.terrainModificationApplied) load+=state.exchange.erosionValue;
            state.particleState.capacity=float(std::clamp(load,0.,double(state.particleState.properties->maxCapacity)));
        }
    }
    if (heights) heightmap->heights=std::move(values);
    else if (layers) layerBasedGrid->getLayers()=std::move(columns);
    else if (voxels) {
        if (binary) {
            for (size_t i=0;i<values.size();++i) {
                const float matter=values[i];
                bool solid=original[i]>0.f;
                if (matter<=0.f) solid=false;
                else if (matter>=binaryUnit) solid=true;
                workingHits[i]=float(((solid?binaryUnit:0.)-matter)*binaryHitScale);
                values[i]=solid?1.f:-1.f;
            }
            *hitStorage=std::move(workingHits);
        }
        voxelGrid->applyModification(values-original);voxelGrid->saveState();
    } else {
        std::vector<std::vector<ErosionPoint>> points(tracks.size());
        for (size_t i=0;i<tracks.size();++i) points[i]=tracks[i].getErosionPoints();
        modifyImplicitTerrain(points);
    }
}

namespace {
void applyNativePoints(ParticleErosion& erosion, TerrainModel* target,
                       const std::vector<std::vector<ErosionPoint>>& points)
{
    if (!target || !std::isfinite(erosion.particleSize) || erosion.particleSize<=0.f)
        throw std::invalid_argument("Native erosion requires a terrain and positive finite radius");
    std::vector<ParticleHistory> tracks(points.size());
    for (size_t i=0;i<points.size();++i) {
        double prefix=0., minimum=0., maximum=0.;
        for (const auto& point : points[i]) {
            if (!finitePosition(point.position) || !std::isfinite(point.erosionValue))
                throw std::invalid_argument("Erosion points must be finite");
            prefix+=point.erosionValue;minimum=std::min(minimum,prefix);maximum=std::max(maximum,prefix);
        }
        const double capacity=maximum-minimum;
        if (!std::isfinite(capacity) || capacity>std::numeric_limits<float>::max())
            throw std::invalid_argument("Erosion point quantities overflow");
        ErosionParticle particle;
        particle.properties=std::make_shared<ParticleProperties>();
        particle.properties->maxCapacity=float(capacity);
        particle.capacity=float(-minimum);
        tracks[i].add(particle);
        for (const auto& point : points[i]) {
            ParticleState state(particle);state.exchange=point;
            state.terrainModificationApplied=point.erosionValue!=0.f;
            tracks[i].add(state);
        }
    }
    TerrainModel* previous=erosion.terrain;
    erosion.balance=ErosionBalance();
    erosion.terrain=target;
    try { erosion.reconcileTerrainExchanges(tracks); }
    catch (...) { erosion.terrain=previous;throw; }
    erosion.terrain=previous;
}
}

void ParticleErosion::modifyHeightmap(const std::vector<std::vector<ErosionPoint>>& points)
{
    applyNativePoints(*this,heightmap,points);
}

void ParticleErosion::modifyVoxelGrid(const std::vector<std::vector<ErosionPoint>>& points)
{
    const auto previous=applyOn;applyOn=DENSITY_VOXELS;
    try { applyNativePoints(*this,voxelGrid,points); }
    catch (...) {applyOn=previous;throw;}
    applyOn=previous;
}

void ParticleErosion::modifyBinaryVoxelGrid(const std::vector<std::vector<ErosionPoint>>& points)
{
    const auto previous=applyOn;applyOn=BINARY_VOXELS;
    try { applyNativePoints(*this,voxelGrid,points); }
    catch (...) {applyOn=previous;throw;}
    applyOn=previous;
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

void ParticleErosion::modifyLayerBased(const std::vector<std::vector<ErosionPoint>>& points)
{
    applyNativePoints(*this,layerBasedGrid,points);
}

std::vector<ParticleHistory> ParticleErosion::process()
{
    if (!terrain || !boundariesTree || quantity < 0 || posAndDirs.size() < size_t(quantity))
        throw std::invalid_argument("Erosion requires a terrain, boundaries and one emission per particle");
    for (float value : {dt, size, matterDensity, flowResponseTime, sedimentRadius, sedimentDensity,
                        airViscosity, waterViscosity, maxSimulationTime, restSpeed})
        if (!std::isfinite(value) || value <= 0.f)
            throw std::invalid_argument("Erosion time, sizes, densities and viscosities must be finite and positive");
    for (float value : {strength, maxCapacityFactor, erosion, deposit, gravity, shearingStressConstantK,
                        shearingRatePower, erosionPowerValue, criticalShearValue, initialCapacity})
        if (!std::isfinite(value) || value < 0.f)
            throw std::invalid_argument("Erosion strengths and shear parameters must be finite and nonnegative");
    balance = ErosionBalance();
    settlingCacheValid=false;
    if (maxTrackingSteps<=0) throw std::invalid_argument("Erosion tracking step limit must be positive");
    settlingAir=clearSettlingSpeed(sedimentRadius,sedimentDensity,1.f,airViscosity,gravity);
    settlingWater=clearSettlingSpeed(sedimentRadius,sedimentDensity,1000.f,waterViscosity,gravity);
    settlingCacheKey={sedimentRadius,sedimentDensity,gravity,airViscosity,waterViscosity};
    settlingCacheValid=true;
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
            } catch (...) {
                if (!failed.exchange(true)) error = std::current_exception();
            }
        }
    });
    if (error) std::rethrow_exception(error);
    auto audit=[&]() {
        balance.initialSediment=balance.finalSediment=balance.exportedSediment=balance.retainedSediment=0.;
    double exchanges=0.;
    for (const auto& track : tracks) {
        if (track.history.empty()) continue;
        balance.initialSediment+=track.history.front().particleState.capacity;
        const float remaining=track.history.back().particleState.capacity;
        balance.finalSediment+=remaining;
        if (track.termination==ParticleHistory::LEFT_DOMAIN) balance.exportedSediment+=remaining;
        else balance.retainedSediment+=remaining;
        for (const auto& state : track.history)
            if (state.terrainModificationApplied) exchanges+=state.exchange.erosionValue;
    }
    balance.requestedTerrainChange=-exchanges;
    balance.sedimentResidual=balance.initialSediment+exchanges-balance.finalSediment;
    };
    audit();
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
            balance.measure=applyOn==BINARY_VOXELS ? ErosionBalance::BINARY_EFFECTIVE_VOLUME : ErosionBalance::DENSITY_SUM;
            const auto values=voxelGrid->getVoxelValues();
            const GridF& hits=binaryHitAccumulator?*binaryHitAccumulator:binaryHits;
            for (size_t i=0;i<values.size();++i) {
                if (applyOn==BINARY_VOXELS) sum+=(values[i]>0.f?double(binaryHitThreshold)/binaryHitScale:0.)-
                    (hits.getDimensions()==values.getDimensions()?hits[i]/double(binaryHitScale):0.);
                else sum+=values[i];
            }
        }
        return sum;
    };
    const double before=measure();
    terrainModifTime = timeIt([&]() { reconcileTerrainExchanges(tracks); });
    audit();
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
