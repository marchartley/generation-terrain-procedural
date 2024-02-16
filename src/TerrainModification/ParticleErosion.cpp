#include "ParticleErosion.h"

#include  "DataStructure/BVH.h"
#include "Graphics/Mesh.h"
#include "TerrainModification/RockErosion.h"

#include "Graphics/DisplayGraphics.h"

Vector3 ErosionParticle::predictNextPos(float dt)
{
    return pos + (velocity + force * dt) * dt;
}

std::vector<Vector3> ErosionParticle::bounce(float dt, SpacePartitioning *boundaries)
{
    int maxNbBounces = 2;
    std::vector<Vector3> allPositions;
    allPositions.reserve(maxNbBounces);
    Vector3 currentPos = this->pos;
    Vector3 currentNextPosition = this->predictNextPos(dt);
    float remainingDistanceToDestination = (currentNextPosition - currentPos).norm();

    size_t intersectingTriangleIndex = -1; // Store the last triangle that was touched, such that we can ignore it on the next bounce
    Vector3 collisionPoint;
    auto collisionPointAndTriangleIndex = boundaries->getIntersectionAndTriangleIndex(currentPos, currentNextPosition, (rolling ? intersectingTriangleIndex : -1));
    std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
    while (maxNbBounces > 0 && remainingDistanceToDestination > 0 && collisionPoint.isValid()) {
        Vector3 triangleNormal = boundaries->triangles[intersectingTriangleIndex].normal;
        Vector3 bounce = (currentNextPosition - currentPos).reflexion(triangleNormal);

        float distanceToCollision = (collisionPoint - currentPos).norm();
        remainingDistanceToDestination -= distanceToCollision;
        currentPos = collisionPoint;
        currentNextPosition = currentPos + bounce.setMag(remainingDistanceToDestination);

        allPositions.push_back(currentNextPosition);

        auto collisionPointAndTriangleIndex = boundaries->getIntersectionAndTriangleIndex(currentPos, currentNextPosition, intersectingTriangleIndex); // (rolling ? intersectingTriangleIndex : -1));
        std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;

        maxNbBounces--;
        if (std::abs(triangleNormal.normalized().dot((currentNextPosition - currentPos).normalized())) < 1e-3) {
            this->rolling = true;
        }
    }
    allPositions.push_back(currentNextPosition);
    return allPositions;
}

void ErosionParticle::addForce(const Vector3 &f)
{
    this->force += f;
}

void ErosionParticle::addVelocity(const Vector3 &v)
{
    this->velocity += v;
}

Vector3 ErosionParticle::eulerIntegration(float dt)
{
    addVelocity(force * dt);
    pos += velocity * dt;
    return this->pos;
}

Vector3 ErosionParticle::verletIntegration(float dt)
{
    // For now, just do Euler
    return eulerIntegration(dt);
}



ParticleErosion::ParticleErosion()
{

}

void ParticleErosion::initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity)
{
    particle.pos = position;
    particle.dir = velocity;
    particle.force = Vector3();
    particle.capacity = initialCapacity;

    particle.properties = std::make_shared<ParticleProperties>();
    particle.properties->density = density;
    particle.properties->maxCapacity = maxCapacity;
    particle.properties->volume = 1.f;
    particle.properties->radius = radius;
    particle.properties->volume = .74 * M_PI * particle.properties->radius * particle.properties->radius;
    particle.properties->mass = particle.properties->density * particle.properties->volume;
}


ParticleHistory ParticleErosion::trackParticlePositions(ErosionParticle &particle)
{
//    std::vector<float> erosionValues;
//    std::vector<Vector3> positions;
//    BSpline tunnel;
    ParticleHistory history;

    particle.capacity = particle.properties->maxCapacity * initialCapacity;

    float flowfieldInfluence = 50.0;
    int maxSteps = 500 / dt; // An estimation of how many step we need
    int steps = maxSteps;
    int maxBounces = (maxCollisions < 0 ? 10000 : maxCollisions);

    if (terrain->checkIsInGround(particle.pos))
            return history; //return {{}, {}, BSpline()};

    history.history.reserve(maxSteps);

    bool continueSimulation = true;
    bool hasBeenAtLeastOnceInside = false;
//    Vector3 lastCollisionPoint(false);
    int lastBouncingTime = 10000;
//    size_t intersectingTriangleIndex = -1;
//    std::vector<Vector3> lastCollisions;
    std::tuple<int, int> firstLastCollisionIndices = {-1, -1};
    while (continueSimulation) {
        float erosionApplied = 0.f;
        float environmentDensity = environmentalDensities.at(particle.pos);
        float particleVolume = .1f;
        float gravityCoefficient = particleVolume * gravity * (matterDensity - environmentDensity) * 0.01f * (matterDensity < 0 ? 0.f : 1.f);
        Vector3 justFlow = flowfieldValues(particle.pos);
        Vector3 justGravity = gravityDefault * gravityCoefficient;

        particle.force *= 0.f;
        particle.addForce(justGravity);
//        particle.addVelocity(justFlow * flowfieldInfluence * dt);
        std::vector<Vector3> positionsWithBounces = particle.bounce(dt, boundariesTree);
        Vector3 nextPos = positionsWithBounces.back();
        bool collisionsOccured = positionsWithBounces.size() > 1;

        if (Vector3::isInBox(nextPos.xy(), -terrainSize.xy(), terrainSize.xy()))
            hasBeenAtLeastOnceInside = true;
//        auto collisionPointAndTriangleIndex = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, nextPos, (particle.rolling ? intersectingTriangleIndex : -1));
//        Vector3 collisionPoint;
//        std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
//        if (particle.rolling || collisionPoint.isValid()) {
        if (particle.rolling || collisionsOccured) {
            erosionApplied = 0.01f;
            lastBouncingTime = 0;
            if (std::get<0>(firstLastCollisionIndices) < 0) std::get<0>(firstLastCollisionIndices) = steps;
            std::get<1>(firstLastCollisionIndices) = steps;

            maxBounces --;

            if (steps < 0 || maxBounces <= 0)
                continueSimulation = false;

            particle.justStartedRolling = false;

            Vector3 normal;
            float vRel;

            if (!particle.rolling && lastBouncingTime < 3 && collisionsOccured) {
                particle.rolling = true;
                particle.justStartedRolling = true;
            } else if (particle.rolling && collisionsOccured) {
                particle.rolling = false;
            }
            vRel = particle.velocity.norm();

//            if (!particle.rolling && lastBouncingTime < 3 && collisionPoint.isValid()) {
//                particle.rolling = true;
//                particle.justStartedRolling = true;
//            } else if (particle.rolling && collisionPoint.isValid()) {
//                particle.rolling = false;
//                lastCollisions.push_back(collisionPoint);
//                if (lastCollisions.size() > 5) {
//                    lastCollisions.erase(lastCollisions.begin());
//                    Vector3 mean;
//                    for (const auto& p : lastCollisions)
//                        mean += p;
//                    mean /= float(lastCollisions.size());
//                    float error = 0.f;
//                    for (const auto& p : lastCollisions)
//                        error += (p - mean).norm2();
//                }

//            }


//            if (particle.rolling && !particle.justStartedRolling) {
//                vRel = particle.velocity.norm();
//            } else {
//                normal = boundariesTree->triangles[intersectingTriangleIndex].normal;
//                float dotNormal = std::abs(particle.velocity.normalized().dot(normal));
//                vRel = particle.velocity.norm() * (1.f - dotNormal); // velocity relative to the surface

//                if (!continueSimulation && normal.z < 0)
//                    continueSimulation = true;
//            }

//            if (particle.justStartedRolling) {
//                particle.velocity = (particle.velocity - (normal * particle.velocity.dot(normal))).setMag(particle.velocity.norm());
//                collisionPointAndTriangleIndex = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, nextPos, intersectingTriangleIndex);
//                std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
//                if (collisionPoint.isValid())
//                    normal = boundariesTree->triangles[intersectingTriangleIndex].normal;
//            }

//            lastBouncingTime = 0;

//            float amountToErode = computeErosionValue(particle, vRel);
//            float amountToDeposit = computeDepositionValue(particle);

////            float amountToErodeIfErosionWasOn = std::pow(std::max(0.f, shear), erosionPowerValue);
////            float amountToDepositIfErosionWasOn = 0.1f;

//            if (amountToErode - amountToDeposit != 0) {
//                particle.capacity += (amountToErode - amountToDeposit);
//                if (nextPos.z >= 0.f) {
////                    erosionValuesAndPositions.push_back({amountToErode - amountToDeposit, nextPos});

//                    erosionValues.push_back(amountToErode - amountToDeposit);
//                    positions.push_back(nextPos);
//                }
//            }
//            if (collisionPoint.isValid()) {
//                // Continue the rock tracing
//                int maxCollisions = 2;
//                do {
//                    maxCollisions --;
//                    if (maxCollisions <= 0 || !continueSimulation) {
//                        continueSimulation = false;
//                        break;
//                    }
//                    Vector3 bounce = particle.velocity.reflexion(normal);
//                    bounce -= normal * bounce.dot(normal) * (1.f - bounciness);
//                    particle.velocity = bounce * bouncingCoefficient;
////                            particle.dir.minMagnitude(.5f);

//                    particle.pos = collisionPoint;
//                    lastCollisionPoint = collisionPoint;
//                    std::tie(collisionPoint, intersectingTriangleIndex) = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, nextPos, intersectingTriangleIndex);
//                } while (collisionPoint.isValid());
//            } else if (particle.rolling){
//                particle.velocity *= bouncingCoefficient;
////                    particle.dir = (particle.dir - (normal * particle.dir.dot(normal))).setMag(particle.dir.norm());
//            }

        } else {

        }
        lastBouncingTime ++;
        steps --;

        particle.velocity = (nextPos - particle.pos) / dt;
        particle.pos = nextPos; //particle.pos + (particle.velocity * dt);
        if (wrapPositions)
            particle.pos = Vector3::wrap(particle.pos, Vector3(0, 0, -100), Vector3(terrain->getSizeX(), terrain->getSizeY(), 1000));
        particle.velocity *= 0.99f;
        if (steps < 0 || (hasBeenAtLeastOnceInside && nextPos.z < -20) || particle.pos.z < -20 || particle.velocity.norm2() < 1e-8 || !continueSimulation || maxBounces <= 0) {
            if (depositFactor > 0.f && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                while (!terrain->checkIsInGround(particle.pos) && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                    particle.pos.z -= .5f;
                }
                erosionApplied -= particle.capacity;
                particle.capacity = 0.f;
            }
            continueSimulation = false;
        }
        history.add(particle);
        if (erosionApplied != 0.f)
            history.history.back().terrainModificationApplied = true;
    }
    return history;
}
//*/


/*
 * The one I'm sure works :

ParticleHistory ParticleErosion::trackParticlePositions(ErosionParticle &particle)
{
    bool rolling = false;
//    std::vector<float> erosionValues;
//    std::vector<Vector3> positions;
//    BSpline tunnel;
//    std::vector<std::pair<float, Vector3>> erosionValuesAndPositions;

    float capacity = particle.properties->maxCapacity * initialCapacity;

    float flowfieldInfluence = 50.0;
    int maxSteps = 500 / dt; // An estimation of how many step we need
    int steps = maxSteps;
    int maxBounces = (maxCollisions < 0 ? 10000 : maxCollisions);

    ParticleHistory history;

    if (terrain->checkIsInGround(particle.pos))
            return history;

    history.history.reserve(maxSteps);

    bool continueSimulation = true;
    bool hasBeenAtLeastOnceInside = false;
    Vector3 lastCollisionPoint(false);
    int lastBouncingTime = 10000;
    size_t intersectingTriangleIndex = -1;
    std::vector<Vector3> lastCollisions;
    std::tuple<int, int> firstLastCollisionIndices = {-1, -1};
    while (continueSimulation) {
        float erosionApplied = 0.f;
        Vector3 nextPos = particle.pos + (particle.dir * dt);
        float environmentDensity = environmentalDensities.at(particle.pos);
        float particleVolume = .1f;
        float gravityCoefficient = particleVolume * gravity * (matterDensity - environmentDensity) * 0.01f * (matterDensity < 0 ? 0.f : 1.f);
        Vector3 justFlow = flowfieldValues(particle.pos);
        Vector3 justGravity = gravityDefault * gravityCoefficient;
        Vector3 flowfield = justFlow + justGravity * dt;

        if (matterDensity < 0)
            particle.dir = (particle.dir + flowfield * flowfieldInfluence) * .5f;
        else
            particle.dir += flowfield * flowfieldInfluence * dt;

        particle.dir.maxMagnitude(maxSpeed);
        if (Vector3::isInBox(nextPos.xy(), -terrainSize.xy(), terrainSize.xy()))
            hasBeenAtLeastOnceInside = true;
        auto collisionPointAndTriangleIndex = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, (rolling ? intersectingTriangleIndex : -1));
        Vector3 collisionPoint;
        std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
        if (rolling || collisionPoint.isValid()) {
            if (std::get<0>(firstLastCollisionIndices) < 0) std::get<0>(firstLastCollisionIndices) = steps;
            std::get<1>(firstLastCollisionIndices) = steps;

            maxBounces --;

            if (steps < 0 || maxBounces <= 0)
                continueSimulation = false;

            bool justStartedRolling = false;

            Vector3 normal;
            float vRel;

            if (!rolling && lastBouncingTime < 3 && collisionPoint.isValid()) {
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
                normal = boundariesTree->triangles[intersectingTriangleIndex].normal;
                float dotNormal = std::abs(particle.dir.normalized().dot(normal));
                vRel = particle.dir.norm() * (1.f - dotNormal); // velocity relative to the surface

                if (!continueSimulation && normal.z < 0)
                    continueSimulation = true;
            }

            if (justStartedRolling) {
                particle.dir = (particle.dir - (normal * particle.dir.dot(normal))).setMag(particle.dir.norm());
                collisionPointAndTriangleIndex = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, intersectingTriangleIndex);
                std::tie(collisionPoint, intersectingTriangleIndex) = collisionPointAndTriangleIndex;
                if (collisionPoint.isValid())
                    normal = boundariesTree->triangles[intersectingTriangleIndex].normal;
            }

            lastBouncingTime = 0;

            float amountToErode = computeErosionValue(particle, vRel);
            float amountToDeposit = computeDepositionValue(particle);

            if (amountToErode - amountToDeposit != 0) {
                erosionApplied = amountToDeposit - amountToErode;
                particle.capacity -= erosionApplied;
                if (nextPos.z >= 0.f) {
//                    erosionValuesAndPositions.push_back({amountToErode - amountToDeposit, nextPos});

//                    erosionValues.push_back(amountToErode - amountToDeposit);
//                    positions.push_back(nextPos);
                }
            }
            if (collisionPoint.isValid()) {
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
                    std::tie(collisionPoint, intersectingTriangleIndex) = boundariesTree->getIntersectionAndTriangleIndex(particle.pos, particle.pos + particle.dir, intersectingTriangleIndex);
                } while (collisionPoint.isValid());
            } else if (rolling){
                particle.dir *= bouncingCoefficient;
//                    particle.dir = (particle.dir - (normal * particle.dir.dot(normal))).setMag(particle.dir.norm());
            }
        } else {

        }
        lastBouncingTime ++;
        steps --;

//        tunnel.points.push_back(particle.pos);

        particle.pos = particle.pos + (particle.dir * dt);
        if (wrapPositions)
            particle.pos = Vector3::wrap(particle.pos, Vector3(0, 0, -100), Vector3(terrain->getSizeX(), terrain->getSizeY(), 1000));
        particle.dir *= 0.99f;
        if (steps < 0 || (hasBeenAtLeastOnceInside && nextPos.z < -20) || particle.pos.z < -20 || particle.dir.norm2() < 1e-4 || !continueSimulation || maxBounces <= 0) {
            if (depositFactor > 0.f && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                while (!terrain->checkIsInGround(particle.pos) && Vector3::isInBox(particle.pos, Vector3(), terrain->getDimensions())) {
                    particle.pos.z -= .5f;
                }
                //Vector3 depositPosition = particle.pos;
                //erosionValues.push_back(-particle.capacity);
                //positions.push_back(depositPosition);
                erosionApplied -= particle.capacity;
                particle.capacity = 0.f;
            }
            continueSimulation = false;
        }
        history.add(particle);
        if (erosionApplied != 0.f)
            history.history.back().terrainModificationApplied = true;
    }
    return history;
}
*/

float ParticleErosion::computeErosionValue(const ErosionParticle& particle, float vRel)
{
    float l = particleSize * 0.001f;
    float theta = vRel / l;
    float shear = shearingStressConstantK * std::pow(theta, shearingRatePower);
    float amountToErode = erosionFactor * std::pow(std::max(0.f, shear)/* - criticalShearValue)*/, erosionPowerValue);
    if (maxCollisions != -1) {
        amountToErode *= 1000000.f;
    }
    amountToErode = std::min(amountToErode, particle.properties->maxCapacity - particle.capacity);
    return amountToErode;
}

float ParticleErosion::computeDepositionValue(const ErosionParticle& particle)
{
    float densitiesRatio = 1.f; //1000.f / matterDensity; // The ground should have a density, I'll assume the value 1000 for now
    //                float u = (2.f / 9.f) * particle.radius * particle.radius * (matterDensity - environmentDensity) * gravity * (1.f - (particle.capacity / particle.maxCapacity));
    float amountToDeposit = std::min(particle.capacity, depositFactor * (densitiesRatio * particle.capacity * (particle.capacity / particle.properties->maxCapacity)));
    return amountToDeposit;
}

std::vector<ErosionPoint> ParticleErosion::computeErosionValuesFromHistory(const ParticleHistory &history, SpacePartitioning &boundaries)
{
    std::vector<ErosionPoint> erosions;
    for (auto& particleState : history.history) {
        auto& particle = particleState.particleState;

        if (particleState.terrainModificationApplied) {
            ErosionPoint point;
            point.position = particle.pos;
            point.erosionValue = computeErosionValue(particle, particle.velocity.norm()) - computeDepositionValue(particle);
            erosions.push_back(point);
        }
    }
    return erosions;
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
    std::vector<ImplicitNaryOperator*> allNary(allErosions.size());
    for (size_t i = 0; i < allNary.size(); i++) {
        allNary[i] = new ImplicitNaryOperator;
    }

    #pragma omp parallel for
    for (size_t i = 0; i < allErosions.size(); i++) {
        const auto& erosionValuesAndPositions = allErosions[i];
        Vector3 previousPos(false);
        for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
            Vector3 pos = erosionValuesAndPositions[iRock].position;
            float val = erosionValuesAndPositions[iRock].erosionValue;
            int size = particleSize;
            val *= 100.f;
            if (iRock == erosionValuesAndPositions.size() - 1) {
                val *= .1f;
            }
//                    if (val <= 0.f) continue; // Fits the paper Paris et al. Terrain Amplification with Implicit 3D Features (2019)
            float dimensions = size;
            if (dimensions < 1.f) continue;
            previousPos = pos;
            auto sphere = dynamic_cast<ImplicitPrimitive*>(ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(dimensions, dimensions, dimensions), std::abs(val)*.1f));
            sphere->material = (val > 0 ? TerrainTypes::AIR : TerrainTypes::DIRT);
            sphere->dimensions = Vector3(dimensions, dimensions, dimensions);
            sphere->supportDimensions = sphere->dimensions;
            sphere->position = pos - sphere->dimensions * .5f;
            allNary[i]->composables.push_back(sphere);
        }
    }

    ImplicitNaryOperator* totalErosion = new ImplicitNaryOperator;
    for (auto& op : allNary)
        if (!op->composables.empty())
            totalErosion->composables.push_back(op);
    if (!totalErosion->composables.empty())
    {
        dynamic_cast<ImplicitNaryOperator*>(implicitTerrain)->addChild(totalErosion);
        implicitTerrain->_cached = false;
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
                        float startZ = std::max(0.f, pos.z - halfHeight);
                        float endZ = pos.z + halfHeight;
                        layerBasedGrid->transformLayer(pos.x + x, pos.y + y, startZ, endZ, (val > 0 ? TerrainTypes::AIR : TerrainTypes::SAND));
                    }
                }
            }
        }
    }
}














//std::tuple<std::vector<BSpline>, int, int, std::vector<std::vector<std::pair<float, Vector3>>>>
std::vector<ParticleHistory> ParticleErosion::process()
{
    matterDensity = std::max(matterDensity, 1e-1f);
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
//        if (pos.z < 0) return 1000.f;
//        return 1.f;
//    };
//    auto flowfieldValues = [&](const Vector3& pos) {
//        return Vector3();
//    };
//    GridV3 gravityfieldValues = GridV3(terrainSize, Vector3(0, 0, -gravity));
//    gravityfieldValues.raiseErrorOnBadCoord = false;
//    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);

    environmentalDensities = voxelGrid->getEnvironmentalDensities(); // Here also
    flowfieldValues = GridV3(terrainSize);
    if (flowType == FLOWFIELD_TYPE::BASIC) {
        Vector3 airDir = Vector3(0, airForce, 0).rotate(0, 0, (airFlowfieldRotation / 180) * PI);
        Vector3 waterDir = Vector3(0, waterForce, 0).rotate(0, 0, (waterFlowfieldRotation / 180) * PI);
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
        flowfieldValues = EnvObject::flowfield.resize(terrainSize).meanSmooth();
        for (auto& v : flowfieldValues)
            v = v.xy();
    }
    flowfieldValues.raiseErrorOnBadCoord = false;
    flowfieldValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

//    Plotter::get()->addImage(flowfieldValues.sliceXY(flowfieldValues.sizeZ / 2));
//    Plotter::get()->show();

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
    std::vector<ErosionParticle> particles(quantity);

//    auto startTime = std::chrono::system_clock::now();
//    float totalCollisionTime = 0.f;
//    float totalOtherTime = 0.f;

    particleSimulationTime = timeIt([&]() {
        #pragma omp parallel for
        for (int i = 0; i < quantity; i++)
        {
            auto [initialPos, initialDir] = posAndDirs[i];

            ErosionParticle& particle = particles[i];
            initializeParticle(particle, initialPos, initialDir, .5f, matterDensity, capacity, maxCapacity);

            auto particleTracking = trackParticlePositions(particle);
            auto erosionValues = particleTracking.getErosionValues();
            auto positions = particleTracking.getErosionPositions();
//            auto tunnel = BSpline(particleTracking.getPositions());
            std::vector<std::pair<float, Vector3>> erosionValuesAndPositions(erosionValues.size());
            for (size_t j = 0; j < erosionValuesAndPositions.size(); j++) {
                erosionValuesAndPositions[j] = {erosionValues[j], positions[j]};
            }
//            nbPos[i] = tunnel.points.size();
//            nbErosions[i] = erosionValuesAndPositions.size();
//            tunnels[i] = tunnel;

//            allErosions[i] = erosionValuesAndPositions;
            tracks[i] = particleTracking;
        }
    });

    int positions = 0, erosions = 0;
    for (int i = 0; i < quantity; i++) {
        positions += tracks[i].history.size();
        erosions += tracks[i].getErosionValues().size();

//        positions += nbPos[i];
//        erosions += nbErosions[i];
    }

    if (!applyTheErosion || strengthValue == 0.f) {
        std::vector<BSpline> tunnels(quantity);
        std::vector<std::vector<std::pair<float, Vector3>>> allErosions(quantity);
        for (size_t i = 0; i < quantity; i++) {
            tunnels[i] = BSpline(tracks[i].getPositions());
            std::vector<Vector3> erosionPositions = tracks[i].getErosionPositions();
            std::vector<float> erosionValues = tracks[i].getErosionValues();
            allErosions[i].resize(erosionValues.size());
            for (size_t j = 0; j < erosionValues.size(); j++) {
                allErosions[i][j] = {erosionValues[j], erosionPositions[j]};
            }
        }
        return tracks;
//        return {tunnels, positions, erosions, allErosions};
    }

    /*std::vector<std::vector<ErosionPoint>> allErosionsAsErosionPoints(allErosions.size());
    for (size_t i = 0; i < allErosionsAsErosionPoints.size(); i++) {
        allErosionsAsErosionPoints[i] = std::vector<ErosionPoint>(allErosions[i].size());
        for (size_t iRock = 0; iRock < allErosions[i].size(); iRock++) {
            auto [val, pos] = allErosions[i][iRock];
            allErosionsAsErosionPoints[i][iRock].erosionValue = val;
            allErosionsAsErosionPoints[i][iRock].position = pos;
        }
    }*/
    std::vector<std::vector<ErosionPoint>> allErosionsAsErosionPoints(quantity);
    for (size_t i = 0; i < quantity; i++) {
        allErosionsAsErosionPoints[i] = this->computeErosionValuesFromHistory(tracks[i], *boundariesTree);
//        allErosionsAsErosionPoints[i] = tracks[i].getErosionPoints();
    }

    terrainModifTime = timeIt([&]() {
        if (asImplicit) {
            modifyImplicitTerrain(allErosionsAsErosionPoints);
        } else if (asVoxels) {
            float sumBefore = voxelGrid->getVoxelValues().sum();
            modifyVoxelGrid(allErosionsAsErosionPoints);
            float sumAfter = voxelGrid->getVoxelValues().sum();
            std::cout << sumBefore << " -> " << sumAfter << " => " << sumAfter - sumBefore << std::endl;
        } else if (asHeightmap) {
            modifyHeightmap(allErosionsAsErosionPoints);
        } else if (asLayers) {
            modifyLayerBased(allErosionsAsErosionPoints);
        }
    });
    return tracks;
}

std::vector<Vector3> ParticleHistory::getPositions()
{
    std::vector<Vector3> positions(this->history.size());
    for (size_t i = 0; i < positions.size(); i++)
        positions[i] = history[i].particleState.pos;
    return positions;
}

std::vector<Vector3> ParticleHistory::getErosionPositions()
{
    std::vector<Vector3> positions;
    for (size_t i = 0; i < history.size(); i++)
        if (history[i].terrainModificationApplied)
            positions.push_back(history[i].particleState.pos);
    return positions;
}
std::vector<float> ParticleHistory::getErosionValues()
{
    std::vector<float> values;
    for (size_t i = 0; i < history.size(); i++)
        if (history[i].terrainModificationApplied && i > 0)
            values.push_back(history[i].particleState.capacity - history[i - 1].particleState.capacity);
    return values;
}

std::vector<ErosionPoint> ParticleHistory::getErosionPoints()
{
    std::vector<float> values = getErosionValues();
    std::vector<Vector3> positions = getErosionPositions();
    std::vector<ErosionPoint> points(values.size());
    for (size_t i = 0; i < values.size(); i++) {
        points[i].position = positions[i];
        points[i].erosionValue = values[i];
    }
    return points;
}

std::vector<float> ParticleHistory::getCapacities()
{
    std::vector<float> capacities(this->history.size());
    for (size_t i = 0; i < capacities.size(); i++)
        capacities[i] = history[i].particleState.capacity;
    return capacities;
}
