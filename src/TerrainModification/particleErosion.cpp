#include "particleErosion.h"

#include  "DataStructure/BVH.h"
#include "Graphics/Mesh.h"

void _initializeParticle(ErosionParticle& particle, Vector3& position, Vector3& velocity, float radius, float density, float initialCapacity, float maxCapacity)
{
    particle.pos = position;
    particle.dir = velocity;
    particle.density = density;
    particle.capacity = initialCapacity;
    particle.maxCapacity = maxCapacity;
    particle.force = Vector3();
    particle.volume = 1.f;
    particle.radius = radius;
    particle.volume = .74 * M_PI * particle.radius * particle.radius;
    particle.mass = particle.density * particle.volume;
}

std::vector<std::vector<Vector3> > erosion()
{
    Vector3 terrainSize(100, 100, 100);
    GridF terrain(terrainSize);
    terrain.iterateParallel([&](const Vector3& p) {
        terrain(p) = std::abs(random_gen::generate_perlin(p.x, p.y)) * (p.z - 100);
    });
    
    // Hydraulic erosion parameters
    float strengthValue = 1.f;
    float maxCapacityFactor = 1.f;
    float erosion = 1.f;
    float deposit = 1.f;
    float maxCapacity = 1.f * strengthValue * maxCapacityFactor;
    float erosionFactor = .01f * strengthValue * erosion;
    float depositFactor = .01f * strengthValue * deposit;
    float initialCapacity = 1.f;
    float gravity = 9.8f;
    float matterDensity = 500.f;
    float shearingStressConstantK = 1.f;
    float shearingRatePower = .5f;
    float criticalShearValue = 0.f;
    float bounciness = 1.f;
    float bouncingCoefficient = 1.f;
    float erosionPowerValue = 1.f;


    float particleSize = 8.f;
    float radius = particleSize * .5f;
    int particleAmount = 1000;

    float flowfieldInfluence = 1.0;

    float dt = .1f;
    int maxSteps = 500 / dt; // An estimation of how many step we need
    int steps = maxSteps;
    int maxBounces = 10000;

    bool wrapPositions = false;

    BVHTree boundariesTree;
    boundariesTree.build(Triangle::vectorsToTriangles(Mesh::applyMarchingCubes(terrain).getTriangles()));

    GridF modifications(terrainSize);

    Vector3 gravityDefault = Vector3(0, 0, -gravity);
    
    GridF environmentalDensities(terrainSize);
    environmentalDensities.iterateParallel([&](const Vector3& p) {
        environmentalDensities(p) = (p.z < terrainSize.z * .3f ? 1.f : 1000.f);
    });
    GridV3 flowfieldValues = GridV3(terrainSize);
    flowfieldValues.raiseErrorOnBadCoord = false;
    flowfieldValues.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;

    std::vector<GridF> submodifications(particleAmount, GridF(modifications.getDimensions()));
    std::vector<std::vector<std::pair<float, Vector3>>> allErosions(particleAmount);

    std::vector<ErosionParticle> particles(particleAmount);

    std::vector<Vector3> initialPositions(particleAmount);
    std::vector<Vector3> initialDirections(particleAmount);
    for (size_t i = 0; i < particleAmount; i++) {
        initialPositions[i] = Vector3::random().xy() + Vector3(0, 0, terrainSize.z);
        initialDirections[i] = Vector3();
    }

    std::vector<BSpline> tunnels(particleAmount);

    auto startTime = std::chrono::system_clock::now();
    float totalCollisionTime = 0.f;

//    #pragma omp parallel for
    for (int i = 0; i < particleAmount; i++)
    {
        BSpline tunnel;
        bool rolling = false;
        std::vector<std::pair<float, Vector3>> erosionValuesAndPositions;
        float capacity = maxCapacity * initialCapacity;

        ErosionParticle& particle = particles[i];
        _initializeParticle(particle, initialPositions[i], initialDirections[i], .5f, matterDensity, capacity, maxCapacity);

        // If particle starts in the ground, just discard it
        if (terrain(particle.pos) > 0)
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
            float gravityCoefficient = particleVolume * gravity * (matterDensity - environmentDensity) * 0.01f;
//            float gravityCoefficient = 0.f;
            Vector3 justFlow = flowfieldValues(particle.pos); // (flowfieldValues.at(pos) + flowfieldValues.at(nextPos)) * .5f;
            Vector3 justGravity = gravityDefault * gravityCoefficient;
//            Vector3 justGravity = (gravityfieldValues.at(particle.pos) * gravityCoefficient); // (gravityfieldValues.at(nextPos) * gravityCoefficient); // .maxMagnitude(2.f);
//            std::cout << "Z=" << particle.pos.z << " -> flow " << justFlow << " -- gravity: " << justGravity << "(" << gravityDefault << " * " << gravityCoefficient << ") density = " << environmentDensity << " (" << environmentalDensities << ")" << std::endl;
            Vector3 flowfield = justFlow + justGravity;
//            if (!rolling)
                particle.dir += flowfield * flowfieldInfluence * dt;

//            particle.dir.maxMagnitude(maxSpeed);
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
//                    std::cout << "CollisionPoint is " << (collisionPoint.isValid() ? "valid" : "NOT valid") << std::endl;
                }

                lastBouncingTime = 0;

                // Using SPH formula
                // Erosion part
                float l = particleSize * 0.001f;
                float theta = vRel / l;
                float shear = shearingStressConstantK * std::pow(theta, shearingRatePower);
                float amountToErode = erosionFactor * std::pow(std::max(0.f, shear) - criticalShearValue, erosionPowerValue); // std::pow(std::max(0.f, shear - criticalShearStress), erosionPowerValue)  * (1.f - resistanceValue * materialImpact);
//                if (rolling)
//                    amountToErode = 0.f;
//                if (maxCollisions != -1) {
//                    amountToErode *= 1000000.f;
//                }
//                std::cout << theta << " " << shear << " " << amountToErode << " " << particle.maxCapacity << " " << particle.capacity << std::endl;
                amountToErode = std::min(amountToErode, particle.maxCapacity - particle.capacity);


//                    // Deposition part
                float densitiesRatio = 1.f; //1000.f / matterDensity; // The ground should have a density, I'll assume the value 1000 for now
//                float u = (2.f / 9.f) * particle.radius * particle.radius * (matterDensity - environmentDensity) * gravity * (1.f - (particle.capacity / particle.maxCapacity));
                float amountToDeposit = std::min(particle.capacity, depositFactor * (densitiesRatio * capacity * (particle.capacity / particle.maxCapacity))); //maxCapacity));
//                    float amountToDeposit = std::min(particle.capacity, depositFactor * u);

                if (amountToErode - amountToDeposit != 0) {
                    particle.capacity += (amountToErode - amountToDeposit);
                    if (nextPos.z >= 1.f) {
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
                particle.pos = Vector3::wrap(particle.pos, Vector3(0, 0, -100), Vector3(terrainSize.xy() + Vector3(0, 0, 1000)));
            particle.dir *= 0.99f;
            if (steps < 0 || (hasBeenAtLeastOnceInside && nextPos.z < -20) || particle.pos.z < -20 || particle.dir.norm2() < 1e-4 || !continueSimulation || maxBounces <= 0) {
                if (depositFactor > 0.f && Vector3::isInBox(particle.pos, Vector3(), terrainSize)) {
                    Vector3 depositPosition = particle.pos;
                    erosionValuesAndPositions.push_back({-particle.capacity, depositPosition});
                }
                continueSimulation = false;
            }
//            if (steps < 0)
//                continueSimulation = false;
        }
        if (std::get<0>(firstLastCollisionIndices) > std::get<1>(firstLastCollisionIndices)) {
            size_t startIndex = maxSteps - std::get<0>(firstLastCollisionIndices);
            size_t endIndex = maxSteps - std::get<1>(firstLastCollisionIndices);
            if (rolling) { // Here we often have errors...
                while (!Vector3::isInBox(tunnel[endIndex], Vector3(), terrainSize)) {
                    endIndex--;
                }
            }
            tunnels[i] = std::vector<Vector3>(tunnel.begin() + startIndex, tunnel.begin() + endIndex);
            if ((tunnels[i].points.front() - tunnels[i].points.back()).norm2() < particleSize * particleSize) {
                erosionValuesAndPositions.clear();
            }
        }
//        nbPos[i] = tunnel.points.size();
//        nbErosions[i] = erosionValuesAndPositions.size();

        allErosions[i] = erosionValuesAndPositions;
    }
    auto endParticleTime = std::chrono::system_clock::now();
/*
    if (!applyTheErosion) {
        auto flatAllErosions = flattenArray(allErosions);
        int positions = 0;
        for (auto nb : nbPos)
            positions += nb;
        int erosions = flatAllErosions.size();
        return {tunnels, positions, erosions, allErosions};
    }


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
//                            float halfHeight = std::abs(val) * std::sqrt(particleSize*particleSize - (x*x + y*y))/particleSize;
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
        #pragma omp parallel for
        for (size_t i = 0; i < allErosions.size(); i++) {
            const auto& erosionValuesAndPositions = allErosions[i];
            Vector3 previousPos(false);
            for (size_t iRock = 0; iRock < erosionValuesAndPositions.size(); iRock++) {
                auto [val, pos] = erosionValuesAndPositions[iRock];
//                val = (val < 0 ? -1 : 1);
                int size = particleSize;
                val *= 100.f;
                if (iRock == erosionValuesAndPositions.size() - 1) {
    //                size *= 2.f;
    //                val *= .5f;
                }
    //            summary += val;

                if (asVoxels) {
    //                std::cout << RockErosion(size, val).createPrecomputedAttackMask(size).sum() - val << " errors" << std::endl;
                    RockErosion(size, val).computeErosionMatrix(submodifications[i], pos - Vector3(.5f, .5f, .5f));
                } else if (asHeightmap) {
                    RockErosion(size, val).computeErosionMatrix2D(submodifications[i], pos);
    //                std::cout << val << " " << submodifications[i].sum() << std::endl;
                } else if (asImplicit) {
                    if (iRock == erosionValuesAndPositions.size() - 1) {
                        val *= .1f;
                    }
//                    if (val <= 0.f) continue; // Fits the paper Paris et al. Terrain Amplification with Implicit 3D Features (2019)
                    float dimensions = size; //std::abs(size * val);
                    if (dimensions < 1.f) continue;
//                    if (previousPos.isValid() && (previousPos - pos).norm2() < 2.f){
//                        std::cout << "Ignored [" << val << "] because dist = " << (previousPos - pos).norm() << "\n";
//                        continue;
//                    }
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
    }
//    std::cout << std::endl;

    for (const auto& sub : submodifications)
        modifications += sub;

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
    } else if (asImplicit) {
        ImplicitNaryOperator* totalErosion = new ImplicitNaryOperator;
        for (auto& op : allNary)
            if (!op->composables.empty())
                totalErosion->composables.push_back(op);
        if (!totalErosion->composables.empty())
        {
    //        totalErosion->composables.push_back(implicitTerrain->copy());
    //        delete implicitTerrain;
    //        implicitTerrain->~ImplicitPatch();
    //        implicitTerrain = new (implicitTerrain) ImplicitNaryOperator;
            dynamic_cast<ImplicitNaryOperator*>(implicitTerrain)->addChild(totalErosion);
            implicitTerrain->_cached = false;
        }
    } else if (asLayers) {
        // ...
    }

    int positions = 0, erosions = 0;
    for (int i = 0; i < particleAmount; i++) {
        positions += nbPos[i];
        erosions += nbErosions[i];
    }
    auto endModificationTime = std::chrono::system_clock::now();

    particleSimulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(endParticleTime - startTime).count();
    terrainModifTime = std::chrono::duration_cast<std::chrono::milliseconds>(endModificationTime - endParticleTime).count();
    return {tunnels, positions, erosions, allErosions};
    */
    return {};
}
