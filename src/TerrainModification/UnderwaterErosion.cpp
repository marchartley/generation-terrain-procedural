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
    : grid(grid), maxRockSize(maxRockSize), rockAmount(rockAmount), maxRockStrength(maxRockStrength)
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

//void updateParticle(Vector3& pos, Vector3& dir)

std::tuple<std::vector<std::vector<Vector3>>, std::vector<std::vector<Vector3>>>
UnderwaterErosion::Apply(Vector3 startingPoint, Vector3 originalDirection, float randomnessFactor, bool fallFromSky,
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
                         float waterForce)
{

    auto startTime = std::chrono::system_clock::now();
    std::vector<std::vector<Vector3>> debugFinishingLines;
    std::vector<std::vector<Vector3>> debugFailingLines;
    std::vector<std::vector<Vector3>> debugLines;
    float starting_distance = pow(grid->getDimensions().maxComp()/2.0, 2); // pow(std::max(grid->sizeX, std::max(grid->sizeY, grid->sizeZ))/2.0, 2);
    starting_distance = sqrt(3 * starting_distance); // same as sqrt(x+y+z)
    starting_distance *= 2.0; // Leave a little bit of gap
    int max_iter = 1000;
    int total_iterations = 0;
    int cpt = 0;

    // Hydraulic erosion parameters
    float maxCapacity = 1.f * this->maxRockStrength * maxCapacityFactor;
    float erosionFactor = .01f * this->maxRockStrength * erosion;
    float depositFactor = .01f * this->maxRockStrength * deposit;
//    std::vector<std::tuple<Vector3, Vector3>> usefulStartingPositions;

    Matrix3<float> initialMapValues = grid->getVoxelValues();
    initialMapValues.raiseErrorOnBadCoord = false;
    initialMapValues.defaultValueOnBadCoord = -1000;
    // Use normals from initial map, it shouldn't change too much...
    Matrix3<Vector3> normals = initialMapValues.gradient(); //initialMapValues.binarize().toDistanceMap().gradient();
    normals.raiseErrorOnBadCoord = false;

    Matrix3<float> modifications(grid->getDimensions());
//    Matrix3<Vector3> flowfieldValues = grid->flowField;

    // Matter information
//    float matterDensity = 1600.0;
    Matrix3<float> particlesSpeeds(grid->getDimensions(), 0.f);

    Matrix3<Vector3> gravityfieldValues = Matrix3<Vector3>(grid->getDimensions(), Vector3(0, 0, -gravity));
    gravityfieldValues.raiseErrorOnBadCoord = false;
    gravityfieldValues.defaultValueOnBadCoord = Vector3(0, 0, -gravity);

    Matrix3<Vector3> flowfieldValues = grid->getFlowfield();
    flowfieldValues.raiseErrorOnBadCoord = false;

    /* BAD PRACTICE, JUST FOR JFIG DEMO*/
    Vector3 airDir = Vector3(0, airForce + .0001f, 0).rotate(0, 0, (airFlowfieldRotation / 180) * PI);
    Vector3 waterDir = Vector3(0, waterForce + .0001f, 0).rotate(0, 0, (waterFlowfieldRotation / 180) * PI);
    for (size_t i = 0; i < flowfieldValues.size(); i++) {
        if (grid->getEnvironmentalDensities().at(i) < 100) { // In the air
            flowfieldValues.at(i) = airDir;
        } else { // In water
            flowfieldValues.at(i) = waterDir;
        }
    }
    /* END BAD */

    std::vector<BSpline> tunnels;
    for (int i = 0; i < this->rockAmount && max_iter > 0; i++)
    {
        BSpline tunnel;
        Matrix3<float> currentMapValues = initialMapValues + modifications;
        currentMapValues.raiseErrorOnBadCoord = false;
        currentMapValues.defaultValueOnBadCoord = -1000;

        float capacity = 0.f;

        float flowfieldInfluence = 1.0;
        cpt ++;
        int steps = 10 * starting_distance; // An estimation of how many step we need
        Vector3 pos;
        Vector3 firstHitPos;
        if (!startingPoint.isValid()) {
            pos = Vector3(random_gen::generate(-1.0, 1.0), random_gen::generate(-1.0, 1.0), random_gen::generate(0.0, 1.0));
            pos.normalize();
            pos *= starting_distance;
        } else {
            pos = startingPoint;
        }
        Vector3 dir = Vector3::random();
        if (originalDirection.isValid())
            dir = originalDirection.normalize();

        if (fallFromSky) {
            pos = Vector3(random_gen::generate(0.f, this->grid->getSizeX()), random_gen::generate(0.f, this->grid->getSizeY()), this->grid->getSizeZ() / 2.f);
            dir = Vector3(0, 0, -.1f);
        }

//        RockErosion rock(random_gen::generate(0.0, this->maxRockSize), random_gen::generate(0.0, this->maxRockStrength) * (random_gen::generate() > .5f ? 1.f : -1.f));
        std::vector<Vector3> coords;

        bool touched = false;
        bool hasBeenAtLeastOnceInside = false;
        bool firstHit = true;
        while (!touched) {
            total_iterations ++;
//            std::cout << "Dir mag = " << dir.norm() << std::endl;
            if (this->grid->contains(pos + dir)) {
                float environmentDensity = this->grid->getEnvironmentalDensities().at(pos + dir);
                float gravityCoefficient = std::max(1.f - (environmentDensity / matterDensity), -1.f); // Keep it between -1 and 1
                Vector3 flowfield = flowfieldValues.at(pos + dir) + gravityfieldValues.at(pos + dir) * gravityCoefficient;
//                dir += flowfield * dirDotGrad * flowfieldFactor * (dist == 0 ? weaknessAgainstFlowfield : 1.0);
                dir += flowfield * flowfieldInfluence;
//                dir.normalize(); // Maybe to remove, who knows...
                hasBeenAtLeastOnceInside = true;
                bool justHit = false;
                if (currentMapValues.at(pos + dir) > 0.0) { // Hit a wall
                    Vector3 normal = normals.at(pos + dir).normalized();
                    float speedRate = (dir.norm() / maxSpeed);
                    float coef = 1.f; // std::max(std::abs(normal.x), std::abs(normal.y));
//                    rock.Apply(this->grid, pos, false, false);
//                    // USING DOT PRODUCT
//                    // Erosion part
                    float amountToErode = std::max(std::min(maxCapacity - capacity, coef * speedRate * (std::abs(dir.normalized().dot(normal))) * erosionFactor * (1.f - std::abs(normal.dot(Vector3(0, 0, 1))))), 0.f);
//                    // Deposition part
                    float amountToDeposit = std::max(std::min(capacity, coef * speedRate * (1 - std::abs(dir.normalized().dot(normal))) * depositFactor), 0.f);
                    // USING SPEED
                    // Erosion part
//                    float amountToErode = std::max(std::min(maxCapacity - capacity, coef * std::min(speedRate, 1.f) * erosionFactor / (materialImpact != 0 ? currentMapValues.at(pos + dir) * materialImpact : 1.f)), 0.f);
                    // Deposition part
//                    float amountToDeposit = std::max(std::min(capacity, (1 - speedRate) * depositFactor * coef), 0.f);
//                    std::cout << amountToErode << " " << amountToDeposit << std::endl;
                    if (!firstHit /*&& (firstHitPos - pos).norm2() > maxRockSize*maxRockSize*/) {
                        if (amountToErode - amountToDeposit != 0) {
                            capacity += (amountToErode - amountToDeposit);
                            int division = 1;
                            RockErosion erosionRock(random_gen::generate(0.0, this->maxRockSize), (amountToErode)/float(division));
                            RockErosion depositRock(random_gen::generate(0.0, this->maxRockSize), (-amountToDeposit)/float(division));
        //                    float previousRockSum = rock.attackMask.sum();
        //                    rock.attackMask *= (amountToErode - amountToDeposit);
        //                    rock.attackMask /= (previousRockSum);
                            // Increase resolution
                            for (int divis = 0; divis < division; divis++) {
                                erosionRock.computeErosionMatrix(modifications, pos + (dir) * (float(divis)/float(division)), false, false);
                                depositRock.computeErosionMatrix(modifications, pos + (dir + Vector3(0, 0, -5.f)) * (float(divis)/float(division)), false, false);
//                                depositRock.computeErosionMatrix(modifications, pos + (dir + gravityfieldValues.at(pos + dir) * gravityCoefficient).normalized() * (float(divis)/float(division)), false, false);
                            }
    //                        rock.computeErosionMatrix(modifications, pos + dir, false, false);
                        }
                        debugFinishingLines.push_back(coords);
                    } else {
                        // Ignore the first collision
                        firstHitPos = pos;
                        if (currentMapValues.at(pos + dir) > 0.0) {
                            firstHit = false;
                            tunnel.points.clear();
                        }
                    }

                    // Continue the rock tracing
//                    touched = true;
                    if (!justHit) {
                        if (std::abs(dir.normalized().dot(normal)) == 1)
                            normal = (normal + Vector3::random(0.1f)).normalized();
                        Vector3 bounce = dir.reflexion(normal) * dir.norm();
                        bounce -= normal * bounce.dot(normal) * (1.f - bounciness);
                        dir = bounce * bouncingCoefficient;
//                        dir = (bounce * (bounciness + .5f)*2.f * bouncingCoefficient + dir * (1 - ((bounciness + .5f)*2.f)))*.5f; // dir.reflexion(normal) * dir.norm() * bouncingCoefficient;
                    }
                    justHit = true;
                    pos += dir;
                    tunnel.points.push_back(pos);
                }
            }
            else {
                dir += Vector3::random(randomnessFactor);
                dir.normalize();
            }
            steps --;
            if (dir.norm2() > maxSpeed*maxSpeed)
                dir = dir.normalize() * maxSpeed;
            else if (dir.norm2() < minSpeed*minSpeed)
                dir = dir.normalize() * minSpeed;
            pos += dir.normalized();
            if (!firstHit) {
                coords.push_back((pos - dir.normalized()));
            }
            if ((hasBeenAtLeastOnceInside && !grid->contains(pos)) || steps < 0 || dir.norm2() < 1e-3) {
                // Failed to hit, is now gone
//                i--;
                max_iter --;
//                if (!grid->contains(pos) && coords.size() > 0)
//                    coords.pop_back();
                if (coords.size() > 1)
                    debugLines.push_back(coords);
//                if (returnEvenLostRocks && hasBeenAtLeastOnceInside) {
//                    debugFailingLines.push_back(coords);
//                }
                if ((steps < 0 || dir.norm2() < 1e-3) && depositFactor > 0.f) {
                    RockErosion rock(random_gen::generate(0.0, this->maxRockSize), -capacity);
                    rock.computeErosionMatrix(modifications, pos - Vector3(0, 0, this->maxRockSize), false);
                }
                break;
            }
        }
        tunnels.push_back(tunnel);
    }
    /*modifications = Matrix3<float>(modifications.getDimensions(), 0.f);
    RockErosion rock(this->maxRockSize, this->maxRockStrength);
    for (BSpline& path : tunnels) {
        float nb_points_on_path = path.length() / (float)(this->maxRockSize / 3.f);
        std::vector<Vector3> coords = path.getPath(nb_points_on_path);
        for (size_t i = 0; i < coords.size() - 1; i++) {
            modifications = rock.computeErosionMatrix(modifications, coords[i]);
            RockErosion speedRock(this->maxRockSize, -(coords[i - 1] - coords[i]).norm());
            particlesSpeeds = speedRock.computeErosionMatrix(particlesSpeeds, coords[i]);
        }
    }*/
//    for (size_t i = 0; i < modifications.size(); i++)
//        if (initialMapValues[i] > 0)
//            modifications[i] = particlesSpeeds[i] * modifications[i] / initialMapValues[i];
    this->grid->applyModification(modifications);
    auto endTime = std::chrono::system_clock::now();
    std::cout << "Simulation en " << std::chrono::duration<float>(endTime - startTime).count() << "s" << std::endl;
    return std::make_tuple(debugLines, debugFailingLines);
//    return std::make_tuple(debugFinishingLines, debugFailingLines);
}


std::vector<Vector3> UnderwaterErosion::CreateTunnel(int numberPoints, bool addingMatter, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    BSpline curve = BSpline(numberPoints); // Random curve
    for (Vector3& coord : curve.points)
        coord = ((coord + Vector3(1.0, 1.0, 1.0)) / 2.0) * grid->getDimensions(); //Vector3(grid->sizeX, grid->sizeY, grid->sizeZ);
    return CreateTunnel(curve, addingMatter, true, applyChanges, startingShape, endingShape);
}
std::vector<Vector3> UnderwaterErosion::CreateTunnel(BSpline path, bool addingMatter, bool usingSpheres, bool applyChanges,
                                                     KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    Matrix3<float> erosionMatrix(grid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
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
        grid->applyModification(erosionMatrix);
    }
//    if (applyChanges)
//        grid->remeshAll();
    return coords;
}

std::vector<std::vector<Vector3> > UnderwaterErosion::CreateMultipleTunnels(std::vector<BSpline> paths, bool addingMatter, bool usingSpheres, bool applyChanges,
                                                                            KarstHolePredefinedShapes startingShape, KarstHolePredefinedShapes endingShape)
{
    Matrix3<float> erosionMatrix(grid->getDimensions()); // this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
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
    grid->applyModification(erosionMatrix);
//    if (applyChanges)
//        grid->remeshAll();
    return allCoords;
}

std::vector<Vector3> UnderwaterErosion::CreateCrack(Vector3 start, Vector3 end, bool applyChanges)
{
    float rx = 3.f, ry = 3.f, rz = 3.f;
    Vector3 ratio(rx, ry, rz);
    Matrix3<int> resizedMap = this->grid->getVoxelValues().resize(grid->getDimensions() / ratio).binarize();
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
    Matrix3<float> erosionMatrix(grid->getDimensions()); //(this->grid->sizeX, this->grid->sizeY, this->grid->sizeZ);
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
        grid->applyModification(erosionMatrix);
    }
//    if (applyChanges)
//        grid->remeshAll();
    return coords;
}

