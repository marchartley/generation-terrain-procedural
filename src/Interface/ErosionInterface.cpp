#include "ErosionInterface.h"

#include "Interface/InterfaceUtils.h"
#include "Interface/TerrainGenerationInterface.h"
#include "EnvObject/EnvObject.h"

#include <chrono>

ErosionInterface::ErosionInterface(QWidget *parent)
    : ActionInterface("erosion", "Rock throwing", "physics", "Erosion on 3D terrain", "erosion.png", parent)
{

}

void ErosionInterface::affectTerrains(std::shared_ptr<Heightmap> heightmap, std::shared_ptr<VoxelGrid> voxelGrid, std::shared_ptr<LayerBasedGrid> layerGrid, std::shared_ptr<ImplicitNaryOperator> implicitPatch)
{
    ActionInterface::affectTerrains(heightmap, voxelGrid, layerGrid, implicitPatch);
    this->erosion = std::make_shared<UnderwaterErosion>(voxelGrid.get(), erosionSize, erosionStrength, erosionQtt);

    const char* vNoShader = "src/Shaders/no_shader.vert";
    const char* fNoShader = "src/Shaders/no_shader.frag";

    this->rocksPathFailure = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->rocksPathFailure.useIndices = false;
    this->rocksPathFailure.shader->setVector("color", std::vector<float>({.7f, .2f, .1f, .5f}));
    this->rocksPathSuccess = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->rocksPathSuccess.useIndices = false;
    this->rocksPathSuccess.shader->setVector("color", std::vector<float>({.1f, .7f, .2f, .5f}));
    this->boundariesMesh = Mesh(std::make_shared<Shader>(vNoShader, fNoShader));
    this->boundariesMesh.useIndices = false;
    this->boundariesMesh.shader->setVector("color", std::vector<float>({.1f, .2f, .7f, .5f}));

    this->erosionProcess = UnderwaterErosion(voxelGrid.get(), this->erosionSize, this->erosionStrength, this->erosionQtt);
    this->erosionProcess.heightmap = heightmap.get();
    this->erosionProcess.implicitTerrain = this->implicitTerrain.get();
    this->erosionProcess.layerBasedGrid = layerGrid.get();

    this->computePredefinedRocksLocations();

    randomObstacles.resize(10);
    for (size_t i = 0; i < randomObstacles.size(); i++) {
        randomObstacles[i] = {Vector3::random(Vector3(100, 100, 100)), random_gen::generate(5.f, 10.f)};
    }

    //    QTimer::singleShot(1000, this, &ErosionInterface::testManyManyErosionParameters);
}

void ErosionInterface::computePredefinedRocksLocations()
{
    if (this->voxelGrid == nullptr) return;
    std::vector<PARTICLE_INITIAL_LOCATION> locs = {SKY, RIVER, RANDOM, RIVER2, UNDERWATER, CENTER_TOP, FROM_X, FROM_BIG_X, EVERYWHERE, JUST_ABOVE_VOXELS, VOLCANO, VOLCANO2, VOLCANO3};

    Vector3 dimensions = voxelGrid->getDimensions();
    for (auto& loc : locs) {
        initialPositionsAndDirections[loc] = std::vector<std::vector<std::pair<Vector3, Vector3>>>(20, std::vector<std::pair<Vector3, Vector3>>(maxParticles));
        auto& poses = initialPositionsAndDirections[loc];
        for (size_t i = 0; i < poses.size(); i++) {
            for (size_t j = 0; j < poses[i].size(); j++) {
                Vector3 position;
                Vector3 direction;
                if (loc == SKY) {
                    position = Vector3::random(Vector3(-30, -30), dimensions.xy() + Vector3(30, 30)) + Vector3(0, 0, dimensions.z);
                    direction = Vector3(0, 0, -1);
                } else if (loc == RIVER) {
                    position = Vector3(15, 50, 50) + Vector3::random(3.f).xy();
                    direction = Vector3(0, 0, 0);
                } else if (loc == RIVER2) {
                    position = Vector3(15, 50, 100) + Vector3::random(25.f).xy();
                    direction = Vector3(0, 0, 0);
                } else if (loc == UNDERWATER) {
                    position = Vector3(5, 0, 0) + Vector3::random(Vector3(0, dimensions.y * .4, dimensions.z * .2f), Vector3(0, dimensions.y *.6, dimensions.z * .5f));
                    direction = (Vector3(1, 0, 0) + Vector3::random(.1f));
                } else if (loc == FROM_X) {
                    position = Vector3::random(Vector3(0, 20, 20), Vector3(0, dimensions.y - 20, dimensions.z - 20)) - Vector3(1, 0, 0);
                    direction = (Vector3(1, 0, 0) + Vector3::random(.5f).xy() + Vector3(0, 0, Vector3::random(.1f).z));
                } else if (loc == FROM_BIG_X) {
                    position = Vector3::random(Vector3(0, 0, 20), Vector3(0, dimensions.y, dimensions.z - 20)) - Vector3(1, 0, 0);
                    direction = (Vector3(1, 0, 0) + Vector3::random(.5f).xy() + Vector3(0, 0, Vector3::random(.1f).z));
                } else if (loc == CENTER_TOP) {
                    position = Vector3::random(Vector3(voxelGrid->getSizeX() * .45f, voxelGrid->getSizeY() * .45f, voxelGrid->getSizeZ() + 2.f), Vector3(voxelGrid->getSizeX() * .55f, voxelGrid->getSizeY()*.55f, voxelGrid->getSizeZ() + 2.f));
                    direction = (Vector3(1, 0, 0) + Vector3::random(.1f));
                } else if (loc == RANDOM) {
                    position = Vector3::random() * 100.f;
                    position.z = std::abs(position.z);
                    direction = -position.normalized();
                    position += voxelGrid->getDimensions().xy() * .5f;
                } else if (loc == EVERYWHERE) {
                    position = Vector3::random(Vector3(), dimensions);
                    direction = Vector3();
                } else if (loc == VOLCANO) {
                    position = Vector3(55, 45, 100) + Vector3::random(8.f).xy();
                    direction = Vector3(0, 0, 0);
                } else if (loc == VOLCANO2) {
                    position = Vector3(50, 50, 100) + Vector3::random(15.f).xy();
                    direction = Vector3(0, 0, 0);
                } else if (loc == VOLCANO3) {
                    position = Vector3(60, 45, 100) + Vector3::random(8.f).xy();
                    direction = Vector3(0, 0, 0);
                }
                poses[i][j] = {position, direction};
            }
        }
    }
//    recomputeAboveVoxelRocksPositions();
}

void ErosionInterface::recomputeAboveVoxelRocksPositions(TerrainModel* terrain)
{
    auto voxels = voxelGrid->getVoxelValues();
    voxels.raiseErrorOnBadCoord = false;
    voxels.defaultValueOnBadCoord = -1;
    auto normals = -voxels.gradient();
//    auto terrainSurface = Mesh().applyMarchingCubes(voxels).vertexArray;
    std::vector<std::pair<Vector3, Vector3>> terrainSurfaceAndNormalInversed;
    Vector3 terrainSize = terrain->getDimensions();
//    GridV3 normals = terrain->getNormals();
    Vector3 down = Vector3(0, 0, -1);
    while (terrainSurfaceAndNormalInversed.size() < maxParticles) {
        for (int x = 0; x < terrainSize.x; x++) {
            for (int y = 0; y < terrainSize.y; y++) {
    //            Vector3 pos = Vector3(x, y) + Vector3::random().xy();
    //            pos.z = terrain->getHeight(pos) + 10.f;
    //            if (random_gen::generate() <= std::pow(pos.z / terrainSize.z, 2.f))
    //                terrainSurfaceAndNormalInversed.push_back({pos, down});
                for (int z = 0; z < terrainSize.z; z++) {
                    Vector3 pos(x + .5f, y + .5f, z + .5f);
                    Vector3 gradient = (z == terrainSize.z - 1 ? Vector3(0, 0, 1) : normals.at(pos).normalized());
                    if (gradient.z > 0 && (!voxels.checkCoord(pos + gradient * 2.f) || (voxels.at(pos) > 0 && voxels.at(pos + gradient * 2.f) < 0))) {
                        terrainSurfaceAndNormalInversed.push_back({pos + gradient * 2.f + Vector3::random(0.5f), -gradient + Vector3::random(0.2f)});
                    }
                }
            }
        }
    }
//    for (int x = 0; x < terrainSize.x; x++) {
//        for (int y = 0; y < terrainSize.y; y++) {
//            std::cout << int(normals(x, y, terrainSize.z - 1).z * 10) * .1f  << " ";
//        }
//        std::cout << std::endl;
//    }
    initialPositionsAndDirections[JUST_ABOVE_VOXELS] = std::vector<std::vector<std::pair<Vector3, Vector3>>>(20); //, std::vector<std::pair<Vector3, Vector3>>(maxParticles));
    auto& poses = initialPositionsAndDirections[JUST_ABOVE_VOXELS];
    for (size_t i = 0; i < poses.size(); i++) {
        std::shuffle(terrainSurfaceAndNormalInversed.begin(), terrainSurfaceAndNormalInversed.end(), random_gen::random_generator);
        poses[i] = std::vector<std::pair<Vector3, Vector3>>(terrainSurfaceAndNormalInversed.begin(), terrainSurfaceAndNormalInversed.begin() + maxParticles);
        /*
        for (size_t j = 0; j < poses[i].size(); j++) {
            Vector3 position;
            Vector3 direction;
            if (terrainSurface.empty()) // No voxel at all
                break;
            position = terrainSurface[int(random_gen::generate(terrainSurface.size()))] + Vector3(0, 0, 1) + Vector3::random(1.f).xy();
            direction = Vector3(0, 0, 0);

            poses[i][j] = {position, direction};
        }
        */
    }
}

std::tuple<float, float, float> ErosionInterface::computeTerrainBoundaries(TerrainModel* terrain, BVHTree* boundariesTree)
{
    float sumGeometry, sumBVH, sumMeshingBoundaries;

    Vector3 terrainDims = terrain->getDimensions();
    std::vector<std::vector<Vector3>> triangles;
    Mesh m;
    Vector3 geomSize = Vector3::min(terrainDims, Vector3(100, 100, 50));
    sumGeometry = timeIt([&]() {
        if (applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN) {
            m = Mesh::applyMarchingCubes(layerGrid->voxelize().resize(geomSize));
        } else {
            m = terrain->getGeometry(geomSize);
        }
    });
    triangles = m.getTriangles();

//            triangles.push_back({terrainDims * Vector3(0, 0, 0.1), terrainDims * Vector3(1, 0, 0.1), terrainDims * Vector3(1, 1, 0.1)});
//            triangles.push_back({terrainDims * Vector3(1, 1, 0.1), terrainDims * Vector3(0, 1, 0.1), terrainDims * Vector3(0, 0, 0.1)});
    if (boundariesTree != nullptr) {
        *boundariesTree = BVHTree();
        sumBVH = timeIt([&]() {
            boundariesTree->build(Triangle::vectorsToTriangles(triangles));
        });
    }
    sumMeshingBoundaries = timeIt([&]() {
        boundariesMesh.fromArray(flattenArray(triangles));
    });

    std::cout << "Boundaries have " << triangles.size() << " triangles." << std::endl;
    return {sumGeometry, sumBVH, sumMeshingBoundaries};
}

void ErosionInterface::ErosionParis2019SeaErosion()
{
    UnderwaterErosion erod(voxelGrid.get(), 0, 0, 0);
    float time = timeIt([&]() {
        erod.ParisSeaErosion();
    });
    std::cout << "Erosion using Paris 2019 (Sea erosion): " << showTime(time) << std::endl;
}

void ErosionInterface::ErosionParis2019InvasionPercolation()
{
    UnderwaterErosion erod(voxelGrid.get(), 0, 0, 0);
    float time = timeIt([&]() {
        erod.ParisInvasionPercolation();
    });
    std::cout << "Erosion using Paris 2019 (Invasion percolation): " << showTime(time) << std::endl;
}

void ErosionInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;
    if (this->displayTrajectories) {
        this->rocksPathSuccess.display(GL_LINES, 3.f);
        this->rocksPathFailure.display(GL_LINES, 3.f);
    }
    if (this->displayBoundaries) {
        this->boundariesMesh.displayWithOutlines(std::vector<float>{0.f, 0.f, 1.f, .4f}, GL_TRIANGLES);
    }
}

void ErosionInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto& parameters = action.at("parameters");
        Vector3 pos = json_to_vec3(parameters.at("position")) + Vector3::random(0.f, 20.f);
        Vector3 dir = json_to_vec3(parameters.at("direction")) + Vector3::random();
        float size = parameters.at("size").get<float>() + random_gen::generate(0.f, 3.f);
        int qtt = parameters.at("quantity").get<int>() + random_gen::generate(0.f, 100.f);
        float strength = parameters.at("strength").get<float>() + random_gen::generate(0.f, 1.f);
        float randomness = parameters.at("randomness").get<float>() + random_gen::generate(0.f, .1f);
        UnderwaterErosion erod(this->voxelGrid.get(), size, strength, qtt);
//        erod.Apply(pos, dir, randomness);
    }
}




void ErosionInterface::throwFromSky()
{
    this->throwFrom(SKY);
}
void ErosionInterface::throwFromCam()
{
    this->throwFrom(JUST_ABOVE_VOXELS);
}
void ErosionInterface::throwFromSide()
{
    this->throwFrom(EVERYWHERE);
}
void ErosionInterface::throwFrom(PARTICLE_INITIAL_LOCATION location)
{
    UnderwaterErosion erod = UnderwaterErosion(voxelGrid.get(), this->erosionSize, this->erosionStrength, this->erosionQtt);
    erod.heightmap = heightmap.get();

    GridF layersHeightmap(layerGrid->getDimensions().x, layerGrid->getDimensions().y);
    for (int x = 0; x < layersHeightmap.sizeX; x++)
        for (int y = 0; y < layersHeightmap.sizeY; y++)
            layersHeightmap.at(x, y) = layerGrid->getHeight(x, y) - 1;
    erod.implicitTerrain = this->implicitTerrain.get();
    erod.layerBasedGrid = layerGrid.get();


    std::vector<BSpline> lastRocksLaunched;
    this->rocksPathSuccess.clear();
    this->rocksPathFailure.clear();
    auto startingTime = std::chrono::system_clock::now();


    TerrainModel *terrain = nullptr;
    BVHTree boundariesTree;
    boundariesTree.useParallel = true;

    GridF densityField;
    GridV3 waterFlowfield = GridV3();
    GridV3 airFlowfield = GridV3();

    int totalPos = 0, totalErosions = 0;
    float sumParticleSimulationTime = 0.f, sumTerrainModifTime = 0.f, sumPreprocess = 0.f, sumGeometry = 0.f, sumBVH = 0.f, sumMeshingBoundaries = 0.f, sumTrajectories = 0.f;
    for (int iteration = 0; iteration < numberOfIterations; iteration++) {
        std::cout << "Iteration " << iteration + 1 << " / " << numberOfIterations << std::endl;
        erod.maxRockSize = this->erosionSize;
        erod.maxRockStrength = (this->erosionStrength * (erosionSize * erosionSize)/(5.f * 5.f)) / 10.f;
        erod.rockAmount = this->erosionQtt;

        int nbPos, nbErosions;
        float particleSimulationTime, terrainModifTime;
        sumPreprocess += timeIt([&]() {

            if (applyOn == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS) {
                terrain = voxelGrid.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP) {
                terrain = heightmap.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN) {
                terrain = implicitTerrain.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN) {
                terrain = layerGrid.get();
            }

    //        auto flowfieldFunction = this->computeFlowfieldFunction();
            if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
                if (this->waterFlowImagePath != "") {
                    waterFlowfield = -GridF::fromImageBW(this->waterFlowImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false).gradient();
                }
                if (this->airFlowImagePath != "") {
                    airFlowfield = -GridF::fromImageBW(this->airFlowImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false).gradient();
                }
            } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::BASIC) {

            } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION) {

            } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS) {

            }


            if (this->densityUsed == UnderwaterErosion::RANDOM_DENSITY) {
                densityField = GridF(voxelGrid->getDimensions());

                FastNoiseLite noise;
                noise.SetFractalType(FastNoiseLite::FractalType_FBm);
                float minDurability = .0f, maxDurability = 1.f;
                for (int x = 0; x < densityField.sizeX; x++) {
                    for (int y = 0; y < densityField.sizeY; y++) {
                        for (int z = 0; z < densityField.sizeZ; z++) {
                            float noiseVal = noise.GetNoise(float(x * 4.f), float(y * 4.f), float(z * 4.f));
                            float alt = std::pow(float(z) / densityField.sizeZ - .5f, 2.f);
                            densityField.at(x, y, z) = noiseVal + alt;

                        }
                    }
                }
                densityField = densityField.normalize();
                for (int x = 0; x < densityField.sizeX; x++) {
                    for (int y = 0; y < densityField.sizeY; y++) {
                        for (int z = 0; z < densityField.sizeZ; z++) {
                            Vector3 pos(x, y, z);
                            for (const auto& [obstaclePos, obstacleRadius] : randomObstacles)
                                if ((obstaclePos.xy() - pos.xy()).norm2() < obstacleRadius * obstacleRadius)
                                    densityField.at(x, y, z) = .1f;
                        }
                    }
                }
                densityField = 1.f - ((1.f - densityField) * (1.f - densityField));
                densityField = (densityField * (maxDurability - minDurability)) + minDurability;
//                voxelGrid->setVoxelValues(densityField - .5f);
                voxelGrid->setVoxelValues(((Matrix3<float>)voxelGrid->getVoxelValues().binarize()) - .9f);
//                std::cout << densityField.displayValues() << std::endl;
            }  else if (this->densityUsed == UnderwaterErosion::LAYERED_DENSITY) {
                densityField = GridF(voxelGrid->getDimensions());
    //            for (auto& v : densityField)
    //                v = std::min(random_gen::generate(5.f), 1.f);

                FastNoiseLite noise;
                noise.SetFractalType(FastNoiseLite::FractalType_FBm);
                float minDurability = 0.f, maxDurability = 1.f;
                for (int x = 0; x < densityField.sizeX; x++) {
                    for (int y = 0; y < densityField.sizeY; y++) {
                        for (int z = 0; z < densityField.sizeZ; z++) {
                            float noiseVal = noise.GetNoise((float) z * 5.f, x*y/100.f);
                            noiseVal = (noiseVal + 1.f) * .5f;
                            densityField.at(x, y, z) = (noiseVal * (maxDurability - minDurability)) + minDurability;
                            if (z > densityField.sizeZ * .8) {
                                float dz = (densityField.sizeZ - z) / (densityField.sizeZ * .2);
                                densityField.at(x, y, z) *= interpolation::wyvill(1.f - dz);
                            } else {
                                densityField.at(x, y, z) *= interpolation::wyvill(1.f - z / (densityField.sizeZ * .8));
                            }
                        }
                    }
                }
                densityField = 1.f - ((1.f - densityField) * (1.f - densityField));
                densityField = (densityField * (maxDurability - minDurability)) + minDurability;

            } else if (this->densityFieldImagePath != "") {
                densityField = GridF::fromImageBW(this->densityFieldImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false);
            }

            if (location == JUST_ABOVE_VOXELS)
                recomputeAboveVoxelRocksPositions(terrain);
            if (location == RIVER) {
                for(auto& [p, d]  : initialPositionsAndDirections[location][iteration % (initialPositionsAndDirections[location].size())]) {
                    p.z = 100; //terrain->getHeight(p.xy()) + 1.f;
                }
            }

            Vector3 terrainDims = terrain->getDimensions();
            std::vector<std::vector<Vector3>> triangles;
            Mesh m;
            Vector3 geomSize = Vector3::min(terrainDims, Vector3(100, 100, 50));
            sumGeometry += timeIt([&]() {
                if (applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN) {
                    m = Mesh::applyMarchingCubes(layerGrid->voxelize().meanSmooth(3, 3, 3, true).resize(geomSize));
                } else {
                    m = terrain->getGeometry(geomSize);
                }
            });
            triangles = m.getTriangles();

        //            triangles.push_back({terrainDims * Vector3(0, 0, 0.1), terrainDims * Vector3(1, 0, 0.1), terrainDims * Vector3(1, 1, 0.1)});
        //            triangles.push_back({terrainDims * Vector3(1, 1, 0.1), terrainDims * Vector3(0, 1, 0.1), terrainDims * Vector3(0, 0, 0.1)});
            boundariesTree = BVHTree();
            sumBVH += timeIt([&]() {
                boundariesTree.build(Triangle::vectorsToTriangles(triangles));
            });
            sumMeshingBoundaries += timeIt([&]() {
                boundariesMesh.fromArray(flattenArray(triangles));
            });

            std::cout << "Boundaries have " << triangles.size() << " triangles." << std::endl;

        });

        if (continuousRotation)
            airFlowfieldRotation += 45.f / 2.f;

        std::vector<std::vector<std::pair<float, Vector3>>> allErosions; // Useless
        std::tie(lastRocksLaunched, nbPos, nbErosions, allErosions) = erod.Apply(this->applyOn,
                                                                    terrain,
                                                                    boundariesTree,
                                                                    particleSimulationTime, terrainModifTime,
                                                                    Vector3(false),
                                                                    Vector3(false),
                                                                    this->rockRandomness,
                                                                    true,
                                                                    gravity,
                                                                    bouncingCoefficient,
                                                                    bounciness,
                                                                    minSpeed,
                                                                    maxSpeed,
                                                                    maxCapacityFactor,
                                                                    erosionFactor,
                                                                    depositFactor,
                                                                    matterDensity, // + .1f,
                                                                    materialImpact,
                                                                    airFlowfieldRotation,
                                                                    waterFlowfieldRotation,
                                                                    airForce,
                                                                    waterForce,
                                                                    dt,
                                                                    shearingStressConstantK,
                                                                    shearingRatePower,
                                                                    erosionPowerValue,
                                                                    criticalShearStress,
                                                                    initialPositionsAndDirections[location][iteration % (initialPositionsAndDirections[location].size())],
                                                                    this->flowfieldUsed,
                                                                    waterFlowfield,
                                                                    airFlowfield,
                                                                    densityUsed,
                                                                    densityField,
                                                                    initialCapacity,
                                                                    selectedSimulationType,
                                                                    wrapParticles,
                                                                    true,
                                                                    particleMaxCollisions
                                                                    );

        totalPos += nbPos;
        totalErosions += nbErosions;
        sumParticleSimulationTime += particleSimulationTime;
        sumTerrainModifTime += terrainModifTime;

        sumTrajectories += timeIt([&]() {
            std::vector<Vector3> asOneVector;
            for (size_t i = 0; i < lastRocksLaunched.size(); i++) {
                auto points = lastRocksLaunched[i].points; // lastRocksLaunched[i].getPath(std::min(100, int(lastRocksLaunched[i].points.size())));
                for (int j = 0; j < int(points.size()) - 1; j++) {
                    if ((points[j] - points[j+1]).norm2() < 20.f * 20.f) {
                        asOneVector.push_back(points[j]);
                        asOneVector.push_back(points[j + 1]);
                    }
                }
            }
            this->rocksPathSuccess.fromArray(asOneVector);
        });
        if (applyOn == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP) {
            voxelGrid->from2DGrid(*heightmap.get());
        } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN) {
            voxelGrid->fromImplicit(implicitTerrain.get());
            implicitTerrain->composables = {ImplicitPrimitive::fromHeightmap(voxelGrid->getVoxelValues())}; // Yeah, I can pass a 3D grid in this function
//            implicitTerrain->_cached = false;
        } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN) {
            layerGrid->reorderLayers();
            voxelGrid->fromLayerBased(*layerGrid, voxelGrid->getSizeZ());
        }
        Q_EMIT this->updated();
        this->computePredefinedRocksLocations();
    }
    auto endingTime = std::chrono::system_clock::now();
    float totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(endingTime - startingTime).count();
    std::cout << "Simulation time for " << numberOfIterations * erosionQtt << " particles: " << showTime(totalTime) << "  (" << totalPos/1000 << "k pos, " << totalErosions/1000 << "k erosions)" << std::endl;
    std::cout << "(" << showTime(sumParticleSimulationTime) << " for particle simulation, " << showTime(sumTerrainModifTime) << " for applying changes, " << showTime(totalTime - (sumParticleSimulationTime + sumTerrainModifTime)) << " for the display) => " << showTime(sumParticleSimulationTime + sumTerrainModifTime) << " for algo." << std::endl;
    std::cout << "More specific: " << showTime(sumGeometry) << " for geometry, " << showTime(sumBVH) << " for BVH, " << showTime(sumMeshingBoundaries) << " for triangles and " << showTime(sumTrajectories) << " for meshing trajectories" << std::endl;

//    UnderwaterErosion erod = UnderwaterErosion(voxelGrid.get(), this->erosionSize, this->erosionStrength, this->erosionQtt);

//    std::vector<BSpline> lastRocksLaunched;
//    for (int iteration = 0; iteration < numberOfIterations; iteration++) {
//        std::cout << "Iteration " << iteration + 1 << " / " << numberOfIterations << std::endl;
//        int nbPos, nbErosions;
//        float particleSimulationTime, terrainModifTime;
//        std::tie(lastRocksLaunched, nbPos, nbErosions) = erod.Apply(this->applyOn, particleSimulationTime, terrainModifTime, pos, dir, this->rockRandomness, false,
//                                                                          gravity,
//                                                                          bouncingCoefficient,
//                                                                          bounciness,
//                                                                          minSpeed,
//                                                                          maxSpeed,
//                                                                          maxCapacityFactor,
//                                                                          erosionFactor,
//                                                                          depositFactor,
//                                                                          matterDensity + 0.1f,
//                                                                          materialImpact,
//                                                                          airFlowfieldRotation,
//                                                                          waterFlowfieldRotation,
//                                                                          airForce,
//                                                                          waterForce,
//                                                                          dt,
//                                                                          shearingStressConstantK,
//                                                                          shearingRatePower,
//                                                                          erosionPowerValue,
//                                                                          criticalShearStress
//                                                                          );
//        std::vector<Vector3> asOneVector;
//        for (size_t i = 0; i < lastRocksLaunched.size(); i++) {
//            auto points = lastRocksLaunched[i].getPath(10);
//            for (size_t j = 0; j < points.size() - 1; j++) {
//                asOneVector.push_back(points[j]);
//                asOneVector.push_back(points[j + 1]);
//            }
//        }
//        this->rocksPathSuccess.fromArray(asOneVector);

//        this->addTerrainAction(nlohmann::json({
//                                                  {"position", vec3_to_json(pos) },
//                                                  {"direction", vec3_to_json(dir) },
//                                                  {"size", erosionSize},
//                                                  {"strength", erosionStrength},
//                                                  {"randomness", rockRandomness},
//                                                  {"quantity", erosionQtt}
//                                              }));

//        Q_EMIT this->updated();
//    }
}




std::map<std::string, float> generateRandomValuesFrom(std::map<std::string, std::tuple<float, float, float>> variables, std::vector<std::string> unlockedVariables = {}, std::vector<std::string> lockedVariables = {}, float t = -1)
{
    if (lockedVariables.empty() && unlockedVariables.size() > 0) {
        for (const auto& var : variables)
            if (!isIn(var.first, unlockedVariables)) lockedVariables.push_back(var.first);
    }
    else if (unlockedVariables.empty() && lockedVariables.size() > 0) {
        for (const auto& var : variables)
            if (!isIn(var.first, lockedVariables)) unlockedVariables.push_back(var.first);
    }
    else if (unlockedVariables.empty() && lockedVariables.empty()) {
        for (const auto& var : variables)
            unlockedVariables.push_back(var.first);
    }
    std::map<std::string, float> results;
    for (const auto& [var, val] : variables) {
        auto [mini, maxi, defaultVal] = val;
        if (isIn(var, unlockedVariables)) {
            if (t == -1) {
                results[var] = random_gen::generate(mini, maxi); // Forget about the step
            } else {
                results[var] = interpolation::inv_linear(t, mini, maxi);
            }
        } else {
            results[var] = defaultVal;
        }
    }
    return results;
}


void ErosionInterface::testManyManyErosionParameters()
{
    /*
    srand(1);

    UnderwaterErosion::EROSION_APPLIED terrainType = this->applyOn; //UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP; // 0 = voxels, 1 = heightmap, 2 = implicit, ...
    if (terrainType == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS)
        viewer->setMapMode(MapMode::VOXEL_MODE);
    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP)
        viewer->setMapMode(MapMode::GRID_MODE);
    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN)
        viewer->setMapMode(MapMode::LAYER_MODE);
    else if (terrainType == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN)
        viewer->setMapMode(MapMode::IMPLICIT_MODE);
    srand(43);

    TerrainGenerationInterface* terrainInterface = static_cast<TerrainGenerationInterface*>(this->viewer->interfaces["terrainGenerationInterface"].get());

    VoxelGrid initialVoxelGrid = *voxelGrid;
    Heightmap initialHeightmap = *heightmap;
    LayerBasedGrid initialLayerGrid = *layerGrid;

//    this->viewer->restoreFromFile("experiments_state2.xml");
//    this->viewer->setStateFileName("experiments_state2.xml");
//    this->viewer->saveStateToFile();

    Vector3 cameraPosition = Vector3(this->viewer->camera()->position());
    Vector3 cameraDirection = Vector3(this->viewer->camera()->viewDirection());

    std::map<std::string, std::tuple<float, float, float>> varyingVariables =
    {
        // Name,            min, max, step
        {"rockRandomness",              {0.f, 1.f, rockRandomness} },
        {"gravity",                     {.5f, 1.5f, gravity} },
        {"bouncingCoefficient",         {.2f, 1.f, bouncingCoefficient} },
        {"bounciness",                  {.2f, 1.f, bounciness} },
        {"minSpeed",                    {0.f, 0.f, minSpeed} }, // Useless
        {"maxSpeed",                    {0.f, 0.f, maxSpeed} }, // Useless
        {"maxCapacityFactor",           {1.f, 10.f, maxCapacityFactor} },
        {"erosionFactor",               {0.f, 3.f, erosionFactor} },
        {"depositFactor",               {0.f, 3.f, depositFactor} },
        {"matterDensity",               {1.f, 2000.f, matterDensity} },
        {"materialImpact",              {0.f, 1.f, materialImpact} },
        {"airFlowfieldRotation",        {0.f, 270.f, airFlowfieldRotation} },
        {"waterFlowfieldRotation",      {0.f, 270.f, waterFlowfieldRotation} },
        {"airForce",                    {0.f, 1.0f, airForce} },
        {"waterForce",                  {0.f, 1.0f, waterForce} },
        {"particleSize",                {2.f, 16.f, erosionSize} },
        {"strength",                    {0.f, .5f, erosionStrength} },
        {"nbParticles",                 {1.f, 50.f, erosionQtt} },
        {"dt",                          {.1f, .1f, dt} },
        {"shearingStressConstantK",     {.8f, 1.2f, shearingStressConstantK} },
        {"shearingRatePower",           {.4f, .6f, shearingRatePower} },
        {"erosionPowerValue",           {.9f, 1.1f, erosionPowerValue} }, // From Wojtan : =1
        {"criticalShearStress",         {.5f, 5.f, criticalShearStress} },
        {"camPos.x",                    {cameraPosition.x, cameraPosition.x, cameraPosition.x} },
        {"camPos.y",                    {cameraPosition.y, cameraPosition.y, cameraPosition.y} },
        {"camPos.z",                    {cameraPosition.z, cameraPosition.z, cameraPosition.z} },
        {"camDir.x",                    {cameraDirection.x, cameraDirection.x, cameraDirection.x} },
        {"camDir.y",                    {cameraDirection.y, cameraDirection.y, cameraDirection.y} },
        {"camDir.z",                    {cameraDirection.z, cameraDirection.z, cameraDirection.z} },
        {"waterLevel",                  {0.f, 1.f, terrainInterface->waterLevel} },
    };

    std::vector<std::string> params;
    for (const auto& var : varyingVariables)
        params.push_back(var.first);

    std::string mainFolder = "";
#ifdef linux
    mainFolder = "/data/erosionsTests/";
#else
    mainFolder = "erosionsTests/";
#endif
    makedir(mainFolder);
    // Create the general CSV
    // Check if already exists
    std::string CSVname = mainFolder + "allData.csv";
    bool exists = checkPathExists(CSVname);
    std::fstream mainCSVfile;
    std::vector<std::map<std::string, float>> alreadyTestedParameters;
    if (exists) {
        mainCSVfile.open(CSVname, std::ios_base::in);
        mainCSVfile.seekg(0);
        mainCSVfile.seekp(0);
        std::string header, lineContent;
        std::vector<std::string> headerValues, sLineValues;
        std::getline(mainCSVfile, header); // get header
        headerValues = split(header, ";");
        while (std::getline(mainCSVfile, lineContent)) {
            sLineValues = split(lineContent, ";");
            std::map<std::string, float> testedParameters;
            for (size_t i = 0; i < sLineValues.size(); i++) {
                float val;
                try {
                    val = std::stof(sLineValues[i]);
                    testedParameters[headerValues[i]] = val;
                } catch (std::exception e) {
                    try {
                        val = std::stof(replaceInString(sLineValues[i], ",", "."));
                        testedParameters[headerValues[i]] = val;
                    } catch (std::exception e2) {

                    }
                 }
            }
            alreadyTestedParameters.push_back(testedParameters);
        }
        mainCSVfile.close();
        mainCSVfile.open(CSVname, std::ios_base::out | std::ios_base::app);

    } else {
        mainCSVfile.open(CSVname, std::ios_base::out | std::ios_base::app);
        for (const auto& param : params)
            mainCSVfile << param << ";";
        mainCSVfile << "folder_name" << std::endl;
    }

    std::vector<BSpline> lastRocksLaunched;
    this->rocksPathSuccess.clear();
    this->rocksPathFailure.clear();

    std::vector<std::vector<std::string>> testedVariables = {
//        { "particleSize" },
//        { "strength" },
//        { "erosionFactor" },
//        { "criticalShearStress" },
//        { "shearingStressConstantK" },
//        { "nbParticles" },
//        { "maxCapacityFactor" },
        { "matterDensity" },
        { "depositFactor" },
//        { "waterLevel" },
        { "waterForce" },
        { "matterDensity" },
    };

    int stopAfterStep = 11; // Set to -1 for infinite tries
    for (int iCombination = 0; iCombination < int(testedVariables.size()); iCombination++) {
        std::cout << "Testing " << join(testedVariables[iCombination], " and ") << std::endl;
        for (int i = 0; i < stopAfterStep || stopAfterStep == -1; i++) {
            // Initialize
            float t = (testedVariables[iCombination].size() == 1 ? float(i)/float(stopAfterStep - 1) : -1);
            auto variables = generateRandomValuesFrom(varyingVariables, testedVariables[iCombination], {}, t);
            bool tested = false;
            for (auto& x : alreadyTestedParameters) {
                tested = true;
                for (auto& [paramName, paramVal] : variables) {
                    if (!startsWith(paramName, "cam") && replaceInString(std::to_string(variables[paramName]), ",", ".") != replaceInString(std::to_string(x[paramName]), ",", ".")) {
                        tested = false;
                    }
                }
                if (tested)
                    break;
            }
            if (tested) {
                continue;
            }
            alreadyTestedParameters.push_back(variables);

            *voxelGrid = initialVoxelGrid;
            *heightmap = initialHeightmap;
            *layerGrid = initialLayerGrid;

            terrainInterface->setWaterLevel(variables["waterLevel"]);

            time_t now = std::time(0);
            tm *gmtm = std::gmtime(&now);
            char s_time[80];
            std::strftime(s_time, 80, "%Y-%m-%d__%H-%M-%S", gmtm);
            std::string subFolderName = join(testedVariables[iCombination], "-") + std::to_string(i) + "__" + std::string(s_time);
            std::string folderName = mainFolder + subFolderName + "/";

            makedir(folderName + "screen/");
            makedir(folderName + "voxels/");
            makedir(folderName + "meshes/");
            makedir(folderName + "particles/");
            makedir(folderName + "heightmap/");

            this->viewer->startRecording(folderName + "screen/");

            std::cout << "Test " << i+1 << " for " << join(testedVariables[iCombination], " and ") << " : \n";
            for (const auto& var : variables)
                std::cout << "\t- \"" << var.first << "\" = " << var.second << "\n";
            std::cout << std::flush;

            GridV3 waterFlowfield = (this->waterFlowImagePath != "" ? -GridF::fromImageBW(this->waterFlowImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false).gradient() : GridV3());
            GridV3 airFlowfield = (this->airFlowImagePath != "" ? -GridF::fromImageBW(this->airFlowImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false).gradient() : GridV3());


            this->erosionProcess.rockAmount = variables["nbParticles"];
            this->erosionProcess.maxRockSize = variables["particleSize"];
            this->erosionProcess.maxRockStrength = variables["strength"];
//            UnderwaterErosion erod = UnderwaterErosion(voxelGrid.get(), variables["particleSize"], variables["strength"], variables["nbParticles"]);
            heightmap->heights.raiseErrorOnBadCoord = false;
//            erod.heightmap = heightmap.get();
            auto startingTime = std::chrono::system_clock::now();
            int totalPos = 0, totalErosions = 0;
            float sumParticleSimulationTime = 0.f, sumTerrainModifTime = 0.f;
            float totalMeshingTime = 0.f;
            int nbIterations = 100;
            for (int iteration = 0; iteration < nbIterations; iteration++) {
                int nbPos, nbErosions;
                float particleSimulationTime, terrainModifTime;


                if (iteration % 10 == 0 and false) {
                    auto startMesh = std::chrono::system_clock::now();
                    if (terrainType == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS) {
                        voxelGrid->saveMap(folderName + "voxels/terrain_" + std::to_string(iteration) + ".data");
                        std::ofstream geomFile; //(folderName + "meshes/terrain_" + std::to_string(iteration) + ".stl");
//                        geomFile << voxelGrid->getGeometry().toSTL();
//                        geomFile.close();
                        geomFile.open(folderName + "meshes/terrain_" + std::to_string(iteration) + ".obj");
                        geomFile << heightmap->getGeometry().toOBJ();
                        geomFile.close();
                    } else if (terrainType == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP) {
                        heightmap->saveMap(folderName + "heightmap/terrain_" + std::to_string(iteration) + ".png");
                        std::ofstream geomFile; // (folderName + "meshes/terrain_" + std::to_string(iteration) + ".stl");
//                        geomFile << heightmap->getGeometry().toSTL();
//                        geomFile.close();
                        geomFile.open(folderName + "meshes/terrain_" + std::to_string(iteration) + ".obj");
                        geomFile << heightmap->getGeometry().toOBJ();
                        geomFile.close();
                    }
                    auto endMesh = std::chrono::system_clock::now();
//                    std::ofstream particleFile(folderName + "particles/flow_" + std::to_string(iteration) + ".txt");
//                    for (const auto& path : lastRocksLaunched) {
//                        for (const auto& p : path) {
//                            particleFile << p.x << " " << p.y << " " << p.z << " ";
//                        }
//                        particleFile << std::endl;
//                    }
//                    particleFile.close();
                    totalMeshingTime += std::chrono::duration_cast<std::chrono::milliseconds>(endMesh - startMesh).count();
                }

//                std::tie(lastRocksLaunched, nbPos, nbErosions) = erod.Apply(terrainType, particleSimulationTime, terrainModifTime,
                std::tie(lastRocksLaunched, nbPos, nbErosions) = this->erosionProcess.Apply(
                            terrainType,
                            particleSimulationTime,
                            terrainModifTime,
                            Vector3(false),
                            Vector3(false),
                            variables["rockRandomness"],
                            true,
                            variables["gravity"],
                            variables["bouncingCoefficient"],
                            variables["bounciness"],
                            variables["minSpeed"],
                            variables["maxSpeed"],
                            variables["maxCapacityFactor"],
                            variables["erosionFactor"],
                            variables["depositFactor"],
                            variables["matterDensity"],
                            variables["materialImpact"],
                            variables["airFlowfieldRotation"],
                            variables["waterFlowfieldRotation"],
                            variables["airForce"],
                            variables["waterForce"],
                            variables["dt"],
                            variables["shearingStressConstantK"],
                            variables["shearingRatePower"],
                            variables["erosionPowerValue"],
                            variables["criticalShearStress"],
                            this->initialPositionsAndDirections[SKY][iteration % (initialPositionsAndDirections[SKY].size())],
                            this->flowfieldUsed,
                            waterFlowfield,
                            airFlowfield
                            );

                totalPos += nbPos;
                totalErosions += nbErosions;
                sumParticleSimulationTime += particleSimulationTime;
                sumTerrainModifTime += terrainModifTime;
//                *heightmap = *erod.heightmap;
                Q_EMIT this->updated();
                qApp->processEvents();
            }
            auto endingTime = std::chrono::system_clock::now();
            float totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(endingTime - startingTime).count();
            std::cout << "Simulation time for " << numberOfIterations * erosionQtt << " particles: " << totalTime << "ms  (" << totalPos/1000 << "k pos, " << totalErosions/1000 << "k erosions)" << std::endl;
            std::cout << "(" << sumParticleSimulationTime << "ms for particle simulation, " << sumTerrainModifTime << "ms for applying changes, " << totalMeshingTime << "ms for meshing and " << totalTime - (sumParticleSimulationTime + sumTerrainModifTime + totalMeshingTime) << "ms for the rest)" << std::endl;

            std::string parametersToCSV = "";
            for (const auto& param : params)
                parametersToCSV += std::to_string(variables[param]) + ";";
            parametersToCSV += subFolderName;
            mainCSVfile << parametersToCSV << std::endl;

            this->viewer->stopRecording();
        }
    }

    mainCSVfile.close();
    std::cout << "Unbelivable but true, it's finished!" << std::endl;*/
}

void ErosionInterface::afterTerrainUpdated()
{

    /*
    if (!this->currentlyModifyingTerrain) {
        GridF layersHeightmap(layerGrid->getDimensions().x, layerGrid->getDimensions().y);
        for (int x = 0; x < layersHeightmap.sizeX; x++)
            for (int y = 0; y < layersHeightmap.sizeY; y++)
                layersHeightmap.at(x, y) = layerGrid->getHeight(x, y) - 1;
        auto implicitHeightmap = ImplicitPrimitive::fromHeightmap(layersHeightmap, "");
        implicitHeightmap->material = TerrainTypes::DIRT;
        implicitHeightmap->position = Vector3();
        this->erosionProcess.implicitTerrain = implicitHeightmap;
    }*/
}

void ErosionInterface::browseWaterFlowFromFile()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Open an image for water flow"), QString::fromStdString("saved_maps/heightmaps/gradients/"));
    std::string filename = q_filename.toStdString();
    if (!q_filename.isEmpty()) {
        this->waterFlowImagePath = q_filename.toStdString();
    }
}

void ErosionInterface::browseAirFlowFromFile()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Open an image for air flow"), QString::fromStdString("saved_maps/heightmaps/gradients/"));
    std::string filename = q_filename.toStdString();
    if (!q_filename.isEmpty()) {
        this->airFlowImagePath = q_filename.toStdString();
    }
}

void ErosionInterface::browseDensityFieldFromFile()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Open an image for density field"), QString::fromStdString("saved_maps/heightmaps/gradients/"));
    std::string filename = q_filename.toStdString();
    if (!q_filename.isEmpty()) {
        this->densityFieldImagePath = q_filename.toStdString();
    }
}

std::function<Vector3 (Vector3)> ErosionInterface::computeFlowfieldFunction()
{
    std::function<Vector3(Vector3)> flowfieldFunction;
    if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE) {
        auto voxelsEnvironmentalDensities = voxelGrid->getEnvironmentalDensities();
        GridV3 waterFlow = -GridF::fromImageBW(this->waterFlowImagePath).resize(Vector3(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f)).flip(true, false, false).gradient();
        for (auto& v : waterFlow)
            v.normalize();
        GridV3 airFlow = -GridF::fromImageBW(this->airFlowImagePath).resize(Vector3(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f)).flip(true, false, false).gradient();
        for (auto& v : airFlow)
            v.normalize();

        flowfieldFunction = [&](const Vector3& pos) {
            return (voxelsEnvironmentalDensities.at(pos) < 500 ? airFlow.at(pos.xy()) : waterFlow.at(pos.xy()));
        };
    } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION) {
        auto fluidSim = voxelGrid->getFlowfield();

        flowfieldFunction = [&](const Vector3& pos) {
            return fluidSim.at(pos);
        };
    } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS) {
        auto fluidSim = EnvObject::flowfield.resize(voxelGrid->getDimensions());

        flowfieldFunction = [&](const Vector3& pos) {
            return fluidSim.at(pos);
        };
    } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::BASIC) {
        flowfieldFunction = nullptr;
    }
    return flowfieldFunction;
}










QLayout *ErosionInterface::createGUI()
{
    this->erosionLayout = new QHBoxLayout();

    FancySlider* rockSizeSlider = new FancySlider(Qt::Horizontal, 0.f, 100.f);
    FancySlider* rockStrengthSlider = new FancySlider(Qt::Horizontal, 0.f, .5f, .01f);
    FancySlider* rockQttSlider = new FancySlider(Qt::Horizontal, 1.f, maxParticles);
    FancySlider* rockRandomnessSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* gravitySlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* bouncingCoefficientSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* bouncinessSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* minSpeedSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* maxSpeedSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* maxCapacityFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 10.f, .01f);
    FancySlider* erosionFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 5.f, .01f);
    FancySlider* depositFactorSlider = new FancySlider(Qt::Horizontal, 0.f, 5.f, .01f);
    FancySlider* matterDensitySlider = new FancySlider(Qt::Horizontal, 0.f, 2000.f, .25f);
    FancySlider* materialImpactSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);

    FancySlider* airFlowfieldRotationSlider = new FancySlider(Qt::Horizontal, 0.f, 360.f, 45.f);
    FancySlider* waterFlowfieldRotationSlider = new FancySlider(Qt::Horizontal, 0.f, 360.f, 45.f);
    FancySlider* airForceSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* waterForceSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);

    FancySlider* dtSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* shearingStressConstantKSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* shearingRatePowerSlider = new FancySlider(Qt::Horizontal, 0.f, 1.f, .01f);
    FancySlider* erosionPowerValueSlider = new FancySlider(Qt::Horizontal, 0.f, 2.f, .01f);
    FancySlider* criticalShearStressSlider = new FancySlider(Qt::Horizontal, 0.f, 5.f, .1f);

    FancySlider* iterationSlider = new FancySlider(Qt::Orientation::Horizontal, 1.f, 500.f);

    FancySlider* initialCapacitySlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, .01f);

    QPushButton* confirmFromRandom = new QPushButton("Random");
    QPushButton* confirmFromSurface = new QPushButton("From surface");
    QPushButton* confirmFromSkyButton = new QPushButton("Rain");
    QPushButton* confirmFromRiverButton = new QPushButton("River");
    QPushButton* confirmFromRiver2Button = new QPushButton("River2");
    QPushButton* confirmFromSideButton = new QPushButton("SideX");
    QPushButton* confirmFromBigSideButton = new QPushButton("AllSideX");
    QPushButton* confirmFromVolcanoButton = new QPushButton("Volcano");
    QPushButton* confirmFromVolcano2Button = new QPushButton("Volcano2");
    QPushButton* confirmFromVolcano3Button = new QPushButton("Volcano3");

    QCheckBox* displayTrajectoriesButton = new QCheckBox("Display path");
    QCheckBox* displayBoundariesButton = new QCheckBox("Display walls");

    QRadioButton* applyOnVoxels = new QRadioButton("on voxels");
    QRadioButton* applyOnHeightmap = new QRadioButton("on heightmap");
    QRadioButton* applyOnImplicit = new QRadioButton("on implicit");
    QRadioButton* applyOnLayers = new QRadioButton("on layers");

    QRadioButton* useBasicFlowfield = new QRadioButton("Basic flowfield");
    QRadioButton* useImageFlowfield = new QRadioButton("Flowfield from image");
    QRadioButton* useSimulatedFlowfield = new QRadioButton("Simulation");
    QRadioButton* useEnvObjFlowfield = new QRadioButton("Env objects");
    QLabel* labWater = new QLabel;
    QLabel* labAir = new QLabel;
    QPushButton* browseWaterFlow = new QPushButton("...");
    QPushButton* browseAirFlow = new QPushButton("...");

    QRadioButton* useRandomDensity = new QRadioButton("Random density");
    QRadioButton* useLayeredDensity = new QRadioButton("Layered");
    QRadioButton* useNativeDensity = new QRadioButton("Native density");
    QRadioButton* useImageDensity = new QRadioButton("Density from image");
    QLabel* labDensityField = new QLabel;
    QPushButton* densityFieldFileChooser = new QPushButton("...");

    QPushButton* lotsOfTestsButton = new QPushButton("Do tests");

    QRadioButton* waterDensity = new QRadioButton("Hydraulic");
    QRadioButton* airDensity = new QRadioButton("Aeolian");

    QRadioButton* particleSizeSmall = new QRadioButton("Small");
    QRadioButton* particleSizeMedium = new QRadioButton("Medium");
    QRadioButton* particleSizeBig = new QRadioButton("Big");

    QCheckBox* continuousRotationButton = new QCheckBox("Continuous rotation");
    QCheckBox* wrapPositionsButton = new QCheckBox("Wrap position");
//    QRadioButton* continuousRotationFalse = new QRadioButton("Without rotation");
//    QRadioButton* continuousRotationTrue = new QRadioButton("With rotation");

    QComboBox* simulationTypeButton = new QComboBox;
    std::vector<FluidSimType> possibleSimTypes = {LBM, FLIP, SPH, STABLE, WARP};
    for (size_t i = 0; i < possibleSimTypes.size(); i++) {
        simulationTypeButton->addItem(QString::fromStdString(stringFromFluidSimType(possibleSimTypes[i])));
    }
//    simulationTypeButton->addItems({"MLB", "FLIP", "SPH", "Stable", "Warp"});
//    simulationTypeButton->addItem("MLB", QVariant(FluidSimType::LBM));
//    simulationTypeButton->addItem("FLIP", QVariant(FluidSimType::FLIP));
//    simulationTypeButton->addItem("SPH", QVariant(FluidSimType::SPH));
//    simulationTypeButton->addItem("Stable", QVariant(FluidSimType::STABLE));
//    simulationTypeButton->addItem("Warp", QVariant(FluidSimType::WARP));

    QRadioButton* noLimitCollisionButton = new QRadioButton("No limit");
    QRadioButton* singleCollisionButton = new QRadioButton("1 collision");
    QRadioButton* twoCollisionButton = new QRadioButton("2 collisions");
    QRadioButton* tenCollisionButton = new QRadioButton("10 collisions");

    erosionLayout->addWidget(createVerticalGroup({
                                                     createMultipleSliderGroup({
//                                                           {"Taille", rockSizeSlider},
                                                           {"Strength", rockStrengthSlider},
                                                           {"Quantity", rockQttSlider},
//                                                           {"gravity", gravitySlider},
                                                           {"bouncing Coefficient", bouncingCoefficientSlider},
                                                           {"bounciness", bouncinessSlider},
//                                                           {"minSpeed", minSpeedSlider},
//                                                           {"maxSpeed", maxSpeedSlider},
                                                           {"max Capacity Factor", maxCapacityFactorSlider},
                                                           {"erosion Factor", erosionFactorSlider},
                                                           {"deposit Factor", depositFactorSlider},
//                                                           {"matter Density", matterDensitySlider},
//                                                           {"material Impact", materialImpactSlider},
                                                           {"air Rotation", airFlowfieldRotationSlider},
                                                           {"water Rotation", waterFlowfieldRotationSlider},
                                                           {"air Force", airForceSlider},
                                                           {"water Force", waterForceSlider},
                                                           {"nb iterations", iterationSlider},
//                                                           {"dt", dtSlider},
//                                                           {"ShearConstantK", shearingStressConstantKSlider},
//                                                           {"ShearRatePower", shearingRatePowerSlider},
//                                                           {"ErosionPower", erosionPowerValueSlider},
//                                                           {"Critical shear stress", criticalShearStressSlider},
                                                           {"Initial capacity", initialCapacitySlider},
                                                       }),
                                                     createHorizontalGroup({
                                                         airDensity, waterDensity
                                                     }),
                                                     createHorizontalGroup({
                                                         particleSizeSmall, particleSizeMedium, particleSizeBig
                                                     }),
                                                     createHorizontalGroup({
                                                         continuousRotationButton, wrapPositionsButton
                                                     }),
                                                     createHorizontalGroup({
                                                         noLimitCollisionButton, singleCollisionButton, twoCollisionButton, tenCollisionButton
                                                     })
                                                 }));
    erosionLayout->addWidget(createVerticalGroup({
                                                     confirmFromSurface,
                                                     confirmFromSkyButton,
                                                     confirmFromRiverButton,
                                                     confirmFromRiver2Button,
                                                     confirmFromRandom,
                                                     confirmFromSideButton,
                                                     confirmFromBigSideButton,
                                                     confirmFromVolcanoButton,
                                                     confirmFromVolcano2Button,
                                                     confirmFromVolcano3Button,

                                                     displayTrajectoriesButton,
                                                     displayBoundariesButton,
                                                     createVerticalGroup({
                                                         applyOnHeightmap,
                                                         applyOnVoxels,
                                                         applyOnImplicit,
                                                         applyOnLayers
                                                     }),
                                                     createVerticalGroup({
                                                         useBasicFlowfield,
                                                         useImageFlowfield,
                                                         useEnvObjFlowfield,
                                                         labWater, browseWaterFlow,
                                                         labAir, browseAirFlow,
                                                         /*createHorizontalGroup({*/useSimulatedFlowfield, simulationTypeButton/*})*/
                                                     }),
//                                                     lotsOfTestsButton,
                                                     createVerticalGroup({
                                                         useImageDensity,
                                                         useLayeredDensity,
                                                         labDensityField,
                                                         densityFieldFileChooser,
                                                         useNativeDensity,
                                                         useRandomDensity
                                                     })
                                                 }));

    QObject::connect(rockSizeSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionSize = newVal; });
    QObject::connect(rockStrengthSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionStrength = newVal; });
    QObject::connect(rockQttSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionQtt = newVal; });
    QObject::connect(rockRandomnessSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->rockRandomness = newVal; });

    QObject::connect(gravitySlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->gravity = newVal; });
    QObject::connect(bouncingCoefficientSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->bouncingCoefficient = newVal; });
    QObject::connect(bouncinessSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->bounciness = newVal; });
    QObject::connect(minSpeedSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->minSpeed = newVal; });
    QObject::connect(maxSpeedSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->maxSpeed = newVal; });
    QObject::connect(maxCapacityFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->maxCapacityFactor = newVal; });
    QObject::connect(erosionFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionFactor = newVal; });
    QObject::connect(depositFactorSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->depositFactor = newVal; });
    QObject::connect(matterDensitySlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->matterDensity = newVal; });
    QObject::connect(materialImpactSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->materialImpact = newVal; });
    QObject::connect(airFlowfieldRotationSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->airFlowfieldRotation = newVal; });
    QObject::connect(waterFlowfieldRotationSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->waterFlowfieldRotation = newVal; });
    QObject::connect(airForceSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->airForce = newVal; });
    QObject::connect(waterForceSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->waterForce = newVal; });
    QObject::connect(iterationSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->numberOfIterations = (int) newVal; });

    QObject::connect(dtSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->dt = newVal; });
    QObject::connect(shearingStressConstantKSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->shearingStressConstantK = newVal; });
    QObject::connect(shearingRatePowerSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->shearingRatePower = newVal; });
    QObject::connect(erosionPowerValueSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->erosionPowerValue = newVal; });
    QObject::connect(criticalShearStressSlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->criticalShearStress = newVal; });
    QObject::connect(initialCapacitySlider, &FancySlider::floatValueChanged, this, [&](float newVal) { this->initialCapacity = newVal; });

    QObject::connect(confirmFromRandom, &QPushButton::pressed, this, [&](){ this->throwFrom(EVERYWHERE); });
    QObject::connect(confirmFromSurface, &QPushButton::pressed, this, [&](){ this->throwFrom(JUST_ABOVE_VOXELS); });
    QObject::connect(confirmFromSkyButton, &QPushButton::pressed, this, [&](){ this->throwFrom(SKY); });
    QObject::connect(confirmFromRiverButton, &QPushButton::pressed, this, [&](){ this->throwFrom(RIVER); });
    QObject::connect(confirmFromRiver2Button, &QPushButton::pressed, this, [&](){ this->throwFrom(RIVER2); });
    QObject::connect(confirmFromSideButton, &QPushButton::pressed, this, [&](){ this->throwFrom(FROM_X); });
    QObject::connect(confirmFromBigSideButton, &QPushButton::pressed, this, [&](){ this->throwFrom(FROM_BIG_X); });
    QObject::connect(confirmFromVolcanoButton, &QPushButton::pressed, this, [&](){ this->throwFrom(VOLCANO); });
    QObject::connect(confirmFromVolcano2Button, &QPushButton::pressed, this, [&](){ this->throwFrom(VOLCANO2); });
    QObject::connect(confirmFromVolcano3Button, &QPushButton::pressed, this, [&](){ this->throwFrom(VOLCANO3); });

    QObject::connect(displayTrajectoriesButton, &QCheckBox::toggled, this, [&](bool checked) { this->displayTrajectories = checked; });
    QObject::connect(displayBoundariesButton, &QCheckBox::toggled, this, [&](bool checked) { this->displayBoundaries = checked; });

    QObject::connect(applyOnVoxels, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS; });
    QObject::connect(applyOnHeightmap, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP; });
    QObject::connect(applyOnImplicit, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN; });
    QObject::connect(applyOnLayers, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN; });

    QObject::connect(useBasicFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::BASIC; });
    QObject::connect(useImageFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE; });
    QObject::connect(useSimulatedFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION; });
    QObject::connect(useEnvObjFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS; });

    QObject::connect(browseWaterFlow, &QPushButton::pressed, this, [=]() { this->browseWaterFlowFromFile(); labWater->setText("Water: " + QString::fromStdString(getFilename(this->waterFlowImagePath)));});
    QObject::connect(browseAirFlow, &QPushButton::pressed, this, [=]() { this->browseAirFlowFromFile(); labAir->setText("Air: " + QString::fromStdString(getFilename(this->airFlowImagePath)));} );

    QObject::connect(useRandomDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::RANDOM_DENSITY; });
    QObject::connect(useLayeredDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::LAYERED_DENSITY; });
    QObject::connect(useNativeDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::NATIVE; });
    QObject::connect(useImageDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::DENSITY_IMAGE; });
    QObject::connect(densityFieldFileChooser, &QPushButton::pressed, this, [=]() { this->browseDensityFieldFromFile(); labAir->setText("File: " + QString::fromStdString(getFilename(this->airFlowImagePath)));} );

    QObject::connect(lotsOfTestsButton, &QPushButton::pressed, this, &ErosionInterface::testManyManyErosionParameters);

    QObject::connect(airDensity, &QRadioButton::pressed, this, [=]() { this->matterDensity = 500.f; });
    QObject::connect(waterDensity, &QRadioButton::pressed, this, [=]() { this->matterDensity = 1662.f; });

    QObject::connect(continuousRotationButton, &QCheckBox::toggled, this, [=](bool checked) { this->continuousRotation = checked; });
    QObject::connect(wrapPositionsButton, &QCheckBox::toggled, this, [=](bool checked) { this->wrapParticles = checked; });

    QObject::connect(particleSizeSmall, &QRadioButton::pressed, this, [=]() { this->erosionSize = 4.f; });
    QObject::connect(particleSizeMedium, &QRadioButton::pressed, this, [=]() { this->erosionSize = 8.f; });
    QObject::connect(particleSizeBig, &QRadioButton::pressed, this, [=]() { this->erosionSize = 12.f; });

    QObject::connect(simulationTypeButton, &QComboBox::currentTextChanged, this, [=](QString text) { this->selectedSimulationType = FluidSimTypeFromString(text.toStdString()); });
//    QObject::connect(simulationTypeButton, &QComboBox::currentIndexChanged, this, [=](int index) { this->selectedSimulationType = simulationTypeButton->item; });

    QObject::connect(noLimitCollisionButton, &QRadioButton::toggled, this, [&](bool checked) { this->particleMaxCollisions = -1; });
    QObject::connect(singleCollisionButton, &QRadioButton::toggled, this, [&](bool checked) { this->particleMaxCollisions = 1; });
    QObject::connect(twoCollisionButton, &QRadioButton::toggled, this, [&](bool checked) { this->particleMaxCollisions = 2; });
    QObject::connect(tenCollisionButton, &QRadioButton::toggled, this, [&](bool checked) { this->particleMaxCollisions = 10; });

    rockSizeSlider->setfValue(this->erosionSize);
    rockStrengthSlider->setfValue(this->erosionStrength);
    rockQttSlider->setfValue(this->erosionQtt);
    rockRandomnessSlider->setfValue(this->rockRandomness);
    gravitySlider->setfValue(this->gravity);
    bouncingCoefficientSlider->setfValue(this->bouncingCoefficient);
    bouncinessSlider->setfValue(this->bounciness);
    minSpeedSlider->setfValue(this->minSpeed);
    maxSpeedSlider->setfValue(this->maxSpeed);
    maxCapacityFactorSlider->setfValue(this->maxCapacityFactor);
    erosionFactorSlider->setfValue(this->erosionFactor);
    depositFactorSlider->setfValue(this->depositFactor);
    matterDensitySlider->setfValue(this->matterDensity);
    materialImpactSlider->setfValue(this->materialImpact);
    airFlowfieldRotationSlider->setfValue(this->airFlowfieldRotation);
    waterFlowfieldRotationSlider->setfValue(this->waterFlowfieldRotation);
    airForceSlider->setfValue(this->airForce);
    waterForceSlider->setfValue(this->waterForce);
    iterationSlider->setfValue(this->numberOfIterations);
    dtSlider->setfValue(this->dt);
    shearingStressConstantKSlider->setfValue(this->shearingStressConstantK);
    shearingRatePowerSlider->setfValue(this->shearingRatePower);
    erosionPowerValueSlider->setfValue(this->erosionPowerValue);
    criticalShearStressSlider->setfValue(this->criticalShearStress);
    initialCapacitySlider->setfValue(this->initialCapacity);

    continuousRotationButton->setChecked(this->continuousRotation);
    wrapPositionsButton->setChecked(this->wrapParticles);

    displayTrajectoriesButton->setChecked(this->displayTrajectories);

    noLimitCollisionButton->setChecked(this->particleMaxCollisions == -1);
    singleCollisionButton->setChecked(this->particleMaxCollisions == 1);
    twoCollisionButton->setChecked(this->particleMaxCollisions == 2);
    tenCollisionButton->setChecked(this->particleMaxCollisions == 10);

    applyOnVoxels->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS);
    applyOnHeightmap->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP);
    applyOnImplicit->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN);
    applyOnLayers->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN);

    labWater->setText("Water: " + QString::fromStdString(getFilename(this->airFlowImagePath)));
    labAir->setText("Air: " + QString::fromStdString(getFilename(this->airFlowImagePath)));

    useBasicFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::BASIC);
    useImageFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE);
    useSimulatedFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION);
    useEnvObjFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_ENVOBJECTS);

    useRandomDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::RANDOM_DENSITY);
    useLayeredDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::LAYERED_DENSITY);
    useNativeDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::NATIVE);
    useImageDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::DENSITY_IMAGE);
    labDensityField->setText("File: " + QString::fromStdString(getFilename(this->densityFieldImagePath)));

    airDensity->setChecked(this->matterDensity < 1000.f);
    waterDensity->setChecked(this->matterDensity >= 1000.f);

    particleSizeSmall->setChecked(this->erosionSize < 5.f);
    particleSizeMedium->setChecked(this->erosionSize >= 5.f && this->erosionSize < 10.f);
    particleSizeBig->setChecked(this->erosionSize >= 10.f);

//    simulationTypeButton->setCurrentIndex(0);
    for (size_t i = 0; i < possibleSimTypes.size(); i++) {
        if (this->selectedSimulationType == possibleSimTypes[i])
            simulationTypeButton->setCurrentIndex(i);
    }

    return erosionLayout;
}

void ErosionInterface::show()
{
    computeTerrainBoundaries(voxelGrid.get(), nullptr);

    this->rocksPathSuccess.show();
    this->rocksPathFailure.show();
    ActionInterface::show();
}

void ErosionInterface::hide()
{
    this->rocksPathSuccess.hide();
    this->rocksPathFailure.hide();
    ActionInterface::hide();
}
