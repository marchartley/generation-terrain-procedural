#include "ErosionInterface.h"

#include "Interface/InterfaceUtils.h"
#include "Interface/TerrainGenerationInterface.h"

#include <chrono>

ErosionInterface::ErosionInterface(QWidget *parent)
    : ActionInterface("rock-throwing", parent)
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

    //    QTimer::singleShot(1000, this, &ErosionInterface::testManyManyErosionParameters);
}

void ErosionInterface::computePredefinedRocksLocations()
{
    if (this->voxelGrid == nullptr) return;
    std::vector<PARTICLE_INITIAL_LOCATION> locs = {SKY, RIVER, RANDOM, RIVER2, UNDERWATER, CENTER_TOP, FROM_X, EVERYWHERE, JUST_ABOVE_VOXELS};

    auto voxels = voxelGrid->getVoxelValues();
    auto terrainSurface = Mesh().applyMarchingCubes(voxels).vertexArray;
    for (auto& loc : locs) {
        initialPositionsAndDirections[loc] = std::vector<std::vector<std::pair<Vector3, Vector3>>>(20, std::vector<std::pair<Vector3, Vector3>>(maxParticles));
        auto& poses = initialPositionsAndDirections[loc];
        for (size_t i = 0; i < poses.size(); i++) {
            for (size_t j = 0; j < poses[i].size(); j++) {
                Vector3 position;
                Vector3 direction;
                if (loc == SKY) {
                    position = Vector3::random(voxelGrid->getDimensions().xy()) + Vector3(0, 0, voxelGrid->getSizeZ());
                    direction = Vector3(0, 0, -1);
                } else if (loc == RIVER) {
                    position = Vector3(5, 0, 0) + Vector3::random(Vector3(0, voxelGrid->getSizeY() * .4, 30), Vector3(0, voxelGrid->getSizeY()*.6, 40));
                    direction = (Vector3(1, 0, 0) + Vector3::random(1.f));
                } else if (loc == RIVER2) {
                    position = Vector3(5, 0, 0) + Vector3::random(Vector3(0, voxelGrid->getSizeY() * .1, 30), Vector3(0, voxelGrid->getSizeY()*.3, 40));
                    direction = (Vector3(1, 0, 0) + Vector3::random(1.f));
                } else if (loc == UNDERWATER) {
                    position = Vector3(5, 0, 0) + Vector3::random(Vector3(0, voxelGrid->getSizeY() * .4, voxelGrid->getSizeZ() * .2f), Vector3(0, voxelGrid->getSizeY()*.6, voxelGrid->getSizeZ() * .5f));
                    direction = (Vector3(1, 0, 0) + Vector3::random(.1f));
                } else if (loc == FROM_X) {
                    position = Vector3(5, 0, 10) + Vector3::random(Vector3(0, 0.f, 0.f), Vector3(0, voxelGrid->getSizeY(), voxelGrid->getSizeZ()));
                    direction = (Vector3(1, 0, 0) + Vector3::random(.1f));
                } else if (loc == CENTER_TOP) {
                    position = Vector3::random(Vector3(voxelGrid->getSizeX() * .45f, voxelGrid->getSizeY() * .45f, voxelGrid->getSizeZ() + 2.f), Vector3(voxelGrid->getSizeX() * .55f, voxelGrid->getSizeY()*.55f, voxelGrid->getSizeZ() + 2.f));
                    direction = (Vector3(1, 0, 0) + Vector3::random(.1f));
                } else if (loc == RANDOM) {
                    position = Vector3::random() * 100.f;
                    position.z = std::abs(position.z);
                    direction = -position.normalized();
                    position += voxelGrid->getDimensions().xy() * .5f;
                } else if (loc == EVERYWHERE) {
                    position = Vector3::random(Vector3(), voxelGrid->getDimensions());
                    direction = Vector3();
                }
                poses[i][j] = {position, direction};
            }
        }
    }
//    recomputeAboveVoxelRocksPositions();
}

void ErosionInterface::recomputeAboveVoxelRocksPositions(TerrainModel* terrain)
{
//    auto voxels = voxelGrid->getVoxelValues();
//    auto terrainSurface = Mesh().applyMarchingCubes(voxels).vertexArray;
    std::vector<std::pair<Vector3, Vector3>> terrainSurfaceAndNormalInversed;
    Vector3 terrainSize = terrain->getDimensions().xy();
//    GridV3 normals = terrain->getNormals();
    Vector3 down = Vector3(0, 0, -1);
    for (int x = 0; x < terrainSize.x; x++) {
        for (int y = 0; y < terrainSize.y; y++) {
            Vector3 pos = Vector3(x, y) + Vector3::random().xy();
            pos.z = terrain->getHeight(pos) + 1.f;
            terrainSurfaceAndNormalInversed.push_back({pos, down});
        }
    }
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

void ErosionInterface::display(Vector3 camPos)
{
    if (!this->isVisible())
        return;
    if (this->displayTrajectories) {
        this->rocksPathSuccess.display(GL_LINES, 3.f);
        this->rocksPathFailure.display(GL_LINES, 3.f);
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

            }


            if (this->densityUsed == UnderwaterErosion::RANDOM_DENSITY) {
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
            } else if (this->densityFieldImagePath != "") {
                densityField = GridF::fromImageBW(this->densityFieldImagePath).resize(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1.f).flip(true, false, false);
            }

            if (applyOn == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS) {
                terrain = voxelGrid.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP) {
                terrain = heightmap.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN) {
                terrain = implicitTerrain.get();
            } else if (applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN) {
                terrain = layerGrid.get();
            }

            if (location == JUST_ABOVE_VOXELS)
                recomputeAboveVoxelRocksPositions(terrain);

            Vector3 terrainDims = terrain->getDimensions();
            std::vector<std::vector<Vector3>> triangles;
            sumGeometry += timeIt([&]() {
                triangles = terrain->getGeometry(/*Vector3(30, 30, 30)*/).getTriangles();
            triangles.push_back({terrainDims * Vector3(0, 0, 0), terrainDims * Vector3(1, 0, 0), terrainDims * Vector3(1, 1, 0)});
            triangles.push_back({terrainDims * Vector3(1, 1, 0), terrainDims * Vector3(0, 1, 0), terrainDims * Vector3(0, 0, 0)});
            });
            boundariesTree = BVHTree();
            sumBVH += timeIt([&]() {
                boundariesTree.build(Triangle::vectorsToTriangles(triangles));
            });
            sumMeshingBoundaries += timeIt([&]() {
                boundariesMesh.fromArray(flattenArray(triangles));
            });
            std::cout << "Boundaries have " << triangles.size() << " triangles." << std::endl;
        });



        std::tie(lastRocksLaunched, nbPos, nbErosions) = erod.Apply(this->applyOn,
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
                                                                    selectedSimulationType
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
    //                if ((points[j] - points[j+1]).norm2() < 20.f) {
                        asOneVector.push_back(points[j]);
                        asOneVector.push_back(points[j + 1]);
    //                }
                }
            }
            this->rocksPathSuccess.fromArray(asOneVector);
            if (this->applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN) {
            }
        });
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
//                        for (const auto& p : path.points) {
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
{/*
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

        flowfieldFunction = [&](Vector3 pos) {
            return (voxelsEnvironmentalDensities.at(pos) < 500 ? airFlow.at(pos.xy()) : waterFlow.at(pos.xy()));
        };
    } else if (this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION) {
        auto fluidSim = voxelGrid->getFlowfield();

        flowfieldFunction = [&](Vector3 pos) {
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

    QPushButton* confirmButton = new QPushButton("Random");
    QPushButton* confirmFromCamButton = new QPushButton("From surface");
    QPushButton* confirmFromSkyButton = new QPushButton("Pluie");

    QCheckBox* displayTrajectoriesButton = new QCheckBox("Display");

    QRadioButton* applyOnVoxels = new QRadioButton("on voxels");
    QRadioButton* applyOnHeightmap = new QRadioButton("on heightmap");
    QRadioButton* applyOnImplicit = new QRadioButton("on implicit");
    QRadioButton* applyOnLayers = new QRadioButton("on layers");

    QRadioButton* useBasicFlowfield = new QRadioButton("Basic flowfield");
    QRadioButton* useImageFlowfield = new QRadioButton("Flowfield from image");
    QRadioButton* useSimulatedFlowfield = new QRadioButton("Simulation");
    QLabel* labWater = new QLabel;
    QLabel* labAir = new QLabel;
    QPushButton* browseWaterFlow = new QPushButton("...");
    QPushButton* browseAirFlow = new QPushButton("...");

    QRadioButton* useRandomDensity = new QRadioButton("Random density");
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

    QComboBox* simulationTypeButton = new QComboBox;
    simulationTypeButton->addItems({"MLB", "FLIP", "SPH", "Stable", "Warp"});
//    simulationTypeButton->addItem("MLB", QVariant(FluidSimType::LBM));
//    simulationTypeButton->addItem("FLIP", QVariant(FluidSimType::FLIP));
//    simulationTypeButton->addItem("SPH", QVariant(FluidSimType::SPH));
//    simulationTypeButton->addItem("Stable", QVariant(FluidSimType::STABLE));
//    simulationTypeButton->addItem("Warp", QVariant(FluidSimType::WARP));

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
                                                     })
                                                 }));
    erosionLayout->addWidget(createVerticalGroup({
                                                     confirmFromCamButton,
                                                     confirmFromSkyButton,
                                                     confirmButton,
                                                     displayTrajectoriesButton,
                                                     createVerticalGroup({
                                                         applyOnHeightmap,
                                                         applyOnVoxels,
                                                         applyOnImplicit,
                                                         applyOnLayers
                                                     }),
                                                     createVerticalGroup({
                                                         useBasicFlowfield,
                                                         useImageFlowfield,
                                                         labWater, browseWaterFlow,
                                                         labAir, browseAirFlow,
                                                         /*createHorizontalGroup({*/useSimulatedFlowfield, simulationTypeButton/*})*/
                                                     }),
//                                                     lotsOfTestsButton,
                                                     createVerticalGroup({
                                                         useImageDensity,
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

    QObject::connect(confirmButton, &QPushButton::pressed, this, &ErosionInterface::throwFromSide);
    QObject::connect(confirmFromCamButton, &QPushButton::pressed, this, &ErosionInterface::throwFromCam);
    QObject::connect(confirmFromSkyButton, &QPushButton::pressed, this, &ErosionInterface::throwFromSky);

    QObject::connect(displayTrajectoriesButton, &QCheckBox::toggled, this, [&](bool checked) { this->displayTrajectories = checked; });

    QObject::connect(applyOnVoxels, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS; });
    QObject::connect(applyOnHeightmap, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP; });
    QObject::connect(applyOnImplicit, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN; });
    QObject::connect(applyOnLayers, &QRadioButton::toggled, this, [&]() { this->applyOn = UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN; });

    QObject::connect(useBasicFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::BASIC; });
    QObject::connect(useImageFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE; });
    QObject::connect(useSimulatedFlowfield, &QRadioButton::toggled, this, [&]() { this->flowfieldUsed = UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION; });

    QObject::connect(browseWaterFlow, &QPushButton::pressed, this, [=]() { this->browseWaterFlowFromFile(); labWater->setText("Water: " + QString::fromStdString(getFilename(this->waterFlowImagePath)));});
    QObject::connect(browseAirFlow, &QPushButton::pressed, this, [=]() { this->browseAirFlowFromFile(); labAir->setText("Air: " + QString::fromStdString(getFilename(this->airFlowImagePath)));} );

    QObject::connect(useRandomDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::RANDOM_DENSITY; });
    QObject::connect(useNativeDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::NATIVE; });
    QObject::connect(useImageDensity, &QRadioButton::toggled, this, [&]() { this->densityUsed = UnderwaterErosion::DENSITY_TYPE::DENSITY_IMAGE; });
    QObject::connect(densityFieldFileChooser, &QPushButton::pressed, this, [=]() { this->browseDensityFieldFromFile(); labAir->setText("File: " + QString::fromStdString(getFilename(this->airFlowImagePath)));} );

    QObject::connect(lotsOfTestsButton, &QPushButton::pressed, this, &ErosionInterface::testManyManyErosionParameters);

    QObject::connect(airDensity, &QRadioButton::pressed, this, [=]() { this->matterDensity = 500.f; });
    QObject::connect(waterDensity, &QRadioButton::pressed, this, [=]() { this->matterDensity = 1662.f; });

    QObject::connect(particleSizeSmall, &QRadioButton::pressed, this, [=]() { this->erosionSize = 4.f; });
    QObject::connect(particleSizeMedium, &QRadioButton::pressed, this, [=]() { this->erosionSize = 8.f; });
    QObject::connect(particleSizeBig, &QRadioButton::pressed, this, [=]() { this->erosionSize = 12.f; });

    QObject::connect(simulationTypeButton, &QComboBox::currentTextChanged, this, [=](QString text) { this->selectedSimulationType = FluidSimTypeFromString(text.toStdString()); });
//    QObject::connect(simulationTypeButton, &QComboBox::currentIndexChanged, this, [=](int index) { this->selectedSimulationType = simulationTypeButton->item; });

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

    displayTrajectoriesButton->setChecked(this->displayTrajectories);

    applyOnVoxels->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::DENSITY_VOXELS);
    applyOnHeightmap->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::HEIGHTMAP);
    applyOnImplicit->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::IMPLICIT_TERRAIN);
    applyOnLayers->setChecked(this->applyOn == UnderwaterErosion::EROSION_APPLIED::LAYER_TERRAIN);

    labWater->setText("Water: " + QString::fromStdString(getFilename(this->airFlowImagePath)));
    labAir->setText("Air: " + QString::fromStdString(getFilename(this->airFlowImagePath)));

    useBasicFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::BASIC);
    useImageFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLOWFIELD_IMAGE);
    useSimulatedFlowfield->setChecked(this->flowfieldUsed == UnderwaterErosion::FLOWFIELD_TYPE::FLUID_SIMULATION);

    useRandomDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::RANDOM_DENSITY);
    useNativeDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::NATIVE);
    useImageDensity->setChecked(this->densityUsed == UnderwaterErosion::DENSITY_TYPE::DENSITY_IMAGE);
    labDensityField->setText("File: " + QString::fromStdString(getFilename(this->densityFieldImagePath)));

    airDensity->setChecked(this->matterDensity < 1000.f);
    waterDensity->setChecked(this->matterDensity >= 1000.f);

    particleSizeSmall->setChecked(this->erosionSize < 5.f);
    particleSizeMedium->setChecked(this->erosionSize >= 5.f && this->erosionSize < 10.f);
    particleSizeBig->setChecked(this->erosionSize >= 10.f);

    simulationTypeButton->setCurrentIndex(0);

    return erosionLayout;
}

void ErosionInterface::show()
{
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
