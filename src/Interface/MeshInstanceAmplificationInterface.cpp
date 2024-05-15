#include "MeshInstanceAmplificationInterface.h"

#include "Interface/CommonInterface.h"

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

#include "Interface/TerrainGenerationInterface.h"

MeshInstanceAmplificationInterface::MeshInstanceAmplificationInterface(QWidget* parent)
    : ActionInterface("meshinstance", "Mesh Instance Amplification", "view", "Amplify the terrain with meshes", "amplification_instances.png", parent)
{
    meshInstancesFile.path = "saved_maps/meshInstances.json";
    meshInstancesFile.onChange([&](std::string content) {
        this->readMeshInstanceFile(content);
    });

    QTimer* hotreloadTimer = new QTimer(this);
    hotreloadTimer->setInterval(500);
    QObject::connect(hotreloadTimer, &QTimer::timeout, this, [&]() {
        meshInstancesFile.check();
    });
    hotreloadTimer->start();
}

void MeshInstanceAmplificationInterface::display(const Vector3& camPos)
{
//    if (!this->isVisible())
//        return;

    /*
    if (this->displayRocks || this->displayCorals) {
        if (displayRocks) {
            for (size_t i = 0; i < rocksIndicesAndPositionAndSize.size(); i++) {
                int iRock;
                Vector3 pos;
                float size;
                std::tie(iRock, pos, size) = rocksIndicesAndPositionAndSize[i];
                possibleRocks[iRock].shader->setVector("instanceOffset", pos);
                possibleRocks[iRock].shader->setFloat("sizeFactor", size);
                possibleRocks[iRock].display();
            }
        }
        if (displayCorals) {
            for (size_t i = 0; i < coralsIndicesAndPositionAndSize.size(); i++) {
                int iCoral;
                Vector3 pos;
                float size;
                std::tie(iCoral, pos, size) = coralsIndicesAndPositionAndSize[i];
                possibleCorals[iCoral].shader->setVector("instanceOffset", pos);
                possibleCorals[iCoral].shader->setFloat("sizeFactor", size);
                possibleCorals[iCoral].display();
            }
        }
    }*/

    if (displayEnvObjects) {

        if (autoUpdateEnvObjLocations) {
            this->regenerateAllTypePositions();
        }
        // allAtOnce.display();
        // return;
        // Dirty, to remove one day :
        float terrainHeightFactor = dynamic_cast<TerrainGenerationInterface*>(viewer->interfaces["terraingeneration"].get())->heightFactor;

        for (auto& meshType : meshesOptions) {
            if (!meshType.displayed) continue;
            for (size_t i = 0; i < meshType.indices.size(); i++) {
                int iMesh = meshType.indices[i];
                Vector3& pos = meshType.positions[i];
                float size = meshType.sizes[i];
                Vector3& orientation = meshType.orientations[i];
                Matrix rotMatrix = Vector3(0, 0, -(orientation.getAngleWith(Vector3(0, 1, 0)))).toRotationMatrix().toHomogeneous();
                rotMatrix[2][2] *= -1.f;
                auto values = rotMatrix.toStdVector();
//                std::tie(iMesh, pos, size) = meshType.indicesAndPositionsAndSizes[i];
                meshType.possibleMeshes[iMesh].shader->setVector("instanceOffset", (pos + Vector3(0, 0, 1.f)) * Vector3(1, 1, terrainHeightFactor) - meshType.requiredTranslation * size);
                meshType.possibleMeshes[iMesh].shader->setFloat("sizeFactor", size);
                meshType.possibleMeshes[iMesh].shader->setMatrix("instanceRotation", values.data(), 4, 4);
//                meshType.possibleMeshes[iMesh].displayWithOutlines(meshType.color);
                meshType.possibleMeshes[iMesh].shader->setVector("color", meshType.color);
                meshType.possibleMeshes[iMesh].display();
            }
        }
    }
}

void MeshInstanceAmplificationInterface::replay(nlohmann::json action)
{

}

void MeshInstanceAmplificationInterface::reloadShaders()
{
    bool verbose = true;

    std::string pathToShaders = "src/Shaders/";

    std::string vRockShader = pathToShaders + "rockShader.vert";
    std::string fRockShader = pathToShaders + "rockShader.frag";
    std::string vTreeShader = pathToShaders + "meshInstancesShader.vert";
    std::string fTreeShader = pathToShaders + "meshInstancesShader.frag";

    // allAtOnce = Mesh(std::make_shared<Shader>(pathToShaders + "no_shader.vert", pathToShaders + "no_shader.frag"));

    std::vector<QString> coralPaths;
    std::vector<QString> rocksPaths;
//    std::vector<QString> algaePaths;

    displayProcessTime("Loading 3D models... ", [&]() {

        QDirIterator itCorals("src/assets/models/coral/", QDir::Files, QDirIterator::Subdirectories);
        std::shared_ptr<Shader> coralsShader = std::make_shared<Shader>(vRockShader, fRockShader);
        while (itCorals.hasNext()) {
            QString dir = itCorals.next();
            coralPaths.push_back(dir);
        }

        size_t nbCorals = coralPaths.size();
        if (this->numberOfLoadedCorals != -1) nbCorals = std::min(nbCorals, (size_t)numberOfLoadedCorals);
        this->possibleCorals = std::vector<Mesh>(nbCorals);

        #pragma omp parallel for
        for (size_t i = 0; i < nbCorals; i++) {
            QString& dir = coralPaths[i];
            // Normalize it and move it upward so the anchor is on the ground
            possibleCorals[i] = Mesh(coralsShader).fromStl(dir.toStdString()).normalize().rotate(deg2rad(180), 0, 0);
        }

        QDirIterator itRocks("src/assets/models/rock/", QDir::Files, QDirIterator::Subdirectories);
        std::shared_ptr<Shader> rocksShader = std::make_shared<Shader>(vRockShader, fRockShader);
        while (itRocks.hasNext()) {
            QString dir = itRocks.next();
            rocksPaths.push_back(dir);
        }
        size_t nbRocks = rocksPaths.size();
        if (this->numberOfLoadedRocks != -1) nbRocks = std::min(nbRocks, (size_t)numberOfLoadedRocks);
        this->possibleRocks = std::vector<Mesh>(nbRocks);
        #pragma omp parallel for
        for (size_t i = 0; i < nbRocks; i++) {
            QString& dir = rocksPaths[i];
            // Normalize it and move it upward so the anchor is on the ground
            possibleRocks[i] = Mesh(rocksShader).fromStl(dir.toStdString()).normalize();
        }

        /*
        meshesOptions.push_back(InstantiationMeshOption("boulder", {3.f, 8.f}, {.2f, 1.f, .5f, 1.f}));
//        meshesOptions.push_back(InstantiationMeshOption("reef", "coral", {1.f, 5.f}, {1.f, .5f, .5f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("coralpolyp", "corals", {5.f, 5.f}, {1.f, .5f, .5f, 1.f}, Vector3(), {5, 20}, 10.f));
        meshesOptions.push_back(InstantiationMeshOption("coralpolypflat", "corals", {5.f, 5.f}, {1.f, .5f, .5f, 1.f}, Vector3(0, 0, .25f), {5, 20}, 10.f));
        meshesOptions.push_back(InstantiationMeshOption("algae", {10.f, 15.f}, {.1f, .5f, .1f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("tree", {20.f, 40.f}, {.1f, 1.f, .1f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("arch", "arche", {20.f, 40.f}, {.8f, .8f, .6f, 1.f}, Vector3(0, 0, .25f)));
        meshesOptions.push_back(InstantiationMeshOption("bigRock", "rocks", {20.f, 40.f}, {.8f, .8f, .6f, 1.f}, Vector3(0, 0, .5f)));
        meshesOptions.push_back(InstantiationMeshOption("smallRock", "rocks", {3.f, 5.f}, {.8f, .8f, .6f, 1.f}, Vector3(0, 0, .25f), {10, 30}, 5.f));
//        meshesOptions.push_back(InstantiationMeshOption("island", {20.f, 40.f}, {.5f, .1f, .5f, 1.f}));

        for (auto& meshType : meshesOptions) {
            meshType.displayed = true;
            QDirIterator it(QString::fromStdString("src/assets/models/" + meshType.folderName + "/"), QDir::Files, QDirIterator::Subdirectories);
            std::shared_ptr<Shader> shader = std::make_shared<Shader>(vTreeShader, fTreeShader);
            std::vector<QString> paths;
            while (it.hasNext()) {
                QString dir = it.next();
                if (endsWith(dir.toStdString(), ".ignore")) continue;
                paths.push_back(dir);
            }
            size_t nbElements = paths.size();
            if (meshType.numberOfLoadedMesh != -1) nbElements = std::min(nbElements, (size_t)meshType.numberOfLoadedMesh);
            meshType.possibleMeshes = std::vector<Mesh>(nbElements);

//            #pragma omp parallel for
            for (size_t i = 0; i < nbElements; i++) {
                QString& dir = paths[i];
                // Normalize it and move it upward so the anchor is on the ground
                meshType.possibleMeshes[i] = Mesh(shader);

                if (dir.endsWith("fbx", Qt::CaseInsensitive)) {
                    meshType.possibleMeshes[i].fromFBX(dir.toStdString());
                } else if (dir.endsWith("stl", Qt::CaseInsensitive)) {
                    meshType.possibleMeshes[i].fromStl(dir.toStdString()).scale(Vector3(1.f, 1.f, -1.f));
                } else if (!dir.endsWith(".ignore")) {
                    std::cerr << "Unable to open file " << dir.toStdString() << std::endl;
                }

                meshType.possibleMeshes[i].normalize().translate(Vector3(0.f, 0.f, (meshType.name == "island" && false ? 0.f : -.5f)));
                meshType.possibleMeshes[i].cullFace = false;
            }
        }*/
    }, verbose);
}

QLayout* MeshInstanceAmplificationInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    QCheckBox* displayCoralsCheckbox = new QCheckBox("Display corals");
    QCheckBox* displayRocksCheckbox = new QCheckBox("Display rocks");

    for (auto& meshType : meshesOptions) {
        CheckboxElement* displayElement = new CheckboxElement("Display " + meshType.name);
        displayElement->setChecked(meshType.displayed);
        displayElement->setOnChecked([&](bool check) { this->setDisplayingType(meshType, check); });
        layout->addWidget(displayElement->get());
    }

//    layout->addWidget(displayCoralsCheckbox);
//    layout->addWidget(displayRocksCheckbox);

    ButtonElement* recomputePositionsButton = new ButtonElement("Recompute positions");
    recomputePositionsButton->setOnClick([&]() { this->afterTerrainUpdated(); }); // Should be modified in the future
    layout->addWidget(recomputePositionsButton->get());

    displayCoralsCheckbox->setChecked(this->displayCorals);
    displayRocksCheckbox->setChecked(this->displayRocks);

    QObject::connect(displayCoralsCheckbox, &QCheckBox::toggled, this, &MeshInstanceAmplificationInterface::setCoralsDisplayed);
    QObject::connect(displayRocksCheckbox, &QCheckBox::toggled, this, &MeshInstanceAmplificationInterface::setRocksDisplayed);

    return layout;
}

std::vector< AABBox > MeshInstanceAmplificationInterface::getAvailablePositionsForMaterial(TerrainTypes target)
{
    std::vector< AABBox > positions;
    auto materialsAndHeightsGrid = layerGrid->getLayers();
    for (int x = 0; x < layerGrid->getSizeX(); x++) {
        for (int y = 0; y < layerGrid->getSizeY(); y++) {
            auto materialsAndHeights = materialsAndHeightsGrid.at(x, y);

            float currentZ = 0.f;
            for (size_t i = 0; i < materialsAndHeights.size(); i++) {
                auto& material = materialsAndHeights[i].first;
                auto& height = materialsAndHeights[i].second;

                if (material == target) {
                    positions.push_back({Vector3(x, y, currentZ), Vector3(x, y, currentZ + height)});
                }
                currentZ += height;
            }
        }
    }
    return positions;
}

std::vector<AABBox> MeshInstanceAmplificationInterface::getCoralAvailablePositions()
{
    auto intervals = getAvailablePositionsForMaterial(TerrainTypes::CORAL);
    std::vector<AABBox> extendedPositions;

    for (auto& interv : intervals) {
        Vector3 start = interv.min();
        Vector3 end = interv.max();
        float minHeight = start.z;
        float maxHeight = end.z;

        extendedPositions.push_back({start.xy() + Vector3(0, 0, minHeight), end});
        extendedPositions.push_back({start.xy() + Vector3(0, 0, maxHeight), end});

        /*for (int z = std::floor(minHeight); z < std::ceil(maxHeight); z++) {
            extendedPositions.push_back({start.xy() + Vector3(0, 0, z), end});
        }*/
    }
    return extendedPositions;
}

std::vector<AABBox> MeshInstanceAmplificationInterface::getRocksAvailablePositions()
{
//    return getAvailablePositionsForMaterial(TerrainTypes::ROCK);
    auto intervals = getAvailablePositionsForMaterial(TerrainTypes::ROCK);
    std::vector<AABBox> extendedPositions;

    for (auto& interv : intervals) {
        Vector3 start = interv.min();
        Vector3 end = interv.max();
        float minHeight = start.z;
        float maxHeight = end.z;

        extendedPositions.push_back({start.xy() + Vector3(0, 0, minHeight), end});
        extendedPositions.push_back({start.xy() + Vector3(0, 0, maxHeight), end});

        /*for (int z = std::floor(minHeight); z < std::ceil(maxHeight); z++) {
            extendedPositions.push_back({start.xy() + Vector3(0, 0, z), end});
        }*/
    }
    return extendedPositions;
}

std::vector<std::tuple<Vector3, float, int> > MeshInstanceAmplificationInterface::getPositionsFor(std::string type)
{
    std::vector<std::tuple<Vector3, float, int>> positionsAndGrowthFactor;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (toUpper(obj->name) != toUpper(type)) continue;
        float growthFactor = obj->computeGrowingState();
        if (auto asPoint = dynamic_cast<EnvPoint*>(obj)) {
            positionsAndGrowthFactor.push_back({Vector3(asPoint->position.x, asPoint->position.y, heightmap->getHeight(asPoint->position)), growthFactor, asPoint->ID});
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(obj)) {
            auto path = asCurve->curve.getPath(60);
            for (auto p : path) {
                p = p + Vector3::random(asCurve->width * .5f);
                positionsAndGrowthFactor.push_back({Vector3(p.x, p.y, heightmap->getHeight(p)), growthFactor, asCurve->ID});
            }
        } else if (auto asArea = dynamic_cast<EnvArea*>(obj)) {
//            float totalArea = asArea->area.computeArea();
            std::vector<Vector3> randomPoints;
            AABBox box = AABBox(asArea->curve.AABBox());
            for (int x = box.min().x; x < box.max().x; x += 5) {
                for (int y = box.min().y; y < box.max().y; y += 5) {
                    Vector3 pos(x, y, 0);
                    if (random_gen::generate_perlin(pos.x * 5.f, pos.y * 5.f) > .5f && asArea->curve.containsXY(pos))
                        randomPoints.push_back(pos);
                }
            }
//            auto randomPoints = asArea->area.randomPointsInside(60);
            for (const auto& p : randomPoints) {
                positionsAndGrowthFactor.push_back({Vector3(p.x, p.y, heightmap->getHeight(p)), growthFactor, asArea->ID});
            }
        }
    }
    return positionsAndGrowthFactor;
}

void MeshInstanceAmplificationInterface::readMeshInstanceFile(const std::string &fileContent)
{
    std::string pathToShaders = "src/Shaders/";
    std::string vTreeShader = pathToShaders + "meshInstancesShader.vert";
    std::string fTreeShader = pathToShaders + "meshInstancesShader.frag";

    meshesOptions.clear();
    auto json = nlohmann::json::parse(toLower(fileContent));
    for (auto& instance : json) {
        std::string name = instance["name"];
        if (startsWith(name, "--")) continue;

        std::string folderName = instance["foldername"];
        if (folderName == "") folderName = name;
        std::vector<std::vector<float>> colors;
        auto jsonColors = instance["colors"];
        for (auto& col : jsonColors) {
            colors.push_back(json_to_color(col));
        }
        Vector3 translation = json_to_vec3(instance["translation"]);
        int minInstances = instance["mininstances"];
        int maxInstances = instance["maxinstances"];
        float minSize = instance["minsize"];
        float maxSize = instance["maxsize"];
        float radius = instance["radius"];

        InstantiationMeshOption meshType(name, folderName, {minSize, maxSize}, colors[0], translation, {minInstances, maxInstances}, radius);
        meshType.displayed = true;

        QDirIterator it(QString::fromStdString("src/assets/models/" + meshType.folderName + "/"), QDir::Files, QDirIterator::Subdirectories);
        std::shared_ptr<Shader> shader = std::make_shared<Shader>(vTreeShader, fTreeShader);
        std::vector<QString> paths;
        while (it.hasNext()) {
            QString dir = it.next();
            if (endsWith(dir.toStdString(), ".ignore")) continue;
            paths.push_back(dir);
        }
        size_t nbElements = paths.size();
        if (meshType.numberOfLoadedMesh != -1) nbElements = std::min(nbElements, (size_t)meshType.numberOfLoadedMesh);
        meshType.possibleMeshes = std::vector<Mesh>(nbElements);

//            #pragma omp parallel for
        for (size_t i = 0; i < nbElements; i++) {
            QString& dir = paths[i];
            // Normalize it and move it upward so the anchor is on the ground
            meshType.possibleMeshes[i] = Mesh(shader);

            if (dir.endsWith("fbx", Qt::CaseInsensitive)) {
                meshType.possibleMeshes[i].fromFBX(dir.toStdString());
            } else if (dir.endsWith("stl", Qt::CaseInsensitive)) {
                meshType.possibleMeshes[i].fromStl(dir.toStdString()).scale(Vector3(1.f, 1.f, -1.f));
            } else if (!dir.endsWith(".ignore")) {
                std::cerr << "Unable to open file " << dir.toStdString() << std::endl;
            }

            meshType.possibleMeshes[i].normalize().translate(Vector3(0.f, 0.f, -.5f));
            meshType.possibleMeshes[i].cullFace = false;
        }

        meshesOptions.push_back(meshType);
    }
}

void MeshInstanceAmplificationInterface::setCoralsDisplayed(bool display)
{
    this->displayCorals = display;
}

void MeshInstanceAmplificationInterface::setRocksDisplayed(bool display)
{
    this->displayRocks = display;
}

void MeshInstanceAmplificationInterface::setDisplayingType(InstantiationMeshOption& options, bool display)
{
    options.displayed = display;
}

void MeshInstanceAmplificationInterface::afterTerrainUpdated()
{
    this->regenerateRocksPositions();
    this->regenerateAllTypePositions();
}

void MeshInstanceAmplificationInterface::regenerateRocksPositions()
{
    auto coralAvailablePositions = (this->possibleCorals.size() > 0 ? this->getCoralAvailablePositions() : std::vector<AABBox>());
    auto rocksAvailablePositions = (this->possibleRocks.size() > 0 ? this->getRocksAvailablePositions() : std::vector<AABBox>());

    std::shuffle(coralAvailablePositions.begin(), coralAvailablePositions.end(), random_gen::random_generator);
    std::shuffle(rocksAvailablePositions.begin(), rocksAvailablePositions.end(), random_gen::random_generator);

    this->rocksIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(std::min(numberOfRocksDisplayed, int(rocksAvailablePositions.size())));
    for (size_t i = 0; i < rocksIndicesAndPositionAndSize.size(); i++) {
        rocksIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleRocks.size())),
                                              Vector3(rocksAvailablePositions[i].min()),
                                                //Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                              random_gen::generate(5.f, 15.f)
                                                     );
    }
    this->coralsIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(std::min(numberOfCoralsDisplayed, int(coralAvailablePositions.size())));
    for (size_t i = 0; i < coralsIndicesAndPositionAndSize.size(); i++) {
        coralsIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleCorals.size())),
                                                Vector3(coralAvailablePositions[i].min()),
                                              // Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                              random_gen::generate(5.f, 15.f)
                                                     );
    }
}

void MeshInstanceAmplificationInterface::regenerateAllTypePositions()
{
    for (auto& meshType : meshesOptions) {
        auto availablePositions = (meshType.possibleMeshes.size() > 0 ? this->getPositionsFor(meshType.name) : std::vector<std::tuple<Vector3, float, int>>());
//        std::shuffle(availablePositions.begin(), availablePositions.end(), random_gen::random_generator);

//        meshType.indicesAndPositionsAndSizes = std::vector<std::tuple<int, Vector3, float>>(); //(std::min(meshType.numberDisplayed, int(availablePositions.size())));
        meshType.clear();
        for (size_t i = 0; i < availablePositions.size(); i++) { //meshType.indicesAndPositionsAndSizes.size(); i++) {
            Vector3 position = std::get<0>(availablePositions[i]);
            float growth = std::get<1>(availablePositions[i]);
            int objIndex = std::get<2>(availablePositions[i]);
//            Vector3 orientation = EnvObject::flowfield(position.xy()).normalized(); //Vector3(0, 0, EnvObject::flowfield(position).getAngleWith(Vector3(1, 0, 0)));

            int nbInstances = int(interpolation::inv_linear((random_gen::generate_perlin(position.x * 1000, position.y * 1000, objIndex) + 1) * .5f, meshType.minMaxInstances.first, meshType.minMaxInstances.second + 1));
            int nbRandomValues = 4;
            std::vector<float> randomVals(nbInstances * nbRandomValues);
            for (int iRand = 0; iRand < randomVals.size(); iRand++) {
                randomVals[iRand] = (random_gen::generate_perlin(position.x * 1000, position.y * 1000, (iRand + objIndex) * 100) + 1) * .5f;
            }
            for (int iInstance = 0; iInstance < nbInstances; iInstance++) {
                size_t randomIdx = iInstance * nbRandomValues;
                Vector3 instancePos = position + (Vector3(randomVals[randomIdx + 0], randomVals[randomIdx + 1], 0.f) - Vector3(.5f, .5f, 0.f)) * meshType.radius;
                if (!Vector3::isInBox(instancePos.xy(), Vector3(.1f, .1f, 0), Vector3(heightmap->getSizeX() - .1f, heightmap->getSizeY() - .1f, 0))) continue;
                instancePos.z = heightmap->getHeight(instancePos.xy());
                Vector3 orientation = EnvObject::flowfield(instancePos.xy()).normalized();
                int meshIndex = int(randomVals[randomIdx + 2] * meshType.possibleMeshes.size());
                float scale = interpolation::inv_linear(randomVals[randomIdx + 3], meshType.minMaxSizes.first * growth, meshType.minMaxSizes.second * growth);
//                instancePos.z -= scale * .15f;
//                meshType.indicesAndPositionsAndSizes.push_back({meshIndex, instancePos, scale});
                meshType.add(meshIndex, instancePos, scale, orientation);
            }

            /*int meshIndex = int(random_gen::generate(0, meshType.possibleMeshes.size()));
            Vector3 position = Vector3(availablePositions[i].first) * (meshType.name == "island" ? Vector3(1, 1, 0) : Vector3(1, 1, 1));
            float scale = random_gen::generate(meshType.minMaxSizes.first, meshType.minMaxSizes.second) * availablePositions[i].second;
            */
//            meshType.indicesAndPositionsAndSizes[i] = {meshIndex, position, scale};
        }
//        std::cout << meshType.name << ": " << availablePositions.size() << " pos available, " << meshType.indicesAndPositionsAndSizes.size() << " instanciated." << std::endl;
    }

    /*allAtOnce.clear();
    for (auto& meshType : meshesOptions) {
        for (int i = 0; i < meshType.indices.size(); i++) {
            auto mesh = meshType.possibleMeshes[meshType.indices[i]];
            float scale = meshType.sizes[i];
            auto& pos = meshType.positions[i];
            allAtOnce.merge(mesh.scale(scale).translate(pos), false);
        }
    }*/
}

void InstantiationMeshOption::clear()
{
    this->indices.clear();
    this->positions.clear();
    this->sizes.clear();
    this->orientations.clear();
}

void InstantiationMeshOption::add(int index, const Vector3 &position, float size, const Vector3 &orientation)
{
    this->indices.push_back(index);
    this->positions.push_back(position);
    this->sizes.push_back(size);
    this->orientations.push_back(orientation);
}
