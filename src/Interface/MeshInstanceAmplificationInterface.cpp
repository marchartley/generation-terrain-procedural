#include "MeshInstanceAmplificationInterface.h"

#include "Interface/CommonInterface.h"

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

MeshInstanceAmplificationInterface::MeshInstanceAmplificationInterface(QWidget* parent)
    : ActionInterface("meshinstance", "Mesh Instance Amplification", "view", "Amplify the terrain with meshes", "amplification_instances.png", parent)
{

}

void MeshInstanceAmplificationInterface::display(const Vector3& camPos)
{
    if (!this->isVisible())
        return;

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
    }

    for (auto& meshType : meshesOptions) {
        if (!meshType.displayed) continue;
        for (size_t i = 0; i < meshType.indicesAndPositionsAndSizes.size(); i++) {
            int iMesh;
            Vector3 pos;
            float size;
            std::tie(iMesh, pos, size) = meshType.indicesAndPositionsAndSizes[i];
            meshType.possibleMeshes[iMesh].shader->setVector("instanceOffset", pos);
            meshType.possibleMeshes[iMesh].shader->setFloat("sizeFactor", size);
            meshType.possibleMeshes[iMesh].displayWithOutlines(meshType.color);
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

        meshesOptions.push_back(InstantiationMeshOption("boulder", {3.f, 8.f}, {.2f, 1.f, .5f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("reef", "coral", {1.f, 5.f}, {1.f, .5f, .5f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("algae", {10.f, 15.f}, {.1f, .5f, .1f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("tree", {20.f, 40.f}, {.1f, 1.f, .1f, 1.f}));
        meshesOptions.push_back(InstantiationMeshOption("island", {20.f, 40.f}, {.5f, .1f, .5f, 1.f}));

        for (auto& meshType : meshesOptions) {
            QDirIterator it(QString::fromStdString("src/assets/models/" + meshType.folderName + "/"), QDir::Files, QDirIterator::Subdirectories);
            std::shared_ptr<Shader> shader = std::make_shared<Shader>(vTreeShader, fTreeShader);
            std::vector<QString> paths;
            while (it.hasNext()) {
                QString dir = it.next();
                paths.push_back(dir);
            }
            size_t nbElements = paths.size();
            if (meshType.numberOfLoadedMesh != -1) nbElements = std::min(nbElements, (size_t)meshType.numberOfLoadedMesh);
            meshType.possibleMeshes = std::vector<Mesh>(nbElements);

            #pragma omp parallel for
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

                meshType.possibleMeshes[i].normalize().translate(Vector3(0.f, 0.f, (meshType.name == "island" && false ? 0.f : -.5f)) + meshType.requiredTranslation);
                meshType.possibleMeshes[i].cullFace = false;
            }
        }
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

std::vector<std::pair<Vector3, float>> MeshInstanceAmplificationInterface::getPositionsFor(std::string type)
{
    std::vector<std::pair<Vector3, float>> positionsAndGrowthFactor;
    for (auto& obj : EnvObject::instantiatedObjects) {
        if (toUpper(obj->name) != toUpper(type)) continue;
        float growthFactor = obj->computeGrowingState();
        if (auto asPoint = dynamic_cast<EnvPoint*>(obj)) {
            positionsAndGrowthFactor.push_back({Vector3(asPoint->position.x, asPoint->position.y, heightmap->getHeight(asPoint->position)), growthFactor});
        } else if (auto asCurve = dynamic_cast<EnvCurve*>(obj)) {
            auto path = asCurve->curve.getPath(60);
            for (auto p : path) {
                p = p + Vector3::random(asCurve->width * .5f);
                positionsAndGrowthFactor.push_back({Vector3(p.x, p.y, heightmap->getHeight(p)), growthFactor});
            }
        } else if (auto asArea = dynamic_cast<EnvArea*>(obj)) {
            auto randomPoints = asArea->area.randomPointsInside(60);
            for (auto p : randomPoints) {
                positionsAndGrowthFactor.push_back({Vector3(p.x, p.y, heightmap->getHeight(p)), growthFactor});
            }
        }
    }
    return positionsAndGrowthFactor;
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
        auto availablePositions = (meshType.possibleMeshes.size() > 0 ? this->getPositionsFor(meshType.name) : std::vector<std::pair<Vector3, float>>());
        std::shuffle(availablePositions.begin(), availablePositions.end(), random_gen::random_generator);

        meshType.indicesAndPositionsAndSizes = std::vector<std::tuple<int, Vector3, float>>(std::min(meshType.numberDisplayed, int(availablePositions.size())));
        for (size_t i = 0; i < meshType.indicesAndPositionsAndSizes.size(); i++) {
            meshType.indicesAndPositionsAndSizes[i] = std::make_tuple<int, Vector3, float>(
                                                  int(random_gen::generate(0, meshType.possibleMeshes.size())),
                                                  Vector3(availablePositions[i].first) * (meshType.name == "island" ? Vector3(1, 1, 0) : Vector3(1, 1, 1)),
                                                    //Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                                  random_gen::generate(meshType.minMaxSizes.first, meshType.minMaxSizes.second) * availablePositions[i].second
                                                         );
        }
//        std::cout << meshType.name << ": " << availablePositions.size() << " pos available, " << meshType.indicesAndPositionsAndSizes.size() << " instanciated." << std::endl;
    }
}
