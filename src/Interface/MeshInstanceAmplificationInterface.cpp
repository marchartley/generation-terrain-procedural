#include "MeshInstanceAmplificationInterface.h"

MeshInstanceAmplificationInterface::MeshInstanceAmplificationInterface(QWidget* parent)
    : ActionInterface("MeshInstanceAmplificationInterface", parent)
{

}

void MeshInstanceAmplificationInterface::display()
{

    if (this->displayRocks || this->displayCorals) {
        // Check if something changed on the terrain :
        /*std::cout << this->previousHistoryIndex << " " << voxelGrid->getCurrentHistoryIndex() << std::flush;
        if (this->previousHistoryIndex != voxelGrid->getCurrentHistoryIndex()) {
            this->previousHistoryIndex = voxelGrid->getCurrentHistoryIndex();
            regenerateRocksPositions();
            std::cout << "Regenerated!" << std::endl;
        } else {
            std::cout << "Nothing" << std::endl;
        }*/

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
        /*
    //                marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values + .5f);

            values.raiseErrorOnBadCoord = false;
            values.defaultValueOnBadCoord = -1.f;
            for (size_t i = 0; i < rocksIndicesAndPositionAndSize.size(); i++) {
                int iRock;
                Vector3 pos;
                float size;
                std::tie(iRock, pos, size) = rocksIndicesAndPositionAndSize[i];
                // Must be on ground, and have water somewhere around
                if (values.at(pos) > 0) {
                    bool isOnBorder = false;
                    for (int dx = -1; dx <= 1 && !isOnBorder; dx++)
                        for (int dy = -1; dy <= 1 && !isOnBorder; dy++)
                            for (int dz = -1; dz <= 1 && !isOnBorder; dz++)
                                if (values.at(pos + Vector3(dx, dy, dz)) < 0) {
                                    isOnBorder = true;
                                }
                    if (isOnBorder) {
                        possibleRocks[iRock].shader->setVector("instanceOffset", pos);
                        possibleRocks[iRock].shader->setFloat("sizeFactor", size);
                        possibleRocks[iRock].display();
                    }
                }
            }
            Matrix3<Vector3> gradientTerrain = values.gradient();
            gradientTerrain.raiseErrorOnBadCoord = false;
            gradientTerrain.defaultValueOnBadCoord = Vector3(0, 0, -1);
            for (size_t i = 0; i < coralsIndicesAndPositionAndSize.size(); i++) {
                int iCoral;
                Vector3 pos;
                float size;
                std::tie(iCoral, pos, size) = coralsIndicesAndPositionAndSize[i];
                // Must be on ground, and have water somewhere around
                if (values.at(pos) > 0) {
                    bool isOnBorder = false;
                    for (int dx = -1; dx <= 1 && !isOnBorder; dx++)
                        for (int dy = -1; dy <= 1 && !isOnBorder; dy++)
                            for (int dz = -1; dz <= 1 && !isOnBorder; dz++)
                                if (values.at(pos + Vector3(dx, dy, dz)) < 0) {
                                    isOnBorder = true;
                                }
                    if (isOnBorder) {
                        Vector3 up = -gradientTerrain.at(pos).normalize();
                        Vector3 fwd = (up == Vector3(0, 1, 0) ? Vector3(1, 0, 0) : up.cross(Vector3(0, 1, 0)));
                        Vector3 right = up.cross(fwd);/*
                        std::vector<std::vector<float>> rotationMatrix({
                                                  {fwd.x, right.x, up.x, 0},
                                                  {fwd.y, right.y, up.y, 0},
                                                  {fwd.z, right.z, up.z, 0},
                                                  {0, 0, 0, 1}
                                              });*/ /*
                        std::vector<std::vector<float>> rotationMatrix({
                                                  {fwd.x, fwd.y, fwd.z, 0},
                                                  {right.x, right.y, right.z, 0},
                                                  {up.x, up.y, up.z, 0},
                                                  {0, 0, 0, 1}
                                              });
                        possibleCorals[iCoral].shader->setVector("instanceOffset", pos);
                        possibleCorals[iCoral].shader->setMatrix("instanceRotation", (std::vector<std::vector<float>>)rotationMatrix);
                        possibleCorals[iCoral].shader->setFloat("sizeFactor", size);
                        possibleCorals[iCoral].display();
                    }
                }
            }
        }*/
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

    if (verbose)
        std::cout << "Loading 3D models... " << std::flush;

    std::vector<QString> coralPaths;
    std::vector<QString> rocksPaths;
//    std::vector<QString> algaePaths;

    auto startTime = std::chrono::system_clock::now();

    QDirIterator itCorals("src/assets/models/corals/", QDir::Files, QDirIterator::Subdirectories);
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

    QDirIterator itRocks("src/assets/models/rocks/", QDir::Files, QDirIterator::Subdirectories);
    std::shared_ptr<Shader> rocksShader = std::make_shared<Shader>(vRockShader, fRockShader);
    while (itRocks.hasNext()) {
        QString dir = itRocks.next();
        rocksPaths.push_back(dir);
    }
    size_t nbRocks = rocksPaths.size();
    if (this->numberOfLoadedRocks != -1) nbRocks = std::min(nbRocks, (size_t)numberOfLoadedRocks);
    this->possibleRocks = std::vector<Mesh>(nbRocks);
#pragma omp parallel for
    for (size_t i = 0; i < nbCorals; i++) {
        QString& dir = rocksPaths[i];
        // Normalize it and move it upward so the anchor is on the ground
        possibleRocks[i] = Mesh(rocksShader).fromStl(dir.toStdString()).normalize();
    }
    auto endTime = std::chrono::system_clock::now();
    if (verbose)
        std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl;

    /*
    {
        auto startTime = std::chrono::system_clock::now();
        QTemporaryDir tempDirCorals;
        QDirIterator itCorals(":/models_3d/corals/", QDir::Files, QDirIterator::Subdirectories);
        this->possibleCorals = std::vector<Mesh>();
        std::shared_ptr<Shader> coralsShader = std::make_shared<Shader>(vRockShader, fRockShader);
        while (itCorals.hasNext()) {
            QString dir = itCorals.next();
            QString newPath = tempDirCorals.path() + QFileInfo(dir).fileName();
            QFile::copy(dir, newPath);
            Mesh m = Mesh(coralsShader);
            // Normalize it and move it upward so the anchor is on the ground
            possibleCorals.push_back(m.fromStl(newPath.toStdString()).normalize().rotate(deg2rad(180), 0, 0)); // .translate(0.f, 0.f, .25f));
        }
        QTemporaryDir tempDirRocks;
        QDirIterator itRocks(":/models_3d/rocks/", QDir::Files, QDirIterator::Subdirectories);
        this->possibleRocks = std::vector<Mesh>();
        std::shared_ptr<Shader> rocksShader = std::make_shared<Shader>(vRockShader, fRockShader);
        while (itRocks.hasNext()) {
            QString dir = itRocks.next();
            QString newPath = tempDirRocks.path() + QFileInfo(dir).fileName();
            QFile::copy(dir, newPath);
            Mesh m = Mesh(rocksShader);
            possibleRocks.push_back(m.fromStl(newPath.toStdString()).normalize());
        }

        auto endTime = std::chrono::system_clock::now();
        std::cout << "Loaded models in " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms (sequential)" << std::endl;
    }
    */
}

QLayout* MeshInstanceAmplificationInterface::createGUI()
{
    QVBoxLayout* layout = new QVBoxLayout();

    QCheckBox* displayCoralsCheckbox = new QCheckBox("Display corals");
    QCheckBox* displayRocksCheckbox = new QCheckBox("Display rocks");

    layout->addWidget(displayCoralsCheckbox);
    layout->addWidget(displayRocksCheckbox);

    displayCoralsCheckbox->setChecked(this->displayCorals);
    displayRocksCheckbox->setChecked(this->displayRocks);

    QObject::connect(displayCoralsCheckbox, &QCheckBox::toggled, this, &MeshInstanceAmplificationInterface::setCoralsDisplayed);
    QObject::connect(displayRocksCheckbox, &QCheckBox::toggled, this, &MeshInstanceAmplificationInterface::setRocksDisplayed);

    return layout;
}

std::vector<std::pair<Vector3, Vector3> > MeshInstanceAmplificationInterface::getAvailablePositionsForMaterial(TerrainTypes target)
{
    std::vector<std::pair<Vector3, Vector3> > positions;
    auto materialsAndHeightsGrid = layerGrid->layers;
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

std::vector<std::pair<Vector3, Vector3> > MeshInstanceAmplificationInterface::getCoralAvailablePositions()
{
    auto intervals = getAvailablePositionsForMaterial(TerrainTypes::CORAL);
    std::vector<std::pair<Vector3, Vector3> > extendedPositions;

    for (auto& interv : intervals) {
        Vector3& start = interv.first;
        Vector3& end = interv.second;
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

std::vector<std::pair<Vector3, Vector3> > MeshInstanceAmplificationInterface::getRocksAvailablePositions()
{
//    return getAvailablePositionsForMaterial(TerrainTypes::ROCK);
    auto intervals = getAvailablePositionsForMaterial(TerrainTypes::ROCK);
    std::vector<std::pair<Vector3, Vector3> > extendedPositions;

    for (auto& interv : intervals) {
        Vector3& start = interv.first;
        Vector3& end = interv.second;
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

void MeshInstanceAmplificationInterface::setCoralsDisplayed(bool display)
{
    this->displayCorals = display;
}

void MeshInstanceAmplificationInterface::setRocksDisplayed(bool display)
{
    this->displayRocks = display;
}

void MeshInstanceAmplificationInterface::afterTerrainUpdated()
{
    this->regenerateRocksPositions();
}

void MeshInstanceAmplificationInterface::regenerateRocksPositions()
{
    auto coralAvailablePositions = this->getCoralAvailablePositions();
    auto rocksAvailablePositions = this->getRocksAvailablePositions();

    std::shuffle(coralAvailablePositions.begin(), coralAvailablePositions.end(), random_gen::random_generator);
    std::shuffle(rocksAvailablePositions.begin(), rocksAvailablePositions.end(), random_gen::random_generator);

    this->rocksIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(std::min(numberOfRocksDisplayed, int(rocksAvailablePositions.size())));
    for (size_t i = 0; i < rocksIndicesAndPositionAndSize.size(); i++) {
        rocksIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleRocks.size())),
                                              Vector3(rocksAvailablePositions[i].first),
                                                //Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                              random_gen::generate(5.f, 15.f)
                                                     );
    }
    this->coralsIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(std::min(numberOfCoralsDisplayed, int(coralAvailablePositions.size())));
    for (size_t i = 0; i < coralsIndicesAndPositionAndSize.size(); i++) {
        coralsIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleCorals.size())),
                                                Vector3(coralAvailablePositions[i].first),
                                              // Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                              random_gen::generate(5.f, 15.f)
                                                     );
    }
}
