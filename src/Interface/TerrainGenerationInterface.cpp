#include "TerrainGenerationInterface.h"

#include "Utils/Utils.h"
#include "Utils/stb_image_write.h"

TerrainGenerationInterface::TerrainGenerationInterface(QWidget *parent) : ActionInterface("terrain_gen", parent)
{
    this->createGUI();
}

void TerrainGenerationInterface::prepareShader()
{
    const char* vShader_mc_voxels = ":/src/Shaders/MarchingCubes.vert";
    const char* gShader_mc_voxels = ":/src/Shaders/MarchingCubes.geom";
    const char* fShader_mc_voxels = ":/src/Shaders/MarchingCubes.frag";

    const char* vShader_grid = ":/src/Shaders/grid.vert";
    const char* fShader_grid = ":/src/Shaders/grid.frag";

    const char* vRockShader = ":/src/Shaders/rockShader.vert";
    const char* fRockShader = ":/src/Shaders/rockShader.frag";

    const char* vParticleShader = ":/src/Shaders/particle.vert";
    const char* fParticleShader = ":/src/Shaders/particle.frag";
    const char* gParticleShader = ":/src/Shaders/particle.geom";

    marchingCubeMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_mc_voxels),
                            true, GL_TRIANGLES);
    Matrix3<float> isoData = this->voxelGrid->getVoxelValues(); //.resize(15, 15, 15);
    std::vector<Vector3> points(isoData.size());
    for (size_t i = 0; i < points.size(); i++) {
        isoData[i] = (std::max(-3.f, std::min(3.f, isoData[i])) / 6.f) + 0.5;
        points[i] = isoData.getCoordAsVector3(i);
    }
    marchingCubeMesh.useIndices = false;
    marchingCubeMesh.fromArray(points);
    marchingCubeMesh.update();
    marchingCubeMesh.shader->use();
    GlobalsGL::f()->glBindVertexArray(marchingCubeMesh.vao);
    marchingCubeMesh.shader->setInt("dataFieldTex", 0);
    marchingCubeMesh.shader->setInt("edgeTableTex", 1);
    marchingCubeMesh.shader->setInt("triTableTex", 2);
    marchingCubeMesh.shader->setFloat("isolevel", 0.f);
    marchingCubeMesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
    //Edge Table texture//
    //This texture store the 256 different configurations of a marching cube.
    //This is a table accessed with a bitfield of the 8 cube edges states
    //(edge cut by isosurface or totally in or out).
    //(cf. MarchingCubes.cpp)

    GlobalsGL::f()->glGenTextures(1, &edgeTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, edgeTableTex);
    //Integer textures must use nearest filtering mode

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 256, 1, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::cubeEdges));

    //Triangle Table texture//
    //This texture store the vertex index list for
    //generating the triangles of each configurations.
    //(cf. MarchingCubes.cpp)

    glGenTextures(1, &triTableTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE2);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, triTableTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, 16, 256, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, &(MarchingCubes::triangleTable));
    marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, voxelGrid->getVoxelValues() / 6.f + .5f);

    QTemporaryDir tempDir;
    QDirIterator it(":/models_3d/", QDir::Files, QDirIterator::Subdirectories);
    this->possibleRocks = std::vector<Mesh>();
    std::shared_ptr<Shader> rocksShader = std::make_shared<Shader>(vRockShader, fRockShader);
    while (it.hasNext()) {
        QString dir = it.next();
        QString newPath = tempDir.path() + QFileInfo(dir).fileName();
        QFile::copy(dir, newPath);
        Mesh m = Mesh(rocksShader);
        possibleRocks.push_back(m.fromStl(newPath.toStdString()).normalize());
    }
    this->regenerateRocksAndParticles();
    this->particlesMesh = Mesh(this->randomParticlesPositions,
                               std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader),
                               true, GL_POINTS);
    particlesMesh.shader->setVector("color", std::vector<float>({1.f, 1.f, 1.f, .2f}));
    this->particlesMesh.shader->setFloat("maxTerrainHeight", voxelGrid->getSizeZ());

    this->heightmapGrid->mesh.shader = std::make_shared<Shader>(vShader_grid, fShader_grid);
}

void TerrainGenerationInterface::regenerateRocksAndParticles()
{
    Matrix3<float> values = this->voxelGrid->getVoxelValues();

    std::vector<Vector3> points(values.size());
    for (size_t i = 0; i < values.size(); i++)
        points[i] = values.getCoordAsVector3(i);
    marchingCubeMesh.fromArray(points);
    this->rocksIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(numberOfRocksDisplayed);
    for (size_t i = 0; i < rocksIndicesAndPositionAndSize.size(); i++) {
        rocksIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleRocks.size())),
                                              Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ)),
                                              random_gen::generate(5.f, 15.f)
                                                     );
    }
    this->randomParticlesPositions = std::vector<Vector3>(numberOfRocksDisplayed * 10);
    for (size_t i = 0; i < randomParticlesPositions.size(); i++) {
        randomParticlesPositions[i] = Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ));
    }
    if (this->particlesMesh.shader != nullptr)
        this->particlesMesh.shader->setFloat("maxTerrainHeight", values.sizeZ);
    this->startingTime = std::chrono::system_clock::now();
}


void TerrainGenerationInterface::display(MapMode mapMode, SmoothingAlgorithm smoothingAlgorithm, bool displayParticles)
{
    if (mapMode == GRID_MODE) {
        if (this->heightmapGrid == nullptr) {
            std::cerr << "No grid to display" << std::endl;
        } else {
            float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();
            this->heightmapGrid->mesh.shader->setFloat("time", time);
            this->heightmapGrid->display(true);
        }
    }
    else if (mapMode == VOXEL_MODE) {
        if (this->voxelGrid == nullptr) {
            std::cerr << "No voxel grid to display" << std::endl;
        } else {
            marchingCubeMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
            Matrix3<float> values = voxelGrid->getVoxelValues();
            // Not the best check, but still pretty good....
            if (marchingCubeMesh.vertexArray.size() != values.size()) {
                regenerateRocksAndParticles();
            }
            marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values / 6.f + .5f);
            marchingCubeMesh.display( GL_POINTS );

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
            if (displayParticles || true) {
//                particlesMesh.fromArray(randomParticlesPositions); // Maybe we shouldn't displace the particles
                float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();
                particlesMesh.shader->setFloat("time", time); // std::chrono::system_clock::now().time_since_epoch().count());
                particlesMesh.display(GL_POINTS);
            }
        }
    }/*
    else if (mapMode == LAYER_MODE) {
        if (this->layerGrid == nullptr) {
            std::cerr << "No layer based grid to display" << std::endl;
        } else {
            this->layerGrid->display();
        }
    }*/
    /*
    */
}

void TerrainGenerationInterface::createTerrainFromNoise(int nx, int ny, int nz, float blockSize, float noise_shifting)
{
    std::shared_ptr<VoxelGrid> tempMap = std::make_shared<VoxelGrid>(nx, ny, nz, blockSize, noise_shifting);
    tempMap->fromIsoData();

    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->heightmapGrid)
        this->heightmapGrid = std::make_shared<Grid>();

    // Okay, super dirty but I don't know how to manage shared_ptr...
    this->voxelGrid->sizeX = tempMap->sizeX;
    this->voxelGrid->sizeY = tempMap->sizeY;
    this->voxelGrid->sizeZ = tempMap->sizeZ;
    this->voxelGrid->blockSize = tempMap->blockSize;
    this->voxelGrid->chunkSize = tempMap->chunkSize;
    this->voxelGrid->noise_shifting = tempMap->noise_shifting;
    this->voxelGrid->initMap();
    this->voxelGrid->tempData = tempMap->tempData;
    this->voxelGrid->fromIsoData();
    this->heightmapGrid->fromVoxelGrid(*voxelGrid);

    this->addTerrainAction(nlohmann::json({
                                              {"from_noise", true},
                                              {"noise_parameters", {
                                                   {"nx", nx},
                                                   {"ny", ny},
                                                   {"nz", nz},
                                                   {"block_size", blockSize},
                                                   {"noise_shifting", noise_shifting}
                                               }}
                                          }));
}

void TerrainGenerationInterface::createTerrainFromFile(std::string filename, std::vector<std::shared_ptr<ActionInterface>> actionInterfaces)
{
    std::string ext = toUpper(getExtention(filename));
    if (!this->heightmapGrid)
        this->heightmapGrid = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    if (ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        heightmapGrid->loadFromHeightmap(filename, 127, 127, 50);
        voxelGrid->from2DGrid(*heightmapGrid);
        voxelGrid->fromIsoData();
    } else if (ext == "JSON") {
        // The JSON file contains the list of actions made on a map
        std::ifstream file(filename);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json json_content = nlohmann::json::parse(content);
        if (!json_content.contains("actions"))
            return;
        for (auto action : json_content.at("actions")) {
            // Let all the interfaces try to replay their actions
            for (auto& possibleAction : actionInterfaces)
                possibleAction->replay(action);
        }
    } else {
        // Then it's our custom voxel grid file
        voxelGrid->retrieveMap(filename);
        voxelGrid->fromIsoData();
        heightmapGrid->fromVoxelGrid(*voxelGrid);
    }
    this->voxelGrid->createMesh();
    this->heightmapGrid->createMesh();

    this->addTerrainAction(nlohmann::json({
                                              {"from_file", filename},
                                              {"from_noise", false}
                                          }));
}

void TerrainGenerationInterface::saveTerrain(std::string filename)
{
    std::string ext = toUpper(getExtention(filename));
    if (ext == "PNG" || ext == "JPG" || ext == "TGA" || ext == "BMP" || ext == "HDR") {
        // To heightmap
        this->heightmapGrid->fromVoxelGrid(*voxelGrid); // Just to be sure to have the last values
        this->heightmapGrid->saveHeightmap(filename);

        /*std::vector<float> toFloatData(heightmapGrid->getSizeX() * heightmapGrid->getSizeY());
        std::vector<int> toIntData(heightmapGrid->getSizeX() * heightmapGrid->getSizeY());

        for (size_t i = 0; i < heightmapGrid->vertices.size(); i++) {
            toFloatData[i] = heightmapGrid->vertices[i].z / heightmapGrid->maxHeight;
            toIntData[i] = toFloatData[i] * 255;
        }
        if (ext == "PNG")
            stbi_write_png(filename.c_str(), heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 1, toIntData.data(), sizeof(int));
        else if (ext == "JPG")
            stbi_write_jpg(filename.c_str(), heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 1, toIntData.data(), 95);
        else if (ext == "BMP")
            stbi_write_bmp(filename.c_str(), heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 1, toIntData.data());
        else if (ext == "TGA")
            stbi_write_tga(filename.c_str(), heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 1, toIntData.data());
        else if (ext == "HDR")
            stbi_write_hdr(filename.c_str(), heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 1, toFloatData.data());*/
    } else if (ext == "JSON") {
        // To JSON file containing the list of actions made on a map
        this->saveAllActions(filename);
    } else {
        // Otherwise it's our custom voxel grid file
        voxelGrid->saveMap(filename);
    }
}


void TerrainGenerationInterface::replay(nlohmann::json action)
{
    if (this->isConcerned(action)) {
        auto parameters = action.at("parameters");
        bool loadFromNoise = parameters.at("from_noise").get<bool>();
        if (loadFromNoise) {
            parameters = parameters.at("noise_parameters");
            int nx = parameters.at("nx").get<int>();
            int ny = parameters.at("ny").get<int>();
            int nz = parameters.at("nz").get<int>();
            int blockSize = parameters.at("block_size").get<float>();
            int noise_shifting = parameters.at("noise_shifting").get<float>();
            this->createTerrainFromNoise(nx, ny, nz, blockSize, noise_shifting);
        } else {
            std::string filename = parameters.at("from_file").get<std::string>();
            this->createTerrainFromFile(filename);
        }
    }
}

void TerrainGenerationInterface::hide()
{
    CustomInteractiveObject::hide();
}

void TerrainGenerationInterface::show()
{
    CustomInteractiveObject::show();
}

QLayout* TerrainGenerationInterface::createGUI()
{
    QLayout* nothing = new QHBoxLayout;
    return nothing;
}
