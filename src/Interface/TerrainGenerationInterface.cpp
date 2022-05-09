#include "TerrainGenerationInterface.h"

#include "Utils/Utils.h"

TerrainGenerationInterface::TerrainGenerationInterface(QWidget *parent) : ActionInterface("terrain_gen", parent)
{
    this->createGUI();
}

void TerrainGenerationInterface::prepareShader()
{
    const char* vShader_mc_voxels = ":/src/Shaders/MarchingCubes_vertex.glsl";
    const char* gShader_mc_voxels = ":/src/Shaders/MarchingCubes_geometry.glsl";
    const char* fShader_mc_voxels = ":/src/Shaders/MarchingCubes_fragment.glsl";

    const char* vRockShader = ":/src/Shaders/rockShader.vert";
    const char* fRockShader = ":/src/Shaders/rockShader.frag";

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
        possibleRocks.push_back(m.fromStl(newPath.toStdString()).normalize().scale(10.f));
        break;
    }
    this->rocksIndicesAndPosition = std::vector<std::tuple<int, Vector3>>();
    for (int i = 0; i < numberOfRocksDisplayed; i++) {
        rocksIndicesAndPosition.push_back(std::make_tuple<int, Vector3>(
                                              int(random_gen::generate(0, possibleRocks.size())),
                                              Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ))
                                              ));
    }
}


void TerrainGenerationInterface::display()
{
    Matrix3<float> values = voxelGrid->getVoxelValues();
    // Not the best check, but still pretty good....
    if (marchingCubeMesh.vertexArray.size() != values.size()) {
        std::vector<Vector3> points(values.size());
        for (size_t i = 0; i < values.size(); i++)
            points[i] = values.getCoordAsVector3(i);
        marchingCubeMesh.fromArray(points);
        this->rocksIndicesAndPosition = std::vector<std::tuple<int, Vector3>>();
        for (int i = 0; i < numberOfRocksDisplayed; i++) {
            rocksIndicesAndPosition.push_back(std::make_tuple<int, Vector3>(
                                                  int(random_gen::generate(0, possibleRocks.size())),
                                                  Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ))
                                                  ));
        }
    }
    marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values / 6.f + .5f);
    marchingCubeMesh.display( GL_POINTS );

    for (size_t i = 0; i < rocksIndicesAndPosition.size(); i++) {
        int iRock;
        Vector3 pos;
        std::tie(iRock, pos) = rocksIndicesAndPosition[i];
        if (this->voxelGrid->getVoxelValue(pos) > 0) {
            possibleRocks[iRock].shader->setVector("instanceOffset", pos);
            possibleRocks[iRock].display();
        }
    }
}

void TerrainGenerationInterface::createTerrainFromNoise(int nx, int ny, int nz, float blockSize, float noise_shifting)
{
    this->voxelGrid = std::make_shared<VoxelGrid>(nx, ny, nz, blockSize, noise_shifting);
    this->voxelGrid->fromIsoData();
    this->heightmapGrid = std::make_shared<Grid>();
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
    const char* vShader_voxels = ":/src/Shaders/voxels_vertex_shader_blinn_phong.glsl";
    const char* fShader_voxels = ":/src/Shaders/voxels_fragment_shader_blinn_phong.glsl";

    std::string ext = toUpper(getExtention(filename));
    if (!this->heightmapGrid)
        this->heightmapGrid = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    if (ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        heightmapGrid->loadFromHeightmap(filename, 127, 127, 255);
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
    this->voxelGrid->createMesh();/*
    this->heightmapGrid->createMesh();
    for(std::shared_ptr<VoxelChunk>& vc : this->voxelGrid->chunks)
        vc->mesh.shader = std::make_shared<Shader>(vShader_voxels, fShader_voxels);
*/
    this->addTerrainAction(nlohmann::json({
                                              {"from_file", filename},
                                              {"from_noise", false}
                                          }));
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
