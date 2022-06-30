#include "TerrainGenerationInterface.h"

#include "Biomes/BiomeInstance.h"
#include "Utils/ConstraintsSolver.h"
#include "Utils/ShapeCurve.h"
#include "Utils/Utils.h"
#include "Utils/Voronoi.h"
#include "Utils/stb_image.h"
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
    const char* gShader_grid = ":/src/Shaders/grid.geom";
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

    this->heightmapMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_grid), true, GL_POINTS);

    Matrix3<float> heightData(this->heightmapGrid->getSizeX(), this->heightmapGrid->getSizeY());
    std::vector<Vector3> positions(heightData.size());
    for (size_t i = 0; i < positions.size(); i++) {
        positions[i] = heightData.getCoordAsVector3(i);
        heightData[i] = heightmapGrid->getHeight(positions[i].x, positions[i].y);
    }
    heightmapMesh.useIndices = false;
    heightmapMesh.fromArray(positions);
    heightmapMesh.update();
    heightmapMesh.shader->use();
    GlobalsGL::f()->glBindVertexArray(heightmapMesh.vao);
    heightmapMesh.shader->setInt("heightmapFieldTex", 3);
    marchingCubeMesh.shader->setInt("heightmapFieldTex", 3);

    float *heightmapData = new float[heightmapGrid->heights.size() * 4];
    Matrix3<Vector3> gradients = heightmapGrid->heights.gradient();
    for (size_t i = 0; i < heightmapGrid->heights.size(); i++) {
        gradients[i] = (gradients[i].normalized() + Vector3(1, 1, 1)) / 2.f;
        heightmapData[i * 4 + 0] = gradients[i].x;
        heightmapData[i * 4 + 1] = gradients[i].y;
        heightmapData[i * 4 + 2] = gradients[i].z;
        heightmapData[i * 4 + 3] = heightmapGrid->heights[i];
    }

    glGenTextures(1, &heightmapFieldTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
//            glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 0,
//            GL_RED, GL_FLOAT, heightmapGrid->heights.data.data());//heightData.data.data());
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 0,
    GL_RGBA, GL_FLOAT, heightmapData);//heightData.data.data());
    delete[] heightmapData;


    glGenTextures(1, &biomeFieldTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE4);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, biomeFieldTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmapGrid->biomeIndices.sizeX, heightmapGrid->biomeIndices.sizeY, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, heightmapGrid->biomeIndices.data.data());
    heightmapMesh.shader->setInt("biomeFieldTex", 4);
    marchingCubeMesh.shader->setInt("biomeFieldTex", 4);

    QTemporaryDir tempTextureDir;
    QDirIterator iTex(":/terrain_textures/", QDir::Files, QDirIterator::Subdirectories);
    Matrix3<Vector3> allColorTextures;
    Matrix3<Vector3> allNormalTextures;
    Matrix3<Vector3> allDisplacementTextures;
    int indexColorTextureClass = 0;
    int indexNormalTextureClass = 0;
    int indexDisplacementTextureClass = 0;
    while (iTex.hasNext()) {
        QString dir = iTex.next();
        QString basename = QFileInfo(dir).baseName();
        if (basename == "color" || basename == "normal" || basename == "displacement") {
            std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
            QString newPath = tempTextureDir.path() + QString::fromStdString(textureClass) + QFileInfo(dir).fileName();
            QFile::copy(dir, newPath);
            int imgW, imgH, c;
            unsigned char *data = stbi_load(newPath.toStdString().c_str(), &imgW, &imgH, &c, 3);
            std::cout << "Loading " << newPath.toStdString().c_str() << " for the class " << textureClass << "..." << std::endl;
            Matrix3<Vector3> map(imgW, imgH);
            for (int x = 0; x < imgW; x++) {
                for (int y = 0; y < imgH; y++) {
                    unsigned char r = data[3 * (x + y * imgW) + 0];
                    unsigned char g = data[3 * (x + y * imgW) + 1];
                    unsigned char b = data[3 * (x + y * imgW) + 2];
                    map.at(x, y) = Vector3(r, g, b);
                }
            }
            map = map.resize(512, 512, 1, RESIZE_MODE::LINEAR);
            stbi_image_free(data);
            if (basename == "color") {
                colorTexturesIndex[textureClass] = indexColorTextureClass++;
                if (allColorTextures.empty()) {
                    allColorTextures = map;
                } else {
                    allColorTextures = allColorTextures.concat(map);
                }
            } else if (basename == "normal") {
                normalTexturesIndex[textureClass] = indexNormalTextureClass++;
                if (allNormalTextures.empty()) {
                    allNormalTextures = map;
                } else {
                    allNormalTextures = allNormalTextures.concat(map);
                }
            } else if (basename == "displacement") {
                displacementTexturesIndex[textureClass] = indexDisplacementTextureClass++;
                if (allDisplacementTextures.empty()) {
                    allDisplacementTextures = map;
                } else {
                    allDisplacementTextures = allDisplacementTextures.concat(map);
                }
            }
        }
    }

    allColorTextures /= 255.f;
    allNormalTextures /= 255.f;
    allDisplacementTextures /= 255.f;
    float *allTexturesColors = new float[allColorTextures.size() * 4];
    for (size_t i = 0; i < allColorTextures.size(); i++) {
        allTexturesColors[4 * i + 0] = allColorTextures[i].x;
        allTexturesColors[4 * i + 1] = allColorTextures[i].y;
        allTexturesColors[4 * i + 2] = allColorTextures[i].z;
        allTexturesColors[4 * i + 3] = 1.f;
    }
    float *allTexturesNormal = new float[allNormalTextures.size() * 4];
    for (size_t i = 0; i < allNormalTextures.size(); i++) {
        allTexturesNormal[4 * i + 0] = allNormalTextures[i].x;
        allTexturesNormal[4 * i + 1] = allNormalTextures[i].y;
        allTexturesNormal[4 * i + 2] = allNormalTextures[i].z;
        allTexturesNormal[4 * i + 3] = 1.f;
    }
    float *allTexturesDisplacement = new float[allDisplacementTextures.size() * 4];
    for (size_t i = 0; i < allDisplacementTextures.size(); i++) {
        allTexturesDisplacement[4 * i + 0] = allDisplacementTextures[i].x;
        allTexturesDisplacement[4 * i + 1] = allDisplacementTextures[i].y;
        allTexturesDisplacement[4 * i + 2] = allDisplacementTextures[i].z;
        allTexturesDisplacement[4 * i + 3] = 1.f;
    }

    glGenTextures(1, &allBiomesColorTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE5);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, allBiomesColorTextures);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allColorTextures.sizeX, allColorTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesColors);
    delete[] allTexturesColors;
    heightmapMesh.shader->setInt("allBiomesColorTextures", 5);
    marchingCubeMesh.shader->setInt("allBiomesColorTextures", 5);

    heightmapMesh.shader->setInt("maxBiomesColorTextures", indexColorTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesColorTextures", indexColorTextureClass);

    glGenTextures(1, &allBiomesNormalTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE6);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, allBiomesNormalTextures);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allNormalTextures.sizeX, allNormalTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesNormal);
    delete[] allTexturesNormal;
    heightmapMesh.shader->setInt("allBiomesNormalTextures", 6);
    marchingCubeMesh.shader->setInt("allBiomesNormalTextures", 6);

    heightmapMesh.shader->setInt("maxBiomesNormalTextures", indexNormalTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesNormalTextures", indexNormalTextureClass);

    glGenTextures(1, &allBiomesDisplacementTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE7);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, allBiomesDisplacementTextures);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allDisplacementTextures.sizeX, allDisplacementTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesDisplacement);
    delete[] allTexturesDisplacement;
    heightmapMesh.shader->setInt("allBiomesDisplacementTextures", 7);
    marchingCubeMesh.shader->setInt("allBiomesDisplacementTextures", 7);

    heightmapMesh.shader->setInt("maxBiomesDisplacementTextures", indexDisplacementTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesDisplacementTextures", indexDisplacementTextureClass);

//    this->heightmapGrid->mesh.shader = std::make_shared<Shader>(vShader_grid, fShader_grid);
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
    /*if (this->heightmapGrid != nullptr) {
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, biomeFieldTex);
        glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmapGrid->biomeIndices.sizeX, heightmapGrid->biomeIndices.sizeY, 0,
        GL_ALPHA_INTEGER_EXT, GL_INT, ((Matrix3<float>)(heightmapGrid->biomeIndices) / 10.f).data.data());
    }*/
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, allBiomesColorTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE6);
    glBindTexture(GL_TEXTURE_2D, allBiomesNormalTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_2D, allBiomesDisplacementTextures);

    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);
//            glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 0,
//            GL_RED, GL_FLOAT, heightmapGrid->heights.data.data());//heightData.data.data());
    float maxHeight = heightmapGrid->heights.max();
    this->heightmapMesh.shader->setFloat("maxHeight", maxHeight);
    float *heightmapData = new float[heightmapGrid->heights.size() * 4];
    Matrix3<Vector3> gradients = heightmapGrid->heights.gradient();
    int maxColorTextureIndex = 0;
    for (auto biome : colorTexturesIndex) {
        maxColorTextureIndex = std::max(maxColorTextureIndex, biome.second + 1);
    }
    int maxNormalTextureIndex = 0;
    for (auto biome : normalTexturesIndex) {
        maxNormalTextureIndex = std::max(maxNormalTextureIndex, biome.second + 1);
    }
    int maxDisplacementTextureIndex = 0;
    for (auto biome : displacementTexturesIndex) {
        maxDisplacementTextureIndex = std::max(maxDisplacementTextureIndex, biome.second + 1);
    }
    int maxBiomeID = 0; //std::max(1, heightmapGrid->biomeIndices.max());
    for (const auto& biomeValues : heightmapGrid->biomeIndices)
        maxBiomeID = std::max(maxBiomeID, (biomeValues.empty() ? 1 : biomeValues.back()));
    Matrix3<std::vector<int>> resizedBiomeIndices = heightmapGrid->biomeIndices; //.resize(heightmapGrid->heights.getDimensions(), RESIZE_MODE::MAX_VAL);
    for (size_t i = 0; i < heightmapGrid->heights.size(); i++) {
//        int biomeID = (resizedBiomeIndices[i].empty() ? 0 : resizedBiomeIndices[i].back());
        float colorTextureOffset = 1.f;
        float normalTextureOffset = 1.f;
        float displacementTextureOffset = 1.f;
        if (resizedBiomeIndices.getDimensions() == heightmapGrid->heights.getDimensions()) {
            auto biome = (!resizedBiomeIndices[i].empty() ? BiomeInstance::instancedBiomes[resizedBiomeIndices[i].back()] : nullptr);
            if (biome != nullptr && !biome->getTextureName().empty()) {
                if (colorTexturesIndex.find(biome->getTextureName()) != colorTexturesIndex.end())
                    colorTextureOffset = colorTexturesIndex[biome->getTextureName()] / (float)(maxColorTextureIndex);
                if (normalTexturesIndex.find(biome->getTextureName()) != normalTexturesIndex.end())
                    normalTextureOffset = normalTexturesIndex[biome->getTextureName()] / (float)(maxNormalTextureIndex);
                if (displacementTexturesIndex.find(biome->getTextureName()) != displacementTexturesIndex.end())
                    displacementTextureOffset = displacementTexturesIndex[biome->getTextureName()] / (float)(maxDisplacementTextureIndex);
            }
        }
//        Vector3 color = HSVtoRGB(biomeID / (float)maxBiomeID, 1, 1);
        gradients[i] = (gradients[i].normalized().abs());// + Vector3(1, 1, 1)) / 2.f;
        heightmapData[i * 4 + 0] = colorTextureOffset;
        heightmapData[i * 4 + 1] = normalTextureOffset; // gradients[i].y;
        heightmapData[i * 4 + 2] = displacementTextureOffset; // gradients[i].z;
        heightmapData[i * 4 + 3] = heightmapGrid->heights[i] / maxHeight;
    }
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmapGrid->getSizeX(), heightmapGrid->getSizeY(), 0,
    GL_RGBA, GL_FLOAT, heightmapData);//heightData.data.data());
    delete[] heightmapData;

    if (mapMode == GRID_MODE) {
        if (this->heightmapGrid == nullptr) {
            std::cerr << "No grid to display" << std::endl;
        } else {
//            float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();

            this->heightmapMesh.display(GL_POINTS);
//            this->heightmapGrid->mesh.shader->setFloat("time", time);
//            this->heightmapGrid->display(true);
        }
    }
    else if (mapMode == VOXEL_MODE) {
        if (this->voxelGrid == nullptr) {
            std::cerr << "No voxel grid to display" << std::endl;
        } else {
            Matrix3<float> values = voxelGrid->getVoxelValues();
            // Not the best check, but still pretty good....
            if (marchingCubeMesh.vertexArray.size() != values.size()) {
                regenerateRocksAndParticles();
            }
            marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values / 6.f + .5f);
            marchingCubeMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
            marchingCubeMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
            marchingCubeMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
            marchingCubeMesh.display( GL_POINTS );
            if (smoothingAlgorithm == SmoothingAlgorithm::NONE) {
                marchingCubeMesh.shader->setBool("displayingIgnoredVoxels", true);
                marchingCubeMesh.shader->setFloat("min_isolevel", -1000.f);
                marchingCubeMesh.shader->setFloat("max_isolevel", this->minIsoLevel/3.f);
                marchingCubeMesh.display( GL_POINTS );
                marchingCubeMesh.shader->setFloat("min_isolevel", this->maxIsoLevel/3.f);
                marchingCubeMesh.shader->setFloat("max_isolevel",  1000.f);
                marchingCubeMesh.display( GL_POINTS );
                marchingCubeMesh.shader->setBool("displayingIgnoredVoxels", false);
            }
/*
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
            }*/
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
        heightmapGrid->loadFromHeightmap(filename, 127, 127, 40);

        voxelGrid->from2DGrid(*heightmapGrid);
        voxelGrid->fromIsoData();

        ConstraintsSolver solver;
        int iArch1 = solver.addItem(new Vector3());
        int iArch2 = solver.addItem(new Vector3());
        int iAlgae1 = solver.addItem(new Vector3());
        int iAlgae2 = solver.addItem(new Vector3());

        solver.addDistanceConstraint(iArch1, iArch2, 30.f);
//        solver.addDistanceConstraint(iArch1, iAlgae1, 30.f);

//        solver.addDistanceConstraint(iArch2, iArch1, 30.f);
        solver.addDistanceConstraint(iAlgae1, iArch1, 30.f);
        solver.addDistanceConstraint(iAlgae2, iArch1, 40.f);
        solver.addDistanceConstraint(iAlgae2, iAlgae1, 70.f);

        solver.addNormalConstraint(iArch1, deg2rad(0), deg2rad(20)); // 0 to 45
        solver.addNormalConstraint(iArch2, deg2rad(20), deg2rad(90)); // 20 to 90
        solver.addNormalConstraint(iAlgae1, deg2rad(0), deg2rad(90)); // 0 to 90

        std::map<int, Vector3> positions = solver.solveWithVoxelGrid(voxelGrid);
        std::vector<Vector3> pos;
        for (auto& tuple : positions) {
            pos.push_back(tuple.second);
            std::cout << "Element " << tuple.first << " is placed at " << tuple.second << "\n";
        }
        std::cout << std::endl;
//        exit(0);

        UnderwaterErosion ue(voxelGrid, 10, 10, 10);
        ue.CreateTunnel(BSpline({pos[0], (pos[0] + pos[1]) / 2 + Vector3(0, 0, 1) * (pos[0] - pos[1]).norm(), pos[1]}), true, true);
        ue.CreateTunnel(BSpline({pos[2], pos[2] + Vector3(0, 0, 50)}), true, true);
        ue.CreateMultipleTunnels({BSpline({pos[3] + Vector3(-5, -5, 0), pos[3] + Vector3(5, 5)}), BSpline({pos[3] + Vector3(-5, 5), pos[3] + Vector3(5, -5)})}, true, true);

    } else if (ext == "JSON") {
        // The JSON file contains the list of actions made on a map
        std::ifstream file(filename);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json json_content = nlohmann::json::parse(content);
        if (!json_content.contains("actions")) {
            this->createTerrainFromBiomes(json_content);
        } else {
            for (auto action : json_content.at("actions")) {
                // Let all the interfaces try to replay their actions
                for (auto& possibleAction : actionInterfaces)
                    possibleAction->replay(action);
            }
        }
    } else if (ext == "DATA") {
        // Then it's our custom voxel grid file
        voxelGrid->retrieveMap(filename);
        voxelGrid->fromIsoData();
        heightmapGrid->fromVoxelGrid(*voxelGrid);
    } else {
        // In any other case, consider that nothing has been done, cancel.
        return;
    }
    this->voxelGrid->createMesh();
    this->heightmapGrid->createMesh();

    this->addTerrainAction(nlohmann::json({
                                              {"from_file", filename},
                                              {"from_noise", false}
                                          }));
}
/*
Vector3 getSurfacePosition(std::shared_ptr<VoxelGrid> grid, Vector3 pos) {
    pos.z = std::max(pos.z, 0.f); // In case of small imprecision
    while (grid->getVoxelValue(pos) > 0) {
        pos += Vector3(0, 0, 1); // Move the position one voxel heigher
        if (!grid->contains(pos)) { // If it gets too high, the whole column must be filled, I guess we should cancel it...
            pos.z = 0;
            break;
        }
    }
    return pos;
}
std::shared_ptr<BiomeInstance> recursivelyCreateBiome(nlohmann::json json_content, Vector3 biomePosition, ShapeCurve area) {
    std::string biomeClass = json_content.at("class").get<std::string>();
    // Should be able to retrieve the parameters of the biome...
    std::shared_ptr<BiomeInstance> instance = std::make_shared<BiomeInstance>(BiomeInstance::fromClass(biomeClass));
    instance->position = biomePosition;
    instance->area = area;
    auto children = json_content.at("children");
    Voronoi diagram(children.size(), area);
    std::vector<BSpline> subarea_borders = diagram.solve();
    for (size_t i = 0; i < children.size(); i++) {
        std::shared_ptr<BiomeInstance> childBiome = recursivelyCreateBiome(children[i], diagram.pointset[i], subarea_borders[i]);
        childBiome->parent = instance;
        instance->instances.push_back(childBiome);
    }
    return instance;
}
Matrix3<float> archeTunnel(BSpline path, float size, float strength, bool addingMatter, std::shared_ptr<VoxelGrid> grid) {
    Matrix3<float> erosionMatrix(grid->getDimensions());
    float nb_points_on_path = path.length() / (size/5.f);
    RockErosion rock(size, strength);
    for (const auto& pos : path.getPath(nb_points_on_path)) {
        erosionMatrix = rock.computeErosionMatrix(erosionMatrix, pos);
    }
    erosionMatrix = erosionMatrix.abs();
    erosionMatrix.toDistanceMap();
    erosionMatrix.normalize();
    for (float& m : erosionMatrix) {
        m = interpolation::linear(m, 0.f, 1.0) * strength * (addingMatter ? 1.f : -1.f);
//        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
    }
    return erosionMatrix;
}*/
void TerrainGenerationInterface::createTerrainFromBiomes(nlohmann::json json_content)
{
    if (!this->heightmapGrid)
        this->heightmapGrid = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    this->heightmapGrid->heights = Matrix3<float>(124, 124, 1, 20.f);
    this->heightmapGrid->maxHeight = 40.f;

    this->voxelGrid->from2DGrid(*this->heightmapGrid);
    this->voxelGrid->fromIsoData();

    this->biomeGenerationNeeded = true;
    this->biomeGenerationModelData = json_content;

#ifdef NONE
    ShapeCurve terrainArea({
                   heightmapGrid->heights.getDimensions() * Vector3(0, 0, 0),
                   heightmapGrid->heights.getDimensions() * Vector3(1, 0, 0),
                   heightmapGrid->heights.getDimensions() * Vector3(1, 1, 0),
                   heightmapGrid->heights.getDimensions() * Vector3(0, 1, 0)
               });
    Vector3 initialSpawn = heightmapGrid->heights.getDimensions() / 2.f + Vector3(10, 0, 0);
    // All the biome hierarchy is created from the json
    // Each biome has a class name, a position and his biome children
    std::shared_ptr<BiomeInstance> rootBiome = std::make_shared<BiomeInstance>(BiomeModel::fromJson(json_content).createInstance(initialSpawn, terrainArea));
    /*std::shared_ptr<BiomeInstance> rootBiome = recursivelyCreateBiome(json_content,
                                                                      initialSpawn,
                                                                      terrainArea);*/

    heightmapGrid->biomeIndices = Matrix3<int>(heightmapGrid->heights.getDimensions(), 0);
    std::vector<std::shared_ptr<BiomeInstance>> biomeQueue;
    std::vector<std::shared_ptr<BiomeInstance>> sortedBiomes;

    biomeQueue.push_back(rootBiome);

    int biomeID = 0;
    while (!biomeQueue.empty()) {
        std::shared_ptr<BiomeInstance> current = biomeQueue.front();
        sortedBiomes.push_back(current);
        biomeQueue.erase(biomeQueue.begin());
        if (!terrainArea.inside(current->position * Vector3(1, 1, 0), true))
            continue;

        biomeQueue.insert(biomeQueue.end(), current->instances.begin(), current->instances.end());
    }
    biomeQueue.push_back(rootBiome);

    // Sort the biomes by level, so that modifications are done top-to-bottom
    std::sort(sortedBiomes.begin(),
              sortedBiomes.end(),
              [](const auto& a, const auto& b) -> bool {
        return a->getLevel() < b->getLevel();
    });

    for (auto& current : sortedBiomes) {
        ShapeCurve area = current->area;
        int level = current->getLevel();
//        area = area.grow(-1 * level); // Shrink the area to be able to see all layers
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        for (int x = AABBoxMin.x; x < AABBoxMax.x; x++) {
            for (int y = AABBoxMin.y; y < AABBoxMax.y; y++) {
                if (area.inside(Vector3(x, y, area.points[0].z)) && heightmapGrid->biomeIndices.checkCoord(Vector3(x, y))) {
                    heightmapGrid->biomeIndices.at(x, y).push_back(biomeID);
                }
            }
        }
        current->instanceID = biomeID;
        BiomeInstance::registerBiome(current);
        biomeID++;

    }
    // First, level the terrain as the biomes design it
    Matrix3<float> heightChange(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1);
    for (auto& current : sortedBiomes) {
        if (!current->depthShape.points.empty()) {
            ShapeCurve area = current->area;
            Vector3 AABBoxMin, AABBoxMax;
            std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
            Matrix3<float> falloff2D(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 0.f);
            for (size_t i = 0; i < falloff2D.size(); i++) {
                Vector3 pos = falloff2D.getCoordAsVector3(i);
                if(AABBoxMin.x <= pos.x && pos.x <= AABBoxMax.x && AABBoxMin.y <= pos.y && pos.y <= AABBoxMax.y) {
                    float dist = area.estimateDistanceFrom(pos);
                    falloff2D[i] = dist < 0 ? 1.f : 0.f;
                }
            }
            float desiredDepth = current->depthShape.getPoint(.5f).y;
            Matrix3<float> newHeight = falloff2D.binarize() * -desiredDepth/5.f;
            heightChange += newHeight;
        }
    }
//    heightChange.meanSmooth(5, 5, 1);
    voxelGrid->add2DHeightModification(heightChange, 1.5f);

    // Now add the primitives on top
    for (auto& current : sortedBiomes) {
        Vector3 pos = current->position;
        Vector3 surfacePos = getSurfacePosition(voxelGrid, pos);
        ShapeCurve area = current->area;
        Vector3 areaSize = area.containingBoxSize() + Vector3(0, 0, heightmapGrid->maxHeight);
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        Matrix3<float> falloff2D(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 1.f);
        Matrix3<float> falloff3D(voxelGrid->getDimensions(), 1.f);
        for (size_t i = 0; i < falloff2D.size(); i++) {
            Vector3 pos = falloff2D.getCoordAsVector3(i);
            if(AABBoxMin.x <= pos.x && pos.x <= AABBoxMax.x && AABBoxMin.y <= pos.y && pos.y <= AABBoxMax.y) {
                float dist = area.estimateDistanceFrom(pos);
                falloff2D[i] = dist < 0 ? 1.f - interpolation::fault_distance(-dist, 0.1f) : 0.f;
            }
            for (int z = 0; z < falloff3D.sizeZ; z++) {
                falloff3D.at(pos.x, pos.y, pos.z) = falloff2D[i];
            }
        }
//        Matrix3<float> subFalloff2D

        Matrix3<float> modifications(voxelGrid->getDimensions(), 0.f);
        Matrix3<float> heightmapModifier(voxelGrid->getSizeX(), voxelGrid->getSizeY(), 1, 0.f);

        bool modif3D = false;
        bool modif2D = false;

        if (current->classname == "main-terrain") {
            // Nothing to do, maybe add some noise on the ground?
//            Matrix3<float> noise
        } else if (current->classname == "colline") {
            // Maybe add a sine wave, or a gaussian height increase
//            heightmapModifier.paste(Matrix3<float>::gaussian(areaSize.x, areaSize.y, 1, 10.f).normalize() * 20.f  * falloff2D, (surfacePos - areaSize * .5f) * Vector3(1, 1, 0));
            modif2D = true;
        } else if (current->classname == "fault") {
            // No idea what to do here, the hole will be drilled in "trench"
        } else if (current->classname == "profondeur") {
            // Maybe add a gaussian height decrease
            Matrix3<float> hole = Matrix3<float>::gaussian(areaSize.x, areaSize.y, areaSize.z, 10.f).normalize() * -1.f;
//            modifications.paste(hole, surfacePos - hole.getDimensions() * .5f);
            modif3D = true;
        } else if (current->classname == "plaine") {
            // Don't think there's anything to do here
        } else if (current->classname == "coral-reef") {
            // No idea
        } else if (current->classname == "arche") {
            // Create an arch from the two "point" children
            Vector3 point1 = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 point4 = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            Vector3 point2 = point1 + (point4 - point1) * .3f + Vector3(0, 0, (point4 - point1).norm() / 2.f); // Midpoint1
            Vector3 point3 = point1 + (point4 - point1) * .6f + Vector3(0, 0, (point4 - point1).norm() / 2.f); // Midpoint2
            modifications = archeTunnel(BSpline({point1, point2, point3, point4}), 5.f, 10.f, true, voxelGrid);
            modif3D = true;
        } else if (current->classname == "algae-area") {
            // Maybe add some big spikes at this place
        } else if (current->classname == "trench") {
            // Dig a tunnel on the surface
            Vector3 start = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 end = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            Vector3 midpoint1 = getSurfacePosition(voxelGrid, start + (end - start) * .4f + (end - start).cross(Vector3(0, 0, 1)) * .1f);
            Vector3 midpoint2 = getSurfacePosition(voxelGrid, start + (end - start) * .6f + (end - start).cross(Vector3(0, 0, 1)) * -.1f);
            modifications = archeTunnel(BSpline({start, midpoint1, midpoint2, end}), 8.f, 3.f, false, voxelGrid);
            modif3D = true;
        } else if (current->classname == "coral") {
            // Add a small wall depending on the "point" children
            Vector3 start = getSurfacePosition(voxelGrid, current->instances[0]->position); // Start
            Vector3 end = getSurfacePosition(voxelGrid, current->instances[1]->position); // End
            modifications = archeTunnel(BSpline({start, end}), 5.f, 3.f, true, voxelGrid);
            modif3D = true;
        } else if (current->classname == "point") {
            // Nothing to do, this is the smallest primitive
        } else {
//            std::cout << "How the fuck did I get here? Class was " << current->classname << std::endl;
        }

        if (modif2D)
            voxelGrid->add2DHeightModification(heightmapModifier, 10.f);
        if (modif3D)
            voxelGrid->applyModification(modifications * falloff3D);
    }

    // Yeah whatever... Just do the translation once again
    this->heightmapGrid->fromVoxelGrid(*voxelGrid);
#endif
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
