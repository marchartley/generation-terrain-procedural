#include "TerrainGenerationInterface.h"

#include "Biomes/BiomeInstance.h"
#include "Utils/ConstraintsSolver.h"
#include "Utils/ShapeCurve.h"
#include "Utils/Utils.h"
#include "Utils/Voronoi.h"
#include "Utils/stb_image.h"
#include "Utils/stb_image_write.h"
#include "DataStructure/Matrix.h"

TerrainGenerationInterface::TerrainGenerationInterface(QWidget *parent) : ActionInterface("terrain_gen", parent)
{
    this->createGUI();
}

void TerrainGenerationInterface::prepareShader()
{
    bool verbose = true;

    colorTexturesIndex = {
        {"algae", 0},
        {"beach_sand", 1},
        {"coral", 2},
        {"coral2", 3},
        {"mountain_rocks", 4},
        {"underwater_ground", 5},
        {"underwater_sand", 6}
    };
    std::map<TerrainTypes, std::string> materialToTexture = {
        {WATER, ""},
        {AIR, ""},
        {CORAL, "coral"},
        {SAND, "beach_sand"},
        {DIRT, "underwater_sand"},
        {ROCK, "underwater_ground"},
        {BEDROCK, "mountain_rocks"},
    };

//    if (verbose)
//        std::cout << "Preparing main shaders..." << std::endl;

    std::string pathToShaders = "src/Shaders/";
    std::string vShader_mc_voxels = pathToShaders + "MarchingCubes.vert";
    std::string gShader_mc_voxels = pathToShaders + "MarchingCubes.geom";
    std::string fShader_mc_voxels = pathToShaders + "MarchingCubes.frag";

    std::string vShader_grid = pathToShaders + "grid.vert";
    std::string gShader_grid = pathToShaders + "grid.geom";
    std::string fShader_grid = pathToShaders + "grid.frag";

    std::string vShader_layers = pathToShaders + "layer_based.vert";
    std::string gShader_layers = pathToShaders + "layer_based.geom";
    std::string fShader_layers = pathToShaders + "MarchingCubes.frag"; // pathToShaders + "layer_based.frag";

    std::string vRockShader = pathToShaders + "rockShader.vert";
    std::string fRockShader = pathToShaders + "rockShader.frag";

    std::string vParticleShader = pathToShaders + "particle.vert";
    std::string fParticleShader = pathToShaders + "particle.frag";
    std::string gParticleShader = pathToShaders + "particle.geom";

    std::string vWaterShader = pathToShaders + "no_shader.vert";
    std::string fWaterShader = pathToShaders + "no_shader.frag";

    waterLevelMesh = Mesh(std::make_shared<Shader>(vWaterShader, fWaterShader));
    waterLevelMesh.shader->setVector("color", std::vector<float>{1., 1., 1., .1});
    waterLevelMesh.cullFace = false;


    layersMesh = Mesh(std::make_shared<Shader>(vShader_layers, fShader_layers, gShader_layers),
                          true, GL_TRIANGLES);
    layersMesh.useIndices = false;

    Matrix3<int> materials; Matrix3<float> matHeights;
    std::tie(materials, matHeights) = layerGrid->getMaterialAndHeightsGrid();
    layersMesh.shader->setTexture3D("matIndicesTex", 1, materials);
    layersMesh.shader->setTexture3D("matHeightsTex", 2, matHeights);

    std::vector<Vector3> layersPoints(materials.size());
    for (size_t i = 0; i < layersPoints.size(); i++) {
        layersPoints[i] = materials.getCoordAsVector3(i);
    }
    layersMesh.fromArray(layersPoints);
    layersMesh.update();


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

    this->regenerateRocksAndParticles();
    this->particlesMesh = Mesh(this->randomParticlesPositions,
                               std::make_shared<Shader>(vParticleShader, fParticleShader, gParticleShader),
                               true, GL_POINTS);
    particlesMesh.shader->setVector("color", std::vector<float>({1.f, 1.f, 1.f, .2f}));
    this->particlesMesh.shader->setFloat("maxTerrainHeight", voxelGrid->getSizeZ());

    this->heightmapMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_grid), true, GL_POINTS);

    Matrix3<float> heightData(this->heightmap->getSizeX(), this->heightmap->getSizeY());
    std::vector<Vector3> positions(heightData.size());
    for (size_t i = 0; i < positions.size(); i++) {
        positions[i] = heightData.getCoordAsVector3(i);
        heightData[i] = heightmap->getHeight(positions[i].x, positions[i].y);
    }
    heightmapMesh.useIndices = false;
    heightmapMesh.fromArray(positions);
    heightmapMesh.update();
    heightmapMesh.shader->use();
    GlobalsGL::f()->glBindVertexArray(heightmapMesh.vao);
    heightmapMesh.shader->setInt("heightmapFieldTex", 3);
    marchingCubeMesh.shader->setInt("heightmapFieldTex", 3);
    layersMesh.shader->setInt("heightmapFieldTex", 3);

    Matrix3<float> heights = heightmap->getHeights();
    float *heightmapData = new float[heights.size() * 4];
    Matrix3<Vector3> gradients = heights.gradient();
    for (size_t i = 0; i < heights.size(); i++) {
        gradients[i] = (gradients[i].normalized() + Vector3(1, 1, 1)) / 2.f;
        heightmapData[i * 4 + 0] = gradients[i].x;
        heightmapData[i * 4 + 1] = gradients[i].y;
        heightmapData[i * 4 + 2] = gradients[i].z;
        heightmapData[i * 4 + 3] = heights[i];
    }

    glGenTextures(1, &heightmapFieldTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
//            glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, heightmap->getSizeX(), heightmap->getSizeY(), 0,
//            GL_RED, GL_FLOAT, heightmap->heights.data.data());//heightData.data.data());
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmap->getSizeX(), heightmap->getSizeY(), 0,
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
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmap->getBiomeIndices().sizeX, heightmap->getBiomeIndices().sizeY, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, heightmap->getBiomeIndices().data.data());
    heightmapMesh.shader->setInt("biomeFieldTex", 4);
    marchingCubeMesh.shader->setInt("biomeFieldTex", 4);
    layersMesh.shader->setInt("biomeFieldTex", 4);

    if (verbose)
        std::cout << "Loading terrain textures... " << std::flush;
//    QTemporaryDir tempTextureDir;
    QDirIterator iTex("src/assets/textures/", QDir::Files, QDirIterator::Subdirectories);
//    QDirIterator iTex(":/terrain_textures/", QDir::Files, QDirIterator::Subdirectories);
    Matrix3<Vector3> allColorTextures;
    Matrix3<Vector3> allNormalTextures;
    Matrix3<Vector3> allDisplacementTextures;
    int indexColorTextureClass = 0;
    int indexNormalTextureClass = 0;
    int indexDisplacementTextureClass = 0;

    std::vector<QString> filesToLoad;
    while (iTex.hasNext()) {
        QString dir = iTex.next();
        QString basename = QFileInfo(dir).baseName();
        if (basename == "color" || basename == "normal") { // || basename == "displacement") {
            if (dir.endsWith(".tiff")) continue; // Ignore TIFF files for now...
            filesToLoad.push_back(dir);
        }
    }
    size_t nbFiles = filesToLoad.size();
    std::vector<Matrix3<Vector3>> texturesMaps(nbFiles);
    int textureSizes = 128;
    auto startTime = std::chrono::system_clock::now();
    #pragma omp parallel for
    for (size_t i = 0; i < nbFiles; i++) {
        QString dir = filesToLoad[i];
        QString basename = QFileInfo(dir).baseName();
            std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
            int imgW, imgH, c;
            unsigned char *data = stbi_load(dir.toStdString().c_str(), &imgW, &imgH, &c, 3);
            texturesMaps[i] = Matrix3<Vector3>(imgW, imgH);
            for (int x = 0; x < imgW; x++) {
                for (int y = 0; y < imgH; y++) {
                    unsigned char r = data[3 * (x + y * imgW) + 0];
                    unsigned char g = data[3 * (x + y * imgW) + 1];
                    unsigned char b = data[3 * (x + y * imgW) + 2];
                    texturesMaps[i].at(x, y) = Vector3(r, g, b);
                }
            }
            texturesMaps[i] = texturesMaps[i].resize(textureSizes, textureSizes, 1, RESIZE_MODE::LINEAR);
            stbi_image_free(data);
    }
    allColorTextures = Matrix3<Vector3>(textureSizes * colorTexturesIndex.size(), textureSizes);
    allNormalTextures = Matrix3<Vector3>(textureSizes * colorTexturesIndex.size(), textureSizes);
    allDisplacementTextures = Matrix3<Vector3>(textureSizes * colorTexturesIndex.size(), textureSizes);
    auto endTempTime = std::chrono::system_clock::now();
    for (size_t i = 0; i < nbFiles; i++) {
        QString dir = filesToLoad[i];
        QString basename = QFileInfo(dir).baseName();
        std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
        Matrix3<Vector3>& map = texturesMaps[i];
        int index = colorTexturesIndex[textureClass];
        if (basename == "color") {
            allColorTextures.paste(map, Vector3(textureSizes * index, 0));
        } else if (basename == "normal") {
            allNormalTextures.paste(map, Vector3(textureSizes * index, 0));
        } else if (basename == "displacement") {
            allDisplacementTextures.paste(map, Vector3(textureSizes * index, 0));
        }
                /*
        if (basename == "color") {
            colorTexturesIndex[textureClass] = indexColorTextureClass++;
            if (allColorTextures.empty()) {
                allColorTextures = map;
            } else {
                allColorTextures = allColorTextures.concat(map);
            }
        } else if (basename == "normal") {
//            normalTexturesIndex[textureClass] = indexNormalTextureClass++;
            if (allNormalTextures.empty()) {
                allNormalTextures = map;
            } else {
                allNormalTextures = allNormalTextures.concat(map);
            }
        } else if (basename == "displacement") {
//            displacementTexturesIndex[textureClass] = indexDisplacementTextureClass++;
            if (allDisplacementTextures.empty()) {
                allDisplacementTextures = map;
            } else {
                allDisplacementTextures = allDisplacementTextures.concat(map);
            }
        }*/
    }
    auto endTime = std::chrono::system_clock::now();
    if (verbose) {
        std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms " << std::endl;
        std::cout << "( " << std::chrono::duration_cast<std::chrono::milliseconds>(endTempTime - startTime).count() << "ms for the parallel and " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - endTempTime).count() << "ms for sequential)" << std::endl;
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
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allColorTextures.sizeX, allColorTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesColors);
    delete[] allTexturesColors;
    heightmapMesh.shader->setInt("allBiomesColorTextures", 5);
    marchingCubeMesh.shader->setInt("allBiomesColorTextures", 5);
    layersMesh.shader->setInt("allBiomesColorTextures", 5);

    heightmapMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);
    layersMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);

    glGenTextures(1, &allBiomesNormalTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE6);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, allBiomesNormalTextures);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allNormalTextures.sizeX, allNormalTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesNormal);
    delete[] allTexturesNormal;
    heightmapMesh.shader->setInt("allBiomesNormalTextures", 6);
    marchingCubeMesh.shader->setInt("allBiomesNormalTextures", 6);
    layersMesh.shader->setInt("allBiomesNormalTextures", 6);

    heightmapMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);
    layersMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);

    glGenTextures(1, &allBiomesDisplacementTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE7);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, allBiomesDisplacementTextures);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, allDisplacementTextures.sizeX, allDisplacementTextures.sizeY, 0,
    GL_RGBA, GL_FLOAT, allTexturesDisplacement);
    delete[] allTexturesDisplacement;
    heightmapMesh.shader->setInt("allBiomesDisplacementTextures", 7);
    marchingCubeMesh.shader->setInt("allBiomesDisplacementTextures", 7);
    layersMesh.shader->setInt("allBiomesDisplacementTextures", 7);

    heightmapMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());//indexDisplacementTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());//indexDisplacementTextureClass);
    layersMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());//indexDisplacementTextureClass);

    /*
     * // TODO one day : put the material data in a texture to be used in the shaders
    std::vector<TerrainTypes> materialTypes = {WATER, AIR, CORAL, SAND, DIRT, ROCK, BEDROCK};
    float* dataForMaterialTypes = new float(materialTypes.size() * 4);
    for (size_t i = 0; i < materialTypes.size(); i++) {
        dataForMaterialTypes[3 * i + 0] = LayerBasedGrid::minDensityFromMaterial(materialTypes[i]);
        dataForMaterialTypes[3 * i + 1] = LayerBasedGrid::maxDensityFromMaterial(materialTypes[i]);
        dataForMaterialTypes[3 * i + 2] = (materialToTexture[materialTypes[i]] == "" ? 0 : 1); // Displayed or not
        dataForMaterialTypes[3 * i + 3] = (materialToTexture[materialTypes[i]] == "" ? 0 : colorTexturesIndex[materialToTexture[materialTypes[i]]]);
    }
    glGenTextures(1, &biomesAndDensitiesTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE8);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, biomesAndDensitiesTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, materialTypes.size(), 1, 0,
    GL_RGBA, GL_FLOAT, dataForMaterialTypes);
    heightmapMesh.shader->setInt("biomesAndDensitiesTex", 8);
    marchingCubeMesh.shader->setInt("biomesAndDensitiesTex", 8);
    layersMesh.shader->setInt("biomesAndDensitiesTex", 8);
    */

    /*for (auto& [lab, ID] : colorTexturesIndex) {
        std::cout << " -> " << lab << " : " << ID << std::endl;
    }
    Plotter plt;
    plt.addImage(allColorTextures);
    plt.exec();*/

    updateDisplayedView(voxelGridOffset, voxelGridScaling);

//    this->heightmap->mesh.shader = std::make_shared<Shader>(vShader_grid, fShader_grid);
    if (verbose)
        std::cout << "Terrain shaders and assets ready." << std::endl;
}

void TerrainGenerationInterface::updateDisplayedView(Vector3 newVoxelGridOffset, float newVoxelGridScaling)
{
    this->voxelGridOffset = newVoxelGridOffset;
    this->voxelGridScaling = newVoxelGridScaling;

//    heightmapMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
//    heightmapMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);
    marchingCubeMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
    marchingCubeMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);

    layersMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
    layersMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);
}

void TerrainGenerationInterface::afterTerrainUpdated()
{

}

void TerrainGenerationInterface::regenerateRocksAndParticles()
{
    this->randomParticlesPositions = std::vector<Vector3>(1000);
    for (size_t i = 0; i < randomParticlesPositions.size(); i++) {
        randomParticlesPositions[i] = Vector3::random(voxelGrid->getDimensions()); //Vector3(random_gen::generate(0, voxelGrid->sizeX), random_gen::generate(0, voxelGrid->sizeY), random_gen::generate(0, voxelGrid->sizeZ));
    }
    if (this->particlesMesh.shader != nullptr)
        this->particlesMesh.shader->setFloat("maxTerrainHeight", voxelGrid->getSizeZ());
    this->startingTime = std::chrono::system_clock::now();
}


void TerrainGenerationInterface::display(MapMode mapMode, SmoothingAlgorithm smoothingAlgorithm, bool displayParticles)
{
    /*if (this->heightmap != nullptr) {
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, biomeFieldTex);
        glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmap->getBiomeIndices().sizeX, heightmap->getBiomeIndices().sizeY, 0,
        GL_ALPHA_INTEGER_EXT, GL_INT, ((Matrix3<float>)(heightmap->getBiomeIndices()) / 10.f).data.data());
    }*/

    GlobalsGL::f()->glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, allBiomesColorTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE6);
    glBindTexture(GL_TEXTURE_2D, allBiomesNormalTextures);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE7);
    glBindTexture(GL_TEXTURE_2D, allBiomesDisplacementTextures);

    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);
//            glTexImage2D( GL_TEXTURE_2D, 0, GL_R32F, heightmap->getSizeX(), heightmap->getSizeY(), 0,
//            GL_RED, GL_FLOAT, heightmap->heights.data.data());//heightData.data.data());
    Matrix3<float> heights = heightmap->getHeights();
    float maxHeight = heights.max();
    this->heightmapMesh.shader->setFloat("maxHeight", maxHeight);
    this->heightmapMesh.shader->setFloat("waterRelativeHeight", waterLevel);
    float *heightmapData = new float[heights.size() * 4];
    Matrix3<Vector3> gradients = heights.gradient();
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
    int maxBiomeID = 0; //std::max(1, heightmap->getBiomeIndices().max());
    for (const auto& biomeValues : heightmap->getBiomeIndices())
        maxBiomeID = std::max(maxBiomeID, (biomeValues.empty() ? 1 : biomeValues.back()));
    Matrix3<std::vector<int>> resizedBiomeIndices = heightmap->getBiomeIndices(); //.resize(heightmap->heights.getDimensions(), RESIZE_MODE::MAX_VAL);
    for (size_t i = 0; i < heights.size(); i++) {
//        int biomeID = (resizedBiomeIndices[i].empty() ? 0 : resizedBiomeIndices[i].back());
        float colorTextureOffset = 1.f;
        float normalTextureOffset = 1.f;
        float displacementTextureOffset = 1.f;
        if (resizedBiomeIndices.getDimensions() == heights.getDimensions()) {
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
        heightmapData[i * 4 + 3] = heights[i] / maxHeight; //1.0; //heightmap->heights[i] / maxHeight;
    }

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmap->getSizeX(), heightmap->getSizeY(), 0,
    GL_RGBA, GL_FLOAT, heightmapData);//heightData.data.data());
    delete[] heightmapData;

    if (mapMode == GRID_MODE) {
        if (this->heightmap == nullptr) {
            std::cerr << "No grid to display" << std::endl;
        } else {
//            float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();

            Matrix3<float> heightData(this->heightmap->getSizeX(), this->heightmap->getSizeY());
            std::vector<Vector3> positions(heightData.size());
            for (size_t i = 0; i < positions.size(); i++) {
                positions[i] = heightData.getCoordAsVector3(i);
                heightData[i] = heightmap->getHeight(positions[i].x, positions[i].y);
            }
            heightmapMesh.useIndices = false;
            heightmapMesh.fromArray(positions);
            heightmapMesh.update();
            this->heightmapMesh.display(GL_POINTS);
//            this->heightmap->mesh.shader->setFloat("time", time);
//            this->heightmap->display(true);
        }
    }
    else if (mapMode == VOXEL_MODE) {
        if (this->voxelGrid == nullptr) {
            std::cerr << "No voxel grid to display" << std::endl;
        } else {
            Matrix3<float> values = voxelGrid->getVoxelValues();
            marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values + .5f);
            marchingCubeMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
            marchingCubeMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
            marchingCubeMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
            marchingCubeMesh.shader->setFloat("waterRelativeHeight", waterLevel);
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
            // Check if something changed on the terrain :
            if (this->voxelsPreviousHistoryIndex != voxelGrid->getCurrentHistoryIndex()) {
                this->voxelsPreviousHistoryIndex = voxelGrid->getCurrentHistoryIndex();
                // Not the best check, but still pretty good....
                if (marchingCubeMesh.vertexArray.size() != values.size()) {
                    regenerateRocksAndParticles();
                }
            }
            if (displayParticles) {
//                particlesMesh.fromArray(randomParticlesPositions); // Maybe we shouldn't displace the particles
                float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();
                particlesMesh.shader->setFloat("time", time); // std::chrono::system_clock::now().time_since_epoch().count());
                particlesMesh.display(GL_POINTS);
            }
        }
    }
    else if (mapMode == LAYER_MODE) {
        if (this->layerGrid == nullptr) {
            std::cerr << "No layer based grid to display" << std::endl;
        } else {
            layersMesh.cullFace = true;
            layersMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
            layersMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
            layersMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
            layersMesh.shader->setFloat("waterRelativeHeight", waterLevel);

            if (this->layersPreviousHistoryIndex != layerGrid->_historyIndex) {
                this->layersPreviousHistoryIndex = layerGrid->_historyIndex;
                Matrix3<int> materials; Matrix3<float> matHeights;
                std::tie(materials, matHeights) = layerGrid->getMaterialAndHeightsGrid();
                layersMesh.shader->setTexture3D("matIndicesTex", 1, materials);
                layersMesh.shader->setTexture3D("matHeightsTex", 2, matHeights);

                std::vector<Vector3> layersPoints(materials.size());
                for (size_t i = 0; i < layersPoints.size(); i++) {
                    layersPoints[i] = materials.getCoordAsVector3(i);
                }
                layersMesh.fromArray(layersPoints);
                layersMesh.update();
            }
            this->layersMesh.display(GL_POINTS);
        }
    }

    // TEST OF GEOMETRY EXTRACTION
//    auto vals = voxelGrid->getVoxelValues();
//    marchingCubeMesh.extractGeometryFromShaders(vals);
}

void TerrainGenerationInterface::displayWaterLevel()
{
    waterLevelMesh.fromArray({
                                Vector3(1.f, 1.f, waterLevel * voxelGrid->getSizeZ()),
                                Vector3(voxelGrid->getSizeX()-1.f, 1.f, waterLevel * voxelGrid->getSizeZ()),
                                Vector3(1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ()),
                                 Vector3(voxelGrid->getSizeX()-1.f, 1.f, waterLevel * voxelGrid->getSizeZ()),
                                Vector3(voxelGrid->getSizeX()-1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ()),
                                Vector3(1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ())
                             });
    waterLevelMesh.display();
}

void TerrainGenerationInterface::createTerrainFromNoise(int nx, int ny, int nz/*, float blockSize*/, float noise_shifting)
{
//    std::shared_ptr<VoxelGrid> tempMap = std::make_shared<VoxelGrid>(nx, ny, nz/*, blockSize*/, noise_shifting);
//    tempMap->fromIsoData();
    VoxelGrid tempMap = VoxelGrid(nx, ny, nz, noise_shifting);

    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>();
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>();

    // Okay, super dirty but I don't know how to manage shared_ptr...
//    this->voxelGrid->sizeX = tempMap->sizeX;
//    this->voxelGrid->sizeY = tempMap->sizeY;
//    this->voxelGrid->sizeZ = tempMap->sizeZ;
//    this->voxelGrid->blockSize = tempMap->blockSize;
//    this->voxelGrid->chunkSize = tempMap->chunkSize;
//    this->voxelGrid->noise_shifting = tempMap->noise_shifting;
//    this->voxelGrid->initMap();
//    this->voxelGrid->tempData = tempMap->tempData;
//    this->voxelGrid->fromIsoData();
    *voxelGrid = tempMap;
    this->heightmap->fromVoxelGrid(*voxelGrid);
    this->layerGrid->fromVoxelGrid(*voxelGrid);

    this->addTerrainAction(nlohmann::json({
                                              {"from_noise", true},
                                              {"noise_parameters", {
                                                   {"nx", nx},
                                                   {"ny", ny},
                                                   {"nz", nz},
//                                                   {"block_size", blockSize},
                                                   {"noise_shifting", noise_shifting}
                                               }}
                                          }));
}

void TerrainGenerationInterface::createTerrainFromFile(std::string filename, std::map<std::string, std::shared_ptr<ActionInterface> > actionInterfaces)
{
    std::string ext = toUpper(getExtention(filename));
    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>();


    // TODO : REMOVE THIS PART VERY SOON !!
    float radius = 20.f;
    Vector3 terrainSize = Vector3(4*radius, 4*radius, 40);
    this->layerGrid = std::make_shared<LayerBasedGrid>(terrainSize.x, terrainSize.y, 1.f);

//    ImplicitPatch* init = new ImplicitPrimitive;
//    ImplicitPatch* p1 = new ImplicitPatch(Vector3(0, 0, -10), Vector3(), Vector3(20, 20, 20), ImplicitPatch::createSphereFunction(0, 20, 20, 5)); p1->name = "sphere1"; p1->densityValue = 0.1f;
//    ImplicitPatch* s1 = new ImplicitPatch(ImplicitPatch::createStack(init, p1, ImplicitPatch::FIXED)); s1->name ="stack1";
//    ImplicitPatch* p2 = new ImplicitPatch(Vector3(20, 20, -10), Vector3(), Vector3(20, 20, 20), ImplicitPatch::createSphereFunction(0, 20, 20, 5)); p2->name = "sphere2"; p2->densityValue = 3.f;
//    ImplicitPatch* s2 = new ImplicitPatch(ImplicitPatch::createStack(s1, p2, ImplicitPatch::ABOVE)); s2->name = "stack2";

    /*
    ImplicitPrimitive* water = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* air = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* coral = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* sand = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* dirt = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* rock = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
    ImplicitPrimitive* bedrock = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);//new ImplicitPrimitive;
*/
    std::vector<TerrainTypes> materials = {WATER, AIR, CORAL, SAND, DIRT, ROCK, BEDROCK};
    std::vector<ImplicitPrimitive*> primitives(materials.size());
    std::vector<ImplicitOperator*> blends(primitives.size() - 1);
    /*
    water->setDimensions(Vector3(20, 20, 20)); water->position = Vector3(0, -10, 0); water->name = "Water"; water->material = TerrainTypes::WATER;
    air->setDimensions(Vector3(20, 20, 20)); air->position = Vector3(0, 10, 0); air->name = "Air"; air->material = TerrainTypes::AIR;
    coral->setDimensions(Vector3(20, 20, 20)); coral->position = Vector3(0, 30, 0); coral->name = "Coral"; coral->material = TerrainTypes::CORAL;
    sand->setDimensions(Vector3(20, 20, 20)); sand->position = Vector3(20, -10, 0); sand->name = "Sand"; sand->material = TerrainTypes::SAND;
    dirt->setDimensions(Vector3(20, 20, 20)); dirt->position = Vector3(20, 10, 0); dirt->name = "Dirt"; dirt->material = TerrainTypes::DIRT;
    rock->setDimensions(Vector3(20, 20, 20)); rock->position = Vector3(20, 30, 0); rock->name = "Rock"; rock->material = TerrainTypes::ROCK;
    bedrock->setDimensions(Vector3(20, 20, 20)); bedrock->position = Vector3(0, -10, 0); bedrock->name = "Bedrock"; bedrock->material = TerrainTypes::BEDROCK;*/
    for (size_t i = 0; i < primitives.size(); i++) {
        primitives[i] = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(ImplicitPatch::Sphere, Vector3(20, 20, 20), 0.f);
        primitives[i]->position = Vector3(-10, -10, 0); //Vector3(-10 + 20 * float(i % (primitives.size() / 2)), 0 + 20 * float(i / (primitives.size() / 2)), 0);
        primitives[i]->setDimensions(Vector3(20, 20, 20));
        primitives[i]->name = stringFromMaterial(materials[i]) + " sphere";
        primitives[i]->material = materials[i];
    }
    for (size_t i = 0; i < blends.size(); i++) {
        auto& blend = blends[i];
        blend = new ImplicitOperator;
        blend->composeFunction = ImplicitPatch::BLEND;
        blend->positionalB = ImplicitPatch::FIXED_POS;
        blend->name = "Blend " + std::to_string(i + 1);
        blend->blendingFactor = 2.f;
        if (i == 0) {
            blend->composableA = primitives[0];
            blend->composableB = primitives[1];
        } else {
            blend->composableA = blends[i - 1];
            blend->composableB = primitives[i + 1];
        }
        blend->updateCache();
    }
//    this->layerGrid->add(blends.back());

    voxelGrid->fromLayerBased(*layerGrid, terrainSize.z);
    voxelGrid->fromIsoData();
    heightmap->fromLayerGrid(*layerGrid);

    Q_EMIT this->terrainUpdated();
    return;

    // END OF SHITTY PART


    if (ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        heightmap->loadFromHeightmap(filename, 128, 128, 40);

        voxelGrid->from2DGrid(*heightmap);
        voxelGrid->fromIsoData();
        /*
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
        */

    } else if (ext == "JSON") {
        // The JSON file contains the list of actions made on a map
        std::ifstream file(filename);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json json_content = nlohmann::json::parse(content);
        if (!json_content.contains("actions")) {
            if (json_content.contains(ImplicitPatch::json_identifier)) {
                this->createTerrainFromImplicitPatches(json_content);
            } else {
                this->createTerrainFromBiomes(json_content);
            }
        } else {
            for (auto action : json_content.at("actions")) {
                // Let all the interfaces try to replay their actions
                for (auto& possibleAction : actionInterfaces)
                    possibleAction.second->replay(action);
            }
        }
    } else /*if (ext == "DATA" || ext.empty())*/ {
        // Then it's our custom voxel grid file
        voxelGrid->retrieveMap(filename);
        voxelGrid->fromIsoData();
        voxelGrid->smoothVoxels();
        voxelGrid->smoothVoxels();
        // Sorry heightmap... You take too long...
        heightmap->fromVoxelGrid(*voxelGrid);

    } /*else {
        // In any other case, consider that nothing has been done, cancel.
        return;
    }*/
//    this->voxelGrid->createMesh();
//    this->heightmap->createMesh();
//    layerGrid->from2DGrid(*heightmap);

/*
    voxelGrid->flowField = Matrix3<Vector3>(voxelGrid->environmentalDensities.getDimensions());
    for (size_t i = 0; i < voxelGrid->environmentalDensities.size(); i++) {
        if (voxelGrid->environmentalDensities.getCoordAsVector3(i).z < voxelGrid->environmentalDensities.sizeZ * waterLevel) {
            voxelGrid->environmentalDensities[i] = 1000.f; // Water density
            voxelGrid->flowField[i] = Vector3(0, -.3f, 0.f); // Water flow
        } else {
            voxelGrid->environmentalDensities[i] = 1.f; // Air density
            voxelGrid->flowField[i] = Vector3(0, .1f, 0.f); // Air flow
        }
    }*/

    this->addTerrainAction(nlohmann::json({
                                              {"from_file", filename},
                                              {"from_noise", false}
                                          }));
    Q_EMIT this->updated();
}

void TerrainGenerationInterface::createTerrainFromBiomes(nlohmann::json json_content)
{
    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    *heightmap = Heightmap(124, 124, 40.f);
    heightmap->raise(Matrix3<float>(heightmap->getSizeX(), heightmap->getSizeY(), 1.f, 20.f));
//    this->heightmap->heights = Matrix3<float>(124, 124, 1, 20.f);
//    this->heightmap->maxHeight = 40.f;

    this->voxelGrid->from2DGrid(*this->heightmap);
    this->voxelGrid->fromIsoData();

    this->biomeGenerationNeeded = true;
    this->biomeGenerationModelData = json_content;

}

void TerrainGenerationInterface::createTerrainFromImplicitPatches(nlohmann::json json_content)
{
    // Dunno...
    // TODO : Link to the PrimitivePatchInterface
}

void TerrainGenerationInterface::saveTerrain(std::string filename)
{
    std::string ext = toUpper(getExtention(filename));
    if (ext == "PNG" || ext == "JPG" || ext == "TGA" || ext == "BMP" || ext == "HDR") {
        // To heightmap
        this->heightmap->fromVoxelGrid(*voxelGrid); // Just to be sure to have the last values
        this->heightmap->saveHeightmap(filename);
    } else if (ext == "JSON") {
        // To JSON file containing the list of actions made on a map
        this->saveAllActions(filename);
    } else {
        // Otherwise it's our custom voxel grid file
        voxelGrid->saveMap(filename);
    }
}

void TerrainGenerationInterface::reloadShaders()
{
    prepareShader();
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
//            int blockSize = parameters.at("block_size").get<float>();
            int noise_shifting = parameters.at("noise_shifting").get<float>();
            this->createTerrainFromNoise(nx, ny, nz/*, blockSize*/, noise_shifting);
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
