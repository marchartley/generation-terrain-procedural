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

    if (verbose)
        std::cout << "Preparing main shaders..." << std::endl;

#ifdef linux
    std::string pathToShaders = "/home/simulateurrsm/Documents/Qt_prog/generation-terrain-procedural/src/Shaders/"; // ":/src/Shaders/"
#else
    std::string pathToShaders = "C:/codes/Qt/generation-terrain-procedural/src/Shaders/"; // ":/src/Shaders/"
#endif
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

    if (verbose)
        std::cout << "Loading 3D models :" << std::endl;

    QTemporaryDir tempDirCorals;
    QDirIterator itCorals(":/models_3d/corals/", QDir::Files, QDirIterator::Subdirectories);
    this->possibleCorals = std::vector<Mesh>();
    std::shared_ptr<Shader> coralsShader = std::make_shared<Shader>(vRockShader, fRockShader);
    while (itCorals.hasNext()) {
        QString dir = itCorals.next();
        if (verbose)
            std::cout << "- " << dir.toStdString() << "... " << std::flush;
        QString newPath = tempDirCorals.path() + QFileInfo(dir).fileName();
        QFile::copy(dir, newPath);
        Mesh m = Mesh(coralsShader);
        // Normalize it and move it upward so the anchor is on the ground
        possibleCorals.push_back(m.fromStl(newPath.toStdString()).normalize().translate(0.f, 0.f, .25f));
        if (verbose)
            std::cout << m.vertexArray.size()/3 << " tris, done." << std::endl;
    }
    QTemporaryDir tempDirRocks;
    QDirIterator itRocks(":/models_3d/rocks/", QDir::Files, QDirIterator::Subdirectories);
    this->possibleRocks = std::vector<Mesh>();
    std::shared_ptr<Shader> rocksShader = std::make_shared<Shader>(vRockShader, fRockShader);
    while (itRocks.hasNext()) {
        QString dir = itRocks.next();
        if (verbose)
            std::cout << "- " << dir.toStdString() << "... " << std::flush;
        QString newPath = tempDirRocks.path() + QFileInfo(dir).fileName();
        QFile::copy(dir, newPath);
        Mesh m = Mesh(rocksShader);
        possibleRocks.push_back(m.fromStl(newPath.toStdString()).normalize());
        if (verbose)
            std::cout << m.vertexArray.size()/3 << " tris, done." << std::endl;
    }
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

    float *heightmapData = new float[heightmap->heights.size() * 4];
    Matrix3<Vector3> gradients = heightmap->heights.gradient();
    for (size_t i = 0; i < heightmap->heights.size(); i++) {
        gradients[i] = (gradients[i].normalized() + Vector3(1, 1, 1)) / 2.f;
        heightmapData[i * 4 + 0] = gradients[i].x;
        heightmapData[i * 4 + 1] = gradients[i].y;
        heightmapData[i * 4 + 2] = gradients[i].z;
        heightmapData[i * 4 + 3] = heightmap->heights[i];
    }

    glGenTextures(1, &heightmapFieldTex);
    GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
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
    glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmap->biomeIndices.sizeX, heightmap->biomeIndices.sizeY, 0,
    GL_ALPHA_INTEGER_EXT, GL_INT, heightmap->biomeIndices.data.data());
    heightmapMesh.shader->setInt("biomeFieldTex", 4);
    marchingCubeMesh.shader->setInt("biomeFieldTex", 4);
    layersMesh.shader->setInt("biomeFieldTex", 4);

    if (verbose)
        std::cout << "Loading terrain textures :" << std::endl;
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
        if (basename == "color" || basename == "normal"/* || basename == "displacement"*/) {
            if (dir.endsWith(".tiff")) continue; // Ignore TIFF files for now...
            std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
            QString newPath = tempTextureDir.path() + QString::fromStdString(textureClass) + QFileInfo(dir).fileName();
            QFile::copy(dir, newPath);
            int imgW, imgH, c;
            unsigned char *data = stbi_load(newPath.toStdString().c_str(), &imgW, &imgH, &c, 3);
            if (verbose)
                std::cout << "- " << dir.toStdString().c_str() << " (" << textureClass << ") ..." << std::flush;
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
            if (verbose)
                std::cout << "Done." << std::endl;
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

    heightmapMesh.shader->setInt("maxBiomesColorTextures", indexColorTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesColorTextures", indexColorTextureClass);
    layersMesh.shader->setInt("maxBiomesColorTextures", indexColorTextureClass);

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

    heightmapMesh.shader->setInt("maxBiomesNormalTextures", indexNormalTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesNormalTextures", indexNormalTextureClass);
    layersMesh.shader->setInt("maxBiomesNormalTextures", indexNormalTextureClass);

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

    heightmapMesh.shader->setInt("maxBiomesDisplacementTextures", indexDisplacementTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesDisplacementTextures", indexDisplacementTextureClass);
    layersMesh.shader->setInt("maxBiomesDisplacementTextures", indexDisplacementTextureClass);

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
    this->coralsIndicesAndPositionAndSize = std::vector<std::tuple<int, Vector3, float>>(numberOfRocksDisplayed);
    for (size_t i = 0; i < coralsIndicesAndPositionAndSize.size(); i++) {
        coralsIndicesAndPositionAndSize[i] = std::make_tuple<int, Vector3, float>(
                                              int(random_gen::generate(0, possibleCorals.size())),
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
    /*if (this->heightmap != nullptr) {
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, biomeFieldTex);
        glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA16I_EXT, heightmap->biomeIndices.sizeX, heightmap->biomeIndices.sizeY, 0,
        GL_ALPHA_INTEGER_EXT, GL_INT, ((Matrix3<float>)(heightmap->biomeIndices) / 10.f).data.data());
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
    float maxHeight = heightmap->heights.max();
    this->heightmapMesh.shader->setFloat("maxHeight", maxHeight);
    this->heightmapMesh.shader->setFloat("waterRelativeHeight", waterLevel);
    float *heightmapData = new float[heightmap->heights.size() * 4];
    Matrix3<Vector3> gradients = heightmap->heights.gradient();
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
    int maxBiomeID = 0; //std::max(1, heightmap->biomeIndices.max());
    for (const auto& biomeValues : heightmap->biomeIndices)
        maxBiomeID = std::max(maxBiomeID, (biomeValues.empty() ? 1 : biomeValues.back()));
    Matrix3<std::vector<int>> resizedBiomeIndices = heightmap->biomeIndices; //.resize(heightmap->heights.getDimensions(), RESIZE_MODE::MAX_VAL);
    for (size_t i = 0; i < heightmap->heights.size(); i++) {
//        int biomeID = (resizedBiomeIndices[i].empty() ? 0 : resizedBiomeIndices[i].back());
        float colorTextureOffset = 1.f;
        float normalTextureOffset = 1.f;
        float displacementTextureOffset = 1.f;
        if (resizedBiomeIndices.getDimensions() == heightmap->heights.getDimensions()) {
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
        heightmapData[i * 4 + 3] = heightmap->heights[i] / maxHeight; //1.0; //heightmap->heights[i] / maxHeight;
    }

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmap->getSizeX(), heightmap->getSizeY(), 0,
    GL_RGBA, GL_FLOAT, heightmapData);//heightData.data.data());
    delete[] heightmapData;

    if (mapMode == GRID_MODE) {
        if (this->heightmap == nullptr) {
            std::cerr << "No grid to display" << std::endl;
        } else {
//            float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();

            this->heightmapMesh.display(GL_POINTS);
//            this->heightmap->mesh.shader->setFloat("time", time);
//            this->heightmap->display(true);
        }
    }
    else if (mapMode == VOXEL_MODE) {
        if (this->voxelGrid == nullptr) {
            std::cerr << "No voxel grid to display" << std::endl;
        } else {
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
            if (this->previousHistoryIndex != voxelGrid->getCurrentHistoryIndex()) {
                this->previousHistoryIndex = voxelGrid->getCurrentHistoryIndex();
                Matrix3<float> values = voxelGrid->getVoxelValues();
                // Not the best check, but still pretty good....
                if (marchingCubeMesh.vertexArray.size() != values.size()) {
                    regenerateRocksAndParticles();
                }
                marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values + .5f);

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
                                                  });*/
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
            layersMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
            layersMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
            layersMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
            layersMesh.shader->setFloat("waterRelativeHeight", waterLevel);
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
            this->layersMesh.display(GL_POINTS);
//            this->layerGrid->display();
        }
    }
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

void TerrainGenerationInterface::createTerrainFromNoise(int nx, int ny, int nz, float blockSize, float noise_shifting)
{
    std::shared_ptr<VoxelGrid> tempMap = std::make_shared<VoxelGrid>(nx, ny, nz, blockSize, noise_shifting);
    tempMap->fromIsoData();

    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->heightmap)
        this->heightmap = std::make_shared<Grid>();
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>();

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
    this->heightmap->fromVoxelGrid(*voxelGrid);
    this->layerGrid->fromVoxelGrid(*voxelGrid);

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

void TerrainGenerationInterface::createTerrainFromFile(std::string filename, std::map<std::string, std::shared_ptr<ActionInterface> > actionInterfaces)
{
    std::string ext = toUpper(getExtention(filename));
    if (!this->heightmap)
        this->heightmap = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>();


    // TODO : REMOVE THIS PART VERY SOON !!

    Vector3 terrainSize = Vector3(40, 40, 40);
    layerGrid->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(terrainSize.x, terrainSize.y);
//    LayerBasedGrid layerGrid(terrainSize.x, terrainSize.y, 1.f);
    voxelGrid->fromLayerBased(*layerGrid, terrainSize.z);
    voxelGrid->fromIsoData();

    // Bedrock as noise function
    FastNoiseLite noise;
    noise = voxelGrid->noise;
    std::function noisyBedrock = [&](Vector3 pos) {
        float noiseValue = noise.GetNoise(pos.x * 2.f, pos.y * 2.f) * 2.f + 2.f;
        return noiseValue;
    };
    ShapeCurve bedrockArea = ShapeCurve({
                                            Vector3(-1, -1),
                                            Vector3(-1, 1),
                                            Vector3(1, 1),
                                            Vector3(1, -1)
                                        }).scale(20.f);
    Patch2D bedrock = Patch2D(terrainSize * .5f, bedrockArea, noisyBedrock);


    std::function noisySandLayer = [&](Vector3 pos) {
        float noiseValue = noise.GetNoise(pos.x * .5f + 500.f, pos.y * .5f + 500.f) * 2.f + 6.f;
        return noiseValue;
    };
    Patch2D sandLayer = Patch2D(terrainSize * .5f, bedrockArea, noisySandLayer);


    std::function noisyRock = [&](Vector3 pos) {
        float noiseValue = noise.GetNoise(pos.x * 2.f + 5000.f, pos.y * 2.f + 5000.f) * 3.f;
        noiseValue += normalizedGaussian(Vector3(5, 5), pos + Vector3(5, 5) * .5f, 2.f) * 5.f;
        return noiseValue;
    };
    std::function antiRockNoise = [&](Vector3 pos) {
        float noiseValue = noise.GetNoise(pos.x * 2.f + 5000.f, pos.y * 2.f + 5000.f) * 3.f;
        noiseValue += (1.f - normalizedGaussian(Vector3(5, 5), pos + Vector3(5, 5) * .5f, 1.f)) * 2.f;
        return noiseValue;
    };
    ShapeCurve rockArea = ShapeCurve({
                                            Vector3(-1, -1),
                                            Vector3(-1, 1),
                                            Vector3(1, 1),
                                            Vector3(1, -1)
                                        }).scale(5.f);


    std::function noisySandBump = [&](Vector3 pos) {
        float noiseValue = noise.GetNoise(pos.x * 2.f + 1000.f, pos.y * 2.f + 1000.f) * 1.f;
        noiseValue += normalizedGaussian(Vector3(30, 30), pos, 15.f) * 5.f;
        return noiseValue;
    };
    ShapeCurve sandArea = ShapeCurve({
                                            Vector3(-1, -1),
                                            Vector3(-1, 1),
                                            Vector3(1, 1),
                                            Vector3(1, -1)
                                        }).scale(20.f);
    Patch2D sand = Patch2D(terrainSize * .5f, sandArea, noisySandBump);




    auto oldVoxels = layerGrid->voxelize(terrainSize.z);
    auto newVoxels = oldVoxels;

    layerGrid->add(bedrock, TerrainTypes::BEDROCK, false);
    newVoxels = layerGrid->voxelize(terrainSize.z);
    voxelGrid->applyModification(newVoxels - oldVoxels);
    oldVoxels = newVoxels;
//    return;

    layerGrid->add(sandLayer, TerrainTypes::SAND, false);
    newVoxels = layerGrid->voxelize(terrainSize.z);
    voxelGrid->applyModification(newVoxels - oldVoxels);
    oldVoxels = newVoxels;

    layerGrid->add(bedrock, TerrainTypes::BEDROCK, false);
    newVoxels = layerGrid->voxelize(terrainSize.z);
    voxelGrid->applyModification(newVoxels - oldVoxels);
    oldVoxels = newVoxels;

    for (int i = 0; i < 10; i++) {
        Vector3 pos = Vector3::random(Vector3(layerGrid->getSizeX(), layerGrid->getSizeY()));
        Patch2D antiRock = Patch2D(pos, rockArea, antiRockNoise);
        Patch2D rock = Patch2D(pos, rockArea, noisyRock);
        layerGrid->add(antiRock, TerrainTypes::WATER, true, .5f);
        layerGrid->add(antiRock, TerrainTypes::WATER, true, .5f);
        layerGrid->add(rock, TerrainTypes::ROCK, true, .5f);
        newVoxels = layerGrid->voxelize(terrainSize.z);
        voxelGrid->applyModification(newVoxels - oldVoxels);
        oldVoxels = newVoxels;
    }
    //    layerGrid->add(sand, TerrainTypes::SAND, true, 1.f);
    layerGrid->add(sandLayer, TerrainTypes::SAND, false);
    newVoxels = layerGrid->voxelize(terrainSize.z);
    voxelGrid->applyModification(newVoxels - oldVoxels);
    oldVoxels = newVoxels;


//    for (int i = 0; i < 100; i++)
//        layerGrid->thermalErosion();
    layerGrid->cleanLayers();
    newVoxels = layerGrid->voxelize(terrainSize.z);
    voxelGrid->applyModification(newVoxels - oldVoxels);
    oldVoxels = newVoxels;

    return;

    // END OF SHITTY PART


    if (ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        heightmap->loadFromHeightmap(filename, 127, 127, 40);

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
            this->createTerrainFromBiomes(json_content);
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
        // Sorry heightmap... You take too long...
//        heightmap->fromVoxelGrid(*voxelGrid);

    } /*else {
        // In any other case, consider that nothing has been done, cancel.
        return;
    }*/
    this->voxelGrid->createMesh();
    this->heightmap->createMesh();


    voxelGrid->flowField = Matrix3<Vector3>(voxelGrid->environmentalDensities.getDimensions());
    for (size_t i = 0; i < voxelGrid->environmentalDensities.size(); i++) {
        if (voxelGrid->environmentalDensities.getCoordAsVector3(i).z < voxelGrid->environmentalDensities.sizeZ * waterLevel) {
            voxelGrid->environmentalDensities[i] = 1000.f; // Water density
            voxelGrid->flowField[i] = Vector3(0, -.3f, 0.f); // Water flow
        } else {
            voxelGrid->environmentalDensities[i] = 1.f; // Air density
            voxelGrid->flowField[i] = Vector3(0, .1f, 0.f); // Air flow
        }
    }

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
    if (!this->heightmap)
        this->heightmap = std::make_shared<Grid>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    this->heightmap->heights = Matrix3<float>(124, 124, 1, 20.f);
    this->heightmap->maxHeight = 40.f;

    this->voxelGrid->from2DGrid(*this->heightmap);
    this->voxelGrid->fromIsoData();

    this->biomeGenerationNeeded = true;
    this->biomeGenerationModelData = json_content;

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
