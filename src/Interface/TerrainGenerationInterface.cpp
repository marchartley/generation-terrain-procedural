#include "TerrainGenerationInterface.h"

#include "Biomes/BiomeInstance.h"
#include "EnvObject/EnvObject.h"
#include "Interface/FancySlider.h"
#include "Interface/InterfaceUtils.h"
#include "Utils/ConstraintsSolver.h"
#include "Utils/ShapeCurve.h"
#include "Utils/Utils.h"
#include "Utils/Voronoi.h"
//#include "Utils/stb_image.h"
//#include "Utils/stb_image_write.h"
#include "DataStructure/Matrix.h"
#include "DataStructure/Image.h"
#include "Graphics/MarchingCubes.h"

TerrainGenerationInterface::TerrainGenerationInterface(QWidget *parent)
    : ActionInterface("terraingeneration", "Terrain generation", "model", "Terrain generation management", "terrain_generation_manager_button.png", parent)
{
//    this->createGUI();
}

void TerrainGenerationInterface::setWaterLevel(float newLevel)
{
    this->waterLevel = newLevel;
    voxelGrid->updateEnvironmentalDensities(newLevel);
    Q_EMIT waterLevelChanged(newLevel);
    Q_EMIT updated();
}

void TerrainGenerationInterface::setAmbiantOcclusion(float newValue)
{
    this->ambiantOcclusionFactor = newValue;
    Q_EMIT updated();
}

void TerrainGenerationInterface::updateDisplayedView(const Vector3& newVoxelGridOffset, float newVoxelGridScaling)
{
    this->voxelGridOffset = newVoxelGridOffset;
    this->voxelGridScaling = newVoxelGridScaling;

//    heightmapMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
//    heightmapMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);
    marchingCubeMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
    marchingCubeMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);

    implicitMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
    implicitMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);

    layersMesh.shader->setVector("subterrainOffset", this->voxelGridOffset);
    layersMesh.shader->setFloat("subterrainScale", this->voxelGridScaling);
}

void TerrainGenerationInterface::afterTerrainUpdated()
{

}

void TerrainGenerationInterface::heightmapToVoxels()
{
    voxelGrid->from2DGrid(*heightmap);
}

void TerrainGenerationInterface::heightmapToLayers()
{
    layerGrid->from2DGrid(*heightmap);
}

void TerrainGenerationInterface::heightmapToImplicit()
{
    implicitTerrain->composables = {ImplicitPrimitive::fromHeightmap(heightmap->getHeights())};
    implicitTerrain->_cached = false;
}

void TerrainGenerationInterface::heightmapToAll()
{
    heightmapToVoxels();
    heightmapToLayers();
    heightmapToImplicit();
}

void TerrainGenerationInterface::voxelsToHeightmap()
{
    heightmap->fromVoxelGrid(*voxelGrid);
}

void TerrainGenerationInterface::voxelsToLayers()
{
    layerGrid->fromVoxelGrid(*voxelGrid);
}

void TerrainGenerationInterface::voxelsToImplicit()
{
    implicitTerrain->composables = {ImplicitPrimitive::fromHeightmap(voxelGrid->getVoxelValues())}; // Yeah, I can pass a 3D grid in this function
    implicitTerrain->_cached = false;
}

void TerrainGenerationInterface::voxelsToAll()
{
    voxelsToHeightmap();
    voxelsToLayers();
//    heightmapToLayers();
    voxelsToImplicit();
}

void TerrainGenerationInterface::layersToVoxels()
{
    voxelGrid->fromLayerBased(*layerGrid);
}

void TerrainGenerationInterface::layersToHeightmap()
{
    heightmap->fromLayerGrid(*layerGrid);
}

void TerrainGenerationInterface::layersToImplicit()
{
    // ..?!
}

void TerrainGenerationInterface::layersToAll()
{
    layersToVoxels();
    layersToHeightmap();
    layersToImplicit();
}

void TerrainGenerationInterface::implicitToVoxels()
{
    voxelGrid->fromImplicit(implicitTerrain.get());
}

void TerrainGenerationInterface::implicitToLayers()
{
    layerGrid->reset();
    layerGrid->add(implicitTerrain.get());
}

void TerrainGenerationInterface::implicitToHeightmap()
{
    // ..?!
}

void TerrainGenerationInterface::implicitToAll()
{
    implicitToHeightmap();
    implicitToVoxels();
    implicitToLayers();
}

void TerrainGenerationInterface::setVisu(MapMode _mapMode, SmoothingAlgorithm _smoothingAlgorithm, bool _displayParticles)
{
    this->mapMode = _mapMode;
    this->smoothingAlgorithm = _smoothingAlgorithm;
//    this->displayParticles = _displayParticles;
}

void TerrainGenerationInterface::displayWaterLevel()
{
    waterLevelMesh.fromArray({
                                Vector3(1.f, 1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor),
                                Vector3(voxelGrid->getSizeX()-1.f, 1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor),
                                Vector3(1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor),
                                 Vector3(voxelGrid->getSizeX()-1.f, 1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor),
                                Vector3(voxelGrid->getSizeX()-1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor),
                                Vector3(1.f, voxelGrid->getSizeY()-1.f, waterLevel * voxelGrid->getSizeZ() * heightFactor)
                             });
    waterLevelMesh.display();
}

void TerrainGenerationInterface::createTerrainFromNoise(int nx, int ny, int nz, bool noise2D, float noiseStrength, float frequency, float lacunarity, float noise_shifting)
{
//    std::cout << noiseStrength << " " << frequency << " " << lacunarity << std::endl;
    GridF values(nx, ny, nz);

    // Create and configure FastNoise object
    if (frequency > 0.01f) {
        FastNoiseLite noise;
        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFrequency(1.f / (float) (values.sizeX * frequency));
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        noise.SetFractalLacunarity(lacunarity);
        noise.SetFractalGain(0.7);
        noise.SetFractalWeightedStrength(0.5);
        noise.SetFractalOctaves(10);

        if (noise2D) {
            for (int x = 0; x < values.sizeX; x++) {
                for (int y = 0; y < values.sizeY; y++) {
                    float grid_height = noise.GetNoise((float)x, (float)y) * noiseStrength * values.sizeZ + this->waterLevel * values.sizeZ;
                    int z = int(std::max(grid_height, 2.f));
                    // Positive values
                    for (int i = 0; i < int(z); i++) {
                        values.at(x, y, i) = .5f;
                    }
                    if (z < values.sizeZ) {
                        values.at(x, y, z) = interpolation::inv_linear(z - int(z), -.5f, .5f);
                        for (int i = z+1; i < values.sizeZ; i++) {
                            values.at(x, y, i) = -.5f;
                        }
                    }
                }
            }
            values = values.meanSmooth();
        } else {
            for (int x = 0; x < values.sizeX; x++) {
                for (int y = 0; y < values.sizeY; y++) {
                    for (int z = 0; z < values.sizeZ; z++) {
                        values.at(x, y, z) = noise.GetNoise((float)x, (float)y, (float)z) * noiseStrength;
                        if (z < this->waterLevel * values.sizeZ && (x > 5))
                            values.at(x, y, z) = std::abs(values.at(x, y, z));
//                        else
//                            values.at(x, y, z) = -std::abs(values.at(x, y, z));
                    }
                }
            }
        }
    } else {
        for (int x = 0; x < values.sizeX; x++) {
            for (int y = 0; y < values.sizeY; y++) {
                for (int z = 0; z < values.sizeZ; z++) {
                    values.at(x, y, z) = (z < this->waterLevel * values.sizeZ ? .1f : -.1f);
                }
            }
        }
    }
    std::cout << frequency << std::endl;

    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();
    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>();
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>();

//    *voxelGrid = tempMap;
//    voxelGrid->_cachedVoxelValues = values;
//    voxelGrid->fromCachedData();
    voxelGrid->setVoxelValues(values);
    this->heightmap->fromVoxelGrid(*voxelGrid);
    this->layerGrid->fromVoxelGrid(*voxelGrid);
    this->initialTerrainValues = values;

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
    Q_EMIT this->updated();
}

void TerrainGenerationInterface::reloadTerrain(std::map<std::string, std::shared_ptr<ActionInterface> > actionInterfaces)
{
    this->createTerrainFromFile(this->lastLoadedMap, actionInterfaces);
    this->setWaterLevel(this->waterLevel);
}

void TerrainGenerationInterface::createTerrainFromFile(std::string filename, std::map<std::string, std::shared_ptr<ActionInterface> > actionInterfaces)
{    
    this->lastLoadedMap = filename;
    std::string ext = toUpper(getExtension(filename));

//    Vector3 terrainSize = Vector3(200, 200, 50); //Vector3(128, 128, 64);
    Vector3 terrainSize = Vector3(200, 100, 50); //Vector3(128, 128, 64);
    if (this->voxelGrid)
        terrainSize = voxelGrid->getDimensions();

    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>(terrainSize.x, terrainSize.y, terrainSize.z);
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>(terrainSize.x, terrainSize.y, terrainSize.z);
    if (!this->layerGrid)
        this->layerGrid = std::make_shared<LayerBasedGrid>(terrainSize.x, terrainSize.y, 0.f);
    if (!this->implicitTerrain)
        this->implicitTerrain = std::make_shared<ImplicitNaryOperator>();

    if (ext == "PGM" || ext == "PNG" || ext == "JPG" || ext == "PNG" || ext == "TGA" || ext == "BMP" || ext == "PSD" || ext == "GIF" || ext == "HDR" || ext == "PIC") {
        // From heightmap
        std::cout << "Heightmap : " << showTime(timeIt([=]() { heightmap->loadFromHeightmap(filename, terrainSize.x, terrainSize.y, terrainSize.z); })) << std::endl;

        std::cout << "Voxels : " << showTime(timeIt([=]() { voxelGrid->from2DGrid(*heightmap); })) << std::endl;
//        layerGrid->fromVoxelGrid(*voxelGrid);
        std::cout << "Layers : " << showTime(timeIt([=]() { layerGrid->from2DGrid(*heightmap); })) << std::endl;

    } else if (ext == "JSON") {
        // The JSON file contains the list of actions made on a map
        std::ifstream file(filename);
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        nlohmann::json json_content = nlohmann::json::parse(content);
        if (!json_content.contains("actions")) {
            if (json_content.contains(ImplicitPatch::json_identifier)) {
                this->createTerrainFromImplicitPatches(json_content);
                this->initialTerrainValues = this->voxelGrid->getVoxelValues();
                return;
            } else {
                this->createTerrainFromBiomes(json_content);
            }
        } else {
            for (const auto &action : json_content.at("actions")) {
                // Let all the interfaces try to replay their actions
                for (auto& possibleAction : actionInterfaces)
                    possibleAction.second->replay(action);
            }
        }
    } else if (ext == "STL") {
        Mesh m;
        m.fromStl(filename);
        FastNoiseLite noise;
        GridF newVoxelValues = (GridF(m.voxelize(voxelGrid->getDimensions() * 2.f)) - .5f).resize(voxelGrid->getDimensions());
//        for (int x = 0; x < newVoxelValues.sizeX; x++) {
//            for (int y = 0; y < newVoxelValues.sizeY; y++) {
//                for (int z = 0; z < 13; z++) {
//                    float noiseVal = abs(noise.GetNoise((float) x, (float) y));
//                    //newVoxelValues(x, y, z) = noiseVal;
//                    newVoxelValues(x, y, z) = noiseVal + .2f;
//                }
//            }
//        }
        GridF swapXY(newVoxelValues.getDimensions().yxz());
        swapXY.iterateParallel([&](int x, int y, int z) {
            swapXY(x, y, z) = newVoxelValues(y, x, z);
        });
//        if (m.isWatertight())
            voxelGrid->setVoxelValues(swapXY.flip(false, true)); //(newVoxelValues);
//        else
//            voxelGrid->setVoxelValues(m.voxelizeSurface(voxelGrid->getDimensions()));
//        voxelGrid->fromCachedData();
        heightmap->fromVoxelGrid(*voxelGrid);
        layerGrid->fromVoxelGrid(*voxelGrid);
    } else if (ext == "FBX") {
        Mesh m;
        if (m.isWatertight())
            voxelGrid->setVoxelValues(m.voxelize(voxelGrid->getDimensions()));
        else
            voxelGrid->setVoxelValues(m.voxelizeSurface(voxelGrid->getDimensions()));
//        voxelGrid->setVoxelValues(m.voxelize(voxelGrid->getDimensions()));
//        voxelGrid->fromCachedData();
        heightmap->fromVoxelGrid(*voxelGrid);
        layerGrid->fromVoxelGrid(*voxelGrid);
    }else /*if (ext == "DATA" || ext.empty())*/ {
        // Then it's our custom voxel grid file
        voxelGrid->retrieveMap(filename);
//        voxelGrid->smoothVoxels();
//        voxelGrid->smoothVoxels();
//        heightmap->fromVoxelGrid(*voxelGrid);
//        layerGrid->fromVoxelGrid(*voxelGrid);
        voxelsToAll();

    } /*else {
        // In any other case, consider that nothing has been done, cancel.
        return;
    }*/
    ImplicitPrimitive* implicitHeightmap = ImplicitPrimitive::fromHeightmap(heightmap->heights, ""); // , dynamic_cast<ImplicitPrimitive*>(implicitTerrain.get()));
    implicitHeightmap->material = TerrainTypes::DIRT;
    this->implicitTerrain->composables = {implicitHeightmap};
    this->implicitTerrain->_cached = false;
    this->addTerrainAction(nlohmann::json({
                                              {"from_file", filename},
                                              {"from_noise", false}
                                          }));
    this->setWaterLevel(this->waterLevel);

    this->initialTerrainValues = this->voxelGrid->getVoxelValues();

    Q_EMIT this->updated();
}

void TerrainGenerationInterface::createTerrainFromBiomes(nlohmann::json json_content)
{
    if (!this->heightmap)
        this->heightmap = std::make_shared<Heightmap>();
    if (!this->voxelGrid)
        this->voxelGrid = std::make_shared<VoxelGrid>();

    *heightmap = Heightmap(124, 124, 40.f);
    heightmap->raise(GridF(heightmap->getSizeX(), heightmap->getSizeY(), 1.f, 20.f));
//    this->heightmap->heights = GridF(124, 124, 1, 20.f);
//    this->heightmap->maxHeight = 40.f;

    this->voxelGrid->from2DGrid(*this->heightmap);
//    this->voxelGrid->fromCachedData();

    this->biomeGenerationNeeded = true;
    this->biomeGenerationModelData = json_content;

}

void TerrainGenerationInterface::createTerrainFromImplicitPatches(nlohmann::json json_content)
{
    std::cout << "Load implicit: " << showTime(timeIt([&]() {
        *this->implicitTerrain = *dynamic_cast<ImplicitNaryOperator*>(ImplicitPatch::fromJson(json_content[ImplicitPatch::json_identifier]));
    })) << std::endl;
/*
    EnvArea* reef = dynamic_cast<EnvArea*>(EnvObject::instantiate("reef"));
    reef->area = ShapeCurve({Vector3(0, 0, 0), Vector3(100, 0, 0), Vector3(100, 100, 0), Vector3(0, 100, 0)});

    EnvCurve* passe = dynamic_cast<EnvCurve*>(EnvObject::instantiate("passe"));
    passe->curve = BSpline({Vector3(10, 10, 0), Vector3(40, 10, 0), Vector3(40, 50, 0), Vector3(70, 70, 0)});

    this->implicitTerrain->addChild(reef->createImplicitPatch());
    this->implicitTerrain->addChild(passe->createImplicitPatch());
*/
/*
    std::cout << "To layers   : " << timeIt([&](){ this->layerGrid->fromImplicit(implicitTerrain.get()); }) << "ms" << std::endl;
    std::cout << "To voxels   : " << timeIt([&](){ this->voxelGrid->fromLayerBased(*layerGrid); }) << "ms" << std::endl;
    std::cout << "To heightmap: " << timeIt([&](){ this->heightmap->fromLayerGrid(*layerGrid); }) << "ms" << std::endl;
    */
//    std::cout << "To heightmap: " << timeIt([&](){ this->heightmap->fromImplicit(implicitTerrain.get()); }) << "ms" << std::endl;

    std::cout << "To voxels: " << showTime(timeIt([&](){ this->voxelGrid->fromImplicit(implicitTerrain.get()); })) << std::endl;
    std::cout << "To layers: " << showTime(timeIt([&](){ this->layerGrid->fromVoxelGrid(*voxelGrid); })) << std::endl;
//    std::cout << "To heightmap: " << timeIt([&](){ this->heightmap->fromImplicit(implicitTerrain.get()); }) << "ms" << std::endl;
    std::cout << "To heightmap: " << showTime(timeIt([&](){ this->heightmap->fromVoxelGrid(*voxelGrid.get()); })) << std::endl;

}

void TerrainGenerationInterface::saveTerrain(std::string filename, Vector3 dimensions)
{
    std::string ext = toUpper(getExtension(filename));
    if (ext == "PNG" || ext == "JPG" || ext == "TGA" || ext == "BMP" || ext == "HDR") {
        // To heightmap
        this->heightmap->fromVoxelGrid(*voxelGrid); // Just to be sure to have the last values
        this->heightmap->saveHeightmap(filename, dimensions);
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
    QLayout* layout = new QHBoxLayout;

    QLabel* heightmapPathLabel = new QLabel(QString::fromStdString(getFilename(this->lastLoadedMap)));
    QPushButton* loadHeightmapButton = new QPushButton("Load");
    QPushButton* reloadButton = new QPushButton("Reload");
    QPushButton* saveHeightmapButton = new QPushButton("Save");
    QDoubleSpinBox* widthEdit = new QDoubleSpinBox(); // (QString::fromStdString(std::to_string(int(voxelGrid->getSizeX()))));
    QDoubleSpinBox* depthEdit = new QDoubleSpinBox(); //(QString::fromStdString(std::to_string(int(voxelGrid->getSizeY()))));
    QDoubleSpinBox* heightEdit = new QDoubleSpinBox(); // (QString::fromStdString(std::to_string(int(voxelGrid->getSizeZ()))));

    FancySlider* noiseStrengthSlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 1.f, 0.01f);
    FancySlider* noiseLacunaritySlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 4.f, 0.01f);
    FancySlider* noiseFrequencySlider = new FancySlider(Qt::Orientation::Horizontal, 0.f, 3.f, 0.01f);
    QPushButton* createFromNoiseButton = new QPushButton("Noise");

    QRadioButton* noise2D = new QRadioButton("2D");
    QRadioButton* noise3D = new QRadioButton("3D");

    layout->addWidget(createVerticalGroup({
                                              createHorizontalGroup({heightmapPathLabel, loadHeightmapButton}),
                                              createHorizontalGroup({reloadButton, saveHeightmapButton}),
                                              createHorizontalGroup({
                                                  widthEdit, new QLabel("x"), depthEdit, new QLabel("x"), heightEdit
                                              }),

                                              createMultipleSliderGroup({
//                                                  {"strength", noiseStrengthSlider},
                                                  {"frequency", noiseFrequencySlider},
                                                  {"lacunarity", noiseLacunaritySlider}
                                              }),
//                                              createFromNoiseButton
                                              createHorizontalGroup({
                                                  noise2D, noise3D
                                              })
                                          }));
    widthEdit->setDecimals(0);
    depthEdit->setDecimals(0);
    heightEdit->setDecimals(0);

    widthEdit->setValue(voxelGrid->getSizeX());
    depthEdit->setValue(voxelGrid->getSizeX());
    heightEdit->setValue(voxelGrid->getSizeX());

    noiseStrengthSlider->setfValue(1.f);
    noiseLacunaritySlider->setfValue(2.f);
    noiseFrequencySlider->setfValue(1.f);

    noise2D->setChecked(false);
    noise3D->setChecked(true);

    QObject::connect(loadHeightmapButton, &QPushButton::pressed, this, [=]() { this->openMapUI(); });
    QObject::connect(reloadButton, &QPushButton::pressed, this, [=]() { this->reloadTerrain(this->actionInterfaces); });
    QObject::connect(saveHeightmapButton, &QPushButton::pressed, this, [=]() { this->saveMapUI(); });
//    QObject::connect(createFromNoiseButton, &QPushButton::pressed, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(widthEdit, &QDoubleSpinBox::editingFinished, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(depthEdit, &QDoubleSpinBox::editingFinished, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(heightEdit, &QDoubleSpinBox::editingFinished, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(noiseFrequencySlider, &FancySlider::floatValueChanged, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(noiseStrengthSlider, &FancySlider::floatValueChanged, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(noiseLacunaritySlider, &FancySlider::floatValueChanged, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), noise2D->isChecked(), noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(noise2D, &QRadioButton::pressed, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), true, noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    QObject::connect(noise3D, &QRadioButton::pressed, this, [=]() { this->createTerrainFromNoise(widthEdit->text().toInt(), depthEdit->text().toInt(), heightEdit->text().toInt(), false, noiseStrengthSlider->getfValue(), noiseFrequencySlider->getfValue(), noiseLacunaritySlider->getfValue()); });
    return layout;
}




void TerrainGenerationInterface::prepareShader(bool reload)
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

    GridI materials; GridF matHeights;
    std::tie(materials, matHeights) = layerGrid->getMaterialAndHeightsGrid();
    layersMesh.shader->setTexture3D("matIndicesTex", 1, materials);
    layersMesh.shader->setTexture3D("matHeightsTex", 2, matHeights);

    std::vector<Vector3> layersPoints(materials.size());
    for (size_t i = 0; i < layersPoints.size(); i++) {
        layersPoints[i] = materials.getCoordAsVector3(i);
    }
    layersMesh.fromArray(layersPoints);
    layersMesh.update();

    implicitMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_mc_voxels),
                            true, GL_TRIANGLES);

    marchingCubeMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_mc_voxels),
                            true, GL_TRIANGLES);
    GridF isoData = this->voxelGrid->getVoxelValues(); //.resize(15, 15, 15);
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
    marchingCubeMesh.shader->setInt("dataChangesFieldTex", 3);
    marchingCubeMesh.shader->setFloat("isolevel", 0.f);
    marchingCubeMesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    marchingCubeMesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    marchingCubeMesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));


    implicitMesh.useIndices = false;
    implicitMesh.fromArray(points);
    implicitMesh.update();
    implicitMesh.shader->use();
    GlobalsGL::f()->glBindVertexArray(implicitMesh.vao);
    implicitMesh.shader->setInt("dataFieldTex", 0);
    implicitMesh.shader->setInt("edgeTableTex", 1);
    implicitMesh.shader->setInt("triTableTex", 2);
    implicitMesh.shader->setFloat("isolevel", 0.f);
    implicitMesh.shader->setVector("vertDecals[0]", Vector3(0.0, 0.0, 0.0));
    implicitMesh.shader->setVector("vertDecals[1]", Vector3(1.0, 0.0, 0.0));
    implicitMesh.shader->setVector("vertDecals[2]", Vector3(1.0, 1.0, 0.0));
    implicitMesh.shader->setVector("vertDecals[3]", Vector3(0.0, 1.0, 0.0));
    implicitMesh.shader->setVector("vertDecals[4]", Vector3(0.0, 0.0, 1.0));
    implicitMesh.shader->setVector("vertDecals[5]", Vector3(1.0, 0.0, 1.0));
    implicitMesh.shader->setVector("vertDecals[6]", Vector3(1.0, 1.0, 1.0));
    implicitMesh.shader->setVector("vertDecals[7]", Vector3(0.0, 1.0, 1.0));
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
    implicitMesh.shader->setTexture3D("dataFieldTex", 0, voxelGrid->getVoxelValues() / 6.f + .5f);


    this->heightmapMesh = Mesh(std::make_shared<Shader>(vShader_mc_voxels, fShader_mc_voxels, gShader_grid), true, GL_POINTS);

    marchingCubeMesh.shader->setTexture3D("dataChangesFieldTex", 0, GridF(voxelGrid->getVoxelValues().getDimensions()));
    heightmapMesh.shader->setTexture3D("dataChangesFieldTex", 0, GridF(voxelGrid->getVoxelValues().getDimensions().x, voxelGrid->getVoxelValues().getDimensions().y, 1));


    GridF heightData(this->heightmap->getSizeX(), this->heightmap->getSizeY());
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
    implicitMesh.shader->setInt("heightmapFieldTex", 3);
    layersMesh.shader->setInt("heightmapFieldTex", 3);

    GridF heights = heightmap->getHeights();
    float *heightmapData = new float[heights.size() * 4];
    GridV3 gradients = heights.gradient();
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
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmap->getSizeX(), heightmap->getSizeY(), 0,
    GL_RGBA, GL_FLOAT, heightmapData);
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
    implicitMesh.shader->setInt("biomeFieldTex", 4);
    layersMesh.shader->setInt("biomeFieldTex", 4);

    if (verbose)
        std::cout << "Loading terrain textures... " << std::flush;
//    QTemporaryDir tempTextureDir;
    QDirIterator iTex("src/assets/textures/", QDir::Files, QDirIterator::Subdirectories);
//    QDirIterator iTex(":/terrain_textures/", QDir::Files, QDirIterator::Subdirectories);
    GridV3 allColorTextures;
    GridV3 allNormalTextures;
    GridV3 allDisplacementTextures;
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
    std::vector<GridV3> texturesMaps(nbFiles);
    int textureSizes = 128;

    float parallelTime = 0;
    float sequentialTime = 0;
    float totalTime = timeIt([&]() {
        parallelTime = timeIt([&]() {
            #pragma omp parallel for
            for (size_t i = 0; i < nbFiles; i++) {
                QString dir = filesToLoad[i];
                QString basename = QFileInfo(dir).baseName();
                    std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
                    int imgW, imgH, c;
                    unsigned char *data = stbi_load(dir.toStdString().c_str(), &imgW, &imgH, &c, 3);
                    texturesMaps[i] = GridV3(imgW, imgH);
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
            allColorTextures = GridV3(textureSizes * colorTexturesIndex.size(), textureSizes);
            allNormalTextures = GridV3(textureSizes * colorTexturesIndex.size(), textureSizes);
            allDisplacementTextures = GridV3(textureSizes * colorTexturesIndex.size(), textureSizes);
        });
        sequentialTime = timeIt([&]() {
            for (size_t i = 0; i < nbFiles; i++) {
                QString dir = filesToLoad[i];
                QString basename = QFileInfo(dir).baseName();
                std::string textureClass = toLower(QFileInfo(QFileInfo(dir).absolutePath()).baseName().toStdString());
                GridV3& map = texturesMaps[i];
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
        });
    });
    if (verbose) {
        std::cout << "Done in " << showTime(totalTime) << std::endl;
        std::cout << "( " << showTime(parallelTime) << " for the parallel and " << showTime(sequentialTime) << " for sequential)" << std::endl;
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
    implicitMesh.shader->setInt("allBiomesColorTextures", 5);
    layersMesh.shader->setInt("allBiomesColorTextures", 5);

    heightmapMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);
    implicitMesh.shader->setInt("maxBiomesColorTextures", colorTexturesIndex.size());//indexColorTextureClass);
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
    implicitMesh.shader->setInt("allBiomesNormalTextures", 6);
    layersMesh.shader->setInt("allBiomesNormalTextures", 6);

    heightmapMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);
    marchingCubeMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);
    implicitMesh.shader->setInt("maxBiomesNormalTextures", colorTexturesIndex.size());//indexNormalTextureClass);
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
    implicitMesh.shader->setInt("allBiomesDisplacementTextures", 7);
    layersMesh.shader->setInt("allBiomesDisplacementTextures", 7);

    heightmapMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());
    marchingCubeMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());
    implicitMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());
    layersMesh.shader->setInt("maxBiomesDisplacementTextures", colorTexturesIndex.size());

    updateDisplayedView(voxelGridOffset, voxelGridScaling);

    if (verbose)
        std::cout << "Terrain shaders and assets ready." << std::endl;
}

GridF getVoxelChanges(std::shared_ptr<VoxelGrid> voxels, GridF initial) {
    return voxels->getVoxelValues() - initial;
}
GridF getHeightmapChanges(std::shared_ptr<VoxelGrid> voxels, GridF initial) {
    auto diff = getVoxelChanges(voxels, initial);

    GridF map(diff.sizeX, diff.sizeY);
    for (int x = 0; x < diff.sizeX; x++) {
        for (int y = 0; y < diff.sizeY; y++) {
            float sum = 0;
            for (int z = 0; z < diff.sizeZ; z++) {
                sum += diff.at(x, y, z);
            }
            if (sum > 0) {
                map.at(x, y) = std::max(sum, 0.f);
            }
            else {
                map.at(x, y) = -std::max(-sum, 0.f);
            }
        }
    }
    return map;
}


void TerrainGenerationInterface::display(const Vector3& camPos)
{
    float maxHeight;
    float initTime;
    float displayTime;
    float meshCreationTime;
    float GLcallTime;
    float realDisplayTime;

    initTime = timeIt([&]() {
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_2D, allBiomesColorTextures);
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE6);
        glBindTexture(GL_TEXTURE_2D, allBiomesNormalTextures);
        GlobalsGL::f()->glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_2D, allBiomesDisplacementTextures);

        GlobalsGL::f()->glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, heightmapFieldTex);

        GridF heights = heightmap->getHeights();
        maxHeight = heights.max();
        float *heightmapData = new float[heights.size() * 4];
        GridV3 gradients = heights.gradient();
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
            gradients[i] = (gradients[i].normalized().abs());// + Vector3(1, 1, 1)) / 2.f;
            heightmapData[i * 4 + 0] = colorTextureOffset;
            heightmapData[i * 4 + 1] = normalTextureOffset; // gradients[i].y;
            heightmapData[i * 4 + 2] = displacementTextureOffset; // gradients[i].z;
            heightmapData[i * 4 + 3] = heights[i] / maxHeight; //1.0; //heightmap->heights[i] / maxHeight;
        }

        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, heightmap->getSizeX(), heightmap->getSizeY(), 0,
        GL_RGBA, GL_FLOAT, heightmapData);//heightData.data.data());
        delete[] heightmapData;
    });

    displayTime = timeIt([&]() {
        if (mapMode == GRID_MODE) {
            if (this->heightmap == nullptr) {
                std::cerr << "No grid to display" << std::endl;
            } else {
    //            float time = std::chrono::duration<float>(std::chrono::system_clock::now() - startingTime).count();

                GridF heightData(this->heightmap->getSizeX(), this->heightmap->getSizeY());
                std::vector<Vector3> positions(heightData.size());
                for (size_t i = 0; i < positions.size(); i++) {
                    positions[i] = heightData.getCoordAsVector3(i);
                    heightData[i] = heightmap->getHeight(positions[i].x, positions[i].y);
                }
                heightmapMesh.shader->setFloat("maxHeight", maxHeight);
                heightmapMesh.shader->setFloat("waterRelativeHeight", waterLevel);
                heightmapMesh.shader->setFloat("ambiantOcclusionFactor", ambiantOcclusionFactor);
                heightmapMesh.shader->setBool("displayAsComparisonTerrain", displayAsComparativeMode);
                heightmapMesh.shader->setFloat("heightFactor", heightFactor);
                if (voxelGrid->getDimensions() == initialTerrainValues.getDimensions())
                    heightmapMesh.shader->setTexture3D("dataChangesFieldTex", 3, getHeightmapChanges(voxelGrid, initialTerrainValues) + 2.f);
                heightmapMesh.fromArray(positions);
                heightmapMesh.update();
                this->heightmapMesh.display(GL_POINTS);
            }
        }
        else if (mapMode == VOXEL_MODE) {
            if (this->voxelGrid == nullptr) {
                std::cerr << "No voxel grid to display" << std::endl;
            } else {
                GridF values;
                meshCreationTime = timeIt([&]() {
                    values = voxelGrid->getVoxelValues().meanSmooth(3, 3, 3, true);//.meanSmooth(3, 3, 3, true);
                    for (int x = 0; x < values.sizeX; x++) {
                        for (int y = 0; y < values.sizeY; y++) {
                            for (int z = 0; z < 2; z++) {
                                values(x, y, z) = std::max(std::abs(values(x, y, z)), .0f);
                            }
                        }
                    }
                    if (marchingCubeMesh.vertexArray.size() != values.size()) {
                        std::vector<Vector3> points(values.size());
                        for (size_t i = 0; i < points.size(); i++) {
                            values[i] = (std::max(-3.f, std::min(3.f, values[i])) / 6.f) + 0.5;
                            points[i] = values.getCoordAsVector3(i);
                        }
                        marchingCubeMesh.useIndices = false;
                        marchingCubeMesh.fromArray(points);
                    }
                });
                GLcallTime = timeIt([&]() {
                    marchingCubeMesh.shader->setTexture3D("dataFieldTex", 0, values + .5f);
        //            marchingCubeMesh.shader->setTexture3D("dataChangesFieldTex", 3, getVoxelChanges(voxelGrid, initialTerrainValues) + 2.f);
                    marchingCubeMesh.shader->setTexture3D("dataChangesFieldTex", 3, (EnvObject::sandDeposit * 10.f + 2.f).resize(values.getDimensions().xy() + Vector3(0, 0, 1)));

                    marchingCubeMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
                    marchingCubeMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
                    marchingCubeMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
                    marchingCubeMesh.shader->setFloat("waterRelativeHeight", waterLevel);
                    marchingCubeMesh.shader->setFloat("ambiantOcclusionFactor", ambiantOcclusionFactor);
                    marchingCubeMesh.shader->setBool("displayAsComparisonTerrain", displayAsComparativeMode);
                    marchingCubeMesh.shader->setFloat("heightFactor", heightFactor);
                });
                realDisplayTime = timeIt([&]() { marchingCubeMesh.display( GL_POINTS ); });
                if (smoothingAlgorithm == SmoothingAlgorithm::NONE) {
                    GLcallTime = timeIt([&]() {
                        marchingCubeMesh.shader->setBool("displayingIgnoredVoxels", true);
                        marchingCubeMesh.shader->setFloat("min_isolevel", -1000.f);
                        marchingCubeMesh.shader->setFloat("max_isolevel", this->minIsoLevel/3.f);
                        marchingCubeMesh.display( GL_POINTS );
                        marchingCubeMesh.shader->setFloat("min_isolevel", this->maxIsoLevel/3.f);
                        marchingCubeMesh.shader->setFloat("max_isolevel",  1000.f);
                    });
                    realDisplayTime += timeIt([&]() { marchingCubeMesh.display( GL_POINTS ); });
                    GLcallTime += timeIt([&]() { marchingCubeMesh.shader->setBool("displayingIgnoredVoxels", false); });
                }
                // Check if something changed on the terrain :
                if (this->voxelsPreviousHistoryIndex != voxelGrid->getCurrentHistoryIndex()) {
                    this->voxelsPreviousHistoryIndex = voxelGrid->getCurrentHistoryIndex();
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
                layersMesh.shader->setFloat("ambiantOcclusionFactor", ambiantOcclusionFactor);
                layersMesh.shader->setBool("displayAsComparisonTerrain", displayAsComparativeMode);
                layersMesh.shader->setFloat("heightFactor", heightFactor);

                //if (this->layersPreviousHistoryIndex != layerGrid->_historyIndex || layerGrid->_historyIndex == -1) {
                    this->layersPreviousHistoryIndex = layerGrid->_historyIndex;
                    GridI materials;
                    GridF matHeights;
                    std::tie(materials, matHeights) = layerGrid->getMaterialAndHeightsGrid();
                    layersMesh.shader->setTexture3D("matIndicesTex", 1, materials);
                    layersMesh.shader->setTexture3D("matHeightsTex", 2, matHeights);

                    std::vector<Vector3> layersPoints(materials.size());
                    for (size_t i = 0; i < layersPoints.size(); i++) {
                        layersPoints[i] = materials.getCoordAsVector3(i);
                    }
                    layersMesh.fromArray(layersPoints);
                    layersMesh.update();
    //            }
                this->layersMesh.display(GL_POINTS);
            }
        } else if (mapMode == IMPLICIT_MODE) {
            if (this->implicitTerrain == nullptr) {
                std::cerr << "No implicit terrain to display" << std::endl;
            } else {
                GridF values;
                values = implicitTerrain->getVoxelized(voxelGrid->getDimensions());
    //            std::cout << values.sum() << std::endl;
                implicitMesh.shader->setTexture3D("dataFieldTex", 0, values + .5f);
                implicitMesh.shader->setBool("useMarchingCubes", smoothingAlgorithm == SmoothingAlgorithm::MARCHING_CUBES);
                implicitMesh.shader->setFloat("min_isolevel", this->minIsoLevel/3.f);
                implicitMesh.shader->setFloat("max_isolevel", this->maxIsoLevel/3.f);
                implicitMesh.shader->setFloat("waterRelativeHeight", waterLevel);
                implicitMesh.shader->setFloat("ambiantOcclusionFactor", ambiantOcclusionFactor);
                implicitMesh.shader->setBool("displayAsComparisonTerrain", displayAsComparativeMode);
                implicitMesh.shader->setFloat("heightFactor", heightFactor);
                implicitMesh.display( GL_POINTS );
                if (smoothingAlgorithm == SmoothingAlgorithm::NONE) {
                    implicitMesh.shader->setBool("displayingIgnoredVoxels", true);
                    implicitMesh.shader->setFloat("min_isolevel", -1000.f);
                    implicitMesh.shader->setFloat("max_isolevel", this->minIsoLevel/3.f);
                    implicitMesh.display( GL_POINTS );
                    implicitMesh.shader->setFloat("min_isolevel", this->maxIsoLevel/3.f);
                    implicitMesh.shader->setFloat("max_isolevel",  1000.f);
                    implicitMesh.display( GL_POINTS );
                    implicitMesh.shader->setBool("displayingIgnoredVoxels", false);
                }
            }
        }
    });
//    std::cout << "Init    : " << showTime(initTime) << "\nDisplay : " << showTime(displayTime) << "\n(mesh: " << showTime(meshCreationTime) << ", GL calls: " << showTime(GLcallTime) << ", display: " << showTime(realDisplayTime) <<  ")" << std::endl;
}



void TerrainGenerationInterface::openMapUI()
{
    QString q_filename = QFileDialog::getOpenFileName(this, QString("Ouvrir une carte"), QString::fromStdString(this->mapSavingFolder));
    this->createTerrainFromFile(q_filename.toStdString(), this->actionInterfaces);
//    this->viewer->setSceneCenter(viewer->voxelGrid->getDimensions() / 2.f);
//    this->terrainGenerationInterface->prepareShader(true);
}

void TerrainGenerationInterface::saveMapUI()
{
    QString q_filename = QFileDialog::getSaveFileName(this, QString("Enregistrer la carte"), QString::fromStdString(this->mapSavingFolder));
    this->saveTerrain(q_filename.toStdString());
}

void TerrainGenerationInterface::reinforceVoxels()
{
    GridF distances = this->voxelGrid->getVoxelValues().binarize().toDistanceMap(false, false);

    for (auto& v : distances)
        v = std::max(1.f, v);
    this->voxelGrid->setVoxelValues(voxelGrid->getVoxelValues() * distances);
}

void TerrainGenerationInterface::saveErosionDepositionTextureMasks(std::string savingFolder, std::string savingName)
{
    std::string terrainFilename = savingFolder + "/" + savingName + "_height.png";
    std::string erosionFilename = savingFolder + "/" + savingName + "_erod.png";
    std::string depositFilename = savingFolder + "/" + savingName + "_depo.png";

    Vector3 dimensions = Vector3(256, 256, 1);

    this->saveTerrain(terrainFilename, dimensions);
    auto diff = getHeightmapChanges(voxelGrid, initialTerrainValues).resize(dimensions) / 2.f;


    GridF erosionMap = diff;
    GridF depositMap = diff;

    for (size_t i = 0; i < diff.size(); i++) {
        erosionMap[i] = std::min(-diff[i], 1.f);
        depositMap[i] = std::min(diff[i], 1.f);
    }
    Image(erosionMap).writeToFile(erosionFilename);
    Image(depositMap).writeToFile(depositFilename);
}

void TerrainGenerationInterface::saveErosionDepositionTextureMasksOnMultiple()
{
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::ExistingFiles);
    QStringList fileNames;
    if (dialog.exec()) {
        fileNames = dialog.selectedFiles();
        for (size_t i = 0; i < fileNames.size(); i++) {
            QString& q_filename = fileNames[i];
            std::string filename = q_filename.toStdString();

            auto path = split(filename, "/");
            std::string basename = split(path.back(), ".")[0];

            path.pop_back();

            std::string folder = join(path, "/") + "/heightmapsAndMasks/";
            this->voxelGrid->setVoxelValues((GridF(Mesh().fromStl(filename).voxelize(voxelGrid->getDimensions())) - .5f).meanSmooth());
            this->saveErosionDepositionTextureMasks(folder, basename);
            std::cout << "Saved " << basename << " (" << (i+1) << "/" << fileNames.size() << ")" << std::endl;
        }
    }
}

void TerrainGenerationInterface::changeDisplayToComparativeMode(bool toComparative)
{
    this->displayAsComparativeMode = toComparative;
    Q_EMIT updated();
}

void TerrainGenerationInterface::setHeightFactor(float newHeightFactor)
{
    this->heightFactor = newHeightFactor;
    Q_EMIT updated();
}
