#include "BiomeInterface.h"
#include "TerrainModification/RockErosion.h"
//#include "Utils/Voronoi.h"
#include "Utils/AdjencySolver.h"

BiomeInterface::BiomeInterface(QWidget* parent)
    : ActionInterface("biome-generation", parent)
{
    createGUI();
}

void BiomeInterface::display(Vector3 camPos)
{
    if (this->isVisible()) {
        auto queue = std::vector<std::shared_ptr<BiomeInstance>>({this->rootBiome});
        while (!queue.empty()) {
            auto& current = queue.front();
            queue.erase(queue.begin());
            if (!isIn(current->classname, allBiomeNames))
                allBiomeNames.push_back(current->classname);
           for (auto& c : current->instances)
               queue.push_back(c);
        }
        for (size_t i = 0; i < selectionPlanes.size(); i++) {
            auto& selectionPlane = selectionPlanes[i];
            if (selectionPlane.shader != nullptr) {
                int classIndex = std::distance(allBiomeNames.begin(), std::find(allBiomeNames.begin(), allBiomeNames.end(), this->selectedBiomes[i]->classname));
                Vector3 color = HSVtoRGB(float(classIndex) / (float)(allBiomeNames.size() - 1), 1.f, 1.f);
                selectionPlane.shader->setVector("color", std::vector<float>{color.x, color.y, color.z, .5f});
            }
            selectionPlane.display();
        }
    }
}


Vector3 getSurfacePosition(std::shared_ptr<VoxelGrid> grid, Vector3 pos, Vector3 offset = Vector3(), float scaling = 1.f) {
//    pos.x = (pos.x - offset.x) * scaling; // clamp(pos.x / scaling - offset.x, 0.f, grid->getDimensions().x - 1);
//    pos.y = (pos.y - offset.y) * scaling; // clamp(pos.y / scaling - offset.y, 0.f, grid->getDimensions().y - 1);
    pos.z = std::max(pos.z, 0.f); // In case of small imprecision
    while (grid->getVoxelValues().at(pos) > 0) {
        pos += Vector3(0, 0, 1); // Move the position one voxel heigher
        if (!grid->contains(pos)) { // If it gets too high, the whole column must be filled, I guess we should cancel it...
            pos.z = 0;
            break;
        }
    }
    return pos;
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
    if (erosionMatrix.max() < 1e-6) return Matrix3<float>(0, 0, 0);
    erosionMatrix.normalize();
    for (float& m : erosionMatrix) {
        m = interpolation::linear(m, 0.f, 1.0) * strength * (addingMatter ? 1.f : -1.f);
//        m = interpolation::quadratic(interpolation::linear(m, 0.f, 5.f)); //(sigmoid(m) - s_0) / (s_1 - s_0);
    }
    return erosionMatrix;
}
Matrix3<float> BiomeInterface::prepareTrench(std::shared_ptr<BiomeInstance> biome) {
    if (biome->getNumberOfPoints() < 2) return Matrix3<float>(0, 0, 0);
    Vector3 start = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(0)->position)); // Start
    Vector3 end = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(1)->position)); // End
    Vector3 midpoint1 = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(start + (end - start) * .4f + (end - start).cross(Vector3(0, 0, 1)) * .1f));
    Vector3 midpoint2 = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(start + (end - start) * .6f + (end - start).cross(Vector3(0, 0, 1)) * -.1f));
    return archeTunnel(BSpline({start, midpoint1, midpoint2, end}), 8.f * this->voxelGridScaleFactor, 3.f, false, voxelGrid);
}
Matrix3<float> BiomeInterface::prepareCoralWall(std::shared_ptr<BiomeInstance> biome) {
    if (biome->getNumberOfPoints() < 2) return Matrix3<float>(0, 0, 0);
    Vector3 start = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(0)->position)); // Start
    Vector3 end = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(1)->position)); // End
    Vector3 gradient = voxelGrid->getVoxelValues().gradient(getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels((start + end) * .5f)));
    Vector3 midpoint = (start + end) * .5f + gradient * (end - start).norm() * .2f;
    return archeTunnel(BSpline({start, midpoint, end}), 5.f * this->voxelGridScaleFactor, 3.f, true, voxelGrid);
}
Matrix3<float> BiomeInterface::prepareArche(std::shared_ptr<BiomeInstance> biome) {
    if (biome->getNumberOfPoints() < 2) return Matrix3<float>(0, 0, 0);
    Vector3 start = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(0)->position)); // Start
    Vector3 end = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->getPointInstance(1)->position)); // End
    Vector3 midpoint1 = start + (end - start) * .3f + Vector3(0, 0, (end - start).norm() / 2.f); // Midpoint1
    Vector3 midpoint2 = start + (end - start) * .6f + Vector3(0, 0, (end - start).norm() / 2.f); // Midpoint2
    return archeTunnel(BSpline({start, midpoint1, midpoint2, end}), 5.f * this->voxelGridScaleFactor, 10.f, true, voxelGrid);
}
Matrix3<float> BiomeInterface::preparePatateCorail(std::shared_ptr<BiomeInstance> biome) {
    float radius = std::min(5.f, biome->area.containingBoxSize().norm() / 2.f) * this->voxelGridScaleFactor;
    Vector3 patatePosition = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(biome->position)) + Vector3(0, 0, radius/10.f);
    return archeTunnel(BSpline({patatePosition, patatePosition + Vector3(0, 0, 0.001f)}), radius, 2.f, true, voxelGrid);
}

Vector3 BiomeInterface::fromHeightmapPosToVoxels(Vector3 pos)
{
    return (pos - this->voxelGridOffsetStart) * this->voxelGridScaleFactor;
}

Vector3 BiomeInterface::fromVoxelsPosToHeightmap(Vector3 pos)
{
    return (pos / this->voxelGridScaleFactor) + this->voxelGridOffsetStart;
}
void BiomeInterface::generateBiomes(std::shared_ptr<BiomeInstance> predefinedBiomeInstance)
{
    BiomeInstance::instancedBiomes.clear();

    Vector3 heightmapDim = Vector3(3*31, 3*31, 1);
//    this->heightmap->heights = Matrix3<float>(heightmapDim, 20.f);
//    this->heightmap->maxHeight = 40.f;
    this->heightmap = std::make_shared<Heightmap>(heightmapDim.x, heightmapDim.y, 40.f);
    heightmap->raise(Matrix3<float>(heightmapDim, 20.f));

    Vector3 voxelGridOffsetEnd = (voxelGridOffsetStart + heightmapDim / voxelGridScaleFactor).floor();
    voxelGridOffsetEnd.z = 1; // Force the Z component to be 1, instead of being rounded to 0 in the division

    std::cout << "From " << voxelGridOffsetStart << " to " << voxelGridOffsetEnd << " with scale = " << voxelGridScaleFactor << std::endl;
    /////// this->voxelGrid->from2DGrid(*this->heightmap, voxelGridOffsetStart, voxelGridOffsetEnd, voxelGridScaleFactor);
//    this->voxelGrid->fromCachedData();

    auto startTime = std::chrono::system_clock::now();
//return;
    ShapeCurve terrainArea({
                   heightmap->getDimensions() * Vector3(0, 0, 0),
                   heightmap->getDimensions() * Vector3(1, 0, 0),
                   heightmap->getDimensions() * Vector3(1, 1, 0),
                   heightmap->getDimensions() * Vector3(0, 1, 0)
               });
    Vector3 initialSpawn = heightmap->getDimensions() / 2.f + Vector3(10, 0, 0);

    // All the biome hierarchy is created from the json
    // Each biome has a class name, a position and his biome children
    if (predefinedBiomeInstance != nullptr) {
        rootBiome = predefinedBiomeInstance;
    } else {
        rootBiome = biomeModel.createInstance(initialSpawn, terrainArea);
        /*for (int _ = 0; _ < 5; _++) {
            rootBiome = biomeModel.createInstance(initialSpawn, terrainArea);
            for (auto& inst: rootBiome->instances) {
                Qt::GlobalColor color = Qt::gray;
                if (inst->classname == "lagon") {
                    color = Qt::cyan;
                } else if (inst->classname == "recif-frangeant") {
                    color = Qt::darkRed;
                } else if (inst->classname == "plage") {
                    color = Qt::yellow;
                } else if (inst->classname == "recif-barriere") {
                    color = Qt::red;
                } else if (inst->classname == "profondeurs") {
                    color = Qt::darkBlue;
                } else if (inst->classname == "ile") {
                    color = Qt::green;
                }
    //            std::cout << inst->classname << std::endl;
                Plotter::getInstance()->addPlot(inst->area.shrink(0.5f).closedPath(), inst->classname, color);
            }
            Plotter::getInstance()->exec();
            Plotter::getInstance()->reset();
        }
        exit(0);*/
    }
    heightmap->getBiomeIndices() = Matrix3<int>(heightmap->getDimensions(), 0);
    std::vector<std::shared_ptr<BiomeInstance>> biomeQueue;
    std::vector<std::shared_ptr<BiomeInstance>> sortedBiomes;

    biomeQueue.push_back(rootBiome);

    int biomeID = 0;
    while (!biomeQueue.empty()) {
        std::shared_ptr<BiomeInstance> current = biomeQueue.front();
        biomeQueue.erase(biomeQueue.begin());
        if (!current->valid) {
            std::cout << "Biome " << current->getInstanceName() << " is not valid" << std::endl;
            continue;
        }
        if (current->classname == "point") {
//            std::cout << "Biome " << current->getInstanceName() << " is not valid" << std::endl;
            continue;
        }
        sortedBiomes.push_back(current);
        biomeQueue.insert(biomeQueue.end(), current->instances.begin(), current->instances.end());
    }
    biomeQueue.push_back(rootBiome);

    // Sort the biomes by level, so that modifications are done top-to-bottom
    std::sort(sortedBiomes.begin(),
              sortedBiomes.end(),
              [](const auto& a, const auto& b) -> bool {
        return a->getLevel() < b->getLevel();
    });

    possibleBiomeInstances.clear();
    selectedBiomeIDs.clear();

    for (auto& current : sortedBiomes) {
        ShapeCurve area = current->area;
        int level = current->getLevel();
//        area = area.grow(1.f); // Grow the area to fill all artefacts
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        BiomeInstance::registerBiomeInstance(current);
        possibleBiomeInstances.push_back(*current);
        std::ostringstream out;
//        std::ostream& out = std::cout;
        out << "Checking for " << current->getInstanceName() << " :" << std::endl;
        bool atLeastOne = false;
        for (int x = AABBoxMin.x; x < AABBoxMax.x; x++) {
            for (int y = AABBoxMin.y; y < AABBoxMax.y; y++) {
                if (area.contains(Vector3(x, y, area.points[0].z)) && heightmap->getBiomeIndices().checkCoord(Vector3(x, y))) {
                    heightmap->getBiomeIndices().at(x, y).push_back(current->instanceID);
                    out << "Found at (" << x << ", " << y << ")" << std::endl;
                    atLeastOne = true;
                }
            }
        }
        if (!atLeastOne) {
            out << "Never seen..." << std::endl;
        }
        out << "Area was\n";
        for (auto& p : current->area.points)
            out << "- " << p << std::endl;
//        current->instanceID = biomeID;
//        biomeID++;
    }

    return;

    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noise.SetFrequency(.01f);
    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalLacunarity(2.0);
    noise.SetFractalGain(0.7);
    noise.SetFractalWeightedStrength(0.5);
    noise.SetFractalOctaves(1);
    heightmap->getBiomeIndices() = heightmap->getBiomeIndices().wrapWithoutInterpolation(Matrix3<Vector3>::fbmNoise2D(noise, heightmap->getBiomeIndices().sizeX, heightmap->getBiomeIndices().sizeY)  * 50.f);

    // First, level the terrain as the biomes design it
    Matrix3<float> heightChange(heightmap->getSizeX(), heightmap->getSizeY(), 1);
    for (auto& current : sortedBiomes) {
        if (!current->depthShape.points.empty()) {
            ShapeCurve area = current->area;
            Vector3 AABBoxMin, AABBoxMax;
            std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
            Matrix3<float> falloff2D(heightmap->getSizeX(), heightmap->getSizeY(), 1, 0.f);
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
    heightChange = heightChange.meanSmooth(15, 15, 1, true);
    this->heightmap->raise(heightChange);
    voxelGrid->add2DHeightModification(heightChange.subset(voxelGridOffsetStart, voxelGridOffsetEnd).resize(voxelGrid->getDimensions().xy() + Vector3(0, 0, 1)) /* * (voxelGrid->getSizeZ() / (float)heightmap->getMaxHeight())*/, 1.5f);

//    Matrix3<float> all2DModifications(heightmapDim);
//    Matrix3<float> all3DModifications(voxelGrid->getDimensions());
//    bool allModif3D = false;
//    bool allModif2D = false;
    // Now add the primitives on top
    for (auto& current : sortedBiomes) {
        Vector3 pos = current->position;
//        Vector3 surfacePos = getSurfacePosition(voxelGrid, pos);
        ShapeCurve area = current->area;
//        Vector3 areaSize = area.containingBoxSize() + Vector3(0, 0, heightmap->getMaxHeight());
        Vector3 AABBoxMin, AABBoxMax;
        std::tie(AABBoxMin, AABBoxMax) = area.AABBox();
        Matrix3<float> falloff2D(heightmap->getSizeX(), heightmap->getSizeY(), 1, 1.f);
        Matrix3<float> falloff3D(voxelGrid->getSizeX(), voxelGrid->getSizeY(), voxelGrid->getSizeZ(), 1.f);
        for (size_t i = 0; i < falloff2D.size(); i++) {
            Vector3 pos = falloff2D.getCoordAsVector3(i);
            if(AABBoxMin.x <= pos.x && pos.x <= AABBoxMax.x && AABBoxMin.y <= pos.y && pos.y <= AABBoxMax.y) {
                float dist = area.estimateDistanceFrom(pos);
                falloff2D[i] = dist < 0 ? 1.f - interpolation::fault_distance(-dist, 0.1f) : 0.f;
            }
            if (falloff3D.checkCoord(fromHeightmapPosToVoxels(pos))) {
                for (int z = 0; z < falloff3D.sizeZ; z++) {
                    falloff3D.at(fromHeightmapPosToVoxels(pos)) = falloff2D[i];
                }
            }
        }
//        Matrix3<float> subFalloff2D

        Matrix3<float> modifications(0, 0, 0);
        Matrix3<float> heightmapModifier(0, 0, 1, 0.f);

        bool modif3D = false;
        bool modif2D = false;

        if (current->classname == "arche") {
            // Create an arch from the two "point" children
            modifications = prepareArche(current);
            modif3D = true;
        } else if (current->classname == "patate-corail") {
            modifications = preparePatateCorail(current);
            modif3D = true;
        } else if (current->classname == "tranchee" || current->classname == "passe-corail") {
            // Dig a tunnel on the surface
            modifications = prepareTrench(current);
            modif3D = true;
        } else if (current->classname == "mur-corail") {
            // Add a small wall depending on the "point" children
            modifications = prepareCoralWall(current);
            modif3D = true;
        } else if (current->classname == "point") {
            // Nothing to do, this is the smallest primitive
        } else {
//            std::cout << "How the fuck did I get here? Class was " << current->classname << std::endl;
        }

//        allModif2D |= modif2D;
//        allModif3D |= modif3D;

        if (modif2D && !heightmapModifier.data.empty()) {

            voxelGrid->add2DHeightModification(heightmapModifier.subset(AABBoxMin.x, AABBoxMax.x, AABBoxMin.y, AABBoxMax.y) * (voxelGrid->getSizeZ() / (float)heightmap->getHeightFactor()), 10.f, AABBoxMin.xy());
        }
        if (modif3D && !modifications.data.empty()) {

            voxelGrid->applyModification(modifications.subset(fromHeightmapPosToVoxels(AABBoxMin).xy(), fromHeightmapPosToVoxels(AABBoxMax).xy()), fromHeightmapPosToVoxels(AABBoxMin).xy());
        }
    }
//    auto endingTime = std::chrono::system_clock::now();
//    std::cout << "Generation took " << std::chrono::duration_cast<std::chrono::milliseconds>(endingTime - startTime).count()/1000.f << "s" << std::endl;

    // Yeah whatever... Just do the translation once again
    voxelGrid->saveState();
//    this->heightmap->fromVoxelGrid(*voxelGrid);
    updateBiomeSelectionGui();
}

void BiomeInterface::randomize()
{

    // If the main biome has already been computed, regenerate it
    if (rootBiome->instances.size() > 0) {
        std::cout << "A" << std::endl;
        this->modifiedBiomeModel = *rootBiome->toBiomeModel();
        std::cout << "B" << std::endl;
        /// TODO : don't modify the original biomeModel ...
        this->biomeModel = *(std::make_shared<BiomeModel>(modifiedBiomeModel)->clone());
    }
    std::cout << "Generating from model :\n" << this->biomeModel.toJson().dump(4) << std::endl;
    this->generateBiomes();
}

void BiomeInterface::replaceBiome(std::shared_ptr<BiomeInstance> biomeToReplace, std::shared_ptr<BiomeInstance> newBiome)
{
    auto previousParent = biomeToReplace->parent;
    Vector3 previousPos = biomeToReplace->position;
    ShapeCurve previousArea = biomeToReplace->area;

    auto it = std::find(previousParent->instances.begin(), previousParent->instances.end(), biomeToReplace);

    // Swap the content of the pointers
    biomeToReplace = newBiome->clone(previousArea);


    biomeToReplace->parent = previousParent;
    biomeToReplace->position = previousPos;
    biomeToReplace->area = previousArea;
    biomeToReplace->updateSubInstances();
    it->swap(biomeToReplace);
//    BiomeInstance::registerBiomeInstance(biomeToReplace); //->instanceID = previousID;

//    biomeToReplace->completeIfNeeded();
}

/*void BiomeInterface::affectVoxelGrid(std::shared_ptr<VoxelGrid> voxelGrid)
{
    this->voxelGrid = voxelGrid;
}
void BiomeInterface::affectHeightmap(std::shared_ptr<Grid> heightmap)
{
    this->heightmap = heightmap;
}*/

void BiomeInterface::replay(nlohmann::json action)
{

}

void BiomeInterface::mouseMoveEvent(QMouseEvent* event)
{
    ActionInterface::mouseMoveEvent(event);
}


void BiomeInterface::setVoxelGridSizeFactor(float newFactor)
{
    newFactor = std::max(newFactor, 1.f);
    Vector3 currentEnd = this->voxelGridOffsetStart + (this->heightmap->getDimensions() / this->voxelGridScaleFactor);
    Vector3 currentCenter = (currentEnd + this->voxelGridOffsetStart) / 2.f;
    Vector3 maxDims = this->heightmap->getDimensions();

    this->voxelGridScaleFactor = newFactor;

    Vector3 newStart = currentCenter - (maxDims / (this->voxelGridScaleFactor * 2.f));
    Vector3 newEnd = currentCenter + (maxDims / (this->voxelGridScaleFactor * 2.f));

    if (newStart.x < 0) {
        newEnd.x -= newStart.x;
        newStart.x = 0;
    }
    if (newStart.y < 0) {
        newEnd.y -= newStart.y;
        newStart.y = 0;
    }
    if (newStart.z < 0) {
        newEnd.z -= newStart.z;
        newStart.z = 0;
    }

    if (newEnd.x >= maxDims.x) {
        newStart.x -= (newEnd.x - maxDims.x);
        newEnd.x = maxDims.x;
    }
    if (newEnd.y >= maxDims.y) {
        newStart.y -= (newEnd.y - maxDims.y);
        newEnd.y = maxDims.y;
    }
    if (newEnd.z >= maxDims.z) {
        newStart.z -= (newEnd.z - maxDims.z);
        newEnd.z = maxDims.z;
    }
    this->voxelGridOffsetStart = newStart.floor();

    Q_EMIT terrainViewModified(this->voxelGridOffsetStart, this->voxelGridScaleFactor);
}

void BiomeInterface::displayAllBiomes()
{
    this->selectedBiomeIDs.clear();
    this->selectedBiomes.clear();
    std::vector<std::shared_ptr<BiomeInstance> > allBiomes = rootBiome->getAllChildrenBreadthFirst();
    for (auto& biome : allBiomes) {
        if (biome->instanceID >= 0 && !biome->isRoot()) {
            selectedBiomeIDs.push_back(biome->instanceID);
            selectedBiomes.push_back(biome);
        }
    }
    updateBiomeSelectionGui();
}

void BiomeInterface::interchangeBiomes()
{
    // Only available if one instance is selected
    if (this->biomeSelectionGui->selectedItems().size() != 1) return;
    BiomeReplacementDialog dialog(this);
//    for (auto& biome : possibleBiomeInstances) {
    for (auto& [id, biome] : BiomeInstance::instancedBiomes) {
        dialog.allAvailableBiomes->addItem(new HierarchicalListWidgetItem(biome->getInstanceName(), biome->instanceID, biome->getLevel(true)));
    }
    int selectedIndex = dynamic_cast<HierarchicalListWidgetItem*>(this->biomeSelectionGui->item(this->biomeSelectionGui->currentRow()))->ID; //this->biomeSelectionGui->currentRow();
    dialog.show();
    int replacementIndex = dialog.exec();
    replacementIndex = tempIndex;
    if (selectedIndex > -1 && replacementIndex > -1)
        this->replaceBiome(BiomeInstance::instancedBiomes[/*this->selectedBiomeIDs[*/selectedIndex/*]*/], BiomeInstance::instancedBiomes/*possibleBiomeInstances*/[replacementIndex]);
}

void BiomeInterface::keyPressEvent(QKeyEvent* event)
{
    /*
    if (this->isVisible()) {
        bool regenMap = false;
        if (event->key() == Qt::Key_Plus) {
            setVoxelGridSizeFactor(this->voxelGridScaleFactor + 1.f);
            regenMap = true;
        }
        else if (event->key() == Qt::Key_Minus) {
            setVoxelGridSizeFactor(this->voxelGridScaleFactor - 1.f);
            regenMap = true;
        }
        else if (event->key() == Qt::Key_Left) {
            this->voxelGridOffsetStart.x = std::max(this->voxelGridOffsetStart.x - this->voxelGrid->getSizeX() / this->voxelGridScaleFactor, 0.f);
            setVoxelGridSizeFactor(this->voxelGridScaleFactor);
            regenMap = true;
        }
        else if (event->key() == Qt::Key_Right) {
            this->voxelGridOffsetStart.x = std::min(this->voxelGridOffsetStart.x + this->voxelGrid->getSizeX() / this->voxelGridScaleFactor, this->heightmap->heights.sizeX - this->voxelGrid->getSizeX() / this->voxelGridScaleFactor);
            setVoxelGridSizeFactor(this->voxelGridScaleFactor);
            regenMap = true;
        }
        else if (event->key() == Qt::Key_Up) {
            this->voxelGridOffsetStart.y = std::min(this->voxelGridOffsetStart.y + this->voxelGrid->getSizeY() / this->voxelGridScaleFactor, this->heightmap->heights.sizeY - this->voxelGrid->getSizeY() / this->voxelGridScaleFactor);
            setVoxelGridSizeFactor(this->voxelGridScaleFactor);
            regenMap = true;
        }
        else if (event->key() == Qt::Key_Down) {
            this->voxelGridOffsetStart.y = std::max(this->voxelGridOffsetStart.y - this->voxelGrid->getSizeY() / this->voxelGridScaleFactor, 0.f);
            setVoxelGridSizeFactor(this->voxelGridScaleFactor);
            regenMap = true;
        } else if (event->key() == Qt::Key_Delete) {
            this->deleteSelectedBiomes();
        }
        if (regenMap)
            this->generateBiomes(rootBiome);
    }
    Q_EMIT updated();
    */
    ActionInterface::keyPressEvent(event);
}
void BiomeInterface::keyReleaseEvent(QKeyEvent* event)
{
    ActionInterface::keyReleaseEvent(event);
}
void BiomeInterface::wheelEvent(QWheelEvent* event)
{
    ActionInterface::wheelEvent(event);
}

void BiomeInterface::mousePressEvent(QMouseEvent *event)
{
    ActionInterface::mousePressEvent(event);
}

void BiomeInterface::mouseDoubleClickOnMapEvent(Vector3 mousePosition, bool mouseInMap, QMouseEvent *event, TerrainModel* model)
{
    if(this->isVisible() && mouseInMap) {
        // Zoom on the area if left click, zoom out with right click
        // We consider that we got the mouse position just before from the mouseClickedOnMapEvent

        // Move the offset so the new center matches the clicked position
        this->voxelGridOffsetStart = mousePosition - (this->heightmap->getDimensions() / (this->voxelGridScaleFactor * 2.f));
        float newFactor = this->voxelGridScaleFactor;
        if (event->button() == Qt::MouseButton::LeftButton) {
            newFactor += 1.f;
        } else if (event->button() == Qt::MouseButton::RightButton) {
            newFactor -= 1.f;
        }
        setVoxelGridSizeFactor(newFactor);
        this->generateBiomes(rootBiome);
    }
}
/*
void BiomeInterface::mouseDoubleClickEvent(QMouseEvent *event)
{
    if(this->isVisible()) {
        // Zoom on the area if left click, zoom out with right click
        // We consider that we got the mouse position just before from the mouseClickedOnMapEvent

        // Move the offset so the new center matches the clicked position
        this->voxelGridOffsetStart = this->tempMousePos - (this->heightmap->heights.getDimensions() / (this->voxelGridScaleFactor * 2.f));

        if (event->button() == Qt::MouseButton::LeftButton) {
            this->voxelGridScaleFactor += 1.f;
        } else if (event->button() == Qt::MouseButton::RightButton) {
            this->voxelGridScaleFactor -= 1.f;
        }
        setVoxelGridSizeFactor(this->voxelGridScaleFactor);
        this->generateBiomes(rootBiome);
    }
    ActionInterface::mouseDoubleClickEvent(event);
}
*/
void BiomeInterface::mouseClickedOnMapEvent(Vector3 mousePosInMap, bool mouseInMap, QMouseEvent* event, TerrainModel* model)
{
    if (this->isVisible()) {
        if (!mouseInMap) return;

        // Get the voxelGrid mouse position and convert it to the heightmap mouse pos
        this->tempMousePos = fromVoxelsPosToHeightmap(model->getTerrainPos(mousePosInMap));
        this->selectedBiomeIDs = this->heightmap->getBiomeIndices().at(model->getTerrainPos(tempMousePos.xy()));
        std::cout << "Biomes : " << (selectedBiomeIDs.empty() ? "None." : "") << std::endl;
        for (auto ID : selectedBiomeIDs) {
            std::cout << "- " << BiomeInstance::instancedBiomes[ID]->getInstanceName() << std::endl;
        }
        updateBiomeSelectionGui();
    }
}

void BiomeInterface::updateSelectionPlaneToFitBiome(int biomeID, int planeIndex, bool callUpdate)
{
    if (biomeID == -1 || BiomeInstance::instancedBiomes[biomeID]->isRoot()) {
        // Remove all triangles
//        selectionPlane.fromArray(std::vector<float>{});
    } else {
        int level = BiomeInstance::instancedBiomes[biomeID]->getLevel();
        ShapeCurve biomeArea = BiomeInstance::instancedBiomes[biomeID]->area;
        biomeArea = biomeArea.shrink(level);
        std::vector<Vector3> upperPoints, lowerPoints;
        float maxHeight = std::numeric_limits<float>::min(),
                minHeight = std::numeric_limits<float>::max();
        for (auto& point : biomeArea.points) {
            Vector3 surfacePoint = getSurfacePosition(voxelGrid, fromHeightmapPosToVoxels(point));
            upperPoints.push_back(surfacePoint);
            lowerPoints.push_back(surfacePoint);
            maxHeight = std::max(maxHeight, surfacePoint.z);
            minHeight = std::min(maxHeight, surfacePoint.z);
        }
        for (size_t i = 0; i < upperPoints.size(); i++) {
            upperPoints[i].z = maxHeight + 5;
            lowerPoints[i].z = minHeight - 5;
        }
        std::vector<Vector3> vertices;
        for (size_t i = 0; i < upperPoints.size(); i++) {
            size_t next_i = (i + 1) % upperPoints.size();
            vertices.insert(vertices.end(), {
                                upperPoints[i],
                                lowerPoints[i],
                                upperPoints[next_i],

                                lowerPoints[i],
                                lowerPoints[next_i],
                                upperPoints[next_i]
                            });
        }
        selectionPlanes[planeIndex].fromArray(vertices);
        selectionPlanes[planeIndex].cullFace = false;
    }
    if (callUpdate)
        Q_EMIT updated();
}

void BiomeInterface::displayUniqueSelection(int selectionIndex)
{
    if (this->selectedBiomeIDs.size() > selectionIndex) {
        this->selectionPlanes.resize(1);
        this->updateSelectionPlaneToFitBiome(this->selectedBiomeIDs[selectionIndex], 0);
    } else {
        this->selectionPlanes.clear();
    }
}

QLayout* BiomeInterface::createGUI()
{
//    if (this->layout != nullptr) return this->layout;

    layout = new QVBoxLayout();
    biomeSelectionGui = new HierarchicalListWidget;

    QPushButton* seeAllBiomesButton = new QPushButton("Tout voir");
    QPushButton* regenerationButton = new QPushButton("Regenerer");
    QPushButton* interchangeBiomeButton = new QPushButton("Changer le biome...");
    QPushButton* randomizeButton = new QPushButton("Randomiser");

    layout->addWidget(biomeSelectionGui);
    layout->addWidget(seeAllBiomesButton);
    layout->addWidget(interchangeBiomeButton);
    layout->addWidget(regenerationButton);
    layout->addWidget(randomizeButton);
    updateBiomeSelectionGui();

    QObject::connect(seeAllBiomesButton, &QPushButton::pressed, this, &BiomeInterface::displayAllBiomes);
    QObject::connect(regenerationButton, &QPushButton::pressed, this, [&]() -> void { this->generateBiomes(rootBiome); });
    QObject::connect(randomizeButton, &QPushButton::pressed, this, &BiomeInterface::randomize);
    QObject::connect(biomeSelectionGui, &QListWidget::currentRowChanged, this, &BiomeInterface::displayUniqueSelection);
    QObject::connect(interchangeBiomeButton, &QPushButton::pressed, this, &BiomeInterface::interchangeBiomes);
    QObject::connect(biomeSelectionGui, &HierarchicalListWidget::itemDoubleClicked, this, [&](QListWidgetItem* item) -> void {
//        int selectedIndex = this->biomeSelectionGui->currentRow();
        auto selectedBiomeItem = dynamic_cast<HierarchicalListWidgetItem*>(item);
        int selectedBiomeID = selectedBiomeItem->ID;
        auto biome = BiomeInstance::instancedBiomes[selectedBiomeID];
        std::cout << "Selected biome : " << biome->getInstanceName() << std::endl;
    });
    QObject::connect(biomeSelectionGui, &HierarchicalListWidget::itemChangedHierarchy, this, [&] (int ID_to_move, int relatedID, HIERARCHY_TYPE relation, QDropEvent* event) {
        std::shared_ptr<BiomeInstance> toMove = BiomeInstance::instancedBiomes[ID_to_move];
        bool createCopy = event != nullptr && event->keyboardModifiers().testFlag(Qt::KeyboardModifier::ControlModifier);
        // If Ctrl is held, create a copy of the instance
        if (createCopy) {
            toMove = toMove->clone(toMove->area, toMove->position);
        }
        std::shared_ptr<BiomeInstance> related = BiomeInstance::instancedBiomes[relatedID];

        if (relation == HIERARCHY_TYPE::SIBLING) {
            // Critical case : trying to make sibling with the root, just consider it as an error and make it a child
            if (related->parent != nullptr)
                related = related->parent;
            relation = HIERARCHY_TYPE::CHILD;
        } else if (relation == HIERARCHY_TYPE::PARENT) {
            std::shared_ptr<BiomeInstance> tmp = toMove;
            toMove = related;
            related = tmp;
            relation = HIERARCHY_TYPE::CHILD;
        }


        std::shared_ptr<BiomeInstance> previousParent = toMove->parent;
        if (!createCopy) {
            // Remove the instance from the parent's tree (and save the index to possibly plug another biome to the parent)
            auto placeToInsert = previousParent->instances.erase(std::find(previousParent->instances.begin(), previousParent->instances.end(), toMove));
            // If we moved a biome in a lower level of his hierarchy...
            if (toMove->getPathToChild(related).size() > 1) {
                // Get the first child that leads to the target and change his parent
                auto child = toMove->getPathToChild(related)[1];
                child->parent = previousParent;
                previousParent->instances.insert(placeToInsert, child);
            }
        }

//        toMove->parent = related;
//        related->instances.push_back(toMove);
        related->addInstance(toMove);

        previousParent->updateSubInstances();
//        this->rootBiome = this->rootBiome->clone(rootBiome->area, rootBiome->position);
        displayAllBiomes();
    });

    this->replaceDialog = new BiomeReplacementDialog(this);
    return layout;
}

void BiomeInterface::hide()
{
    for (auto& selectionPlane : selectionPlanes)
        selectionPlane.hide();
    CustomInteractiveObject::hide();
}

void BiomeInterface::show()
{
    for (auto& selectionPlane : selectionPlanes)
        selectionPlane.show();
    CustomInteractiveObject::show();
}

void BiomeInterface::addTunnel(KarstHole &hole)
{
    if (hole.path.points.empty()) return;

    Vector3 startPoint = this->fromVoxelsPosToHeightmap(hole.path.points.front());
//    Vector3 endPoint = this->fromVoxelsPosToHeightmap(hole.path.points.back());

    if (this->heightmap->getBiomeIndices().checkCoord(startPoint.xy()) && !this->heightmap->getBiomeIndices().at(startPoint.xy()).empty()) {
        std::shared_ptr<BiomeInstance> startingBiome = BiomeInstance::instancedBiomes[this->heightmap->getBiomeIndices().at(startPoint.xy()).back()];
    //    std::shared_ptr<BiomeInstance> endingBiome = BiomeInstance::instancedBiomes[this->heightmap->getBiomeIndices().at(endPoint.xy()).back()];

        std::shared_ptr<BiomeModel> archeModel = nullptr;
        for (auto& [id, biomeInstance] : BiomeInstance::instancedBiomes) {
            Q_UNUSED(id);
            if (biomeInstance->classname == "arche" && biomeInstance->model != nullptr) {
                archeModel = biomeInstance->model->clone();
                break;
            }
        }
        if (archeModel == nullptr) {
            std::cerr << "Unfortunately, no current Arche has a model, so until I define a static model for it, the regeneration of the model will not include this arche." << std::endl;
        }
        std::shared_ptr<BiomeInstance> newArche = std::make_shared<BiomeInstance>(BiomeInstance::fromClass("arche"));
        newArche->model = archeModel;
        startingBiome->addInstance(newArche);
        newArche->completeIfNeeded();
    }
    Q_EMIT updated();
}

void BiomeInterface::setBindings()
{

}

void BiomeInterface::updateBiomeSelectionGui()
{
//    qDeleteAll(this->biomeSelectionGui->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly));
    biomeSelectionGui->clear();

    std::vector<int> filteredIDs;

    for (size_t i = 0; i < this->selectedBiomeIDs.size(); i++) {
        int biomeID = this->selectedBiomeIDs[i];
        if (BiomeInstance::instancedBiomes.find(biomeID) == BiomeInstance::instancedBiomes.end())
            continue; // Biome not registered... This shouldn't occur.
        auto biome = BiomeInstance::instancedBiomes[biomeID];
//        std::cout << "Adding " << biome->getInstanceName() << " #" << biome->instanceID << " (depth : " << biome->getLevel(true) << ")" << std::endl;
        if (biome != nullptr && !biome->isRoot()) {
            filteredIDs.push_back(biome->instanceID);
            biomeSelectionGui->addItem(new HierarchicalListWidgetItem(biome->getInstanceName(), biome->instanceID, biome->getLevel(true)));
        }
    }
    this->selectedBiomeIDs = filteredIDs;

    this->selectionPlanes.resize(this->selectedBiomeIDs.size());
    for (size_t i = 0; i < this->selectedBiomeIDs.size(); i++)
        this->updateSelectionPlaneToFitBiome(this->selectedBiomeIDs[i], i, false);
    Q_EMIT updated();
}

void BiomeInterface::deleteSelectedBiomes()
{
    for(auto selection : this->biomeSelectionGui->selectedItems()) {
        auto biomeSelection = dynamic_cast<HierarchicalListWidgetItem*>(selection);
        if (biomeSelection != nullptr) {
            this->deleteBiomeFromID(biomeSelection->ID);
        }
    }
    this->updateBiomeSelectionGui();
}

void BiomeInterface::deleteBiomeFromID(int ID)
{
    auto biome = BiomeInstance::instancedBiomes[ID];
    if (biome != nullptr) {

        // Remove from parent
        biome->parent->instances.erase(std::find(biome->parent->instances.begin(), biome->parent->instances.end(), biome));
        biome->parent->updateSubInstances();
        biome->parent = nullptr;

        for (auto& child : biome->instances)
            deleteBiomeFromID(child->instanceID);

    }
    // Maybe even remove it completely ?
    BiomeInstance::instancedBiomes.erase(ID);
}

BiomeReplacementDialog::BiomeReplacementDialog(BiomeInterface* caller)
    : QDialog(), caller(caller)
{
    QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
    allAvailableBiomes = new HierarchicalListWidget(this);
    cancelButton = new QPushButton("Annuler", this);
    validButton = new QPushButton("Confirmer", this);

    vBoxLayout->addWidget(allAvailableBiomes);
    vBoxLayout->addWidget(cancelButton);
    vBoxLayout->addWidget(validButton);

    QObject::connect(cancelButton, &QPushButton::pressed, this, &BiomeReplacementDialog::cancel);
    QObject::connect(validButton, &QPushButton::pressed, this, &BiomeReplacementDialog::confirm);

    setLayout(vBoxLayout);
    setSizeGripEnabled(true);
}

void BiomeReplacementDialog::open()
{
    QDialog::open();
}

void BiomeReplacementDialog::cancel()
{
    selectedBiomeIndex = -1;
    setResult(-1);
    caller->tempIndex = -1;
    this->close();
}

void BiomeReplacementDialog::confirm()
{
    if (allAvailableBiomes->currentRow() >= 0) {
        selectedBiomeIndex = dynamic_cast<HierarchicalListWidgetItem*>(allAvailableBiomes->item(allAvailableBiomes->currentRow()))->ID; //allAvailableBiomes->currentRow();
        setResult(selectedBiomeIndex);
        caller->tempIndex = selectedBiomeIndex;
        this->close();
    } else {
        cancel();
    }
}
