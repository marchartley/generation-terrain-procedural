#include "TerrainGen/LayerBasedGrid.h"

//#include "TerrainModification/UnderwaterErosion.h"

//#include "Graphics/Shader.h"

//#include "Utils/FastNoiseLit.h"
//#include "Utils/Globals.h"

std::map<TerrainTypes, std::pair<float, float>> LayerBasedGrid::materialLimits = {
    {AIR, {-10.f, -8.f} },
    {CURRENT_MIDDLE, {-8.f, -7.f}},
    {CURRENT_TOP, {-7.f, -6.f}},
    {CURRENT_BOTTOM, {-6.f, -5.f}},
    {WATER, {-5.f, 0.f}},
    {CORAL, {0.f, .2f}},
    {SAND, {.2f, .5f}},
    {DIRT, {.5f, 1.f}},
    {ROCK, {1.f, 2.f}},
    {BEDROCK, {2.f, 3.f}},
};

std::vector<TerrainTypes> LayerBasedGrid::invisibleLayers = {
    AIR,
    WATER,
    CURRENT_MIDDLE,
    CURRENT_TOP,
    CURRENT_BOTTOM
};

std::vector<TerrainTypes> LayerBasedGrid::instanciableLayers = {
    BEDROCK,
    CORAL
};

LayerBasedGrid::LayerBasedGrid() : LayerBasedGrid(10, 10, 1.f)
{

}
LayerBasedGrid::LayerBasedGrid(int nx, int ny, float nz)
{
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(nx, ny, 1, {{TerrainTypes::BEDROCK, nz}}); // Initial terrain = 1 layer of bedrock
    this->previousState = layers;

    this->transformationRules = {
        {{{CURRENT_BOTTOM, 1.f}, {BEDROCK, .5f}},       {{CURRENT_BOTTOM, 1.f}, {SAND, .5f}}},
        {{{CURRENT_BOTTOM, 1.f}},                       {{WATER, 1.f}}},
        {{{CURRENT_MIDDLE, 1.f}},                       {{WATER, 1.f}}},
        {{{CURRENT_TOP, 1.f}},                          {{WATER, 1.f}}},
    };
}

TerrainTypes LayerBasedGrid::getValue(Vector3 pos)
{
    return getValue(pos.x, pos.y, pos.z);
}

TerrainTypes LayerBasedGrid::getValue(float x, float y, float z)
{
    if (!this->layers.checkCoord(x, y))
        return TerrainTypes::WATER;
    float currentHeight = 0;
    for (auto materialTuple : this->layers.at(x, y)) {
        if (currentHeight <= z && z < currentHeight + materialTuple.second)
            return materialTuple.first;
        currentHeight += materialTuple.second;
    }
    // No material found at this position, let's assume z < 0 is ground and z > max is water
    return (z < 0 ? TerrainTypes::DIRT : TerrainTypes::WATER);
}

void LayerBasedGrid::addLayer(Vector3 position, float height, TerrainTypes material)
{
    auto layers = this->layers.at(position);
    auto newLayers = layers;
    if (layers.empty()) {
        // Nothing here yet
        newLayers.push_back({material, height});
    } else if (layers.back().first == material){
        // Merge layers
        newLayers.back().second += height;
    } else {
        // Add layer
        newLayers.push_back({material, height});
    }
    this->layers.at(position) = newLayers;
}

/*bool checkOrder(std::vector<std::pair<TerrainTypes, float>>& stack) {
    for (size_t i = 0; i < stack.size() - 1; i++) {
        float densA = LayerBasedGrid::densityFromMaterial(stack[i].first);
        float densB = LayerBasedGrid::densityFromMaterial(stack[i + 1].first);
        if (densA > 0.f && densB > 0.f) {
            if (densA < densB)
                return false;
        }
    }
    return true;
}*/
void LayerBasedGrid::reorderLayers()
{
    size_t nbStacks = layers.size();
    #pragma omp parallel for
    for (size_t iStack = 0; iStack < nbStacks; iStack++) {
        auto& stack = layers[iStack];
        bool reorderingNeeded = true;
        while (reorderingNeeded) {
            reorderingNeeded = false;
            for (size_t i = 0; i < stack.size() - 1; i++) {
//                if (densA > 0.f && densB > 0.f) { // If it's AIR or WATER, don't do anything, I guess
                if ((stack[i].first == TerrainTypes::SAND || stack[i + 1].first == TerrainTypes::SAND) && (isIn(stack[i].first, LayerBasedGrid::invisibleLayers) || isIn(stack[i+1].first, LayerBasedGrid::invisibleLayers))) {
                    float densA = LayerBasedGrid::densityFromMaterial(stack[i].first);
                    float densB = LayerBasedGrid::densityFromMaterial(stack[i + 1].first);
                    if (densA < densB) {
                        auto tmp = stack[i];
                        stack[i] = stack[i + 1];
                        stack[i + 1] = tmp;
                        reorderingNeeded = true;
                    }
                }
                if (stack[i].first == stack[i + 1].first) {
                    stack[i].second += stack[i + 1].second;
                    stack.erase(stack.begin() + i + 1);
                    reorderingNeeded = true;
                }
            }
        }
    }
    this->_historyIndex++;
    this->cleanLayers();
}

float LayerBasedGrid::getSizeZ()
{
    if (this->layers.empty()) return 0;
    float maxHeight = 0.f;
    for (const auto& stack : layers) {
        if (stack.empty()) continue;
        float height = 0;
        for (size_t i = 0; i < stack.size() - 1; i++) {
            height += stack[i].second;
        }
        if (!isIn(stack.back().first, LayerBasedGrid::invisibleLayers)) { // != TerrainTypes::AIR && stack.back().first != TerrainTypes::WATER) {
            height += stack.back().second;
        }
        maxHeight = std::max(maxHeight, height);
    }
    return maxHeight;
//    for (auto materialPair : this->layers[0])
//        height += materialPair.second;
    //    return height;
}

float LayerBasedGrid::getHeight(float x, float y)
{
    float currentHeight = 0.f;
    bool startCounting = false;
    auto& materials = this->layers.at(x, y);
    for (int i = materials.size() - 1; i >= 0; i--) {
        auto [mat, height] = materials[i];
        if (!startCounting)
            startCounting = !isIn(mat, LayerBasedGrid::invisibleLayers); //(mat != AIR && mat != WATER);
        if (startCounting) {
            currentHeight += height;
        }
    }
    return currentHeight;
}

//float LayerBasedGrid::getHeight(Vector3 pos)
//{
//    return this->getHeight(pos.x, pos.y);
//}

void LayerBasedGrid::from2DGrid(Heightmap grid)
{
    // The materials for now will only be dirt and water
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(grid.getSizeX(), grid.getSizeY(), 1);
    for (int x = 0; x < grid.getSizeX(); x++) {
        for (int y = 0; y < grid.getSizeY(); y++) {
            this->layers.at(x, y) = {
                { TerrainTypes::SAND, grid.getHeight(x, y) },
                { TerrainTypes::WATER, grid.getMaxHeight() - grid.getHeight(x, y) }
            };
        }
    }
}

void LayerBasedGrid::fromVoxelGrid(VoxelGrid& voxelGrid)
{
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(voxelGrid.getSizeX(), voxelGrid.getSizeY(), 1);
    Matrix3<float> voxelValues = voxelGrid.getVoxelValues();
    for (int x = 0; x < voxelGrid.getSizeX(); x++) {
        for (int y = 0; y < voxelGrid.getSizeY(); y++) {
            for (int z = 0; z < voxelGrid.getSizeZ(); z++) {
                float voxelValue = voxelValues.at(x, y, z);
                TerrainTypes material = LayerBasedGrid::materialFromDensity(voxelValue);
                if (this->layers.at(x, y).empty()) {
                    this->layers.at(x, y).push_back({material, 1.f});
                } else {
                    if (this->layers.at(x, y).back().first == material)
                        this->layers.at(x, y).back().second += 1.f;
                    else
                        this->layers.at(x, y).push_back({material, 1.f});
                }
            }
        }
    }
}

VoxelGrid LayerBasedGrid::toVoxelGrid()
{
    VoxelGrid vox;
    vox.fromLayerBased(*this);
    return vox;
}

Matrix3<float> LayerBasedGrid::voxelize(int fixedHeight, float kernelSize)
{
    float maxHeight = (fixedHeight == -1 ? this->getSizeZ() : fixedHeight);
    Matrix3<float> values = Matrix3<float>(this->getSizeX(), this->getSizeY(), maxHeight, LayerBasedGrid::densityFromMaterial(TerrainTypes::WATER));

    int sizeX = values.sizeX;
    int sizeY = values.sizeY;
    int sizeZ = maxHeight;
//    auto start = std::chrono::system_clock::now();
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            for (int z = 0; z < sizeZ; z++) {
                auto materialAndVolume = this->getKernel(Vector3(x +.5f, y + .5f, z +.5f), kernelSize);
                float solidMaterial = 0.f;
                float airMaterial = 0.f;
                float totalDensity = 0.f;
                for (auto& [mat, vol] : materialAndVolume) {
                    if (isIn(mat, LayerBasedGrid::invisibleLayers)) { //mat == AIR || mat == WATER) {
                        airMaterial += vol;
                    } else if (isIn(mat, LayerBasedGrid::instanciableLayers)) { // For rocks and corals
                        airMaterial += vol;
                    } else {
                        solidMaterial += vol;
                    }
                    totalDensity += vol * LayerBasedGrid::densityFromMaterial(mat);
                }
//                float val = (2.f * solidMaterial / (solidMaterial + airMaterial)) - 1.f;
//                values.at(x, y, z) = (val > 0.f ? LayerBasedGrid::densityFromMaterial(SAND) : LayerBasedGrid::densityFromMaterial(AIR));
                values.at(x, y, z) = (totalDensity / (airMaterial + solidMaterial)); // Total density / total volume
            }
            /*
            auto& stack = this->layers.at(x, y);
            float currentHeight = 0.f;
            for (size_t i = 0; i < stack.size(); i++) {
                float density = LayerBasedGrid::densityFromMaterial(stack[i].first);
                float height = stack[i].second;

                for (int z = std::floor(currentHeight); z < std::min(std::floor(currentHeight + height), maxHeight); z++) {
                    values.at(x, y, z) = density;
                }
                currentHeight += height;
            }*/
        }
    }
//    auto end = std::chrono::system_clock::now();
//    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(end -start).count() << "ms" << std::endl;
    return values;
}

std::map<TerrainTypes, float> LayerBasedGrid::getKernel(Vector3 pos, float kernelSize)
{
    std::map<TerrainTypes, float> finalDensities;
    float halfK = kernelSize * .5f;
    float minX = std::max(pos.x - halfK, 0.f);
    float maxX = std::min(pos.x + halfK, float(this->getSizeX()));
    float minY = std::max(pos.y - halfK, 0.f);
    float maxY = std::min(pos.y + halfK, float(this->getSizeY()));
    float minZ = std::max(pos.z - halfK, 0.f);
    float maxZ = pos.z + halfK; // No way to know the Z max in a fast way

    for (int x = std::floor(minX); x < std::ceil(maxX); x++) {
        for (int y = std::floor(minY); y < std::ceil(maxY); y++) {
            float widthX = 1.f; //std::min(x+1.f, maxX) - std::max(x+0.f, minX);
            float widthY = 1.f; // std::min(y+1.f, maxY) - std::max(y+0.f, minY);
            float currentHeight = 0.f;
            int currentStack = 0;
            while (currentHeight < maxZ && currentStack < this->layers.at(x, y).size()) {
                auto [blockMaterial, blockHeight] = this->layers.at(x, y)[currentStack];
                float start = clamp(currentHeight, minZ, maxZ);
                float end   = clamp(currentHeight + blockHeight, minZ, maxZ);
                float volume = (end - start) * (widthX * widthY);
                currentHeight += blockHeight;
                if (finalDensities.count(blockMaterial) == 0)
                    finalDensities[blockMaterial] = 0.f;
                finalDensities[blockMaterial] += volume;
                currentStack++;
            }
            if (finalDensities.count(AIR) == 0)
                finalDensities[AIR] = 0.f;
            finalDensities[AIR] += std::max(0.f, maxZ - currentHeight) * (widthX * widthY);
        }
    }
    return finalDensities;
}

std::pair<TerrainTypes, float> LayerBasedGrid::getMaterialAndHeight(Vector3 pos)
{
    if (!this->layers.checkCoord(pos.xy()))
        return {TerrainTypes::AIR, 0.f};
    float currentHeight = 0;
    for (auto materialTuple : this->layers.at(pos.xy())) {
        if (currentHeight <= pos.z && pos.z < currentHeight + materialTuple.second)
            return materialTuple;
        currentHeight += materialTuple.second;
    }
    // No material found at this position, let's assume z < 0 is ground and z > max is water
//    return (pos.z < 0 ? {TerrainTypes::DIRT, 0.f} : {TerrainTypes::WATER, 0.f});
    return {TerrainTypes::WATER, 0.f};
}

Vector3 LayerBasedGrid::getFirstIntersectingStack(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    if (!minPos.isValid()) minPos = Vector3();
    if (!maxPos.isValid()) maxPos = this->getDimensions();

    Vector3 currPos = origin;
//    auto values = this->getVoxelValues();
//    values.raiseErrorOnBadCoord = false;
    float distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
    float distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    // Continue while we are in the grid or we are heading towards the grid
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {
        float isoval = LayerBasedGrid::densityFromMaterial(this->getMaterialAndHeight(currPos).first);
        if (isoval > 0.0) {
            return currPos;
        }
        currPos += dir;
        distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
        distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    }
    return Vector3(false);
}

Vector3 LayerBasedGrid::getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    return this->getFirstIntersectingStack(origin, dir, minPos, maxPos);
}

std::pair<Matrix3<int>, Matrix3<float>> LayerBasedGrid::getMaterialAndHeightsGrid()
{
    if (currentHistoryIndex != _historyIndex)
    {
        int maxStackSize = 0;
        for (auto& cell : layers)
            maxStackSize = std::max(maxStackSize, int(cell.size()));

        Matrix3<int> materials(this->getSizeX(), this->getSizeY(), maxStackSize);
        Matrix3<float> heights(materials.getDimensions());

        for (size_t i = 0; i < layers.size(); i++) {
            for (size_t iStack = 0; iStack < layers[i].size(); iStack++) {
                Vector3 index = layers.getCoordAsVector3(i) + Vector3(0, 0, iStack);
                materials.at(index) = (int) layers[i][iStack].first;
                heights.at(index) = (float) layers[i][iStack].second;
            }
        }
        this->_cachedMaterialAndHeights = {materials, heights};
        currentHistoryIndex = _historyIndex;
    }
    return this->_cachedMaterialAndHeights;
}

void LayerBasedGrid::thermalErosion()
{
    this->cleanLayers();
    for (int x = 0; x < layers.sizeX; x++) {
        for (int y = 0; y < layers.sizeY; y++) {
            auto& stack = layers.at(x, y);
            float currentHeight = 0.f;
            for (size_t i = 0; i < stack.size(); i++) {
                TerrainTypes type = stack[i].first;
                float height = stack[i].second;
                float higherPos = currentHeight + height;
                if (type == TerrainTypes::SAND) {
                    float maxMatterPerNeighbor = height / 8.f;
                    // If WATER or AIR is around, transfer some SAND
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            float sharableMatter = maxMatterPerNeighbor;
                            Vector3 nPos = Vector3(x + dx, y + dy);
                            if (!layers.checkCoord(nPos)) continue;

                            float currentNHeight = 0.f;
                            for (size_t neighbor = 0; neighbor < layers.at(nPos).size() && currentNHeight < higherPos; neighbor++) {
                                TerrainTypes nType = layers.at(nPos)[neighbor].first;
                                float nHeight = layers.at(nPos)[neighbor].second;
                                float nDensity = LayerBasedGrid::densityFromMaterial(nType);
                                float higherNPos = currentNHeight + nHeight;

                                if (nDensity <= 0.f) {
                                    // If it's WATER or AIR, we can share layers
                                    float minHeight = std::clamp(currentNHeight, currentHeight, higherPos);
                                    float maxHeight = std::clamp(higherNPos, currentHeight, higherPos);

                                    float sharedMatter = std::min(sharableMatter, maxHeight - minHeight);
                                    if (sharedMatter > 0.01f) {
                                        sharableMatter -= sharedMatter;

                                        layers.at(nPos)[neighbor].second -= sharedMatter;
                                        layers.at(nPos).insert(layers.at(nPos).begin() + neighbor, {type, sharedMatter});

                                        higherPos -= sharedMatter;
                                        stack[i].second -= sharedMatter;
                                        stack.insert(stack.begin() + i + 1, {nType, sharedMatter});
                                    }
                                }

                                currentNHeight = higherNPos;
                            }

                            if (currentNHeight < higherPos) {
                                float sharedMatter = std::min(sharableMatter, higherPos - currentNHeight);
                                if (sharedMatter > 0.01f) {
                                    sharableMatter -= sharedMatter;

                                    layers.at(nPos).insert(layers.at(nPos).end(), {type, sharedMatter});

//                                    higherPos -= sharedMatter;
                                    stack[i].second -= sharedMatter;
                                }
                            }
                        }
                    }
                }
                currentHeight = higherPos;
            }
        }
    }
    this->reorderLayers();
}

void LayerBasedGrid::cleanLayers(float minLayerHeight)
{
    size_t nbStacks = layers.size();
    for (size_t iStack = 0; iStack < nbStacks; iStack++) {
        auto& stack = layers[iStack];
        while (stack.size() > 0 && isIn(stack.back().first, LayerBasedGrid::invisibleLayers)) //(stack.back().first == WATER || stack.back().first == AIR))
            stack.pop_back();
        for (size_t i = 0; i < stack.size(); i++) {
            if (stack[i].second < minLayerHeight) {
                stack.erase(stack.begin() + i);
            }
        }
    }
}
/*
void LayerBasedGrid::add(Patch2D patch, TerrainTypes material, bool applyDistanceFalloff, float distancePower)
{
//    auto [minPos, maxPos] = patch.shape.AABBox();
    Vector3 minPos = patch.position;
    Vector3 maxPos = patch.dimension;
    minPos += patch.position;
    maxPos += patch.position;
    int minX = std::max(0, int(minPos.x));
    int maxX = std::min(this->getSizeX(), int(maxPos.x));
    int minY = std::max(0, int(minPos.y));
    int maxY = std::min(this->getSizeY(), int(maxPos.y));
    Matrix3<float> distanceMap = Matrix3<float>(maxX - minX, maxY - minY, 1, 1.f);
    if (applyDistanceFalloff) {
        distanceMap = distanceMap.toDistanceMap(true);
        distanceMap.normalize();
        for (auto& val : distanceMap)
            val = std::pow(val, distancePower);
    }
    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            float height = patch.getMaxHeight(Vector3(x, y) - patch.position);
            float strength = distanceMap.at(x - minX, y - minY);
            this->addLayer(Vector3(x, y), height * strength, material);
        }
    }
//    this->reorderLayers();
    _historyIndex ++;
    return;
}

void LayerBasedGrid::add(Patch3D patch, TerrainTypes material, bool applyDistanceFalloff, float distancePower)
{
    Vector3 minPos = patch.position;
    Vector3 maxPos = patch.dimension;
    minPos += patch.position;
    maxPos += patch.position;
    int minX = std::max(0, int(minPos.x));
    int maxX = std::min(this->getSizeX(), int(maxPos.x));
    int minY = std::max(0, int(minPos.y));
    int maxY = std::min(this->getSizeY(), int(maxPos.y));
    float minZ = std::max(0.f, minPos.z);
    float maxZ = float(maxPos.z); // std::min(this->getSizeZ(), float(maxPos.z));
    Matrix3<float> distanceMap = Matrix3<float>(maxX - minX, maxY - minY, 1, 1.f);
    if (applyDistanceFalloff) {
        distanceMap = distanceMap.toDistanceMap(true);
        distanceMap.normalize();
        for (auto& val : distanceMap)
            val = std::pow(val, distancePower);
    }
    float zResolution = .10f;
    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            float z = minZ;
            while ( z < maxZ) {
                float eval = patch.evaluate(Vector3(x, y, z) - patch.position);
                if (eval != 0.f) {
                    float layerDens = patch.get(Vector3(x, y, z) - patch.position);
                    if (layerDens > 0) {
                        int a = 0;
                    }
                    this->addLayer(Vector3(x, y), zResolution, LayerBasedGrid::materialFromDensity(layerDens));
                } else {
                    this->addLayer(Vector3(x, y), zResolution, TerrainTypes::AIR);
                }
                z += zResolution;
            }
        }
    }
//    this->reorderLayers();
    _historyIndex ++;
    return;
}*/

void LayerBasedGrid::add(ImplicitPatch* patch)
{
    auto AABBox = patch->getSupportBBox();
    Vector3 minPos = AABBox.first;
    Vector3 maxPos = AABBox.second;

    int minX = std::max(0, int(minPos.x));
    int maxX = std::min(int(this->getSizeX()), int(maxPos.x));
    int minY = std::max(0, int(minPos.y));
    int maxY = std::min(int(this->getSizeY()), int(maxPos.y));
    float minZ = std::max(0.f, minPos.z);
    float maxZ = float(maxPos.z);
    std::cout << "Evaluation from " << Vector3(minX, minY, minZ) << " to " << Vector3(maxX, maxY, maxZ) << std::endl;
    float zResolution = 1.f;
    int nbEvaluations = 0;
    auto start = std::chrono::system_clock::now();
    #pragma omp parallel for collapse(2)
    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            float z = minZ;
            while ( z < maxZ) {
                nbEvaluations++;
                auto [totalValue, initialDensityAndProportion] = patch->getMaterialsAndTotalEvaluation(Vector3(x, y, z)); //patch->getDensityAtPosition(Vector3(x, y, z));
                auto densityAndProportion = initialDensityAndProportion;
                if (totalValue < ImplicitPatch::isovalue) { // Not enough material, just add air
                    this->addLayer(Vector3(x, y), zResolution, TerrainTypes::AIR);
                } else { // Otherwise, get the material to display and add it to the terrain
                    // Just testing rules
                    for (size_t i = 0; i < transformationRules.size(); i++) {
                        auto [input, output] = transformationRules[i];
//                        bool transformPossible = true;
                        float maxTransform = 10000.f;
                        for (auto [inMaterial, inDose] : input) {
                            maxTransform = std::min(maxTransform, densityAndProportion[inMaterial] / inDose);
                        }
                        if (maxTransform > 1e-3) { //(transformPossible) {
                            for (auto [inMaterial, inDose] : input) {
                                densityAndProportion[inMaterial] -= inDose * maxTransform;
                            }
                            for (auto [outMaterial, outDose] : output) {
                                densityAndProportion[outMaterial] += outDose * maxTransform;
                            }
                        }
                    }
                    // End testing rules

                    // If changes has been made with rules, bring back the ratio
                    float newTotalValue = 0.f;
                    for (const auto& [mat, val] : densityAndProportion) // Get new material total
                        newTotalValue += val;
                    for (auto& [mat, val] : densityAndProportion) // Adjust to get back the initial material quantity
                        val *= totalValue / newTotalValue;

                    // Find the most used material
                    TerrainTypes selectedMaterial;
                    float maxValue = 0.f;
                    for (const auto& [material, value] : densityAndProportion) {
                        if (value > maxValue) {
                            selectedMaterial = material;
                            maxValue = value;
                        }
                    }
//                    if (totalValue >= .5f)
                    this->addLayer(Vector3(x, y), zResolution, selectedMaterial); //LayerBasedGrid::materialFromDensity(selectedDensity));
                }
//                else
//                    this->addLayer(Vector3(x, y), zResolution, TerrainTypes::AIR);
                /*
                float eval = patch->evaluate(Vector3(x, y, z));
                if (eval >= 0.5f) {
//                    float layerDens = patch->get(Vector3(x, y, z) - patch->position);
                    this->addLayer(Vector3(x, y), zResolution, material); // LayerBasedGrid::materialFromDensity(layerDens));
                } else {
                    this->addLayer(Vector3(x, y), zResolution, TerrainTypes::AIR);
                }*/
                z += zResolution;
            }
        }
    }
    auto end = std::chrono::system_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Evaluation time: " << totalTime << "ms for " << nbEvaluations << " points (mean = " << totalTime/(float)nbEvaluations << "ms/eval)" << std::endl;
    this->cleanLayers();
//    this->reorderLayers();
    _historyIndex ++;
    return;
}

Mesh LayerBasedGrid::getGeometry()
{
    std::vector<Vector3> vertices;

//    auto layersAndHeights = this->getMaterialAndHeightsGrid();
//    auto layers = layersAndHeights.first;
//    auto heights = layersAndHeights.second;

    for (int x = 0; x < this->getSizeX(); x++) {
        for (int y = 0; y < this->getSizeY(); y++) {
            auto layers = this->layers.at(x, y);
            float currentHeight = 0.f;

            for (size_t i = 0; i < layers.size(); i++) {
                float height = layers[i].second; //heights.at(x, y, i);
                TerrainTypes material = layers[i].first; //(TerrainTypes) layers.at(x, y, i);
                if (!isIn(material, LayerBasedGrid::invisibleLayers)) { //material != TerrainTypes::AIR && material != TerrainTypes::WATER) {
                    auto cube = CubeMesh::cubesVertices;

                    for (auto& vert : cube) {
                        vert = vert * Vector3(1.f, 1.f, height) + Vector3(x, y, currentHeight);
                    }
                    vertices.insert(vertices.end(), cube.begin(), cube.end());
                }
                currentHeight += height;

            }
        }
    }

    Mesh m;
    m.useIndices = false;
    m.fromArray(vertices);
    return m;
}

TerrainTypes LayerBasedGrid::materialFromDensity(float density)
{
    for (auto materialPair : LayerBasedGrid::materialLimits) {
        float valMin, valMax;
        std::tie(valMin, valMax) = materialPair.second;
        if (valMin <= density && density < valMax)
            return materialPair.first;
    }
    // No material found at this density, let's assume if dens < 0 it's water and z > 0 it's ground
    return (density > 0 ? (TerrainTypes)(TerrainTypes::LAST - 1): TerrainTypes::WATER);
}

float LayerBasedGrid::densityFromMaterial(TerrainTypes material)
{
    float valMin, valMax;
    std::tie(valMin, valMax) = LayerBasedGrid::materialLimits[material];
    return (valMin + valMax) / 2.f;
}

float LayerBasedGrid::minDensityFromMaterial(TerrainTypes material)
{
    float valMin, valMax;
    std::tie(valMin, valMax) = LayerBasedGrid::materialLimits[material];
    return valMin;
}

float LayerBasedGrid::maxDensityFromMaterial(TerrainTypes material)
{
    float valMin, valMax;
    std::tie(valMin, valMax) = LayerBasedGrid::materialLimits[material];
    return valMax;
}

