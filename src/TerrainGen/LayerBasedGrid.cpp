#include "TerrainGen/LayerBasedGrid.h"

//#include "TerrainModification/UnderwaterErosion.h"

//#include "Graphics/Shader.h"

//#include "Utils/FastNoiseLit.h"
//#include "Utils/Globals.h"
#include "Graphics/CubeMesh.h"

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

std::set<TerrainTypes> LayerBasedGrid::invisibleLayers = {
    AIR,
    WATER,
    CURRENT_MIDDLE,
    CURRENT_TOP,
    CURRENT_BOTTOM
};

std::set<TerrainTypes> LayerBasedGrid::instanciableLayers = {
    BEDROCK,
    CORAL
};

LayerBasedGrid::LayerBasedGrid() : LayerBasedGrid(10, 10, 1.f)
{

}
LayerBasedGrid::LayerBasedGrid(int nx, int ny, float nz)
{
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(nx, ny, 1, {{TerrainTypes::SAND, nz}}); // Initial terrain = 1 layer of bedrock
    this->previousState = layers;
}

TerrainTypes LayerBasedGrid::getValue(const Vector3& pos)
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

void LayerBasedGrid::addLayer(const Vector3& position, float height, TerrainTypes material)
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
            for (int i = 0; i < int(stack.size() - 1); i++) {
//                if (densA > 0.f && densB > 0.f) { // If it's AIR or WATER, don't do anything, I guess
                if ((isIn(stack[i].first, std::vector<TerrainTypes>{TerrainTypes::SAND, TerrainTypes::DIRT}) ||
                     isIn(stack[i+1].first, std::vector<TerrainTypes>{TerrainTypes::SAND, TerrainTypes::DIRT})) &&
                        (isIn(stack[i].first, LayerBasedGrid::invisibleLayers) || isIn(stack[i+1].first, LayerBasedGrid::invisibleLayers))) {
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

float LayerBasedGrid::getSizeZ() const
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

//float LayerBasedGrid::getHeight(const Vector3& pos)
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
    this->cleanLayers();
    this->_historyIndex = -1;
}

void LayerBasedGrid::fromVoxelGrid(VoxelGrid& voxelGrid)
{
    int sizeX = voxelGrid.getSizeX();
    int sizeY = voxelGrid.getSizeY();
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(sizeX, sizeY, 1);
    GridF voxelValues = voxelGrid.getVoxelValues();
#pragma omp parallel for collapse(2)
    for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
            float prevVoxelValue = -1;  // Assume a voxel value of 0 at z = -1 for simplicity
            for (int z = 0; z < voxelGrid.getSizeZ(); z++) {
                float currentVoxelValue = voxelValues(x, y, z);
                TerrainTypes material = currentVoxelValue < 0 ? TerrainTypes::AIR : TerrainTypes::BEDROCK; //LayerBasedGrid::materialFromDensity(voxelValue);
                float height = 1.f;
                if ((prevVoxelValue < 0 && currentVoxelValue >= 0) || (prevVoxelValue >= 0 && currentVoxelValue < 0)) {
                    height = 1.0f - std::abs(prevVoxelValue) / (std::abs(prevVoxelValue) + std::abs(currentVoxelValue));
                }
                layers.at(x, y).push_back({material, height});
                prevVoxelValue = currentVoxelValue;
                /*
                if (this->layers.at(x, y).empty()) {
                    this->layers.at(x, y).push_back({material, 1.f});
                } else {
                    if (this->layers.at(x, y).back().first == material)
                        this->layers.at(x, y).back().second += 1.f;
                    else
                        this->layers.at(x, y).push_back({material, 1.f});
                }*/
            }
        }
    }
    cleanLayers();

    /*return; // TODO: Remove this return very soon
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(voxelGrid.getSizeX(), voxelGrid.getSizeY(), 1);
    GridF voxelValues = voxelGrid.getVoxelValues();
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
    }*/
}

void LayerBasedGrid::fromImplicit(ImplicitPatch *implicitTerrain)
{
    this->reset();
    this->add(implicitTerrain);
}

VoxelGrid LayerBasedGrid::toVoxelGrid()
{
    VoxelGrid vox;
    vox.fromLayerBased(*this);
    return vox;
}

GridF LayerBasedGrid::voxelize(int fixedHeight, float kernelSize)
{
    float maxHeight = (fixedHeight == -1 ? this->getSizeZ() : fixedHeight);
    GridF values = GridF(this->getSizeX(), this->getSizeY(), maxHeight, LayerBasedGrid::densityFromMaterial(TerrainTypes::WATER));

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

GridF LayerBasedGrid::getVoxelized(const Vector3& dimensions, const Vector3& scale)
{
    return voxelize();
}

std::map<TerrainTypes, float> LayerBasedGrid::getKernel(const Vector3& pos, float kernelSize)
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
            size_t currentStack = 0;
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

std::pair<TerrainTypes, float> LayerBasedGrid::getMaterialAndHeight(const Vector3& pos)
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

Vector3 LayerBasedGrid::getFirstIntersectingStack(const Vector3& origin, const Vector3& dir, const Vector3& minPos, const Vector3& maxPos)
{
    AABBox myAABBox = AABBox(Vector3::max(minPos, Vector3()), Vector3::max(maxPos, this->getDimensions()));
//    if (!minPos.isValid()) minPos = Vector3();
//    if (!maxPos.isValid()) maxPos = this->getDimensions();

    Vector3 currPos = origin;
//    auto values = this->getVoxelValues();
//    values.raiseErrorOnBadCoord = false;
    float distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, myAABBox.min(), myAABBox.max());
    float distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, myAABBox.min(), myAABBox.max());
    // Continue while we are in the grid or we are heading towards the grid
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {
        if (Vector3::isInBox(currPos, myAABBox.min(), myAABBox.max())) {
            float isoval = LayerBasedGrid::densityFromMaterial(this->getMaterialAndHeight(currPos).first);
            if (isoval > 0.0) {
                return currPos;
            }
        }
        currPos += dir;
        distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, myAABBox.min(), myAABBox.max());
        distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, myAABBox.min(), myAABBox.max());
    }
    return Vector3(false);
}

Vector3 LayerBasedGrid::getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos, const Vector3& maxPos)
{
    return this->getFirstIntersectingStack(origin, dir, minPos, maxPos);
}

std::pair<GridI, GridF> LayerBasedGrid::getMaterialAndHeightsGrid()
{
//    int x = 0;
    if (true || currentHistoryIndex != _historyIndex)
    {
        int maxStackSize = 0;
        for (auto& cell : layers)
            maxStackSize = std::max(maxStackSize, int(cell.size()));

        GridI materials(this->getSizeX(), this->getSizeY(), maxStackSize);
        GridF heights(materials.getDimensions());

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

GridV3 LayerBasedGrid::getNormals()
{
    GridV3 normals = -this->voxelize().gradient();
    for (auto& n : normals)
        n.normalize();
    normals.raiseErrorOnBadCoord = false;
    normals.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    return normals;
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

void LayerBasedGrid::cleanLayer(int x, int y, float minLayerHeight)
{
    auto& stack = layers.at(x, y);
    while (stack.size() > 0 && isIn(stack.back().first, LayerBasedGrid::invisibleLayers)) //(stack.back().first == WATER || stack.back().first == AIR))
        stack.pop_back();
    for (size_t i = 0; i < stack.size(); i++) {
        if (stack[i].second <= minLayerHeight) {
            stack.erase(stack.begin() + i);
        }
    }
    for (int i = stack.size() - 2; i >= 0; i--) {
        if (stack[i].first == stack[i+1].first) {
            stack[i].second += stack[i+1].second;
            stack.erase(stack.begin() + i + 1);
//            stack[i].second = 0;
//            stack[i+1].second = 0;
        }
    }
}

void LayerBasedGrid::cleanLayers(float minLayerHeight)
{
    size_t nbStacks = layers.size();
    for (size_t iStack = 0; iStack < nbStacks; iStack++) {
        auto& stack = layers[iStack];
        Vector3 pos = layers.getCoordAsVector3(iStack);
        cleanLayer(pos.x, pos.y, minLayerHeight);
        /*while (stack.size() > 0 && isIn(stack.back().first, LayerBasedGrid::invisibleLayers)) //(stack.back().first == WATER || stack.back().first == AIR))
            stack.pop_back();
        for (size_t i = 0; i < stack.size(); i++) {
            if (stack[i].second < minLayerHeight) {
                stack.erase(stack.begin() + i);
            }
        }
        for (int i = stack.size() - 2; i >= 0; i--) {
            if (stack[i].first == stack[i+1].first) {
//                stack[i].second += stack[i+1].second;
                stack.erase(stack.begin() + i + 1);
            }
        }*/
    }
}

LayerBasedGrid* LayerBasedGrid::transformLayer(int x, int y, float startZ, float endZ, TerrainTypes material)
{
//    material = AIR;
    auto initialColumn = this->layers.at(x, y);
    auto zLevels = std::vector<std::pair<TerrainTypes, std::pair<float, float>>>();
    float currentHeight = 0.f;
    for (size_t i = 0; i < initialColumn.size(); i++) {
        auto [mat, height] = initialColumn[i];

        float startHeight = currentHeight;
        float endHeight = currentHeight + height;

        if (startZ > endHeight || endZ < startHeight) {
            zLevels.push_back({mat, {startHeight, endHeight}}); // Nothing to do on this layer
        } else {
            // Under new layer
            float startUnder = startHeight;
            float endUnder = std::min(endHeight, startZ); // Min is not really useful since we checked that endHeight >= startZ.
            if (startUnder != endUnder)
                zLevels.push_back({mat, {startUnder, endUnder}});

            // Inside new layer
            float startIn = endUnder;
            float endIn = std::min(endHeight, endZ);
            if (startIn != endIn)
                zLevels.push_back({material, {startIn, endIn}});

            // Above new layer
            float startAbove = endIn;
            float endAbove = endHeight;
            if (startAbove != endAbove)
                zLevels.push_back({mat, {startAbove, endAbove}});
        }
        currentHeight += height;
        startZ = currentHeight;
    }
    if (currentHeight < endZ) {
        if (currentHeight < startZ) {
            zLevels.push_back({AIR, {currentHeight, startZ}});
            currentHeight = startZ;
        }
        zLevels.push_back({material, {currentHeight, endZ}});
    }
    this->layers.at(x, y).resize(zLevels.size());
    for (size_t i = 0; i < zLevels.size(); i++) {
        this->layers.at(x, y)[i] = {zLevels[i].first, (zLevels[i].second.second - zLevels[i].second.first)};
    }
    this->cleanLayer(x, y);
    return this;


    /*auto initialColumn = this->layers.at(x, y);
    auto& column = this->layers.at(x, y);
    float currentHeight = 0.f;
    int stackIndex = 0;
    for (auto& [mat, height] : initialColumn) {
        if (currentHeight > endZ) break;
        if (currentHeight + height >= startZ && mat != material) {
            if (currentHeight < startZ && currentHeight + height > endZ) {
                // Cut the layer in 3 : before and after and middle (with material change)
                column[stackIndex].second = startZ - currentHeight;
                column.insert(column.begin() + stackIndex + 1, {material, endZ - startZ});
                column.insert(column.begin() + stackIndex + 2, {column[stackIndex].first, (currentHeight + height - endZ)});
                stackIndex += 2; // We added 2 layers
            } else if (currentHeight >= startZ && currentHeight + height <= endZ) {
                // The whole layer material is changed
                column[stackIndex].first = material;
            } else if (currentHeight + height < endZ) {
                // The top part of the layer is changed
                column[stackIndex].second = startZ - currentHeight;
                column.insert(column.begin() + stackIndex + 1, {material, (currentHeight + height) - startZ});
                stackIndex ++;
            } else {
                // The bottom part of the layer is changed
                column.insert(column.begin() + stackIndex + 1, {column[stackIndex].first, endZ});
                column[stackIndex].second = endZ - currentHeight;
                stackIndex++;
            }
        }
        currentHeight += height;
        stackIndex ++;
    }
    return this;*/
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
    GridF distanceMap = GridF(maxX - minX, maxY - minY, 1, 1.f);
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
    GridF distanceMap = GridF(maxX - minX, maxY - minY, 1, 1.f);
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
    Vector3 minPos = (AABBox.min() / this->scaling) - this->translation;
    Vector3 maxPos = (AABBox.max() / this->scaling) - this->translation;

    int minX = std::max(0, int(minPos.x));
    int maxX = std::min(int(this->getSizeX()), int(maxPos.x));
    int minY = std::max(0, int(minPos.y));
    int maxY = std::min(int(this->getSizeY()), int(maxPos.y));
    float minZ = std::max(0.f, minPos.z);
    float maxZ = std::min(float(maxPos.z), 40.f);
    float zResolution = 1.f;
    int nbEvaluations = 0;
    #pragma omp parallel for collapse(2)
    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            float z = minZ;
//            if (x == 50 && y == 50) {
//                int a = 0;
//            }
            while ( z < maxZ) {
                nbEvaluations++;
                auto [totalValue, initialDensityAndProportion] = patch->getMaterialsAndTotalEvaluation((Vector3(x, y, z) + this->translation) * this->scaling); //patch->getDensityAtPosition(Vector3(x, y, z));
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
                    if (maxValue < ImplicitPatch::isovalue)
                        selectedMaterial = AIR;
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
    this->cleanLayers();
//    this->reorderLayers();
    _historyIndex ++;
    return;
}

Mesh LayerBasedGrid::getGeometry(const Vector3& dimensions)
{
    std::vector<Vector3> vertices;

    Vector3 originalDimensions = Vector3(this->getSizeX(), this->getSizeY(), 1);
    Vector3 finalDimensions = dimensions;
    if (!dimensions.isValid())
        finalDimensions = originalDimensions;
    finalDimensions.z = 1; // Force the Z to 1

    auto copiedLayers = this->layers.resizeNearest(finalDimensions);

//    auto layersAndHeights = this->getMaterialAndHeightsGrid();
//    auto layers = layersAndHeights.first;
//    auto heights = layersAndHeights.second;

    for (int x = 0; x < finalDimensions.x; x++) {
        for (int y = 0; y < finalDimensions.y; y++) {
            auto layers = copiedLayers.at(x, y);
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
    m.scale(originalDimensions / finalDimensions);
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

bool LayerBasedGrid::checkIsInGround(const Vector3& position)
{
    return !(isIn(this->getValue(position), LayerBasedGrid::invisibleLayers));
}

