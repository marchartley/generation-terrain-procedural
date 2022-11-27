#include "TerrainGen/LayerBasedGrid.h"

#include "TerrainModification/UnderwaterErosion.h"

#include "Graphics/Shader.h"

#include "Utils/FastNoiseLit.h"
#include "Utils/Globals.h"

std::map<TerrainTypes, std::pair<float, float>> LayerBasedGrid::materialLimits = {
    {TerrainTypes::AIR, {-10.f, -1.f} },
    {TerrainTypes::WATER, {-1.f, 0.f}},
    {TerrainTypes::SAND, {0.f, .5f}},
    {TerrainTypes::DIRT, {.5f, 1.f}},
    {TerrainTypes::ROCK, {1.f, 2.f}},
    {TerrainTypes::BEDROCK, {2.f, 3.f}},
};

LayerBasedGrid::LayerBasedGrid() : LayerBasedGrid(10, 10, 1.f)
{

}
LayerBasedGrid::LayerBasedGrid(int nx, int ny, float nz)
{
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(nx, ny, 1, {{TerrainTypes::BEDROCK, nz}}); // Initial terrain = 1 layer of bedrock
    this->previousState = layers;
}

void LayerBasedGrid::createMesh() {
}

void LayerBasedGrid::display() {
    this->mesh.display();
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
    for (auto& stack : layers) {
        bool reorderingNeeded = true;
        while (reorderingNeeded) {
            reorderingNeeded = false;
            for (size_t i = 0; i < stack.size() - 1; i++) {
                float densA = LayerBasedGrid::densityFromMaterial(stack[i].first);
                float densB = LayerBasedGrid::densityFromMaterial(stack[i + 1].first);
//                if (densA > 0.f && densB > 0.f) { // If it's AIR or WATER, don't do anything, I guess
                if ((stack[i].first == TerrainTypes::SAND || stack[i + 1].first == TerrainTypes::SAND) && (densA <= 0.f || densB < 0.f)) {
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
                }
            }
        }
    }
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
        if (stack.back().first != TerrainTypes::AIR && stack.back().first != TerrainTypes::WATER) {
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
            startCounting = (mat != AIR && mat != WATER);
        if (startCounting) {
            currentHeight += height;
        }
    }
    return currentHeight;
}

float LayerBasedGrid::getHeight(Vector3 pos)
{
    return this->getHeight(pos.x, pos.y);
}

void LayerBasedGrid::from2DGrid(Grid grid)
{
    // The materials for now will only be dirt and water
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(grid.getSizeX(), grid.getSizeY(), 1);
    for (int x = 0; x < grid.getSizeX(); x++) {
        for (int y = 0; y < grid.getSizeY(); y++) {
            this->layers.at(x, y) = {
                { TerrainTypes::SAND, grid.getHeight(x, y) },
                { TerrainTypes::WATER, grid.maxHeight - grid.getHeight(x, y) }
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

    for (int x = 0; x < values.sizeX; x++) {
        for (int y = 0; y < values.sizeY; y++) {
            for (int z = 0; z < maxHeight; z++) {
                auto materialAndVolume = this->getKernel(Vector3(x, y, z), kernelSize);
                float solidMaterial = 0.f;
                float airMaterial = 0.f;
                for (auto& [mat, vol] : materialAndVolume) {
                    if (mat == AIR || mat == WATER)
                        airMaterial += vol;
                    else
                        solidMaterial += vol;
                }
                float val = (2.f * solidMaterial / (solidMaterial + airMaterial)) - 1.f;
                values.at(x, y, z) = (val > 0.f ? LayerBasedGrid::densityFromMaterial(SAND) : LayerBasedGrid::densityFromMaterial(AIR));
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
    float distanceToGrid = Vector3::signedDistanceToBoundaries(currPos, minPos, maxPos);
    float distanceToGridDT = Vector3::signedDistanceToBoundaries(currPos + dir, minPos, maxPos);
    // Continue while we are in the grid or we are heading towards the grid
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {
        float isoval = LayerBasedGrid::densityFromMaterial(this->getMaterialAndHeight(currPos).first);
        if (isoval > 0.0) {
            return currPos;
        }
        currPos += dir;
        distanceToGrid = Vector3::signedDistanceToBoundaries(currPos, minPos, maxPos);
        distanceToGridDT = Vector3::signedDistanceToBoundaries(currPos + dir, minPos, maxPos);
    }
    return Vector3(false);
}

Vector3 LayerBasedGrid::getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    return this->getFirstIntersectingStack(origin, dir, minPos, maxPos);
}

std::pair<Matrix3<int>, Matrix3<float> > LayerBasedGrid::getMaterialAndHeightsGrid()
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
    return {materials, heights};
}

void LayerBasedGrid::thermalErosion()
{
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
    for (auto& stack : layers) {
        for (size_t i = 0; i < stack.size(); i++) {
            if (stack[i].second < minLayerHeight) {
                stack.erase(stack.begin() + i);
            }
        }
    }
}

void LayerBasedGrid::add(Patch2D patch, TerrainTypes material, bool applyDistanceFalloff, float distancePower)
{
//    auto [minPos, maxPos] = patch.shape.AABBox();
    Vector3 minPos = patch.boundMin;
    Vector3 maxPos = patch.boundMax;
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
    return;
}

void LayerBasedGrid::add(Patch3D patch, TerrainTypes material, bool applyDistanceFalloff, float distancePower)
{
    Vector3 minPos = patch.boundMin;
    Vector3 maxPos = patch.boundMax;
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
    return;
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

Patch::Patch()
    : Patch(Vector3(0, 0, 0), Vector3(-1, -1, -1), Vector3(1, 1, 1), [](Vector3 _pos) {return 0.f; })
{

}

Patch::Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, evalFunction, [&](Vector3 pos) { return this->evaluate(pos); }, densityValue)
{
    this->compositionFunction = [&](Vector3 pos) {
        float eval = this->evaluate(pos);
        return eval;
    };
}

Patch::Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, evalFunction, compositionFunction, [&](Vector3 _pos) { return 1.f; }, densityValue)
{
    this->alphaDistanceFunction = [&](Vector3 _pos) {
        float distPower = 2.f;
        float minDistance = 0.8f;
        float distance = -Vector3::signedDistanceToBoundaries(_pos, this->boundMin, this->boundMax);
        float maxDistance = ((this->boundMax - this->boundMin) * .5f).maxComp();
        float ratio = 1.f - std::clamp(distance / maxDistance, 0.f, 1.f);
        return 1.f - (ratio <= minDistance ? 0.f : std::pow((ratio - minDistance) / (1.f - minDistance), distPower));
    };
}

Patch::Patch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : position(pos), boundMin(boundMin), boundMax(boundMax), evalFunction(evalFunction), densityValue(densityValue), compositionFunction(compositionFunction), alphaDistanceFunction(alphaDistanceFunction)
{
    this->_cachedHeights = Matrix3<float>(boundMax - boundMin, -1000.f);
    this->_cachedHeights.raiseErrorOnBadCoord = false;
    this->_cachedHeights.defaultValueOnBadCoord = 0.f;
}

float Patch::getMaxHeight(Vector3 position)
{
    return 0.f;
}

float Patch::evaluate(Vector3 position)
{
    return 0.f;
}

//Patch::~Patch()
//{

//}


Patch2D::Patch2D()
    : Patch2D(Vector3(0, 0, 0), Vector3(-1, -1, 0), Vector3(1, 1, 0), [](Vector3 _pos) {return 0.f; })
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, heightFunction, densityValue)
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, heightFunction, compositionFunction, densityValue)
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, heightFunction, compositionFunction, alphaDistanceFunction, densityValue)
{

}

float Patch2D::getMaxHeight(Vector3 position)
{
    Vector3 p = position; // - this->position;
    if (boundMin.x > p.x || p.x > boundMax.x || boundMin.y > p.y || p.y > boundMax.y/* || boundMin.z > p.z || p.z > boundMax.z*/)
        return 0.f;
//    if (this->shape.grow(.1f).translate(Vector3(.001f, 0.f)).inside(p))
    return this->evalFunction(p.xy());
//    return 0.f;
}

float Patch2D::evaluate(Vector3 position)
{
    Vector3 p = position; // - this->position;
//    if (!shape.inside(p.xy())) return 0.f;

    if (boundMin.x > p.x || p.x > boundMax.x || boundMin.y > p.y || p.y > boundMax.y/* || boundMin.z > p.z || p.z > boundMax.z*/)
        return 0.f;

    float val = this->getMaxHeight(p.xy());
    return (val > p.z ? 1.f : 0.f);
}

Patch3D::Patch3D()
    : Patch3D(Vector3(0, 0, 0), Vector3(-1, -1, -1), Vector3(1, 1, 1), [](Vector3 _pos) {return 0.f; })
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, evalFunction, densityValue)
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, evalFunction, compositionFunction, densityValue)
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : Patch(pos, boundMin, boundMax, evalFunction, compositionFunction, alphaDistanceFunction, densityValue)
{
}

Patch3D::Patch3D(const Patch3D &copy)
{
    this->position = copy.position;
    this->boundMin = copy.boundMin;
    this->boundMax = copy.boundMax;
    this->evalFunction = copy.evalFunction;
    this->compositionFunction = copy.compositionFunction;
    this->alphaDistanceFunction = copy.alphaDistanceFunction;
    this->densityValue = copy.densityValue;
    float test = compositionFunction(Vector3(10,10,10));
    std::cout << test << std::endl;
    this->_cachedHeights = Matrix3<float>(boundMax - boundMin, -1000.f);
    this->_cachedHeights.raiseErrorOnBadCoord = false;
    this->_cachedHeights.defaultValueOnBadCoord = 0.f;
}

float Patch3D::getMaxHeight(Vector3 position)
{
    if (this->_cachedHeights.at(position.xy()) < -1.f) {
        float resolution = 0.1;
        float maxHeight = this->boundMax.z;
        while (maxHeight >= boundMin.z && evaluate(Vector3(position.x, position.y, maxHeight)) <= 0.f)
            maxHeight -= resolution;
        this->_cachedHeights.at(position.xy()) = maxHeight;
        return maxHeight;
    }
    return this->_cachedHeights.at(position.xy());
}

float Patch3D::evaluate(Vector3 position)
{
    Vector3 p = position; // - this->position;
    if (boundMin.x > p.x || p.x > boundMax.x || boundMin.y > p.y || p.y > boundMax.y || boundMin.z > p.z || p.z > boundMax.z)
        return 0.f;
    return (this->evalFunction(p) > 0.f ? 1.f : 0.f);
}

Patch3D Patch3D::stack(Patch3D *P1, Patch3D *P2)
{
    Vector3 P1center = P1->position;
    Vector3 P2center = P2->position;

    Vector3 boundMin = Vector3::min(P1->boundMin + P1center, P2->boundMin + P2center);
    Vector3 boundMax = Vector3::max(P1->boundMax + P1center, P2->boundMax + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->boundMax.z - P1->boundMin.z) + (P2->boundMax.z - P2->boundMin.z);

    Vector3 newCenter = boundMin;

    boundMax -= newCenter;
    boundMin -= newCenter;

    std::function<float(Vector3)> evalFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float zOffset = P1->getMaxHeight(pos1);
        pos2.z -= zOffset;
        float val1 = P1->evaluate(pos1);
        float val2 = P2->evaluate(pos2);
        return (val1 + val2 > 0.f ? 1.f : 0.f);
    };

    std::function<float(Vector3)> compositionFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float zOffset = P1->getMaxHeight(pos1);
        pos2.z -= zOffset;
        float val1 = P1->get(pos1);
        float val2 = P2->get(pos2);
        return (val1 > 0.f ? val1 : val2);
    };

    std::function<float(Vector3)> alphaDistanceFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float zOffset = P1->getMaxHeight(pos1);
        pos2.z -= zOffset;
        return std::clamp(std::max(P1->getAlphaValue(pos), P2->getAlphaValue(pos)), 0.f, 1.f);
    };

    return Patch3D(newCenter, boundMin, boundMax, evalFunction, compositionFunction, alphaDistanceFunction);
}

Patch3D Patch3D::replace(Patch3D *P1, Patch3D *P2)
{
    Vector3 P1center = P1->position;
    Vector3 P2center = P2->position;

    Vector3 boundMin = Vector3::min(P1->boundMin + P1center, P2->boundMin + P2center);
    Vector3 boundMax = Vector3::max(P1->boundMax + P1center, P2->boundMax + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->boundMax.z - P1->boundMin.z) + (P2->boundMax.z - P2->boundMin.z);

    Vector3 newCenter = boundMin;

    boundMax -= newCenter;
    boundMin -= newCenter;

    std::function<float(Vector3)> evalFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float val1 = P1->evaluate(pos1);
        float val2 = P2->evaluate(pos2);
        return (val1 + val2 > 0.f ? 1.f : 0.f);
    };

    std::function<float(Vector3)> compositionFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float val1 = P1->get(pos1);
        float val2 = P2->get(pos2);
        float check_val2 = val2;

        float alpha2 = P2->getAlphaValue(pos2);

        return (std::abs(check_val2) > 0.f ? val1 * (1.f - alpha2) + val2 * alpha2 : val1);
    };
    std::function<float(Vector3)> alphaDistanceFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float zOffset = P1->getMaxHeight(pos1);
        pos2.z -= zOffset;
        return std::clamp(std::max(P1->getAlphaValue(pos1), P2->getAlphaValue(pos2)), 0.f, 1.f);
    };
    return Patch3D(newCenter, boundMin, boundMax, evalFunction, compositionFunction, alphaDistanceFunction);
}

Patch3D Patch3D::blend(Patch3D *P1, Patch3D *P2)
{
    Vector3 P1center = P1->position;
    Vector3 P2center = P2->position;

    Vector3 boundMin = Vector3::min(P1->boundMin + P1center, P2->boundMin + P2center);
    Vector3 boundMax = Vector3::max(P1->boundMax + P1center, P2->boundMax + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->boundMax.z - P1->boundMin.z) + (P2->boundMax.z - P2->boundMin.z);

    Vector3 newCenter = boundMin;

    boundMax -= newCenter;
    boundMin -= newCenter;

    std::function<float(Vector3)> evalFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float val1 = P1->evaluate(pos1);
        float val2 = P2->evaluate(pos2);
        return (val1 + val2 > 0.f ? 1.f : 0.f);
    };

    std::function<float(Vector3)> compositionFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float alpha1 = P1->getAlphaValue(pos1);
        float alpha2 = P2->getAlphaValue(pos2);
        float val1 = P1->get(pos1);
        float val2 = P2->get(pos2);
        if (val1 <= 0.f) {
            alpha1 = 0.f;
            alpha2 = 1.f;
        } else if (val2 <= 0.f) {
            alpha1 = 1.f;
            alpha2 = 0.f;
        }
        return (alpha1 * val1 + alpha2 * val2) / (alpha1 + alpha2);
//        float check_val1 = val1;
//        float check_val2 = val2;

//        return (check_val1 > 0.f && check_val2 > 0.f ? (val1 + val2) * .5f : (check_val1 > 0.f ? val1 : val2));
    };
    std::function<float(Vector3)> alphaDistanceFunction = [P1, P2, P1center, P2center, newCenter](Vector3 pos) {
        Vector3 pos1 = pos - (P1center - newCenter);
        Vector3 pos2 = pos - (P2center - newCenter);
        float zOffset = P1->getMaxHeight(pos1);
        pos2.z -= zOffset;
        return std::clamp(std::max(P1->getAlphaValue(pos1), P2->getAlphaValue(pos2)), 0.f, 1.f);
    };
    return Patch3D(newCenter, boundMin, boundMax, evalFunction, compositionFunction, alphaDistanceFunction);
}
