#include "TerrainGen/TerrainSampling.h"
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
    return getValue(pos.x(), pos.y(), pos.z());
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
                        ((LayerBasedGrid::invisibleLayers.count(stack[i].first) != 0) || (LayerBasedGrid::invisibleLayers.count(stack[i+1].first) != 0))) {
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
        if ((LayerBasedGrid::invisibleLayers.count(stack.back().first) == 0)) { // != TerrainTypes::AIR && stack.back().first != TerrainTypes::WATER) {
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
            startCounting = (LayerBasedGrid::invisibleLayers.count(mat) == 0); //(mat != AIR && mat != WATER);
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

void LayerBasedGrid::from2DGrid(const Heightmap& grid)
{
    layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(grid.getSizeX(), grid.getSizeY(), 1);
    for (int x = 0; x < grid.getSizeX(); ++x) {
        for (int y = 0; y < grid.getSizeY(); ++y) {
            const float height = std::max(0.f, grid.heights(x, y));
            if (height > 0.f) layers(x, y).push_back({SAND, height});
        }
    }
    scaling = grid.scaling;
    translation = grid.translation;
    previousState = layers;
    ++_historyIndex;
}

void LayerBasedGrid::fromVoxelGrid(VoxelGrid& voxelGrid)
{
    const auto values = voxelGrid.getVoxelValues();
    layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(values.sizeX, values.sizeY, 1);
    for (int x = 0; x < int(values.sizeX); ++x) {
        for (int y = 0; y < int(values.sizeY); ++y) {
            if (values.sizeZ == 0) continue;
            float start = 0.f;
            float previous = values(x, y, 0);
            bool solid = previous > 0.f;
            for (int z = 1; z <= int(values.sizeZ); ++z) {
                const float current = z < int(values.sizeZ) ? values(x, y, z) : TerrainSampling::aboveColumn(values, x, y);
                const bool nextSolid = current > 0.f;
                if (nextSolid != solid) {
                    const float boundary = z - 1.f + previous / (previous - current);
                    if (boundary > start)
                        layers(x, y).push_back({solid ? BEDROCK : AIR, boundary - start});
                    start = boundary;
                    solid = nextSolid;
                }
                previous = current;
            }
        }
    }
    cleanLayers();
    scaling = voxelGrid.scaling;
    translation = voxelGrid.translation;
    previousState = layers;
    ++_historyIndex;
}

void LayerBasedGrid::fromImplicit(ImplicitPatch* implicitTerrain)
{
    if (!implicitTerrain) throw std::invalid_argument("Missing implicit terrain");
    const auto bounds = implicitTerrain->getBBox();
    const int nx = std::max(1, int(std::ceil(bounds.max().x())));
    const int ny = std::max(1, int(std::ceil(bounds.max().y())));
    const float top = std::max(1.f, std::ceil(bounds.max().z()) + 1.f);
    layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(nx, ny, 1);
    const ImplicitPatch* materialSource = implicitTerrain;
    while (typeid(*materialSource) == typeid(ImplicitNaryOperator)) {
        const auto* group = static_cast<const ImplicitNaryOperator*>(materialSource);
        if (group->composables.size() != 1) break;
        materialSource = group->composables.front();
    }
    auto materialAt = [&](int x, int y, float z) {
        const Vector3 p(x, y, z);
        if (implicitTerrain->getSignedEvaluation(p) <= 0.f) return AIR;
        TerrainTypes selected = SAND;
        float greatest = 0.f;
        for (const auto& [material, amount] : materialSource->getMaterials(p)) {
            if ((LayerBasedGrid::invisibleLayers.count(material) == 0) && amount > greatest) {
                selected = material;
                greatest = amount;
            }
        }
        return selected;
    };
    // Refine every detected interface instead of rounding its height to a cell.
    const float step = ImplicitPatch::zResolution;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            float start = 0.f;
            float previousZ = 0.f;
            auto previousMaterial = materialAt(x, y, 0.f);
            for (int i = 1; i <= int(std::ceil(top / step)); ++i) {
                const float z = std::min(top, i * step);
                const auto material = materialAt(x, y, z);
                if (material != previousMaterial) {
                    float lower = previousZ, upper = z;
                    for (int refinement = 0; refinement < 14; ++refinement) {
                        const float middle = (lower + upper) * .5f;
                        if (materialAt(x, y, middle) == previousMaterial) lower = middle;
                        else upper = middle;
                    }
                    const float boundary = (lower + upper) * .5f;
                    if (boundary > start)
                        layers(x, y).push_back({previousMaterial, boundary - start});
                    start = boundary;
                    previousMaterial = material;
                }
                previousZ = z;
            }
            if (top > start) layers(x, y).push_back({previousMaterial, top - start});
        }
    }
    cleanLayers();
    scaling = implicitTerrain->scaling;
    translation = implicitTerrain->translation;
    previousState = layers;
    ++_historyIndex;
}

VoxelGrid LayerBasedGrid::toVoxelGrid()
{
    VoxelGrid vox;
    vox.fromLayerBased(*this);
    return vox;
}

float LayerBasedGrid::sampleSignedField(const Vector3& pos) const
{
    if (layers.empty() || pos.x() < 0.f || pos.y() < 0.f ||
        pos.x() > getSizeX() - 1.f || pos.y() > getSizeY() - 1.f) return -1.f;
    auto columnField = [&](int x, int y) {
        const auto& stack = layers(x, y);
        if (stack.empty()) return -std::max(1.f, pos.z());
        float height = 0.f;
        float distance = std::numeric_limits<float>::infinity();
        bool inside = pos.z() < 0.f && (LayerBasedGrid::invisibleLayers.count(stack.front().first) == 0);
        for (size_t i = 0; i < stack.size(); ++i) {
            const auto [material, thickness] = stack[i];
            const bool solid = (LayerBasedGrid::invisibleLayers.count(material) == 0);
            if (height <= pos.z() && pos.z() < height + thickness) inside = solid;
            height += thickness;
            const bool nextSolid = i + 1 < stack.size() && (LayerBasedGrid::invisibleLayers.count(stack[i + 1].first) == 0);
            if (solid != nextSolid) distance = std::min(distance, std::abs(pos.z() - height));
        }
        if (!std::isfinite(distance)) return -1.f;
        return inside ? distance : -distance;
    };
    const int x0 = int(std::floor(pos.x())), y0 = int(std::floor(pos.y()));
    const int x1 = std::min(x0 + 1, int(getSizeX()) - 1), y1 = std::min(y0 + 1, int(getSizeY()) - 1);
    const float tx = pos.x() - x0, ty = pos.y() - y0;
    return (1.f - ty) * ((1.f - tx) * columnField(x0, y0) + tx * columnField(x1, y0)) +
                  ty * ((1.f - tx) * columnField(x0, y1) + tx * columnField(x1, y1));
}

namespace {
GridF sampleLayerField(const LayerBasedGrid& layers, const Vector3& extent, const Vector3& scale)
{
    GridF values(TerrainSampling::sampleDimensions(extent, scale), -1.f);
    // Each edge crossing constrains the ratio of its endpoint values. Propagate
    // those ratios through consecutive crossings so narrow caves don't drift
    // when the interval representation is sampled back onto a regular grid.
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < int(values.sizeX); ++x) {
        for (int y = 0; y < int(values.sizeY); ++y) {
            for (int z = 0; z < int(values.sizeZ); ++z)
                values(x, y, z) = layers.sampleSignedField(Vector3(x, y, z) / scale);
            int chainStart = 0;
            for (int z = 1; z < int(values.sizeZ); ++z) {
                float& previous = values(x, y, z - 1);
                float& current = values(x, y, z);
                if (previous * current >= 0.f) {
                    chainStart = z;
                    continue;
                }
                float lower = (z - 1.f) / scale.z(), upper = z / scale.z();
                const bool solidBelow = previous > 0.f;
                for (int i = 0; i < 22; ++i) {
                    const float middle = (lower + upper) * .5f;
                    const float field = layers.sampleSignedField(Vector3(x / scale.x(), y / scale.y(), middle));
                    if ((field > 0.f) == solidBelow) lower = middle;
                    else upper = middle;
                }
                const float fraction = std::clamp((lower + upper) * .5f * scale.z() - (z - 1.f), 1e-6f, 1.f - 1e-6f);
                current = -previous * (1.f - fraction) / fraction;
                if (std::abs(current) > 1000.f) {
                    const float divisor = std::abs(current);
                    for (int i = chainStart; i <= z; ++i) values(x, y, i) /= divisor;
                }
            }
        }
    }
    values.defaultValueOnBadCoord = -1.f;
    return values;
}
}

GridF LayerBasedGrid::voxelize(int fixedHeight, float kernelSize) const
{
    const float depth = fixedHeight < 0 ? std::max(2.f, std::ceil(getSizeZ()) + 2.f) : fixedHeight;
    if (kernelSize <= 0.f)
        return sampledField(Vector3(getSizeX(), getSizeY(), depth), Vector3(1, 1, 1));
    GridF values(TerrainSampling::sampleDimensions(Vector3(getSizeX(), getSizeY(), depth), Vector3(1, 1, 1)), -1.f);
    values.iterateParallel([&](int x, int y, int z) {
        if (kernelSize <= 0.f) {
            values(x, y, z) = sampleSignedField(Vector3(x, y, z));
        } else {
            float solidVolume = 0.f;
            for (const auto& [material, volume] : getKernel(Vector3(x + .5f, y + .5f, z), kernelSize))
                if ((LayerBasedGrid::invisibleLayers.count(material) == 0)) solidVolume += volume;
            values(x, y, z) = 2.f * solidVolume / (kernelSize * kernelSize * kernelSize) - 1.f;
        }
    });
    values.defaultValueOnBadCoord = -1.f;
    return values;
}

GridF LayerBasedGrid::sampledField(const Vector3& extent, const Vector3& scale) const
{
    const auto counts = TerrainSampling::sampleDimensions(extent, scale);
    // Public layer stacks are mutable, so history counters alone cannot detect
    // edits. Compare the intervals, not the much larger sampled volume.
    if (_sampledField.getDimensions() != counts || _sampledScale != scale ||
        _fieldLayersSnapshot.getDimensions() != layers.getDimensions() || _fieldLayersSnapshot.data != layers.data) {
        _sampledField = sampleLayerField(*this, extent, scale);
        _fieldLayersSnapshot = layers;
        _sampledScale = scale;
    }
    return _sampledField;
}

GridF LayerBasedGrid::getVoxelized(const Vector3& dimensions, const Vector3& scale)
{
    const auto extent = dimensions.isValid() ? dimensions :
        Vector3(getSizeX(), getSizeY(), std::max(2.f, std::ceil(getSizeZ()) + 2.f));
    return sampledField(extent, scale);
}

std::map<TerrainTypes, float> LayerBasedGrid::getKernel(const Vector3& pos, float kernelSize) const
{
    if (!std::isfinite(kernelSize) || kernelSize <= 0.f)
        throw std::invalid_argument("Layer kernel requires a positive finite size");
    const float half = kernelSize * .5f;
    const Vector3 lower = pos - Vector3(half, half, half);
    const Vector3 upper = pos + Vector3(half, half, half);
    std::map<TerrainTypes, float> volumes{{AIR, kernelSize * kernelSize * kernelSize}};
    const int minX = std::max(0, int(std::floor(lower.x())));
    const int maxX = std::min(int(getSizeX()), int(std::ceil(upper.x())));
    const int minY = std::max(0, int(std::floor(lower.y())));
    const int maxY = std::min(int(getSizeY()), int(std::ceil(upper.y())));
    for (int x = minX; x < maxX; ++x) {
        for (int y = minY; y < maxY; ++y) {
            const float area = (std::min(x + 1.f, upper.x()) - std::max(float(x), lower.x())) *
                               (std::min(y + 1.f, upper.y()) - std::max(float(y), lower.y()));
            auto addVolume = [&](TerrainTypes material, float bottom, float top) {
                const float volume = std::max(0.f, std::min(top, upper.z()) - std::max(bottom, lower.z())) * area;
                if (material != AIR) {
                    volumes[material] += volume;
                    volumes[AIR] -= volume;
                }
            };
            const auto& stack = layers(x, y);
            if (!stack.empty() && lower.z() < 0.f) addVolume(stack.front().first, lower.z(), 0.f);
            float height = 0.f;
            for (auto [material, thickness] : stack) {
                addVolume(material, height, height + thickness);
                height += thickness;
            }
        }
    }
    return volumes;
}

std::pair<TerrainTypes, float> LayerBasedGrid::getMaterialAndHeight(const Vector3& pos)
{
    if (!this->layers.checkCoord(pos.xy()))
        return {TerrainTypes::AIR, 0.f};
    float currentHeight = 0;
    for (auto materialTuple : this->layers.at(pos.xy())) {
        if (currentHeight <= pos.z() && pos.z() < currentHeight + materialTuple.second)
            return materialTuple;
        currentHeight += materialTuple.second;
    }
    // No material found at this position, let's assume z < 0 is ground and z > max is water
//    return (pos.z() < 0 ? {TerrainTypes::DIRT, 0.f} : {TerrainTypes::WATER, 0.f});
    return {TerrainTypes::WATER, 0.f};
}

Vector3 LayerBasedGrid::getFirstIntersectingStack(const Vector3& origin, const Vector3& dir, const Vector3& minPos, const Vector3& maxPos)
{
    AABBox myAABBox = AABBox(Vector3::max(minPos, Vector3()), Vector3::max(maxPos, this->getDimensions()));
//    if (!minPos.isValid()) minPos = Vector3();
//    if (!maxPos.isValid()) maxPos = this->getDimensions();

    Vector3 currPos = origin;
    if (!Vector3::isInBox(currPos, myAABBox.min(), myAABBox.max()))
        currPos = Collision::intersectionRayAABBox(origin, dir, myAABBox);
    myAABBox.expand({myAABBox.min() - Vector3(1, 1, 1), myAABBox.max() + Vector3(1, 1, 1)});
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
    if (currPos.z() <= 0 && Vector3::isInBox(currPos.xy(), myAABBox.min().xy(), myAABBox.max().xy()))
        return currPos.xy(); // There is no ground here, still want to detect the collision... I guess...
    return Vector3::invalid;
}

Vector3 LayerBasedGrid::getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos, const Vector3& maxPos)
{
    return this->getFirstIntersectingStack(origin, dir, minPos, maxPos);
}

std::pair<GridI, GridF> LayerBasedGrid::getMaterialAndHeightsGrid()
{
//    int x = 0;
    if (_materialLayersSnapshot.getDimensions() != layers.getDimensions() ||
        _materialLayersSnapshot.data != layers.data)
    {
        int maxStackSize = 0;
        for (auto& cell : layers)
            maxStackSize = std::max(maxStackSize, int(cell.size()));

        GridI materials(this->getSizeX(), this->getSizeY(), maxStackSize);
        GridF heights(materials.getDimensions());
        _cachedLayerTops = GridF(materials.getDimensions());

        for (size_t i = 0; i < layers.size(); i++) {
            float top = 0.f;
            for (size_t iStack = 0; iStack < layers[i].size(); iStack++) {
                Vector3 index = layers.getCoordAsVector3(i) + Vector3(0, 0, iStack);
                materials.at(index) = (int) layers[i][iStack].first;
                heights.at(index) = (float) layers[i][iStack].second;
                top += layers[i][iStack].second;
                _cachedLayerTops.at(index) = top;
            }
        }
        this->_cachedMaterialAndHeights = {materials, heights};
        _materialLayersSnapshot = layers;
        currentHistoryIndex = _historyIndex;
    }
    return this->_cachedMaterialAndHeights;
}

GridF LayerBasedGrid::getLayerTopHeights()
{
    getMaterialAndHeightsGrid();
    return _cachedLayerTops;
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
    auto& stack = layers(x, y);
    stack.erase(std::remove_if(stack.begin(), stack.end(), [=](const auto& layer) {
        return !std::isfinite(layer.second) || layer.second <= minLayerHeight;
    }), stack.end());
    for (size_t i = 1; i < stack.size();) {
        if (stack[i - 1].first == stack[i].first) {
            stack[i - 1].second += stack[i].second;
            stack.erase(stack.begin() + i);
        } else ++i;
    }
    while (!stack.empty() && (LayerBasedGrid::invisibleLayers.count(stack.back().first) != 0)) stack.pop_back();
}

void LayerBasedGrid::cleanLayers(float minLayerHeight)
{
    size_t nbStacks = layers.size();
    for (size_t iStack = 0; iStack < nbStacks; iStack++) {
        auto& stack = layers[iStack];
        Vector3 pos = layers.getCoordAsVector3(iStack);
        cleanLayer(pos.x(), pos.y(), minLayerHeight);
        /*while (stack.size() > 0 && (LayerBasedGrid::invisibleLayers.count(stack.back().first) != 0)) //(stack.back().first == WATER || stack.back().first == AIR))
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

LayerBasedGrid& LayerBasedGrid::transformLayer(int x, int y, float startZ, float endZ, TerrainTypes material)
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
    return *this;


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
    return *this;*/
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
//    AABBox.expand({Vector3(0, 0, 0), Vector3(100, 100, 0)});
    Vector3 minPos = (AABBox.min() / this->scaling) - this->translation;
    Vector3 maxPos = (AABBox.max() / this->scaling) - this->translation;

    int minX = std::max(0, int(minPos.x()));
    int maxX = std::min(int(this->getSizeX()), int(maxPos.x()));
    int minY = std::max(0, int(minPos.y()));
    int maxY = std::min(int(this->getSizeY()), int(maxPos.y()));
    float minZ = std::max(0.f, minPos.z());
    float maxZ = std::min(float(maxPos.z()), 40.f);
    float zResolution = 1.f;
    int nbEvaluations = 0;
    if (auto as2Dnary = dynamic_cast<Implicit2DNary*>(patch)) {
        as2Dnary->update();
        as2Dnary->updateCache();
        as2Dnary->reevaluateAll();
    }
//    #pragma omp parallel for collapse(2)
    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            float z = minZ;
            /*if (x == 59 && y == 53) {
                int a = 0;
            }*/
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
    auto values = voxelize();
    const Vector3 originalDimensions = values.getDimensions();
    const auto finalDimensions = dimensions.isValid() ? dimensions : originalDimensions;
    values = values.resize(finalDimensions);
    auto mesh = Mesh::applyMarchingCubes(values);
    mesh.scale(originalDimensions / finalDimensions);
    return mesh;
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
    return !((LayerBasedGrid::invisibleLayers.count(this->getValue(position)) != 0));
}

