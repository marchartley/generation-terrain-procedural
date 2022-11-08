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
    //: sizeX(nx), sizeY(ny), sizeZ(nz)
{
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(nx, ny, 1, {{TerrainTypes::BEDROCK, nz}}); // Initial terrain = 1 layer of bedrock
}

void LayerBasedGrid::createMesh() {
/*
//    this->computeGroups();

    std::vector<Vector3> voxelVertices;
    std::vector<Vector3> colors;
    this->applyToVoxels([&](std::shared_ptr<Voxel> v) -> void {
        if ((bool)*v) {
            // Add the vertices to the global mesh
            std::vector<Vector3> vertice = v->getMeshVertices();
            voxelVertices.insert(voxelVertices.end(), vertice.begin(), vertice.end());
            // Add the colors to each vertex
            int X = 6; // Start with 6 faces

            for(auto& n : v->neighbors)
                if (n.second && (bool)*n.second)
                    X--;    // Remove a face per neighbor
            X *= 6; // Multiply the number of face by the 6 vertex that defines it (2 triangles)
            for (int x = 0; x < X; x++) {
                colors.push_back(Vector3(random_gen::generate(), random_gen::generate(), random_gen::generate())); //(v->isOnGround ? 1.0 : 0.0), (v->isOnGround ? 0.0 : 1.0), 1.0));
//                        colors->push_back(HSVtoRGB((voxels[i][j][k]->group/((float)Voxel::voxelGroups.size()+1)), 1.0, 1.0));
            }
        }
    });

    this->mesh.colorsArray = colors;
    this->mesh.fromArray(voxelVertices);
    this->mesh.update();*/
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

void LayerBasedGrid::from2DGrid(Grid grid)
{
    // The materials for now will only be dirt and water
    this->layers = Matrix3<std::vector<std::pair<TerrainTypes, float>>>(grid.getSizeX(), grid.getSizeY(), 1);
    for (int x = 0; x < grid.getSizeX(); x++) {
        for (int y = 0; y < grid.getSizeY(); y++) {
            this->layers.at(x, y) = {
                { TerrainTypes::DIRT, grid.getHeight(x, y) },
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

Matrix3<float> LayerBasedGrid::voxelize(int fixedHeight)
{
    float maxHeight = (fixedHeight == -1 ? this->getSizeZ() : fixedHeight);
    Matrix3<float> values = Matrix3<float>(this->getSizeX(), this->getSizeY(), maxHeight, LayerBasedGrid::densityFromMaterial(TerrainTypes::WATER));

    for (int x = 0; x < values.sizeX; x++) {
        for (int y = 0; y < values.sizeY; y++) {
            auto& stack = this->layers.at(x, y);
            float currentHeight = 0.f;
            for (size_t i = 0; i < stack.size(); i++) {
                float density = LayerBasedGrid::densityFromMaterial(stack[i].first);
                float height = stack[i].second;

                for (int z = std::floor(currentHeight); z < std::min(std::floor(currentHeight + height), maxHeight); z++) {
                    values.at(x, y, z) = density;
                }
                currentHeight += height;
            }
        }
    }
    return values;
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
                            for (size_t neighbor = 0; neighbor < layers.at(nPos).size(); neighbor++) {
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
                        }
                    }
                }
                currentHeight = higherPos;
            }
        }
    }
    this->cleanLayers();
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

void LayerBasedGrid::add(Patch2D &patch, TerrainTypes material, bool applyDistanceFalloff, float distancePower)
{
    auto [minPos, maxPos] = patch.shape.AABBox();
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
            float height = patch.getHeight(Vector3(x, y));
            float strength = distanceMap.at(x - minX, y - minY);
            this->addLayer(Vector3(x, y), height * strength, material);
        }
    }
    this->reorderLayers();
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
    return (density > 0 ? TerrainTypes::DIRT : TerrainTypes::WATER);
}

float LayerBasedGrid::densityFromMaterial(TerrainTypes material)
{
    float valMin, valMax;
    std::tie(valMin, valMax) = LayerBasedGrid::materialLimits[material];
    return (valMin + valMax) / 2.f;
}

Patch2D::Patch2D()
    : Patch2D(Vector3(0, 0, 0), ShapeCurve({Vector3(-1, -1, 0), Vector3(-1, 1, 0), Vector3(1, 1, 0), Vector3(1, -1, 0)}), [](Vector3 _pos) {return 0.f; })
{

}

Patch2D::Patch2D(Vector3 pos, ShapeCurve shape, std::function<float (Vector3)> heightFunction)
    : position(pos), shape(shape), heightFunction(heightFunction)
{

}

float Patch2D::getHeight(Vector3 position)
{
    Vector3 p = position - this->position;
//    if (this->shape.grow(.1f).translate(Vector3(.001f, 0.f)).inside(p))
        return this->heightFunction(p);
//    return 0.f;
}
