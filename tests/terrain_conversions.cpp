#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "Utils/stb_image.h"
#include "Utils/stb_image_write.h"
#include "Utils/stb_image_resize.h"
#include "TerrainGen/LayerBasedGrid.h"
#include "TerrainGen/TerrainSampling.h"
#include <iostream>
#include <memory>
#include <stdexcept>

static int checks = 0;
static void require(bool condition, const std::string& message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
static void near(float actual, float expected, const std::string& message, float tolerance = 1e-4f)
{
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
            message + ": expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}
static void sameGrid(const GridF& a, const GridF& b, const std::string& message)
{
    require(a.getDimensions() == b.getDimensions(), message + " dimensions");
    for (size_t i = 0; i < a.size(); ++i) near(a[i], b[i], message);
}
static void heightmaps()
{
    for (float elevation : {0.f, .25f, 1.f, 5.f, 5.25f, 55.75f}) {
        Heightmap source;
        source.heights = GridF(6, 5, 1, elevation);
        source.scaling = Vector3(2, 3, 4);
        source.translation = Vector3(1, 2, 3);
        VoxelGrid voxels(6, 5, 10);
        voxels.from2DGrid(source);
        Heightmap result;
        result.fromVoxelGrid(voxels);
        sameGrid(result.heights, source.heights, "heightmap/voxel round trip");
        require(result.scaling == source.scaling && result.translation == source.translation, "terrain transform copied");
        LayerBasedGrid layers;
        layers.from2DGrid(source);
        voxels.fromLayerBased(layers);
        result.fromVoxelGrid(voxels);
        sameGrid(result.heights, source.heights, "heightmap/layers/voxel round trip");
        std::unique_ptr<ImplicitPrimitive> implicit(ImplicitPrimitive::fromHeightmap(source.heights));
        std::unique_ptr<ImplicitPatch> decodedHeight(ImplicitPatch::fromJson(implicit->toJson()));
        sameGrid(decodedHeight->getVoxelized(), implicit->getVoxelized(), "heightmap snapshot JSON round trip");
        layers.fromImplicit(implicit.get());
        result.fromLayerGrid(layers);
        sameGrid(result.heights, source.heights, "heightmap/implicit/layers round trip");
    }
    Heightmap source;
    source.heights = GridF(6, 5, 1);
    source.heights.iterate([&](int x, int y, int) { source.heights(x, y) = .25f + x * .35f + y * .17f; });
    VoxelGrid voxels;
    voxels.from2DGrid(source);
    Heightmap result;
    result.fromVoxelGrid(voxels);
    sameGrid(result.heights, source.heights, "sloping heightmap round trip");
    source.scaling = Vector3(2, 3, 1);
    source.translation = Vector3(4, 5, 0);
    voxels.from2DGrid(source, Vector3(1, 1), Vector3(5, 4), 2.f);
    near(voxels.getWorldPos(Vector3(2, 2)).x(), source.getWorldPos(Vector3(2, 2)).x(), "subsection world X");
    near(voxels.getWorldPos(Vector3(2, 2)).y(), source.getWorldPos(Vector3(2, 2)).y(), "subsection world Y");
    auto dense = source.getVoxelized(Vector3(6, 5, 8), Vector3(2, 2, 2));
    require(dense.getDimensions() == Vector3(12, 10, 16), "heightmap requested sampling density");
    near(dense(4, 4, 2), source.heights(2, 2) - 1.f, "heightmap sample coordinates");
}
static void columns()
{
    for (const auto& column : std::vector<std::vector<float>>{
            {1}, {0}, {-1}, {1, 1, 1, 1}, {1, -1, 1, -1}, {-1, -1, 1, 1},
            {0, 1, 0, -1}, {4, -1, 1, -4}, {4, -1, .1f, -4}}) {
        GridF field(1, 1, column.size());
        for (size_t z = 0; z < column.size(); ++z) field(0, 0, z) = column[z];
        VoxelGrid voxels;
        voxels.setVoxelValues(field);
        LayerBasedGrid layers;
        layers.fromVoxelGrid(voxels);
        near(layers.getHeight(0, 0), TerrainSampling::topHeight(field, 0, 0), "column top surface");
        const auto reconstructed = layers.voxelize();
        float height = 0.f;
        for (const auto& [material, thickness] : layers.getLayers()(0, 0)) {
            require(thickness > 0.f && std::isfinite(thickness), "valid interval thickness");
            height += thickness;
            if (height < float(field.sizeZ - 1)) {
                const int z = int(std::floor(height));
                near(field(0, 0, z) * (1.f - (height - z)) + field(0, 0, z + 1) * (height - z),
                     0.f, "layer interface follows scalar zero crossing");
                const float first = reconstructed(0, 0, z), second = reconstructed(0, 0, z + 1);
                near(z + first / (first - second), height, "reconstructed cave interface");
            }
        }
        std::unique_ptr<ImplicitPrimitive> primitive(ImplicitPrimitive::fromVoxelField(field));
        ImplicitNaryOperator container;
        container.addChild(primitive.release());
        VoxelGrid restored;
        restored.fromImplicit(&container);
        sameGrid(restored.getVoxelValues(), field, "voxel/implicit round trip including singleton Z");
        auto encoded = container.composables.front()->toJson();
        std::unique_ptr<ImplicitPatch> decoded(ImplicitPatch::fromJson(encoded));
        sameGrid(decoded->getVoxelized(), field, "converted field JSON round trip");
        container.deleteAllChildren();
    }
    VoxelGrid voxels;
    GridF alternating(1, 1, 4);
    alternating(0, 0, 0) = 1; alternating(0, 0, 1) = -1;
    alternating(0, 0, 2) = 1; alternating(0, 0, 3) = -1;
    voxels.setVoxelValues(alternating);
    LayerBasedGrid layers;
    layers.fromVoxelGrid(voxels);
    const auto& stack = layers.getLayers()(0, 0);
    require(stack.size() == 3, "cave intervals retained");
    near(stack[0].second, .5f, "first interface");
    near(stack[1].second, 1.f, "cave height");
    near(stack[2].second, 1.f, "roof thickness");
}
static void reconstruction()
{
    LayerBasedGrid layers(4, 4, 1.f);
    GridF reference;
    for (auto material : {SAND, ROCK, BEDROCK, CORAL}) {
        for (auto& stack : layers.getLayers()) stack = {{material, 5.25f}};
        const auto values = layers.voxelize();
        if (reference.empty()) reference = values;
        else sameGrid(values, reference, "material-independent geometry");
        near(TerrainSampling::topHeight(values, 1, 1), 5.25f, "fractional layer top");
    }
    for (auto& stack : layers.getLayers()) stack = {{SAND, 1.f}};
    for (float width : {.5f, 1.f, 1.5f, 2.5f}) {
        for (auto center : {Vector3(1.1f, 1.4f, 1.f), Vector3(0, 0, 0)}) {
            float sum = 0.f;
            for (auto [material, volume] : layers.getKernel(center, width)) {
                require(volume >= -1e-5f, "nonnegative kernel overlap");
                sum += volume;
            }
            near(sum, width * width * width, "kernel volume conservation");
        }
    }
    near(layers.getKernel(Vector3(1.5f, 1.5f, 1.f), 1.f)[SAND], .5f, "half-filled kernel");
    near(layers.voxelize(-1, 1.f)(1, 1, 1), 0.f, "box reconstruction half occupancy surface");
    auto dense = layers.getVoxelized(Vector3(4, 4, 4), Vector3(2, 2, 2));
    require(dense.getDimensions() == Vector3(8, 8, 8), "layer requested sampling density");
    near(dense(2, 2, 2), 0.f, "layer sample coordinates");
    Heightmap source;
    source.heights = GridF(3, 2, 1, 4.f);
    layers.from2DGrid(source);
    layers.getLayers()(0, 0).clear();
    layers.reset();
    require(layers.getSizeX() == 3 && layers.getSizeY() == 2, "reset keeps converted dimensions");
    near(layers.getHeight(0, 0), 4.f, "reset keeps converted baseline");
    layers.getLayers()(0, 0) = {{SAND, 0}, {AIR, 0}, {SAND, 2}, {SAND, 0}, {WATER, 0}};
    layers.cleanLayers();
    require(layers.getLayers()(0, 0).size() == 1, "all consecutive zero layers removed");
}
static void layerCacheInvalidation()
{
    LayerBasedGrid layers(4, 3, 2.f);
    const auto first = layers.getVoxelized(Vector3(4, 3, 8));
    sameGrid(first, layers.getVoxelized(Vector3(4, 3, 8)), "unchanged layer field cache");
    const auto history = layers._historyIndex;
    layers.getLayers()(1, 1) = {{ROCK, 4.25f}};
    require(layers._historyIndex == history, "test covers edits without history updates");
    auto changed = layers.getVoxelized(Vector3(4, 3, 8));
    near(TerrainSampling::topHeight(changed, 1, 1), 4.25f, "direct edits invalidate layer field cache");
    layers.getLayers()(1, 1) = {{ROCK, 1.25f}, {AIR, 2.f}, {SAND, 1.f}};
    const auto maps = layers.getMaterialAndHeightsGrid();
    const auto tops = layers.getLayerTopHeights();
    near(maps.second(1, 1, 1), 2.f, "thickness API retained");
    near(tops(1, 1, 0), 1.25f, "first prefix height");
    near(tops(1, 1, 1), 3.25f, "air prefix height");
    near(tops(1, 1, 2), 4.25f, "upper prefix height");
    layers.getLayers()(1, 1)[0].second = 2.f;
    near(layers.getLayerTopHeights()(1, 1, 2), 5.f, "prefix cache detects direct thickness edit");
    auto dense = layers.getVoxelized(Vector3(4, 3, 8), Vector3(2, 2, 2));
    require(dense.getDimensions() == Vector3(8, 6, 16), "layer cache sampling scale key");
    auto normal = layers.getVoxelized(Vector3(4, 3, 8));
    require(normal.getDimensions() == Vector3(4, 3, 8), "layer cache restores normal resolution");
}

static void implicitSampling()
{
    ImplicitPrimitive source;
    source.position = Vector3(0, 0, 8);
    source.setDimensions(Vector3(4, 3, 70));
    source.setSupportDimensions(Vector3(4, 3, 70));
    source.material = DIRT;
    source.evalFunction = [](const Vector3& p) { return std::clamp(.5f + std::min(p.z(), 62.25f - p.z()), 0.f, 1.f); };
    LayerBasedGrid layers;
    layers.fromImplicit(&source);
    require(layers.getSizeX() == 4 && layers.getSizeY() == 3, "implicit replaces destination dimensions");
    const auto& stack = layers.getLayers()(1, 1);
    require(stack.size() == 2 && stack[0].first == AIR && stack[1].first == DIRT, "implicit preserves material and empty space below");
    near(stack[0].second, 8.f, "translated implicit bottom");
    near(layers.getHeight(1, 1), 70.25f, "implicit height above old 40-unit cap");
    Heightmap heightmap;
    heightmap.fromImplicit(&source);
    near(heightmap.getHeight(1, 1), 70.25f, "implicit heightmap projection");
    auto a = source.getVoxelized(Vector3(4, 3, 80));
    auto b = source.getVoxelized(Vector3(2, 3, 20));
    require(b.getDimensions() == Vector3(2, 3, 20), "implicit cache size key");
    auto c = source.getVoxelized(Vector3(2, 3, 20), Vector3(2, 2, 2));
    require(c.getDimensions() == Vector3(4, 6, 40), "implicit cache scale key");
    near(c(2, 2, 19), source.getSignedEvaluation(Vector3(1, 1, 9.5f)), "implicit sampling coordinate");
    auto d = source.getVoxelized(Vector3(4, 3, 80));
    sameGrid(a, d, "implicit cache recomputes expanded domain");
}
int main()
{
    try {
        heightmaps(); columns(); reconstruction(); implicitSampling(); layerCacheInvalidation();
        std::cout << "PASS: " << checks << " terrain conversion checks\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
