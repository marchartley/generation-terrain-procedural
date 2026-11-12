#include "TerrainGen/TerrainSampling.h"
#include "ImplicitPatch.h"
#include "Curves/Contour.h"
//#include "Utils/stb_image.h"
#include "Utils/Utils.h"

#include "TerrainGen/VoxelGrid.h"
#include "TerrainGen/Heightmap.h"
#include "TerrainGen/LayerBasedGrid.h"


#include "serialization/Serializer.h"
#include <memory>
#include <filesystem>
#include <random>
#include <atomic>
#include <exception>
#include <mutex>
#include <QSaveFile>

float ImplicitPatch::isovalue = .5f;
const float ImplicitPatch::zResolution = .1f;
int ImplicitPatch::currentMaxIndex = -1;
std::string ImplicitPatch::json_identifier = "implicitPatches";


ImplicitPatch::ImplicitPatch()
{
    this->setIndex();

    this->_cachedMaxHeight = GridF(1, 1, 1, 0.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;
    this->_cachedMinHeight = GridF(1, 1, 1, 0.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

namespace {
thread_local bool snapshotSerialization = false;
struct SnapshotScope {
    bool previous = snapshotSerialization;
    SnapshotScope() { snapshotSerialization = true; }
    ~SnapshotScope() { snapshotSerialization = previous; }
};
struct FieldHash {
    size_t value = 0;
    template<class T> void add(const T& input) {
        const size_t h = std::hash<T>{}(input);
        value ^= h + size_t(0x9e3779b9) + (value << 6) + (value >> 2);
    }
    void vector(const Vector3& v) { add(v.isValid()); add(v.x()); add(v.y()); add(v.z()); }
};
size_t primitiveShapeSignature(const ImplicitPrimitive& primitive)
{
    FieldHash hash;
    hash.add(int(primitive.predefinedShape));
    hash.add(dynamic_cast<const Implicit2DPrimitive*>(&primitive) != nullptr);
    hash.vector(primitive.dimensions);
    for (float parameter : primitive.parametersProvided) hash.add(parameter);
    for (const auto& point : primitive.optionalCurve) hash.vector(point);
    hash.add(primitive.sampledVolume);
    hash.add(primitive.heightmapFilename);
    hash.vector(primitive.cachedHeightmap.getDimensions());
    for (float sample : primitive.cachedHeightmap) hash.add(sample);
    return hash.value;
}
void prepareImplicitField(ImplicitPatch* node)
{
    if (auto* primitive = dynamic_cast<ImplicitPrimitive*>(node)) {
        const bool generated = primitive->evaluationFunctionType() == primitive->_generatedFunctionType;
        if (generated && primitive->predefinedShape != ImplicitPatch::None &&
            primitiveShapeSignature(*primitive) != primitive->_generatedShapeSignature) primitive->update();
    }
    if (auto* group = dynamic_cast<ImplicitNaryOperator*>(node))
        for (auto* child : group->composables) if (child) prepareImplicitField(child);
}
size_t implicitFieldSignature(const ImplicitPatch& node, bool& cacheable)
{
    FieldHash hash;
    hash.add(typeid(node).hash_code());
    hash.add(ImplicitPatch::isovalue);
    for (auto material : LayerBasedGrid::invisibleLayers) hash.add(int(material));
    hash.add(node.mirrored);
    hash.vector(node.scaling); hash.vector(node.translation);
    if (const auto* primitive = dynamic_cast<const ImplicitPrimitive*>(&node)) {
        hash.add(primitiveShapeSignature(*primitive));
        hash.vector(primitive->position); hash.vector(primitive->supportDimensions);
        hash.add(int(primitive->material));
        hash.vector(primitive->sampledMaterialColumns.getDimensions());
        for (const auto& column : primitive->sampledMaterialColumns)
            for (const auto& [material, thickness] : column) { hash.add(int(material)); hash.add(thickness); }
        hash.add(primitive->evaluationFunctionType());
        if (primitive->evaluationFunctionType() != primitive->_generatedFunctionType) cacheable = false;
    }
    if (const auto* binary = dynamic_cast<const ImplicitBinaryOperator*>(&node)) {
        hash.add(int(binary->composeFunction)); hash.add(int(binary->positionalB));
        hash.add(binary->blendingFactor); hash.add(binary->withIntersectionOnB);
    }
    if (const auto* unary = dynamic_cast<const ImplicitUnaryOperator*>(&node)) {
        hash.add(unary->getTransformCount()); hash.vector(unary->_noise);
        hash.vector(unary->_translation); hash.vector(unary->_rotation); hash.vector(unary->_scale);
        hash.add(unary->noiseFunctionType());
        if (!unary->hasOnlySerializableTransforms() || unary->noiseFunctionType() != unary->_generatedNoiseType)
            cacheable = false;
    }
    if (const auto* group = dynamic_cast<const ImplicitNaryOperator*>(&node)) {
        hash.add(group->composables.size());
        for (const auto* child : group->composables) {
            hash.add(reinterpret_cast<size_t>(child));
            if (child) hash.add(implicitFieldSignature(*child, cacheable));
        }
    }
    return hash.value;
}
void readImplicitTerrainFrame(ImplicitPatch& patch, const nlohmann::json& content)
{
    const Vector3 scale = content.value("terrainScaling", patch.scaling);
    const Vector3 translation = content.value("terrainTranslation", patch.translation);
    auto finite = [](const Vector3& value) {
        return value.isValid() && std::isfinite(value.x()) && std::isfinite(value.y()) && std::isfinite(value.z());
    };
    if (!finite(scale) || scale.x() == 0.f || scale.y() == 0.f || scale.z() == 0.f || !finite(translation))
        throw std::invalid_argument("Implicit terrain frame requires finite translation and finite nonzero scale");
    patch.scaling = scale;
    patch.translation = translation;
}
float implicitColumnBoundary(const ImplicitPatch& patch, const Vector3& position, bool top)
{
    const auto bounds = patch.getSupportBBox();
    if (!bounds.contains(Vector3(position.x(), position.y(), bounds.center().z()))) return 0.f;
    const float start = top ? bounds.max().z() : bounds.min().z();
    const float finish = top ? bounds.min().z() : bounds.max().z();
    if (!std::isfinite(start) || !std::isfinite(finish) || finish == start) return 0.f;
    auto inside = [&](float z) { return patch.getSignedEvaluation(Vector3(position.x(), position.y(), z)) > 0.f; };
    if (inside(start)) return start;
    float previous = start;
    const double count = std::ceil(std::abs(double(finish) - start) / ImplicitPatch::zResolution);
    for (double i = 1; i <= count; ++i) {
        const float z = top ? std::max(finish, float(start - i * ImplicitPatch::zResolution))
                            : std::min(finish, float(start + i * ImplicitPatch::zResolution));
        if (inside(z)) {
            float outsideZ = previous, insideZ = z;
            for (int refinement = 0; refinement < 18; ++refinement) {
                const float middle = outsideZ + (insideZ - outsideZ) * .5f;
                if (inside(middle)) insideZ = middle;
                else outsideZ = middle;
            }
            return outsideZ + (insideZ - outsideZ) * .5f;
        }
        previous = z;
    }
    return 0.f; // Existing empty-column convention; never a cache-validity sentinel.
}
}

float ImplicitPatch::getMaxHeight(const Vector3& pos)
{
    // Positioning is invoked by parallel Z samples. Do not lazily write a
    // shared integer-XY cache from these continuous-coordinate queries.
    return implicitColumnBoundary(*this, pos, true);
}

float ImplicitPatch::getMinHeight(const Vector3& pos)
{
    return implicitColumnBoundary(*this, pos, false);
}

float ImplicitPatch::getMinimalHeight(AABBox BBox)
{
    return this->getMinimalHeight(BBox.min(), BBox.max());
}

float ImplicitPatch::getMaximalHeight(AABBox BBox)
{
    return this->getMaximalHeight(BBox.min(), BBox.max());
}

float ImplicitPatch::getMinimalHeight(const Vector3& minBox, const Vector3& maxBox)
{
    auto BBox = this->getBBox();
    if (minBox.isValid()) BBox.mini = Vector3::max(BBox.min(), minBox);
    if (maxBox.isValid()) BBox.maxi = Vector3::min(BBox.max(), maxBox);
//    if (!minBox.isValid()) minBox = BBox.min();
//    if (!maxBox.isValid()) maxBox = BBox.max();

    float minHeight = std::numeric_limits<float>::infinity();
    for (double x = std::ceil(BBox.min().x()); x < BBox.max().x(); x += 1.) {
        for (double y = std::ceil(BBox.min().y()); y < BBox.max().y(); y += 1.) {
            minHeight = std::min(minHeight, this->getMinHeight(Vector3(x, y)));
        }
    }
    return std::isfinite(minHeight) ? minHeight : 0.f;
}

float ImplicitPatch::getMaximalHeight(const Vector3& minBox, const Vector3& maxBox)
{
    auto BBox = this->getBBox();
    if (minBox.isValid()) BBox.mini = Vector3::max(BBox.min(), minBox);
    if (maxBox.isValid()) BBox.maxi = Vector3::min(BBox.max(), maxBox);
//    if (!minBox.isValid()) minBox = BBox.min();
//    if (!maxBox.isValid()) maxBox = BBox.max();

    float maxHeight = -std::numeric_limits<float>::infinity();
    for (double x = std::ceil(BBox.min().x()); x < BBox.max().x(); x += 1.) {
        for (double y = std::ceil(BBox.min().y()); y < BBox.max().y(); y += 1.) {
            maxHeight = std::max(maxHeight, this->getMaxHeight(Vector3(x, y)));
        }
    }
    return std::isfinite(maxHeight) ? maxHeight : 0.f;
}

Vector3 ImplicitPatch::getSamplingDimensions() const
{
    if (typeid(*this) == typeid(ImplicitNaryOperator)) {
        const auto* group = static_cast<const ImplicitNaryOperator*>(this);
        if (group->composables.size() == 1 && group->composables.front()) return group->composables.front()->getSamplingDimensions();
        if (group->composables.empty()) return Vector3(1, 1, 2);
    }
    if (const auto* primitive = dynamic_cast<const ImplicitPrimitive*>(this)) {
        if (primitive->predefinedShape == ImplicitHeightmap && primitive->sampledVolume)
            return Vector3::max((primitive->position + primitive->cachedHeightmap.getDimensions()).roundedUp(), Vector3(1, 1, 1));
    }
    return Vector3::max(getBBox().max().roundedUp() + Vector3(0, 0, 2), Vector3(1, 1, 2));
}

float ImplicitPatch::getSignedEvaluation(const Vector3& pos) const
{
    if (const auto* heightTerrain = dynamic_cast<const Implicit2DNary*>(this))
        return heightTerrain->computeHeight(pos) - pos.z();
    if (const auto* heightPrimitive = dynamic_cast<const Implicit2DPrimitive*>(this)) {
        if (!LayerBasedGrid::invisibleLayers.count(heightPrimitive->material)) {
            const float height = heightPrimitive->evaluate(pos.xy());
            return height > 0.f ? height - pos.z() : -1.f;
        }
    }
    if (const auto* unary = dynamic_cast<const ImplicitUnaryOperator*>(this)) {
        if (unary->composableA() && unary->_noise.x() == 0.f && unary->noiseAt(pos) == 0.f)
            return unary->composableA()->getSignedEvaluation(unary->inverseTransform(pos));
    }
    // A one-child container must not threshold away the child's exterior field.
    if (typeid(*this) == typeid(ImplicitNaryOperator)) {
        const auto* group = static_cast<const ImplicitNaryOperator*>(this);
        if (group->composables.size() == 1 && group->composables.front()) return group->composables.front()->getSignedEvaluation(pos);
    }
    if (const auto* primitive = dynamic_cast<const ImplicitPrimitive*>(this)) {
        if (primitive->hasGeneratedField() && primitive->predefinedShape == ImplicitHeightmap && !primitive->cachedHeightmap.empty() &&
            (!primitive->sampledMaterialColumns.empty() || LayerBasedGrid::invisibleLayers.count(primitive->material) == 0)) {
            const Vector3 local = pos - primitive->position;
            if (primitive->sampledVolume || primitive->cachedHeightmap.sizeZ > 1) {
                const float value = primitive->cachedHeightmap.interpolate(local);
                if (value > 0.f && !primitive->sampledMaterialColumns.empty()) {
                    const auto materials = primitive->getMaterials(pos);
                    if (!materials.empty() && LayerBasedGrid::invisibleLayers.count(materials.begin()->first)) return -value;
                }
                return value;
            }
            if (local.x() < 0.f || local.y() < 0.f || local.x() > primitive->dimensions.x() - 1.f || local.y() > primitive->dimensions.y() - 1.f)
                return -1.f;
            return TerrainSampling::heightAt(primitive->cachedHeightmap, local.x(), local.y()) - local.z();
        }
    }
    float value = 0.f;
    for (const auto& [material, amount] : getMaterials(pos))
        value += (LayerBasedGrid::invisibleLayers.count(material) != 0) ? -amount : amount;
    return value - isovalue;
}

std::pair<float, std::map<TerrainTypes, float> > ImplicitPatch::getMaterialsAndTotalEvaluation(const Vector3& pos) const
{
    return this->getBinaryMaterialsAndTotalEvaluation(pos);
    /*std::map<TerrainTypes, float> materials = this->getMaterials(pos);
    float totalValue = 0.f;

    for (const auto& [mat, val] : materials)
        totalValue += val;

    return {totalValue, materials};*/
}

std::pair<float, std::map<TerrainTypes, float> > ImplicitPatch::getBinaryMaterialsAndTotalEvaluation(const Vector3& pos) const
{
    std::map<TerrainTypes, float> materials = this->getMaterials(pos);
    float totalValue = 0.f;

    for (const auto& [mat, val] : materials)
        totalValue += val * ((LayerBasedGrid::invisibleLayers.count(mat) != 0) ? -1.f : 1.f);

    return {std::abs(totalValue), std::move(materials)};
}

Vector3 ImplicitPatch::getDimensions() const
{
    auto AABBox = this->getBBox();
    return AABBox.dimensions();
}

Vector3 ImplicitPatch::getSupportDimensions() const
{
    auto AABBox = this->getSupportBBox();
    return AABBox.dimensions();
}

GridV3 ImplicitPatch::getNormals()
{
    GridV3 normals = -this->getVoxelized().gradient();
    for (auto& n : normals)
        n.normalize();
    normals.raiseErrorOnBadCoord = false;
    normals.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    return normals;
}

Vector3 ImplicitPatch::getNormal(const Vector3& pos) const
{
    float delta = 1e-3; // Not too big, not too small...

    float minusX = this->getSignedEvaluation(pos - Vector3(delta, 0, 0));
    float plusX = this->getSignedEvaluation(pos + Vector3(delta, 0, 0));
    float minusY = this->getSignedEvaluation(pos - Vector3(0, delta, 0));
    float plusY = this->getSignedEvaluation(pos + Vector3(0, delta, 0));
    float minusZ = this->getSignedEvaluation(pos - Vector3(0, 0, delta));
    float plusZ = this->getSignedEvaluation(pos + Vector3(0, 0, delta));

    return Vector3(minusX - plusX, minusY - plusY, minusZ - plusZ).normalize();
}

void ImplicitPatch::setIndex(int newIndex)
{
    static std::mutex indexMutex;
    std::lock_guard<std::mutex> lock(indexMutex);
    if (newIndex < -1) throw std::invalid_argument("Implicit node index must be nonnegative or automatically assigned");
    if (newIndex == -1) {
        if (currentMaxIndex == std::numeric_limits<int>::max()) throw std::overflow_error("Implicit node indices exhausted");
        newIndex = currentMaxIndex + 1;
    }
    index = newIndex;
    currentMaxIndex = std::max(index, currentMaxIndex);
}

ImplicitPatch *ImplicitPatch::fromJson(nlohmann::json content)
{
    namespace fs = std::filesystem;
    static thread_local unsigned depth = 0;
    static thread_local std::vector<fs::path> includes;
    if (depth >= 256) throw std::invalid_argument("Implicit JSON nesting exceeds 256 nodes");
    struct DepthScope { unsigned& value; explicit DepthScope(unsigned& v) : value(v) { ++value; } ~DepthScope() { --value; } } scope(depth);
    if (!content.is_object()) throw std::invalid_argument("Implicit JSON node must be an object");
    if (content.contains(json_identifier)) {
        if (content.value("schemaVersion", 1) != 1) throw std::invalid_argument("Unsupported implicit JSON schema version");
        return fromJson(content.at(json_identifier));
    }
    std::unique_ptr<ImplicitPatch> result;
    if (content.contains("type")) {
        const std::string type = toUpper(content.at("type"));
        if (type == "PRIMITIVE") result.reset(ImplicitPrimitive::fromJson(content));
        else if (type == "COMPOSE" || type == "BINARY") result.reset(ImplicitBinaryOperator::fromJson(content));
        else if (type == "UNARY") result.reset(ImplicitUnaryOperator::fromJson(content));
        else if (type == "NARY") result.reset(ImplicitNaryOperator::fromJson(content));
        else throw std::invalid_argument("Unknown implicit JSON node type: " + type);
    } else if (content.contains("file")) {
        fs::path path(content.at("file").get<std::string>());
        if (path.empty()) throw std::invalid_argument("Empty implicit include filename");
        if (path.is_relative() && !includes.empty()) {
            const auto relative = includes.back().parent_path() / path;
            if (fs::exists(relative)) path = relative;
        }
        if (!fs::exists(path)) throw std::runtime_error("Implicit file not found: " + path.string());
        path = fs::canonical(path);
        if (std::find(includes.begin(), includes.end(), path) != includes.end())
            throw std::invalid_argument("Recursive implicit file include: " + path.string());
        includes.push_back(path);
        struct IncludeScope { std::vector<fs::path>& paths; ~IncludeScope() { paths.pop_back(); } } includeScope{includes};
        const bool jsonFile = toUpper(path.extension().string()) == ".JSON";
        if (jsonFile) {
            std::ifstream input(path);
            if (!input) throw std::runtime_error("Cannot read implicit file: " + path.string());
            result.reset(fromJson(nlohmann::json::parse(input)));
        } else {
            result.reset(ImplicitPrimitive::fromHeightmap(path.string(), content.value("dimensions", Vector3::invalid)));
        }
        if (content.contains("material")) {
            const auto material = materialFromString(content.at("material"));
            std::vector<ImplicitPatch*> nodes{result.get()};
            while (!nodes.empty()) {
                auto* node = nodes.back(); nodes.pop_back();
                if (auto* primitive = dynamic_cast<ImplicitPrimitive*>(node)) {
                    primitive->material = material;
                    primitive->sampledMaterialColumns = ImplicitPrimitive::MaterialColumns();
                }
                if (auto* group = dynamic_cast<ImplicitNaryOperator*>(node))
                    for (auto* child : group->composables) if (child) nodes.push_back(child);
            }
        }
        const auto originalBounds = result->getBBox();
        const Vector3 desiredPosition = content.value("position", originalBounds.min());
        const Vector3 desiredDimensions = content.value("dimensions", originalBounds.dimensions());
        if (jsonFile && desiredDimensions != originalBounds.dimensions()) {
            if (originalBounds.dimensions().minComp() <= 0.f || desiredDimensions.minComp() <= 0.f)
                throw std::invalid_argument("Cannot resize a degenerate implicit include");
            auto wrapper = std::make_unique<ImplicitUnaryOperator>();
            wrapper->scaling = result->scaling; wrapper->translation = result->translation;
            wrapper->addChild(result.release());
            wrapper->scale(desiredDimensions / originalBounds.dimensions());
            result = std::move(wrapper);
        }
        if (desiredPosition != result->getBBox().min()) {
            if (auto* primitive = dynamic_cast<ImplicitPrimitive*>(result.get())) primitive->position += desiredPosition - result->getBBox().min();
            else {
                const Vector3 shift = desiredPosition - result->getBBox().min();
                auto wrapper = std::make_unique<ImplicitTranslation>();
                wrapper->scaling = result->scaling; wrapper->translation = result->translation;
                wrapper->addChild(result.release());
                wrapper->translate(shift);
                result = std::move(wrapper);
            }
        }
        result->used_json_filename = path.string();
        result->referenceProperties = content;
    } else throw std::invalid_argument("Implicit JSON node requires a type or file");
    readImplicitTerrainFrame(*result, content);
    result->updateCache();
    return result.release();
}

void ImplicitPatch::saveMap(const std::string& filename)
{
    const std::string payload = nlohmann::json({{"schemaVersion", 1}, {json_identifier, toJson()}}).dump(1, '\t');
    QSaveFile file(QString::fromStdString(filename));
    if (!file.open(QIODevice::WriteOnly) || file.write(payload.data(), qint64(payload.size())) != qint64(payload.size()) || !file.commit())
        throw std::runtime_error("Cannot atomically save implicit terrain: " + filename);
}

void ImplicitPatch::retrieveMap(const std::string& filename)
{
    std::unique_ptr<ImplicitPatch> decoded(fromJson({{"file", filename}}));
    if (typeid(*this) == typeid(ImplicitNaryOperator) || typeid(*this) == typeid(Implicit2DNary)) {
        auto* target = static_cast<ImplicitNaryOperator*>(this);
        if (typeid(*decoded) == typeid(*this)) *target = *static_cast<ImplicitNaryOperator*>(decoded.get());
        else {
            ImplicitNaryOperator replacement;
            replacement.addChild(decoded.release());
            *target = replacement;
        }
    } else if (auto* target = dynamic_cast<ImplicitPrimitive*>(this)) {
        auto* source = dynamic_cast<ImplicitPrimitive*>(decoded.get());
        if (!source || typeid(*target) != typeid(*source)) throw std::invalid_argument("Implicit load cannot change the receiver's primitive type");
        auto* oldParent = parent;
        const int oldIndex = index;
        *target = *source;
        parent = oldParent;
        index = oldIndex;
        target->update();
    } else throw std::invalid_argument("Load transformed terrain through an implicit root container");
    // Loading a whole map creates an editable snapshot rather than an include
    // that would serialize back to the map's own filename.
    used_json_filename.clear();
    referenceProperties = nlohmann::json::object();
    updateCache();
}

void ImplicitPatch::updateCache()
{
    this->_cached = false;
    Vector3 dim = Vector3::max(Vector3(1, 1), this->getDimensions().roundedUp());

    this->_cachedMaxHeight = GridF(dim.x(), dim.y(), 1, -1.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;

    this->_cachedMinHeight = GridF(dim.x(), dim.y(), 1, -1.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

bool ImplicitPatch::checkIsInGround(const Vector3& position)
{
    return getSignedEvaluation(position) > 0.f;
}

Mesh ImplicitPatch::getGeometry(const Vector3& dimensions)
{
    const Vector3 extent = dimensions.isValid() ? dimensions : getSamplingDimensions();
    VoxelGrid voxels;
    voxels.setVoxelValues(getVoxelized(extent));
    return voxels.getGeometry();
}

Vector3 ImplicitPatch::getIntersection(const Vector3& origin, const Vector3& dir, const Vector3& minPos, const Vector3& maxPos)
{
    // Density is not a distance bound. Walk a bounded ray interval, then refine
    // the first detected air/solid transition using the rendering field.
    auto finite = [](const Vector3& v) {
        return v.isValid() && std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
    };
    if (!finite(origin) || !finite(dir) || dir.norm2() == 0.f) return Vector3::invalid;
    const Vector3 direction = dir.normalized();
    auto bounds = getSupportBBox();
    if (!finite(bounds.min()) || !finite(bounds.max())) return Vector3::invalid;
    if (minPos.isValid()) bounds.mini = Vector3::max(bounds.min(), minPos);
    if (maxPos.isValid()) bounds.maxi = Vector3::min(bounds.max(), maxPos);
    if ((bounds.max() - bounds.min()).minComp() < 0.f) return Vector3::invalid;
    float entry = 0.f, exit = std::numeric_limits<float>::infinity();
    for (int axis = 0; axis < 3; ++axis) {
        if (direction[axis] == 0.f) {
            if (origin[axis] < bounds.min()[axis] || origin[axis] > bounds.max()[axis])
                return Vector3::invalid;
            continue;
        }
        float first = (bounds.min()[axis] - origin[axis]) / direction[axis];
        float last = (bounds.max()[axis] - origin[axis]) / direction[axis];
        if (first > last) std::swap(first, last);
        entry = std::max(entry, first);
        exit = std::min(exit, last);
        if (entry > exit) return Vector3::invalid;
    }
    if (!std::isfinite(exit)) return Vector3::invalid;
    float previous = entry;
    if (getSignedEvaluation(origin + direction * entry) > 0.f)
        return origin + direction * entry;
    // Use double loop counters so floating-point addition cannot stall at
    // large coordinates. Features thinner than zResolution remain sampled.
    const double count = std::ceil((double(exit) - entry) / zResolution);
    for (double i = 1; i <= count; ++i) {
        const float t = float(std::min(double(exit), double(entry) + i * zResolution));
        if (getSignedEvaluation(origin + direction * t) > 0.f) {
            float lower = previous, upper = t;
            for (int refinement = 0; refinement < 20; ++refinement) {
                const float middle = lower + (upper - lower) * .5f;
                if (getSignedEvaluation(origin + direction * middle) > 0.f) upper = middle;
                else lower = middle;
            }
            return origin + direction * (lower + (upper - lower) * .5f);
        }
        previous = t;
    }
    return Vector3::invalid;
}


float ImplicitPatch::getHeight(float x, float y)
{
    return this->getHeight(Vector3(x, y, 0));
}

float ImplicitPatch::getHeight(const Vector3& pos)
{
    return this->getMaxHeight(pos);
}

bool ImplicitPatch::contains(const Vector3& v)
{
    auto BBox = this->getSupportBBox();
    return BBox.contains(v);
}

bool ImplicitPatch::contains(float x, float y, float z)
{
    return this->contains(Vector3(x, y, z));
}

Vector3 ImplicitPatch::getGlobalPositionOf(const Vector3& posInsidePatch)
{
    Vector3 position = posInsidePatch + getBBox().min();
    for (ImplicitPatch* child = this; child->parent; child = child->parent) {
        auto* ancestor = child->parent;
        if (auto* unary = dynamic_cast<ImplicitUnaryOperator*>(ancestor)) position = unary->applyTransform(position);
        else if (auto* binary = dynamic_cast<ImplicitBinaryOperator*>(ancestor)) {
            if (binary->composableB() == child) {
                const float zero = binary->getEvaluationPositionForComposableB(position.xy()).z();
                const float one = binary->getEvaluationPositionForComposableB(position.xy() + Vector3(0, 0, 1)).z();
                if (!std::isfinite(one - zero) || std::abs(one - zero) < 1e-6f)
                    throw std::domain_error("Noninvertible implicit positional operation");
                position.z() = (position.z() - zero) / (one - zero);
            }
        }
    }
    return position;
}

Vector3 ImplicitPatch::getLocalPositionOf(const Vector3& globalPosition)
{
    Vector3 position = globalPosition;
    std::vector<std::pair<ImplicitNaryOperator*, ImplicitPatch*>> ancestry;
    for (ImplicitPatch* child = this; child->parent; child = child->parent) ancestry.emplace_back(child->parent, child);
    for (auto it = ancestry.rbegin(); it != ancestry.rend(); ++it) {
        if (auto* unary = dynamic_cast<ImplicitUnaryOperator*>(it->first)) position = unary->inverseTransform(position);
        else if (auto* binary = dynamic_cast<ImplicitBinaryOperator*>(it->first))
            if (binary->composableB() == it->second) position = binary->getEvaluationPositionForComposableB(position);
    }
    return position - getBBox().min();
}

size_t ImplicitPatch::getCurrentHistoryIndex() const
{
    bool cacheable = true;
    const size_t signature = implicitFieldSignature(*this, cacheable);
    static std::atomic<size_t> customRevision{0};
    return cacheable ? signature : ++customRevision;
}

GridF ImplicitPatch::getVoxelized(const Vector3& dimensions, const Vector3& scale)
{
    prepareImplicitField(this);
    if (!std::isfinite(isovalue)) throw std::invalid_argument("Implicit isovalue must be finite");
    const Vector3 extent = dimensions.isValid() ? dimensions : getSamplingDimensions();
    const auto counts = TerrainSampling::sampleDimensions(extent, scale);
    bool cacheable = true;
    const size_t signature = implicitFieldSignature(*this, cacheable);
    if (_cached && cacheable && _cachedFieldSignature == signature &&
        _cachedVoxelized.getDimensions() == counts && _cachedVoxelScale == scale) return _cachedVoxelized;
    _cached = false;
    _cachedVoxelized = GridF(counts, -1.f);
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    _cachedVoxelized.iterateParallel([&](int x, int y, int z) {
        if (failed.load(std::memory_order_relaxed)) return;
        try {
            const float value = getSignedEvaluation(Vector3(x, y, z) / scale);
            if (!std::isfinite(value)) throw std::domain_error("Implicit field returned a nonfinite sample");
            _cachedVoxelized(x, y, z) = value;
        } catch (...) {
            if (!failed.exchange(true)) error = std::current_exception();
        }
    });
    if (error) std::rethrow_exception(error);
    _cachedVoxelized.defaultValueOnBadCoord = -1.f;
    _cachedVoxelScale = scale;
    _cachedFieldSignature = signature;
    _cached = cacheable;
    return _cachedVoxelized;
}

ImplicitPrimitive *ImplicitPatch::createPredefinedShape(PredefinedShapes shape, const Vector3& dimensions, float additionalParam, CatmullRomSpline parametricCurve, bool in2D)
{
    std::unique_ptr<ImplicitPrimitive> primitive;
    if (in2D)
        primitive.reset(new Implicit2DPrimitive());
    else
        primitive.reset(new ImplicitPrimitive());
    primitive->predefinedShape = shape;
    primitive->evalFunction = ImplicitPatch::createPredefinedShapeFunction(shape, dimensions, additionalParam, parametricCurve, in2D);
    primitive->optionalCurve = parametricCurve;
    primitive->parametersProvided = {additionalParam};
    primitive->dimensions = dimensions;
    primitive->_generatedFunctionType = primitive->evaluationFunctionType();
    primitive->_generatedShapeSignature = primitiveShapeSignature(*primitive);
    return primitive.release();
}
ImplicitPatch *ImplicitPatch::createIdentity()
{
    ImplicitPrimitive* primitive = new ImplicitPrimitive();
    primitive->predefinedShape = ImplicitPatch::PredefinedShapes::None;
    primitive->evalFunction = [](const Vector3&) { return 0.f; };
    primitive->optionalCurve = {};
    primitive->parametersProvided = {0.f};
    primitive->_generatedFunctionType = primitive->evaluationFunctionType();
    primitive->_generatedShapeSignature = primitiveShapeSignature(*primitive);
    return primitive;
}

std::function<float (const Vector3&)> ImplicitPatch::createPredefinedShapeFunction(PredefinedShapes shape, const Vector3& dimensions, float additionalParam, CatmullRomSpline parametricCurve, bool in2D)
{
    if (shape == None) return [](const Vector3&) { return 0.f; };
    if (shape == ImplicitHeightmap) throw std::invalid_argument("Use fromHeightmap or fromVoxelField for sampled primitives");
    if (!dimensions.isValid() || !std::isfinite(dimensions.x()) || !std::isfinite(dimensions.y()) ||
        !std::isfinite(dimensions.z()) || dimensions.minComp() < 0.f || !std::isfinite(additionalParam))
        throw std::invalid_argument("Implicit shapes require finite nonnegative dimensions and finite parameters");
    if (dimensions.minComp() == 0.f) return [](const Vector3&) { return 0.f; };
    if ((shape == Gaussian || shape == ParametricTunnel || shape == MountainChain) && additionalParam <= 0.f)
        return [](const Vector3&) { return 0.f; };
    if (shape == ParametricTunnel || shape == MountainChain || shape == Polygon || shape == DistanceMap) {
        const size_t required = (shape == Polygon || shape == DistanceMap) ? 3 : 2;
        if (parametricCurve.size() < required) throw std::invalid_argument("Implicit curve has too few control points");
        for (const auto& point : parametricCurve)
            if (!point.isValid() || !std::isfinite(point.x()) || !std::isfinite(point.y()) || !std::isfinite(point.z()))
                throw std::invalid_argument("Implicit curve has a nonfinite control point");
    }
    std::function<float(const Vector3&)> func;
    switch(shape) {
    case Sphere:
        func = ImplicitPatch::createSphereFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Block:
        func = ImplicitPatch::createBlockFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Gaussian:
        func = ImplicitPatch::createGaussianFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Rock:
        func = ImplicitPatch::createRockFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Mountain:
        func = ImplicitPatch::createMountainFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Dune:
        func = ImplicitPatch::createDuneFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Basin:
        func = ImplicitPatch::createBasinFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Cave:
        func = ImplicitPatch::createCaveFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Arch:
        func = ImplicitPatch::createArchFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Noise2D:
        func = ImplicitPatch::createNoise2DFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Cylinder:
        func = ImplicitPatch::createCylinderFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case MountainChain:
        func = ImplicitPatch::createMountainChainFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), parametricCurve, in2D);
        break;
    case Polygon:
        func = ImplicitPatch::createPolygonFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), parametricCurve, in2D);
        break;
    case DistanceMap:
        func = ImplicitPatch::createDistanceMapFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), parametricCurve, in2D);
        break;
    case ImplicitHeightmap:
        throw std::invalid_argument("Sampled primitive requires sample data");
    case ParametricTunnel:
        func = ImplicitPatch::createParametricTunnelFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), parametricCurve, in2D);
        break;
    case None:
        func = ImplicitPatch::createIdentityFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    case Ripple:
        func = ImplicitPatch::createRippleFunction(additionalParam, dimensions.x(), dimensions.y(), dimensions.z(), in2D);
        break;
    }
    if (!func) throw std::invalid_argument("Unknown implicit shape enumeration");
    return func;
}




void ImplicitPrimitive::setEvaluationFunction(std::function<float(const Vector3&)> function)
{
    if (!function) throw std::invalid_argument("Empty implicit evaluation callback");
    evalFunction = std::move(function);
    _generatedFunctionType = 0;
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;
}

ImplicitPrimitive::ImplicitPrimitive()
{
    this->evalFunction = [](const Vector3&) -> float { return 0.f; };
    this->parametersProvided = {-1.f}; // Legacy identity parameter.
    position = Vector3();
    dimensions = Vector3(1, 1, 1);
    _generatedFunctionType = evalFunction.target_type().hash_code();
    _generatedShapeSignature = primitiveShapeSignature(*this);
}

float ImplicitPrimitive::evaluate(const Vector3& pos) const
{
    auto [minPos, maxPos] = this->getBBox();
    auto [minSupportPos, maxSupportPos] = this->getSupportBBox();
    if (!Vector3::isInBox(pos, minSupportPos, maxSupportPos))
        return 0.f;
    Vector3 supportWidth = maxSupportPos - minSupportPos;
    Vector3 width = maxPos - minPos;

    float evaluation = 0.f;

    Vector3 evalPos = pos - this->position;
    if (Vector3::isInBox(pos, minPos, maxPos)) {
        evaluation = this->evalFunction(evalPos);
    } else {
        if (!mirrored) {
            evaluation = this->evalFunction(evalPos);
        } else {
            float distanceToBBox = Vector3::signedDistanceToBoundaries(pos, minPos, maxPos);
            float distanceFactor = std::max(1.f, (supportWidth - width).maxComp() * .5f); // Don't get a 0 here
            distanceToBBox /= distanceFactor; // Make it depending on the support area
            float distanceFalloff = interpolation::wyvill(std::clamp(distanceToBBox, 0.f, 1.f));

            evalPos = evalPos.abs();
            evalPos.x() = (evalPos.x() < width.x() ? evalPos.x() : 2.f * width.x() - evalPos.x());
            evalPos.y() = (evalPos.y() < width.y() ? evalPos.y() : 2.f * width.y() - evalPos.y());
            evalPos.z() = (evalPos.z() < width.z() ? evalPos.z() : 2.f * width.z() - evalPos.z());
            evaluation = std::clamp(this->evalFunction(evalPos), 0.f, 1.f) * distanceFalloff * .5f;
        }
    }

    if (!std::isfinite(evaluation)) return 0.f;
    evaluation = std::clamp(evaluation, 0.f, 1.f);
    return evaluation;
}

std::map<TerrainTypes, float> ImplicitPrimitive::getMaterials(const Vector3& pos) const
{
    TerrainTypes selected = material;
    if (!sampledMaterialColumns.empty()) {
        const Vector3 local = pos - position;
        selected = AIR;
        if (local.x() >= 0.f && local.y() >= 0.f && local.z() >= 0.f &&
            local.x() < sampledMaterialColumns.sizeX && local.y() < sampledMaterialColumns.sizeY) {
            float top = 0.f;
            for (const auto& [label, thickness] : sampledMaterialColumns(int(local.x()), int(local.y()))) {
                top += thickness;
                if (local.z() < top) { selected = label; break; }
            }
        }
    }
    return {{selected, evaluate(pos)}};
}

AABBox ImplicitPrimitive::getSupportBBox() const
{
    Vector3 margin = (this->supportDimensions - this->dimensions) * .5f;
    if (!this->supportDimensions.isValid() || this->supportDimensions == Vector3())
        margin = this->dimensions * 1.f; // Nice big margin
    return {this->position - margin, this->position + this->dimensions + margin};
}

AABBox ImplicitPrimitive::getBBox() const
{
    return {this->position, this->position + this->dimensions};
}

void ImplicitPrimitive::update()
{
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;


    if (this->predefinedShape == ImplicitHeightmap) {
        for (float value : cachedHeightmap) if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite implicit sample");
        if (cachedHeightmap.empty() && heightmapFilename.empty()) throw std::invalid_argument("Sampled primitive requires samples or a filename");
        if (this->cachedHeightmap.size() == 0) {
            this->cachedHeightmap = Heightmap(this->heightmapFilename, this->dimensions.x(), this->dimensions.y()).heights;
        }
        if (sampledVolume || this->cachedHeightmap.sizeZ > 1) {
            float maxMagnitude = 1.f;
            for (float value : cachedHeightmap) maxMagnitude = std::max(maxMagnitude, std::abs(value));
            sampledFieldScale = .5f / maxMagnitude;
            cachedHeightmap.defaultValueOnBadCoord = -maxMagnitude;
            this->evalFunction = [this](const Vector3& pos) -> float {
                return isovalue + cachedHeightmap.interpolate(pos) * sampledFieldScale;
            };
        } else {
            const float range = std::max(1.f, cachedHeightmap.max() + 2.f);
            sampledFieldScale = .5f / range;
            this->evalFunction = [this](const Vector3& pos) -> float {
                return isovalue + (TerrainSampling::heightAt(cachedHeightmap, pos.x(), pos.y()) - pos.z()) * sampledFieldScale;
            };
        }
    } else {
        this->evalFunction = ImplicitPatch::createPredefinedShapeFunction(this->predefinedShape, this->dimensions, this->parametersProvided.at(0), this->optionalCurve, dynamic_cast<Implicit2DPrimitive*>(this) != nullptr);
    }
    _generatedFunctionType = evalFunction.target_type().hash_code();
    _generatedShapeSignature = primitiveShapeSignature(*this);
}

std::string ImplicitPrimitive::toString() const
{
    std::string name = this->name;
    if (name == "") {
        name = stringFromPredefinedShape(this->predefinedShape);
    }
    return name + ": Primitive " + name + " #" + std::to_string(this->index);
}

nlohmann::json ImplicitPrimitive::toJson()
{
    if (!hasGeneratedField() && !snapshotSerialization) throw std::runtime_error("Sample a custom evaluation callback before serializing it");
    nlohmann::json content;
    if (!this->used_json_filename.empty() && !snapshotSerialization) {
        content = referenceProperties;
        content["file"] = this->used_json_filename;
        content["position"] = this->position;
        content["dimensions"] = this->dimensions;
        content["material"] = stringFromMaterial(this->material);
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = patchOffset;
        content["scale"] = patchRescale); */
    } else {
        content["name"] = this->name;
        content["type"] = "primitive";
        content["heightfield"] = dynamic_cast<Implicit2DPrimitive*>(this) != nullptr;
        content["position"] = this->position;
        content["dimensions"] = this->dimensions;
        content["supportDimensions"] = this->supportDimensions;
        content["densityValue"] = LayerBasedGrid::densityFromMaterial(this->material);
        content["material"] = stringFromMaterial(this->material);
        content["index"] = this->index;
        content["shape"] = stringFromPredefinedShape(this->predefinedShape);
        content["sigmaValue"] = this->parametersProvided[0];
        content["optionalCurve"] = this->optionalCurve;
        content["mirrored"] = this->mirrored;
        content["filename"] = this->heightmapFilename;
        if (predefinedShape == ImplicitHeightmap && !cachedHeightmap.empty()) {
            content["sampledVolume"] = sampledVolume;
            content["samples"] = cachedHeightmap;
            if (!sampledMaterialColumns.empty()) {
                nlohmann::json columns = nlohmann::json::array();
                for (const auto& column : sampledMaterialColumns) {
                    nlohmann::json intervals = nlohmann::json::array();
                    for (const auto& [label, thickness] : column)
                        intervals.push_back({{"material", stringFromMaterial(label)}, {"thickness", thickness}});
                    columns.push_back(std::move(intervals));
                }
                content["materialColumns"] = {{"sizeX", sampledMaterialColumns.sizeX}, {"sizeY", sampledMaterialColumns.sizeY}, {"columns", columns}};
            }
        }
    }

    content["terrainScaling"] = scaling;
    content["terrainTranslation"] = translation;
    return content;
}

ImplicitPatch *ImplicitPrimitive::fromJson(nlohmann::json content)
{
    ImplicitPrimitive* patch = content.value("heightfield", false) ? new Implicit2DPrimitive : new ImplicitPrimitive;
    std::unique_ptr<ImplicitPrimitive> owner(patch);
    patch->name = content["name"];
    patch->position = content["position"];
    patch->setDimensions(content["dimensions"]);
    patch->supportDimensions = content["supportDimensions"];
    patch->material = LayerBasedGrid::materialFromDensity(content["densityValue"]);
    if (content.contains("material"))
        patch->material = materialFromString(content["material"]);
    if (content.contains("optionalCurve"))
        patch->optionalCurve = content["optionalCurve"];
    patch->setIndex(content["index"]); // Advance beyond IDs recorded by older sessions.
    patch->setIndex(); // This imported node is a new live instance, not an alias.
    patch->predefinedShape = predefinedShapeFromString(content["shape"]);
    patch->parametersProvided = {content["sigmaValue"]};
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    if (content.contains("filename"))
        patch->heightmapFilename = content["filename"];
    patch->sampledVolume = content.value("sampledVolume", false);
    if (content.contains("samples")) patch->cachedHeightmap = content["samples"].get<GridF>();
    if (content.contains("materialColumns")) {
        const auto& saved = content.at("materialColumns");
        const size_t nx = saved.at("sizeX"), ny = saved.at("sizeY");
        if (!patch->sampledVolume || patch->predefinedShape != ImplicitHeightmap ||
            nx != patch->cachedHeightmap.sizeX || ny != patch->cachedHeightmap.sizeY ||
            !saved.at("columns").is_array() || saved.at("columns").size() != nx * ny)
            throw std::invalid_argument("Material columns do not match the sampled field");
        patch->sampledMaterialColumns = MaterialColumns(nx, ny, 1);
        for (size_t i = 0; i < nx * ny; ++i) {
            float total = 0.f;
            if (!saved.at("columns").at(i).is_array()) throw std::invalid_argument("Material column must be an interval array");
            for (const auto& interval : saved.at("columns").at(i)) {
                const float thickness = interval.at("thickness");
                total += thickness;
                if (!std::isfinite(thickness) || thickness <= 0.f || !std::isfinite(total)) throw std::invalid_argument("Invalid material interval thickness");
                patch->sampledMaterialColumns.data[i].emplace_back(materialFromString(interval.at("material")), thickness);
            }
        }
    }
    patch->update();
    readImplicitTerrainFrame(*patch, content);
    return owner.release();
}

void ImplicitPrimitive::setDimensions(const Vector3& newDimensions)
{
    this->dimensions = newDimensions;
    if (predefinedShape != None && predefinedShape != ImplicitHeightmap) update();
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;
    this->updateCache();
}

void ImplicitPrimitive::setSupportDimensions(const Vector3& newSupportDimensions)
{
    this->supportDimensions = newSupportDimensions;
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;
}

std::vector<ImplicitPatch *> ImplicitPrimitive::findAll(PredefinedShapes shape)
{
    if (this->contains(shape)) {
        return {this};
    }
    return {};
}

ImplicitPatch *ImplicitPrimitive::copy() const
{
    ImplicitPrimitive* copy = dynamic_cast<const Implicit2DPrimitive*>(this)
        ? static_cast<ImplicitPrimitive*>(new Implicit2DPrimitive(*static_cast<const Implicit2DPrimitive*>(this)))
        : new ImplicitPrimitive(*this);
    copy->parent = nullptr;
    copy->setIndex();
    if (hasGeneratedField() && copy->predefinedShape != None) copy->update();
    return copy;
}

ImplicitPrimitive *ImplicitPrimitive::fromHeightmap(const std::string& filename, const Vector3& dimensions, ImplicitPrimitive *prim)
{
    if (!checkPathExists(filename)) {
        throw std::runtime_error("Error: Impossible to load '" + filename + "', file not found");
    }
    std::string ext = toUpper(getExtension(filename));
    if (isIn(ext, std::set<std::string>{"JPG", "PNG", "TGA", "BMP", "PSD", "GIF", "HDR", "PIC"})) {
        int imgW, imgH, nbChannels;
        unsigned char *data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
        if (data == NULL)
        {
            throw std::runtime_error("Error : Impossible to load heightmap at " + filename + "\n" +
                                    "Either file is not found, or type is incorrect. Available file types are : \n" +
                                    "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC");
        }
        GridF map(imgW, imgH);
        for (int x = 0; x < imgW; x++) {
            for (int y = 0; y < imgH; y++) {
                float value = (float)data[x + y * imgW];
                map.at(x, y) = value;
            }
        }
        stbi_image_free(data);

        if (dimensions.isValid()) {
            map = map.resize(dimensions.x(), dimensions.y(), 1);
            map = map.normalize() * dimensions.z();
        }
        return ImplicitPrimitive::fromHeightmap(map, filename,  prim);
    } else {
        VoxelGrid voxels = VoxelGrid();
        voxels.retrieveMap(filename);
        GridF map = voxels.getVoxelValues();
        auto* result = ImplicitPrimitive::fromVoxelField(map, prim);
        result->heightmapFilename = filename;
        return result;
    }
}

ImplicitPrimitive *ImplicitPrimitive::fromHeightmap(GridF heightmap, std::string filename, ImplicitPrimitive *prim)
{
    if (heightmap.empty()) throw std::invalid_argument("Cannot create an implicit terrain from an empty heightmap");
    for (float value : heightmap) if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite implicit height sample");

    if (heightmap.sizeZ > 1) return fromVoxelField(heightmap, prim);
    if (prim == nullptr)
        prim = new ImplicitPrimitive;
    prim->sampledVolume = false;
    prim->sampledMaterialColumns = MaterialColumns();
    prim->used_json_filename.clear();
    prim->referenceProperties = nlohmann::json::object();
    if (heightmap.height() > 1)
        prim->setDimensions(heightmap.getDimensions());
    else
        prim->setDimensions(heightmap.getDimensions().xy() + Vector3(0, 0, heightmap.max()));
    prim->position = Vector3(0, 0, 0);
//    prim->evalFunction = ImplicitPrimitive::convert2DfunctionTo3Dfunction([=](const Vector3& pos) -> float {
//        return heightmap.data[heightmap.getIndex(pos.xy())];
//    });

    prim->name = "Heightmap";
    prim->heightmapFilename = filename;
    prim->predefinedShape = ImplicitHeightmap;
    prim->cachedHeightmap = heightmap;
    prim->cachedHeightmap.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    prim->cachedHeightmap.raiseErrorOnBadCoord = false;
    prim->material = ROCK;
    prim->update();
    return prim;
}

ImplicitPrimitive* ImplicitPrimitive::fromVoxelField(const GridF& values, ImplicitPrimitive* prim)
{
    if (values.empty()) throw std::invalid_argument("Cannot create an implicit terrain from an empty field");
    for (float value : values) if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite implicit volume sample");
    if (!prim) prim = new ImplicitPrimitive;
    prim->sampledVolume = true;
    prim->sampledMaterialColumns = MaterialColumns();
    prim->heightmapFilename.clear();
    prim->used_json_filename.clear();
    prim->referenceProperties = nlohmann::json::object();
    prim->cachedHeightmap = values;
    prim->cachedHeightmap.raiseErrorOnBadCoord = false;
    prim->cachedHeightmap.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    prim->position = Vector3();
    prim->setDimensions(values.getDimensions());
    prim->setSupportDimensions(values.getDimensions());
    prim->predefinedShape = ImplicitHeightmap;
    prim->material = ROCK;
    prim->name = "Sampled terrain";
    prim->update();
    return prim;
}

ImplicitPrimitive* ImplicitPrimitive::fromLayerGrid(LayerBasedGrid& layers, ImplicitPrimitive* prim)
{
    for (const auto& column : layers.getLayers()) {
        float total = 0.f;
        for (const auto& [label, thickness] : column) {
            total += thickness;
            if (label < AIR || label >= LAST || !std::isfinite(thickness) || thickness <= 0.f || !std::isfinite(total))
                throw std::invalid_argument("Invalid layered snapshot interval");
        }
    }
    std::unique_ptr<ImplicitPrimitive> owned;
    if (!prim) { owned = std::make_unique<ImplicitPrimitive>(); prim = owned.get(); }
    fromVoxelField(layers.voxelize(), prim);
    prim->sampledMaterialColumns = layers.getLayers();
    prim->scaling = layers.scaling;
    prim->translation = layers.translation;
    prim->name = "Layered terrain snapshot";
    return owned ? owned.release() : prim;
}







ImplicitBinaryOperator::ImplicitBinaryOperator()
{
    this->composables = std::vector<ImplicitPatch*>(2);
    composeFunction = STACK;
    positionalB = FIXED_POS;
}

float ImplicitBinaryOperator::evaluate(const Vector3& pos) const
{
    return getMaterialsAndTotalEvaluation(pos).first;
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterials(const Vector3& pos) const
{
    // Get materials on A, get materials on B
    // Take the operator final evaluation
    // Use the ratio "mat on A / final eval" and "mat on B / final eval"
    // to determine final materials
    auto [contributionA, materialsAndEvalsForA] = this->getMaterialsAndTotalEvaluationA(pos);
    auto [contributionB, materialsAndEvalsForB] = this->getMaterialsAndTotalEvaluationB(pos);
    float operationEvaluation = this->evaluateFromAandB(contributionA, contributionB);

    std::map<TerrainTypes, float> result;


    if (operationEvaluation == 0) {
        return {};
    }
    if (composeFunction == DIFFERENCE && !withIntersectionOnB) {
        for (const auto& [material, amount] : materialsAndEvalsForA)
            result[material] = contributionA > 0.f ? amount * operationEvaluation / contributionA : 0.f;
        return result;
    }
    float mu_A = contributionA;
    float mu_B = contributionB;
    float rA = mu_A / (mu_A + mu_B);
    float rB = mu_B / (mu_A + mu_B);
    float tempA = rA * mu_A;
    float tempB = rB * mu_B;
    float globalCoef = operationEvaluation / (tempA + tempB);
    float ratioA = rA * globalCoef;
    float ratioB = rB * globalCoef;
    if (composeFunction == NONE && !withIntersectionOnB) {
        return materialsAndEvalsForA;
    }
    if (this->composeFunction == REPLACE || this->composeFunction == ONE_SIDE_BLEND) {
        if (contributionB > ImplicitPatch::isovalue) {
            ratioA = 0.f;
            ratioB = (contributionB > 0 ? (operationEvaluation / contributionB) : 0.f);
        } else {
            ratioA = (contributionA > 0 ? (operationEvaluation / contributionA) : 0.f);
            ratioB = 0.f;
        }
    }

    // Normally, we should get : ratioA + ratioB = 1.0

    for (const auto& [mat, val] : materialsAndEvalsForA) {
        if (val > 0) {
            if (result.count(mat) == 0)
                result[mat] = 0.f;
            result[mat] += val * ratioA;
        }
    }
    for (const auto& [mat, val] : materialsAndEvalsForB) {
        if (val > 0) {
            if (result.count(mat) == 0)
                result[mat] = 0.f;
            result[mat] += val * ratioB;
        }
    }

    return result;
}

float ImplicitBinaryOperator::evaluateFromAandB(float evalA, float evalB) const
{
    if (!std::isfinite(blendingFactor) || blendingFactor <= 0.f)
        throw std::invalid_argument("Implicit blend exponent must be finite and positive");
    if (!std::isfinite(evalA) || !std::isfinite(evalB) || !std::isfinite(isovalue))
        throw std::invalid_argument("Implicit composition requires finite field values and threshold");
    if (withIntersectionOnB || composeFunction == INTERSECTION) return std::min(evalA, evalB);
    float evaluation = 0.f;
    if (this->composeFunction == STACK) {
        evaluation = std::max(evalA, evalB);
    } else if (this->composeFunction == REPLACE) {
        evaluation = (evalB > ImplicitPatch::isovalue ? evalB : evalA);
    } else if (composeFunction == DIFFERENCE) {
        evaluation = std::min(evalA, std::clamp(2.f * isovalue - evalB, 0.f, 1.f));
    } else if (this->composeFunction == BLEND || this->composeFunction == ONE_SIDE_BLEND) {
        if (composeFunction == ONE_SIDE_BLEND) {
            float weight = std::clamp(evalB / std::max(isovalue, 1e-6f), 0.f, 1.f);
            weight = weight * weight * (3.f - 2.f * weight);
            evalB *= weight;
        }
        evaluation = std::pow(std::pow(evalA, this->blendingFactor) + std::pow(evalB, this->blendingFactor), 1.f/this->blendingFactor);
    } else if (composeFunction == NONE) {
        evaluation = evalA;
    } else {
        throw std::invalid_argument("Unknown implicit composition function");
    }
    evaluation = std::clamp(evaluation, 0.f, 1.f);
    return evaluation;
}

float ImplicitBinaryOperator::evaluateA(const Vector3& pos) const
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);

    return composableA() ? composableA()->evaluate(evaluationPosA) : 0.f;
}

float ImplicitBinaryOperator::evaluateB(const Vector3& pos) const
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);

    return composableB() ? composableB()->evaluate(evaluationPosB) : 0.f;
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterialsA(const Vector3& pos) const
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return composableA() ? composableA()->getMaterials(evaluationPosA) : std::map<TerrainTypes, float>{};
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterialsB(const Vector3& pos) const
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return composableB() ? composableB()->getMaterials(evaluationPosB) : std::map<TerrainTypes, float>{};
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitBinaryOperator::getMaterialsAndTotalEvaluationA(const Vector3& pos) const
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return composableA() ? composableA()->getMaterialsAndTotalEvaluation(evaluationPosA) : std::pair<float, std::map<TerrainTypes, float>>{0.f, {}};
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitBinaryOperator::getMaterialsAndTotalEvaluationB(const Vector3& pos) const
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return composableB() ? composableB()->getMaterialsAndTotalEvaluation(evaluationPosB) : std::pair<float, std::map<TerrainTypes, float>>{0.f, {}};
}

bool ImplicitBinaryOperator::contains(PredefinedShapes shape)
{
    return (composableA() && composableA()->contains(shape)) || (composableB() && composableB()->contains(shape));
}

std::vector<ImplicitPatch *> ImplicitBinaryOperator::findAll(PredefinedShapes shape)
{
    return vectorMerge(composableA() ? composableA()->findAll(shape) : std::vector<ImplicitPatch*>{},
                       composableB() ? composableB()->findAll(shape) : std::vector<ImplicitPatch*>{});
}

namespace {
AABBox implicitBinaryBounds(const ImplicitBinaryOperator& operation, bool support)
{
    const auto* a = operation.composableA();
    const auto* b = operation.composableB();
    if (!a && !b) return AABBox(Vector3(), Vector3());
    auto box = [&](const ImplicitPatch* node) { return support ? node->getSupportBBox() : node->getBBox(); };
    if (!a) return box(b);
    if (!b || operation.withIntersectionOnB || operation.composeFunction == ImplicitPatch::INTERSECTION ||
        operation.composeFunction == ImplicitPatch::DIFFERENCE || operation.composeFunction == ImplicitPatch::NONE) return box(a);
    const auto boundsA = box(a);
    auto boundsB = box(b);
    if (operation.positionalB != ImplicitPatch::FIXED_POS) {
        const auto heights = a->getSupportBBox();
        float low = std::min(0.f, heights.min().z());
        float high = std::max(0.f, heights.max().z());
        if (operation.positionalB == ImplicitPatch::INSIDE_TOP) {
            low -= b->getDimensions().z();
            high -= b->getDimensions().z();
        }
        boundsB.mini.z() += low;
        boundsB.maxi.z() += high;
    }
    return AABBox(Vector3::min(boundsA.min(), boundsB.min()), Vector3::max(boundsA.max(), boundsB.max()));
}
}

AABBox ImplicitBinaryOperator::getSupportBBox() const
{
    return implicitBinaryBounds(*this, true);
}

AABBox ImplicitBinaryOperator::getBBox() const
{
    return implicitBinaryBounds(*this, false);
}

void ImplicitBinaryOperator::update()
{
    // Guess we have nothing to do...
}

std::string ImplicitBinaryOperator::toString() const
{
    return this->name + ": Operation between " + (composableA() ? "#" + std::to_string(composableA()->index) : "undefined") + " and " + (composableB() ? "#" + std::to_string(composableB()->index) : "undefined");
}

nlohmann::json ImplicitBinaryOperator::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty() && !snapshotSerialization) {
        content = referenceProperties;
        content["file"] = this->used_json_filename;
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = patchOffset;
        content["scale"] = patchRescale); */
    } else {
        content["type"] = "compose";
        content["name"] = this->name;
        content["operator"] = stringFromCompositionOperation(this->composeFunction);
        content["blendingFactor"] = this->blendingFactor;
        content["index"] = this->index;
        content["positionalB"] = stringFromPositionalLabel(this->positionalB);
        content["useIntersection"] = this->withIntersectionOnB;
        if (this->composableA())
            content["composableA"] = this->composableA()->toJson();
        if (this->composableB())
            content["composableB"] = this->composableB()->toJson();
        content["optionalCurve"] = this->optionalCurve;
        content["mirrored"] = this->mirrored;
    }

    content["terrainScaling"] = scaling;
    content["terrainTranslation"] = translation;
    return content;
}

ImplicitPatch *ImplicitBinaryOperator::fromJson(nlohmann::json content)
{
    ImplicitBinaryOperator* patch = new ImplicitBinaryOperator;
    std::unique_ptr<ImplicitBinaryOperator> owner(patch);
    patch->addChild(ImplicitPatch::fromJson(content["composableA"]), 0); //patch->composableA() = ImplicitPatch::fromJson(content["composableA"]);
    patch->addChild(ImplicitPatch::fromJson(content["composableB"]), 1); //patch->composableB() = ImplicitPatch::fromJson(content["composableB"]);
    patch->name = content["name"];
    patch->composeFunction = compositionOperationFromString(content["operator"]);
    patch->blendingFactor = content["blendingFactor"];
    patch->positionalB = positionalLabelFromString(content["positionalB"]);
    patch->setIndex(content["index"]); // Advance beyond IDs recorded by older sessions.
    patch->setIndex(); // This imported node is a new live instance, not an alias.
    if (content.contains("useIntersection"))
        patch->withIntersectionOnB = content["useIntersection"];
    if (content.contains("optionalCurve"))
        patch->optionalCurve = content["optionalCurve"];
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    patch->update();
    readImplicitTerrainFrame(*patch, content);
    return owner.release();
}

void ImplicitBinaryOperator::updateCache()
{
    ImplicitPatch::updateCache();
    if (this->composableA() != nullptr)
        this->composableA()->updateCache();
    if (this->composableB() != nullptr)
        this->composableB()->updateCache();
}

void ImplicitBinaryOperator::swapAB()
{
    std::swap(composables[0], composables[1]);
    for (auto* child : composables) if (child) child->updateCache();
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;
}

void ImplicitBinaryOperator::deleteAllChildren()
{
    for (auto& child : composables) {
        delete child;
    }
    this->composables.assign(2, nullptr);
}

Vector3 ImplicitBinaryOperator::getEvaluationPositionForComposableA(const Vector3& pos) const
{
    return pos; // Nothing to do
}

Vector3 ImplicitBinaryOperator::getEvaluationPositionForComposableB(const Vector3& pos) const
{
    if (!composableA() || !composableB()) return pos;
    // Get the correct evaluation position for the B component
    float offsetB = 0.f;
    Vector3 positionB = pos;
    if (this->positionalB == ABOVE) {
        offsetB = this->composableA()->getMaxHeight(positionB);
    } else if (this->positionalB == INSIDE_TOP) {
        offsetB = this->composableA()->getMaxHeight(positionB) - this->composableB()->getDimensions().z();
    } else if (this->positionalB == INSIDE_BOTTOM) {
        offsetB = this->composableA()->getMinHeight(positionB);
    } else if (this->positionalB == FIXED_POS) {
        offsetB = 0.f;
    } else if (this->positionalB == SMOOTH_ABOVE) {
        // Smooth the supporting height, rather than multiplying world Z by
        // a possibly zero factor and turning B into an unbounded column.
        const int weights[3] = {1, 2, 1};
        float sum = 0.f;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                sum += weights[dx + 1] * weights[dy + 1] * composableA()->getMaxHeight(pos.xy() + Vector3(dx, dy));
        offsetB = sum / 16.f;
    } else {
        std::cerr << "Wrong position label" << std::endl;
    }
    positionB.z() -= offsetB;
    return positionB;
}

ImplicitPatch *ImplicitBinaryOperator::copy() const
{
    ImplicitBinaryOperator* copy = new ImplicitBinaryOperator(*this);
    return copy;
}

void ImplicitBinaryOperator::addChild(ImplicitPatch *newChild, int index)
{
    if (index != 0 && index != 1)
        throw std::out_of_range("On a binary operator, only the child[0] or child[1] can be affected. In the function, the index given is " + std::to_string(index));

    ImplicitNaryOperator::addChild(newChild, index);
}

ImplicitPatch *ImplicitBinaryOperator::composableA() const
{
    return this->composables[0];
}
ImplicitPatch* ImplicitBinaryOperator::composableB() const
{
    return this->composables[1];
}





ImplicitUnaryOperator::ImplicitUnaryOperator()
{
    this->composables = std::vector<ImplicitPatch*>(1);
    this->noiseFunction = [](const Vector3&) { return 0.f; };
    _generatedNoiseType = noiseFunction.target_type().hash_code();
}
float ImplicitUnaryOperator::evaluate(const Vector3& pos) const
{
    return getMaterialsAndTotalEvaluation(pos).first;
}

void ImplicitUnaryOperator::addTranslation(const Vector3 &translation) {
    this->_translation += translation;
    auto translationTransform = [translation](const Vector3& point) {
        return point + translation;
    };
    auto inverseTranslationTransform = [translation](const Vector3& point) {
        return point - translation;
    };
    addTransformation(translationTransform, inverseTranslationTransform);
    operationDescriptions.back() = {{"kind", "translation"}, {"value", translation}};
}

void ImplicitUnaryOperator::addRotation(const Vector3& rotationAngles) {
    this->_rotation += rotationAngles;
    Matrix rotationMatrix = rotationAngles.toRotationMatrix();
    Matrix inverseRotationMatrix = rotationMatrix.inverse();
    auto rotationTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + Vector3::fromMatrix(rotationMatrix.product((point - center).toMatrix()));
    };
    auto inverseRotationTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + Vector3::fromMatrix(inverseRotationMatrix.product((point - center).toMatrix()));
    };
    addTransformation(rotationTransform, inverseRotationTransform);
    operationDescriptions.back() = {{"kind", "rotation"}, {"value", rotationAngles}};
}

void ImplicitUnaryOperator::addScaling(const Vector3 &scaleFactors) {
    if (!scaleFactors.isValid() || !std::isfinite(scaleFactors.x()) ||
        !std::isfinite(scaleFactors.y()) || !std::isfinite(scaleFactors.z()) ||
        scaleFactors.x() == 0.f || scaleFactors.y() == 0.f || scaleFactors.z() == 0.f)
        throw std::invalid_argument("Implicit scaling requires finite nonzero factors");
    this->_scale *= scaleFactors;
    auto scalingTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + (point - center) * scaleFactors; //Vector3(point.x() * scaleFactors.x, point.y() * scaleFactors.y, point.z() * scaleFactors.z);
    };
    auto inverseScalingTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + (point - center) / scaleFactors; // Vector3(point.x() / scaleFactors.x, point.y() / scaleFactors.y, point.z() / scaleFactors.z);
    };
    addTransformation(scalingTransform, inverseScalingTransform);
    operationDescriptions.back() = {{"kind", "scaling"}, {"value", scaleFactors}};
}

void ImplicitUnaryOperator::addDisplacementField(const GridV3 &displacementField) {
    if (displacementField.empty()) throw std::invalid_argument("Empty implicit displacement field");
    auto field = std::make_shared<GridV3>(displacementField);
    for (const auto& v : *field)
        if (!v.isValid() || !std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z()))
            throw std::invalid_argument("Nonfinite implicit displacement");
    field->raiseErrorOnBadCoord = false;
    field->returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    addDisplacementFunction([field](const Vector3& point) { return field->interpolate(point); },
                            Vector3::min(field->data), Vector3::max(field->data));
    operationDescriptions.back() = {{"kind", "displacement"}, {"field", *field}};
}

void ImplicitUnaryOperator::addDisplacementFunction(TransformFunction displacement, const Vector3& minimum, const Vector3& maximum)
{
    auto forward = [displacement](const Vector3& point) { return point + displacement(point); };
    auto inverse = [displacement](const Vector3& target) {
        Vector3 point = target;
        for (int iteration = 0; iteration < 100; ++iteration) {
            Vector3 next = target - displacement(point);
            if ((next - point).length() <= 1e-5f) return next;
            point = next;
        }
        throw std::runtime_error("Implicit displacement inverse did not converge");
    };
    addTransformation(forward, inverse);
    boundTransformations.back() = [minimum, maximum](const AABBox& box) {
        return AABBox(box.min() + minimum, box.max() + maximum);
    };
}

Vector3 ImplicitUnaryOperator::applyTransform(const Vector3& pos) const
{
    Vector3 finalPos = pos;
    for (int i = 0; i < transformations.size(); i++)
        finalPos = transformations[i].first(finalPos);
    return finalPos;
}

Vector3 ImplicitUnaryOperator::inverseTransform(const Vector3& pos) const
{
    Vector3 finalPos = pos;
    for (int i = transformations.size() - 1; i >= 0; i--)
        finalPos = transformations[i].second(finalPos);
    return finalPos;
}

std::map<TerrainTypes, float> ImplicitUnaryOperator::getMaterials(const Vector3& pos) const
{
    if (!composableA()) return {};
    auto [value, materials] = composableA()->getMaterialsAndTotalEvaluation(inverseTransform(pos));
    if (value > 0.f) {
        const float noisy = std::clamp(value + noiseFunction(pos), 0.f, 1.f);
        for (auto& [material, amount] : materials) amount *= noisy / value;
    }
    return materials;
}

AABBox ImplicitUnaryOperator::transformBounds(AABBox bounds) const
{
    for (size_t i = 0; i < transformations.size(); ++i) {
        if (boundTransformations[i]) bounds = boundTransformations[i](bounds);
        else {
            auto vertices = Vector3::getAABBoxVertices(bounds.min(), bounds.max());
            for (auto& vertex : vertices) vertex = transformations[i].first(vertex);
            bounds = AABBox(Vector3::min(vertices), Vector3::max(vertices));
        }
    }
    return bounds;
}

AABBox ImplicitUnaryOperator::getSupportBBox() const
{
    return composableA() ? transformBounds(composableA()->getSupportBBox()) : AABBox(Vector3(), Vector3());
}

AABBox ImplicitUnaryOperator::getBBox() const
{
    return composableA() ? transformBounds(composableA()->getBBox()) : AABBox(Vector3(), Vector3());
}

std::string ImplicitUnaryOperator::toString() const
{
    return this->name + ": Unary operation on " + (composableA() ? "#" + std::to_string(composableA()->index) : "undefined");
}

nlohmann::json ImplicitUnaryOperator::toJson()
{
    if (noiseFunctionType() != _generatedNoiseType && !snapshotSerialization)
        throw std::runtime_error("Sample a custom noise callback before serializing it");
    nlohmann::json content;
    if (!this->used_json_filename.empty() && !snapshotSerialization) {
        content = referenceProperties;
        content["file"] = this->used_json_filename;
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = patchOffset;
        content["scale"] = patchRescale); */
    } else {
        content["type"] = "unary";
        content["name"] = this->name;
        content["index"] = this->index;
        if (this->composableA())
            content["composableA"] = this->composableA()->toJson();
        content["optionalCurve"] = this->optionalCurve;
        for (const auto& operation : operationDescriptions)
            if (operation.is_null()) throw std::runtime_error("Cannot serialize an arbitrary implicit transform callback");
        content["operations"] = operationDescriptions;
        content["noise"] = _noise;
        content["spreadFactor"] = this->_spreadingFactor;
        content["mirrored"] = this->mirrored;
    }

    content["terrainScaling"] = scaling;
    content["terrainTranslation"] = translation;
    return content;
}

ImplicitPatch *ImplicitUnaryOperator::fromJson(nlohmann::json content)
{
    ImplicitUnaryOperator* patch; // = new ImplicitUnaryOperator;
    Vector3 rotation = Vector3(0, 0, 0);
    Vector3 translation = Vector3(0, 0, 0);
    Vector3 scale = Vector3(1, 1, 1);
    Vector3 noise = Vector3(0, 0, 0);
    Vector3 distortion = Vector3(0, 0, 0);
    if (content.contains("translation") && content["translation"] != Vector3()) {
        patch = new ImplicitTranslation();
        translation = content["translation"];
    } else if (content.contains("scale") && content["scale"] != Vector3(1, 1, 1)) {
        patch =  new ImplicitScaling();
        scale = content["scale"];
        if (scale == Vector3())
            scale = Vector3(1.f, 1.f, 1.f);
    } else if (content.contains("rotation") && content["rotation"] != Vector3()) {
        patch =  new ImplicitRotation();
        rotation = content["rotation"];
    } else if (content.contains("noise") && content["noise"] != Vector3()) {
        patch = new ImplicitNoise();
        noise = content["noise"];
    } else if (content.contains("distortion") && content["distortion"] != Vector3()) {
        patch = new ImplicitWraping();
        distortion = content["distortion"];
    } else if (content.contains("spreadFactor") && content["spreadFactor"] != 0) {
        patch = new ImplicitSpread();
        // Apply spreading only after the child has been attached.
    } else {
        patch = new ImplicitUnaryOperator();
    }
    std::unique_ptr<ImplicitUnaryOperator> owner(patch);
    patch->addChild(ImplicitPatch::fromJson(content["composableA"]));
    translation = content.value("translation", Vector3());
    rotation = content.value("rotation", Vector3());
    scale = content.value("scale", Vector3(1, 1, 1));
    distortion = content.value("distortion", Vector3());
    patch->name = content["name"];
    patch->setIndex(content["index"]); // Advance beyond IDs recorded by older sessions.
    patch->setIndex(); // This imported node is a new live instance, not an alias.

    if (content.contains("operations")) {
        for (const auto& operation : content.at("operations")) {
            const std::string kind = operation.at("kind");
            if (kind == "translation") patch->translate(operation.at("value").get<Vector3>());
            else if (kind == "rotation") patch->rotate(operation.at("value").get<Vector3>());
            else if (kind == "scaling") patch->scale(operation.at("value").get<Vector3>());
            else if (kind == "displacement") patch->addDisplacementField(operation.at("field").get<GridV3>());
            else if (kind == "wrap") patch->addWrapFunction(operation.at("field").get<GridV3>());
            else if (kind == "spread") patch->spread(operation.at("value"));
            else if (kind == "wavelets") patch->addWavelets();
            else if (kind == "randomWrap") {
                const Vector3 value = operation.at("value");
                patch->addRandomWrap(value.x(), value.y(), value.z());
            } else throw std::invalid_argument("Unknown serialized implicit transformation: " + kind);
        }
    } else {
        if (translation != Vector3()) patch->translate(translation);
        if (scale != Vector3(1, 1, 1)) patch->scale(scale);
        if (rotation != Vector3()) patch->rotate(rotation);
        if (distortion != Vector3()) patch->addRandomWrap(distortion.x(), distortion.y(), distortion.z());
        if (content.contains("spreadFactor") && content["spreadFactor"] != 0) patch->spread(content["spreadFactor"]);
    }
    if (content.contains("noise")) noise = content["noise"];
    if (noise != Vector3()) patch->addRandomNoise(noise.x(), noise.y(), noise.z());

    if (content.contains("optionalCurve"))
        patch->optionalCurve = content["optionalCurve"];
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];

    patch->update();
    readImplicitTerrainFrame(*patch, content);
    return owner.release();
}

ImplicitPatch* ImplicitUnaryOperator::copy() const
{
    SnapshotScope snapshot;
    auto* result = ImplicitUnaryOperator::fromJson(const_cast<ImplicitUnaryOperator*>(this)->toJson());
    auto* unary = static_cast<ImplicitUnaryOperator*>(result);
    unary->deleteAllChildren();
    if (composableA()) unary->addChild(composableA()->copy());
    if (noiseFunctionType() != _generatedNoiseType) unary->setNoiseFunction(noiseFunction);
    result->parent = nullptr;
    std::vector<ImplicitPatch*> nodes{result};
    while (!nodes.empty()) {
        auto* node = nodes.back(); nodes.pop_back();
        node->setIndex();
        if (auto* group = dynamic_cast<ImplicitNaryOperator*>(node))
            for (auto* child : group->composables) if (child) nodes.push_back(child);
    }
    return result;
}

bool ImplicitUnaryOperator::contains(PredefinedShapes shape)
{
    return composableA() && composableA()->contains(shape);
}

std::vector<ImplicitPatch *> ImplicitUnaryOperator::findAll(PredefinedShapes shape)
{
    return composableA() ? composableA()->findAll(shape) : std::vector<ImplicitPatch*>{};
}

void ImplicitUnaryOperator::addChild(ImplicitPatch *newChild, int index)
{

    if (index != 0)
        throw std::out_of_range("On a unary operator, only the child[0] can be affected. In the function, the index given is " + std::to_string(index));

    ImplicitNaryOperator::addChild(newChild, index);
}

void ImplicitUnaryOperator::deleteAllChildren()
{
    for (auto& child : composables) {
        delete child;
    }
    this->composables.assign(1, nullptr);
}

ImplicitPatch* ImplicitUnaryOperator::composableA() const
{
    return this->composables[0];
}

void ImplicitUnaryOperator::translate(const Vector3& translation)
{
    addTranslation(translation);
}

void ImplicitUnaryOperator::rotate(const Vector3& angles)
{
    addRotation(angles);
}

void ImplicitUnaryOperator::scale(const Vector3& factor)
{
    addScaling(factor);
}

void ImplicitUnaryOperator::addRandomNoise(float amplitude, float period, float offset)
{
    if (!std::isfinite(amplitude) || !std::isfinite(period) || !std::isfinite(offset))
        throw std::invalid_argument("Implicit noise requires finite parameters");
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;

    this->_noise = Vector3(amplitude, period, offset);
    this->noiseFunction = [=](const Vector3& pos) -> float {
        FastNoiseLite noise;
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseVal = noise.GetNoise(pos.x() * period + offset, pos.y() * period + offset, pos.z() * period + offset) * 1.f;
        return noiseVal * amplitude;
    };
    _generatedNoiseType = noiseFunction.target_type().hash_code();
}

void ImplicitUnaryOperator::addRandomWrap(float amplitude, float period, float offset)
{
    if (!std::isfinite(amplitude) || !std::isfinite(period) || !std::isfinite(offset))
        throw std::invalid_argument("Nonfinite implicit noise displacement");
    _distortion = Vector3(amplitude, period, offset);
    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    auto displacement = [noise, amplitude, period, offset](const Vector3& p) mutable {
        return Vector3(1, 1, 1) * (amplitude * noise.GetNoise(p.x() * period + offset,
                                                          p.y() * period + offset, p.z() * period + offset));
    };
    const Vector3 margin( std::abs(amplitude), std::abs(amplitude), std::abs(amplitude));
    addDisplacementFunction(displacement, -margin, margin);
    operationDescriptions.back() = {{"kind", "randomWrap"}, {"value", _distortion}};
}

void ImplicitUnaryOperator::addWrapFunction(GridV3 field)
{
    if (field.empty() || !composableA()) throw std::invalid_argument("Wrapping requires a field and child");
    const AABBox domain = composableA()->getSupportBBox();
    for (const auto& value : field)
        if (!value.isValid() || !std::isfinite(value.x()) || !std::isfinite(value.y()) || !std::isfinite(value.z()))
            throw std::invalid_argument("Nonfinite implicit wrap displacement");
    if (domain.dimensions().minComp() <= 0.f) throw std::invalid_argument("Degenerate wrap domain");
    field.raiseErrorOnBadCoord = false;
    field.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    const Vector3 minimum = Vector3::min(field.data) * domain.dimensions();
    const Vector3 maximum = Vector3::max(field.data) * domain.dimensions();
    addDisplacementFunction([field, domain](const Vector3& point) {
        const Vector3 uv = Vector3::min(Vector3(1, 1, 1), Vector3::max(Vector3(), ((point - domain.min()) / domain.dimensions())));
        return field.interpolate(uv * (field.getDimensions() - Vector3(1, 1, 1))) * domain.dimensions();
    }, minimum, maximum);
    operationDescriptions.back() = {{"kind", "wrap"}, {"field", field}};
}

void ImplicitUnaryOperator::spread(float factor)
{
    if (!composableA() || !std::isfinite(factor) || std::abs(factor) > 10.f)
        throw std::invalid_argument("Spreading requires a child and a finite factor in [-10,10]");
    const auto box = getBBox();
    if (box.dimensions().z() <= 0.f) throw std::invalid_argument("Degenerate spreading height");
    // Exponential radial stretch stays invertible for both positive and
    // negative factors. Z is unchanged, so the inverse is analytic.
    auto stretch = [box, factor](const Vector3& point, bool inverse) {
        const float h = std::clamp((point.z() - box.min().z()) / box.dimensions().z(), 0.f, 1.f);
        const float scale = std::exp((inverse ? -factor : factor) * (1.f - 2.f * h));
        const auto center = box.center();
        return Vector3(center.x() + (point.x() - center.x()) * scale,
                       center.y() + (point.y() - center.y()) * scale, point.z());
    };
    addTransformation([stretch](const Vector3& p) { return stretch(p, false); },
                      [stretch](const Vector3& p) { return stretch(p, true); });
    const float maximumScale = std::exp(std::abs(factor));
    boundTransformations.back() = [box, maximumScale](const AABBox& input) {
        const auto center = box.center();
        const float rx = std::max(std::abs(input.min().x() - center.x()), std::abs(input.max().x() - center.x())) * maximumScale;
        const float ry = std::max(std::abs(input.min().y() - center.y()), std::abs(input.max().y() - center.y())) * maximumScale;
        return AABBox(Vector3(center.x() - rx, center.y() - ry, input.min().z()),
                      Vector3(center.x() + rx, center.y() + ry, input.max().z()));
    };
    _spreadingFactor += factor;
    operationDescriptions.back() = {{"kind", "spread"}, {"value", factor}};
}

void ImplicitUnaryOperator::addWavelets()
{
    addTransformation([](const Vector3& p) { return p + Vector3(0, 0, .5f * std::cos(p.x())); },
                      [](const Vector3& p) { return p - Vector3(0, 0, .5f * std::cos(p.x())); });
    boundTransformations.back() = [](const AABBox& box) {
        return AABBox(box.min() - Vector3(0, 0, .5f), box.max() + Vector3(0, 0, .5f));
    };
    operationDescriptions.back() = {{"kind", "wavelets"}};
}



std::function<float (const Vector3&)> ImplicitPatch::createSphereFunction(float sigma, float width, float depth, float height, bool in2D)
{
    const Vector3 radius = Vector3(width, depth, height) * .5f;
    if (in2D) return [radius](const Vector3& pos) {
        const Vector3 q = (pos.xy() - radius.xy()) / Vector3(radius.x(), radius.y(), 1.f);
        const float squared = q.norm2();
        return squared > 1.f ? 0.f : radius.z() * (1.f + std::sqrt(1.f - squared));
    };
    return [radius, sigma](const Vector3& pos) {
        return std::max(0.f, (1.f - ((pos - radius) / radius).norm()) * sigma);
    };
}

std::function<float (const Vector3&)> ImplicitPatch::createBlockFunction(float sigma, float width, float depth, float height, bool in2D)
{
    if (in2D) {
        return [=](const Vector3& pos) {
            return height;
        };
    } else {
        return [sigma, width, depth, height] (const Vector3& pos) {
            bool onlyUseVerticalDistance = (sigma > 10.f); // Hidden cheat code
            Vector3 minPos = Vector3(-width * 0.f, -depth * 0.f, -height * 0.f);
            Vector3 maxPos = Vector3(width * 1.f, depth * 1.f, height * 1.f);
            Vector3 normalizedPos = (pos - minPos) / (maxPos - minPos);
            float distanceToBBox = 0.f;
            if (onlyUseVerticalDistance)
                distanceToBBox = Vector3::signedDistanceToBoundaries(normalizedPos, Vector3(-1000, -1000, 0), Vector3(1000, 1000, 1)); // /*(pos / maxPos)*/pos, minPos, maxPos /*Vector3(1, 1, 1)*/);
            else
                distanceToBBox = Vector3::signedDistanceToBoundaries(normalizedPos, Vector3(0, 0, 0), Vector3(1, 1, 1)); // /*(pos / maxPos)*/pos, minPos, maxPos /*Vector3(1, 1, 1)*/);
            float distanceFalloff = 1.f - (distanceToBBox + 0.5f);
            float evaluation = std::clamp(distanceFalloff, 0.f, 1.f);
            return evaluation;
        };
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createGaussianFunction(float sigma, float width, float depth, float height, bool in2D)
{
    std::function<float (const Vector3&)> func2D = [width, depth, sigma, height](const Vector3& pos) { return normalizedGaussian(Vector3(width, depth, 0), pos.xy(), sigma) * height; };
    if (in2D) return func2D;
    else return ImplicitPatch::convert2DfunctionTo3Dfunction(func2D);
}

std::function<float (const Vector3&)> ImplicitPatch::createCylinderFunction(float sigma, float width, float depth, float height, bool in2D)
{
    const Vector3 radius(width * .5f, depth * .5f, 1.f);
    if (in2D) return [radius, height](const Vector3& pos) {
        return ((pos.xy() - radius.xy()) / radius).norm2() <= 1.f ? height : 0.f;
    };
    return [radius, height](const Vector3& pos) {
        const float radial = 1.f - .5f * ((pos.xy() - radius.xy()) / radius).norm();
        const float cap = .5f + std::min(pos.z(), height - pos.z()) / height;
        return std::clamp(std::min(radial, cap), 0.f, 1.f);
    };
}

std::function<float (const Vector3&)> ImplicitPatch::createRockFunction(float sigma, float width, float depth, float height, bool in2D)
{
    auto sphereFunction = ImplicitPatch::createSphereFunction(sigma, width, depth, height, in2D);
    std::function<float (const Vector3&)> rockFunction;
    if (in2D) {
        rockFunction = [=](const Vector3& pos) {
            float noiseOffset = sigma * 1000.f;
            float noiseValue = (random_gen::generate_perlin(pos.x() * (100.f / width) + noiseOffset, pos.y() * (100.f / width) + noiseOffset, pos.z() * (100.f / width) + noiseOffset) + 1.f) * width * .5f;
            return sphereFunction(pos) + noiseValue;
        };
    } else {
        rockFunction = [=] (const Vector3& pos) {
            /*FastNoiseLite noise;
    //        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
            noise.SetFractalType(FastNoiseLite::FractalType_FBm);*/
            float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
            float noiseValue = (random_gen::generate_perlin(pos.x() * (100.f / width) + noiseOffset, pos.y() * (100.f / width) + noiseOffset, pos.z() * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
            return sphereFunction(pos + Vector3(0, 0, 0/*width * .25f*/)) - (noiseValue * .5f);
        };
    }
    return rockFunction;
}

std::function<float (const Vector3&)> ImplicitPatch::createMountainFunction(float sigma, float width, float depth, float height, bool in2D)
{
    std::function<float (const Vector3&)> mountainFunction = [=](const Vector3& pos) {
        Vector3 translatedPos = pos - (Vector3(width, depth, 0) * .5f);
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(width, depth, 1.f) * .5f);
//        return normalizedPos.norm();
        float functionsHeight = std::clamp((1.f - normalizedPos.norm()), 0.f, 1.f) * height;
        return functionsHeight;
    };
    if (in2D) {
        return mountainFunction;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(mountainFunction);
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createDuneFunction(float sigma, float width, float depth, float height, bool in2D)
{
    const float crest = sigma > 0.f ? std::clamp(sigma / (1.f + sigma), .1f, .9f) : .65f;
    auto profile = [=](const Vector3& pos) {
        const float x = pos.x() / width, y = 2.f * pos.y() / depth - 1.f;
        if (x < 0.f || x > 1.f || std::abs(y) > 1.f) return 0.f;
        const float ridge = x < crest ? x / crest : (1.f - x) / (1.f - crest);
        return height * ridge * (1.f - y * y);
    };
    return in2D ? std::function<float(const Vector3&)>(profile) : convert2DfunctionTo3Dfunction(profile);
}

std::function<float (const Vector3&)> ImplicitPatch::createBasinFunction(float sigma, float width, float depth, float height, bool in2D)
{
    auto profile = [=](const Vector3& pos) {
        const Vector3 q = (pos.xy() - Vector3(width, depth) * .5f) / Vector3(width * .5f, depth * .5f, 1.f);
        const float radius = q.norm();
        return radius > 1.f ? 0.f : height * std::pow(radius, sigma > 0.f ? sigma : 2.f);
    };
    return in2D ? std::function<float(const Vector3&)>(profile) : convert2DfunctionTo3Dfunction(profile);
}

std::function<float (const Vector3&)> ImplicitPatch::createCaveFunction(float sigma, float width, float depth, float height, bool in2D)
{
    if (in2D) {
        return [=](const Vector3& pos) {
            return 0.f; // Nothing to do in 2D...
        };
    } else {
        std::function<float (const Vector3&)> caveFunc = [=] (const Vector3& pos) {
            Vector3 start = Vector3(width * .5f, depth * 1.f, height * 1.f);
            Vector3 p1 = Vector3(width * .5f, depth * .7f, height * .3f);
            Vector3 p2 = Vector3(width * .5f, depth * .4f, height * .2f);
            Vector3 end = Vector3(width * .5f, depth * .2f, height * .2f);
            CatmullRomSpline curve = CatmullRomSpline({start, p1, p2, end});
            float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed
            return 1.f - (curve.estimateDistanceFrom(pos) / (width * .5f));
        };
        return caveFunc;
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createArchFunction(float sigma, float width, float depth, float height, bool in2D)
{
    std::function<float (const Vector3&)> archFunc;

    Vector3 start = Vector3(width * .5f, depth * 0.f, height * 0.f);
    Vector3 p1 = Vector3(width * .5f, depth * .3f, height * 1.f);
    Vector3 p2 = Vector3(width * .5f, depth * .7f, height * 1.f);
    Vector3 end = Vector3(width * .5f, depth * 1.f, height * 0.f);
    CatmullRomSpline curve = CatmullRomSpline({start, p1, p2, end});
    CatmullRomSpline curve2D = CatmullRomSpline({start.xy(), p1.xy(), p2.xy(), end.xy()});
    float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed

    if (in2D) {
        archFunc = [=](const Vector3& pos) {
            float closestTime = curve2D.estimateClosestTime(pos.xy());
            float dist = (curve2D.getPoint(closestTime) - pos.xy()).norm2();
            float w = (width * width * .25f);
            if (dist > w) return 0.f;
            float h = std::sqrt(w - dist);
            h += curve.getPoint(closestTime).z();
            return h;
        };
    } else {
        archFunc = [=] (const Vector3& pos) {
            return 1.f - (curve.estimateDistanceFrom(pos) / (width * .5f));
        };
    }
    return archFunc;
}

std::function<float (const Vector3&)> ImplicitPatch::createNoise2DFunction(float sigma, float width, float depth, float height, bool in2D)
{
    std::function<float (const Vector3&)> noiseFunc = [sigma, width, depth, height](const Vector3& pos) -> float {
        FastNoiseLite noise;
    //    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseValue = (noise.GetNoise(pos.x() * (100.f * sigma / width) + noiseOffset, pos.y() * (100.f * sigma / depth) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
        return noiseValue * height;
    };
    if (in2D) {
        return noiseFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(noiseFunc);
    }
    //float noiseValue = (noise.GetNoise(pos.x() * (100.f / width) + noiseOffset, pos.y() * (100.f / width) + noiseOffset, pos.z() * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
}

std::function<float (const Vector3&)> ImplicitPatch::createMountainChainFunction(float sigma, float width, float depth, float height, CatmullRomSpline path, bool in2D)
{
    for (auto& point : path) point.z() = 0.f;
    std::function<float (const Vector3&)> chainFunc = [=] (const Vector3& pos) -> float {
        float closestTime = path.estimateClosestTime(pos.xy());
        Vector3 closestPoint = path.getPoint(closestTime);
//        Vector3 vertical(0, 0, 1);
//        Vector3 islandDirection = vertical.cross(path.getDirection(closestTime)).normalize();

        Vector3 translatedPos = pos - closestPoint;
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(sigma, sigma, 1.f) * .5f);
        float distance = normalizedPos.norm();
        float functionsHeight = std::clamp((1.f - distance), 0.f, 1.f) * sigma;
        return functionsHeight;
    };
    if (in2D) {
        return chainFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(chainFunc);
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createPolygonFunction(float sigma, float width, float depth, float height, CatmullRomSpline path, bool in2D)
{
    Contour polygon(path);
    std::function<float (const Vector3&)> polygonFunc = [=] (const Vector3& pos) -> float {
        return (polygon.containsXY(pos.xy(), false) ? height : 0.f);
    };

    if (in2D) {
        return polygonFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(polygonFunc);
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createDistanceMapFunction(float sigma, float width, float depth, float height, CatmullRomSpline path, bool in2D)
{
    std::mt19937 generator(static_cast<unsigned>(std::hash<float>{}((path.getPoint(0) + path.getPoint(.1f)).divergence())));
    std::uniform_real_distribution<float> uniform(0.f, 1.f);
    int nbPoints = 50;
    int nbPaths = 50;
    Contour polygon(path);
    for (auto& p : *polygon.curve)
        p.z() = 0;
    AABBox bbox = polygon.curve->getAABBox();
    Vector3 size = bbox.dimensions();
    polygon.translate(-bbox.min()); //(-(bbox.min()));
    // auto points = ShapeCurve(path).randomPointsInside(50);
    int maxTries = 1000 * nbPoints;
    std::vector<Vector3> points;
    while(points.size() < nbPoints && maxTries > 0) {
        Vector3 newPoint = Vector3(uniform(generator) * size.x(), uniform(generator) * size.y());
        if (polygon.containsXY(newPoint)) points.push_back(newPoint);
        else maxTries--;
    }
    if (points.size() == 0) return [](const Vector3&) {return 0.f;};

    GridF heights(size.x(), size.y(), 1);

    std::vector<CatmullRomSpline> randomPaths(nbPaths);
    for (int i = 0; i < randomPaths.size(); i++) {
        auto& path = randomPaths[i];
        int i0 = std::uniform_int_distribution<int>(0, int(points.size()) - 1)(generator);
        int i1 = std::uniform_int_distribution<int>(0, int(points.size()) - 1)(generator);
        // while (i0 == i1)
        // i1 = std::uniform_int_distribution<int>(0, int(points.size()) - 1)(generator);
        path = CatmullRomSpline({points[i0], points[i1]}).resamplePoints(4);
        float length = (path[-1] - path[0]).norm();
        for (auto& p : path) {
            p += Vector3(2.f * uniform(generator) - 1.f, 2.f * uniform(generator) - 1.f) * (length / 10.f);
        }
    }
    heights.iterateParallel([&] (const Vector3& pos) {
        for (int i = 0; i < randomPaths.size(); i++) {
            auto& path = randomPaths[i];
            auto closestPoint = toPolyline(path).estimateClosestPos(pos);
            heights(pos) += std::exp(-(pos - closestPoint).norm2() * i);
        }
        heights(pos) += random_gen::generate_fbm(pos.x(), pos.y()) * (.1f * std::abs(heights(pos)));
    });

    heights.normalize();
    heights.iterateParallel([&] (size_t i) {
        heights[i] = std::exp(heights[i] * 3.f);
    });
    heights = heights.meanSmooth(5, 5, 1).normalized();
    // heights = heights.gaussianSmooth(2.f).normalized();
    heights.iterateParallel([&] (const Vector3& pos) {
        // heights(pos) -= 2.f;
        if (!polygon.containsXY(pos)) heights(pos) = 0;
    });
    heights = heights.meanSmooth(5, 5, 1).normalize();
    // heights = heights.gaussianSmooth(2.f).normalize();
    heights.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    heights.defaultValueOnBadCoord = 0.f;

    Vector3 offsets = Vector3(width - size.x(), depth - size.y()) * .5f;

    std::function<float (const Vector3&)> distFunc = [=] (const Vector3& pos) -> float {
        float finalHeight = sigma * std::clamp(heights(pos.xy() - offsets), 0.f, 1.f);
        return finalHeight;
    };

    if (in2D) {
        return distFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(distFunc);
    }

    /*
    ShapeCurve polygon(path.points);
    Vector3 center = polygon.centroid().xy();
    std::vector<Vector3> randomPoints = polygon.randomPointsInside(10);
    std::function<float (const Vector3&)> distFunc = [=] (const Vector3& pos) -> float {
        if (!polygon.containsXY(pos.xy(), false))
            return 0.f;
        Vector3 closest = polygon.estimateClosestPos(pos.xy()).xy();
        Vector3 closestCenter;
        float closestDist = std::numeric_limits<float>::max();
        for (auto& p : randomPoints) {
            float dist = (p - pos).norm2();
            if (dist < closestDist) {
                closestDist = dist;
                closestCenter = p;
            }
        }
        float dist = (closest - pos.xy()).norm() / (closestCenter - closest).norm(); // Take the closest point and use only the "xy" components!
        // float dist = (closest - pos.xy()).norm() / (center - closest).norm(); // Take the closest point and use only the "xy" components!
        return height * std::clamp(dist, 0.f, 1.f);
    };
    if (in2D) {
        return distFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(distFunc);
    }*/
}

std::function<float (const Vector3&)> ImplicitPatch::createParametricTunnelFunction(float sigma, float width, float depth, float height, CatmullRomSpline path, bool in2D)
{
    float epsilon = 1e-1;

    if (in2D) {
        auto flatPath = path;
        for (auto & p : flatPath)
            p = p.xy();
        return [=](const Vector3& pos) {
            float d = flatPath.estimateSqrDistanceFrom(pos.xy());
            float r = (sigma * sigma);
            if (d > r) return 0.f;
            return std::sqrt(r - d);
        };
    } else {
        return [=] (const Vector3& pos) -> float {
            float closestT = path.estimateClosestTime(pos);
            auto closestPoint = path.getPoint(closestT);
            float distance = (pos - closestPoint).norm();
            return 1.f - (distance / (sigma));
        };
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createRippleFunction(float sigma, float width, float depth, float height, bool in2D)
{
    std::function<float (const Vector3&)> rippleFunc = [=] (const Vector3& pos) -> float {
        if(!Vector3::isInBox(pos, Vector3(), Vector3(width, depth, height)))
            return 0.f;
        float x = 2 * pos.x() / width;
        float y = 2.f * pos.y() / depth - 1;
        float maxHeight = (x < 1.f ? x : 1.f/(x*x));
        float damping = std::max(1 - y*y, 0.f);
//        std::cout << x << " " << y << " " << maxHeight << " " << damping << std::endl;
        return height * maxHeight * damping;
    };
    if (in2D) {
        return rippleFunc;
    } else {
        return ImplicitPatch::convert2DfunctionTo3Dfunction(rippleFunc);
    }
}

std::function<float (const Vector3&)> ImplicitPatch::createIdentityFunction(float sigma, float width, float depth, float height, bool in2D)
{
    return [](const Vector3&) -> float { return 0.f; };
}

std::function<float (const Vector3&)> ImplicitPatch::convert2DfunctionTo3Dfunction(std::function<float (const Vector3&)> func)
{
    std::function<float (const Vector3&)> _3Dfunction = [func](const Vector3& pos) {
        float height = func(pos);
        if (!std::isfinite(height) || height <= 0.f) return 0.f;
        float inverse_isovalue = .5f - (pos.z() / height); //height * .5f - pos.z;
        float isovalue = 1.f - std::abs(inverse_isovalue);
        return std::max(0.f, isovalue);
    };
    return _3Dfunction;
}

ImplicitNaryOperator* ImplicitPatch::getParent() const
{
    return this->parent;
}

void ImplicitPatch::setParent(ImplicitNaryOperator *newParent)
{
    this->parent = newParent;
}

ImplicitPatch::CompositionFunction compositionOperationFromString(const std::string& _name)
{
    std::string name = toUpper(_name);
    if (name == "NONE")
        return ImplicitPatch::CompositionFunction::NONE;
    else if (name == "STACK")
        return ImplicitPatch::CompositionFunction::STACK;
    else if (name == "BLEND")
        return ImplicitPatch::CompositionFunction::BLEND;
    else if (name == "REPLACE")
        return ImplicitPatch::CompositionFunction::REPLACE;
    else if (name == "ONE_SIDE_BLEND")
        return ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND;
    else if (name == "DIFFERENCE" || name == "SUBTRACT")
        return ImplicitPatch::DIFFERENCE;
    else if (name == "INTERSECTION")
        return ImplicitPatch::INTERSECTION;

    // Error
    throw std::invalid_argument("Unknown implicit composition: " + name);
}

std::string stringFromCompositionOperation(ImplicitPatch::CompositionFunction operation)
{
    if (operation == ImplicitPatch::CompositionFunction::NONE)
        return "None";
    else if (operation == ImplicitPatch::CompositionFunction::STACK)
        return "Stack";
    else if (operation == ImplicitPatch::CompositionFunction::BLEND)
        return "Blend";
    else if (operation == ImplicitPatch::CompositionFunction::REPLACE)
        return "Replace";
    else if (operation == ImplicitPatch::CompositionFunction::ONE_SIDE_BLEND)
        return "One_side_blend";
    else if (operation == ImplicitPatch::DIFFERENCE)
        return "Difference";
    else if (operation == ImplicitPatch::INTERSECTION)
        return "Intersection";

    // Impossible
    std::cerr << "Unknown operation given." << std::endl;
    return "";
}

ImplicitPatch::PredefinedShapes predefinedShapeFromString(const std::string& _name)
{
    std::string name = simplify(_name);
    if (name == "sphere")
        return ImplicitPatch::PredefinedShapes::Sphere;
    else if (name == "block")
        return ImplicitPatch::PredefinedShapes::Block;
    else if (name == "gaussian")
        return ImplicitPatch::PredefinedShapes::Gaussian;
    else if (name == "cylinder")
        return ImplicitPatch::PredefinedShapes::Cylinder;
    else if (name == "rock")
        return ImplicitPatch::PredefinedShapes::Rock;
    else if (name == "mountain")
        return ImplicitPatch::PredefinedShapes::Mountain;
    else if (name == "dune")
        return ImplicitPatch::PredefinedShapes::Dune;
    else if (name == "basin")
        return ImplicitPatch::PredefinedShapes::Basin;
    else if (name == "cave")
        return ImplicitPatch::PredefinedShapes::Cave;
    else if (name == "arch")
        return ImplicitPatch::PredefinedShapes::Arch;
    else if (name == "noise2d")
        return ImplicitPatch::PredefinedShapes::Noise2D;
    else if (name == "mountainchain")
        return ImplicitPatch::PredefinedShapes::MountainChain;
    else if (name == "polygon")
        return ImplicitPatch::PredefinedShapes::Polygon;
    else if (name == "implicitheightmap")
        return ImplicitPatch::PredefinedShapes::ImplicitHeightmap;
    else if (name == "parametrictunnel")
        return ImplicitPatch::PredefinedShapes::ParametricTunnel;
    else if (name == "none")
        return ImplicitPatch::PredefinedShapes::None;
    else if (name == "ripple")
        return ImplicitPatch::PredefinedShapes::Ripple;
    else if (name == "distancemap")
        return ImplicitPatch::PredefinedShapes::DistanceMap;

    // Error
    throw std::invalid_argument("Unknown implicit shape: " + name);
}

std::string stringFromPredefinedShape(ImplicitPatch::PredefinedShapes shape)
{
    if (shape == ImplicitPatch::PredefinedShapes::Sphere)
        return "Sphere";
    else if (shape == ImplicitPatch::PredefinedShapes::Block)
        return "Block";
    else if (shape == ImplicitPatch::PredefinedShapes::Gaussian)
        return "Gaussian";
    else if (shape == ImplicitPatch::PredefinedShapes::Cylinder)
        return "Cylinder";
    else if (shape == ImplicitPatch::PredefinedShapes::Rock)
        return "Rock";
    else if (shape == ImplicitPatch::PredefinedShapes::Mountain)
        return "Mountain";
    else if (shape == ImplicitPatch::PredefinedShapes::Dune)
        return "Dune";
    else if (shape == ImplicitPatch::PredefinedShapes::Basin)
        return "Basin";
    else if (shape == ImplicitPatch::PredefinedShapes::Cave)
        return "Cave";
    else if (shape == ImplicitPatch::PredefinedShapes::Arch)
        return "Arch";
    else if (shape == ImplicitPatch::PredefinedShapes::Noise2D)
        return "Noise2D";
    else if (shape == ImplicitPatch::PredefinedShapes::MountainChain)
        return "MountainChain";
    else if (shape == ImplicitPatch::PredefinedShapes::Polygon)
        return "Polygon";
    else if (shape == ImplicitPatch::PredefinedShapes::ImplicitHeightmap)
        return "Implicit_Heightmap";
    else if (shape == ImplicitPatch::PredefinedShapes::ParametricTunnel)
        return "Parametric_Tunnel";
    else if (shape == ImplicitPatch::PredefinedShapes::None)
        return "None";
    else if (shape == ImplicitPatch::PredefinedShapes::Ripple)
        return "Ripple";
    else if (shape == ImplicitPatch::PredefinedShapes::DistanceMap)
        return "DistanceMap";

    // Impossible
    std::cerr << "Unknown predefined shape given." << std::endl;
    return "";
}

ImplicitPatch::PositionalLabel positionalLabelFromString(const std::string& _name)
{
    std::string name = toUpper(_name);

    if (name == "ABOVE")
        return ImplicitPatch::PositionalLabel::ABOVE;
    else if (name == "INSIDE_TOP")
        return ImplicitPatch::PositionalLabel::INSIDE_TOP;
    else if (name == "INSIDE_BOTTOM")
        return ImplicitPatch::PositionalLabel::INSIDE_BOTTOM;
    else if (name == "FIXED_POS")
        return ImplicitPatch::PositionalLabel::FIXED_POS;
    else if (name == "SMOOTH_ABOVE")
        return ImplicitPatch::PositionalLabel::SMOOTH_ABOVE;

    // Error
    throw std::invalid_argument("Unknown implicit positional label: " + name);
}

std::string stringFromPositionalLabel(ImplicitPatch::PositionalLabel label)
{
    if (label == ImplicitPatch::PositionalLabel::ABOVE)
        return "ABOVE";
    else if (label == ImplicitPatch::PositionalLabel::INSIDE_TOP)
        return "INSIDE_TOP";
    else if (label == ImplicitPatch::PositionalLabel::INSIDE_BOTTOM)
        return "INSIDE_BOTTOM";
    else if (label == ImplicitPatch::PositionalLabel::FIXED_POS)
        return "FIXED_POS";
    else if (label == ImplicitPatch::PositionalLabel::SMOOTH_ABOVE)
        return "SMOOTH_ABOVE";

    // Impossible
    std::cerr << "Unknown positional label given." << std::endl;
    return "";
}

TerrainTypes materialFromString(const std::string& _name)
{
    std::string name = toUpper(_name);

    if (name == "AIR")
        return AIR;
    else if (name == "CURRENT_MIDDLE")
        return CURRENT_MIDDLE;
    else if (name == "CURRENT_TOP")
        return CURRENT_TOP;
    else if (name == "CURRENT_BOTTOM")
        return CURRENT_BOTTOM;
    else if (name == "WATER")
        return WATER;
    else if (name == "CORAL")
        return CORAL;
    else if (name == "SAND")
        return SAND;
    else if (name == "DIRT")
        return DIRT;
    else if (name == "ROCK")
        return ROCK;
    else if (name == "BEDROCK")
        return BEDROCK;

    // Impossible
    throw std::invalid_argument("Unknown implicit material: " + name);
}

std::string stringFromMaterial(TerrainTypes material)
{
    if (material == AIR)
        return "AIR";
    else if (material == CURRENT_MIDDLE)
        return "CURRENT_MIDDLE";
    else if (material == CURRENT_TOP)
        return "CURRENT_TOP";
    else if (material == CURRENT_BOTTOM)
        return "CURRENT_BOTTOM";
    else if (material == WATER)
        return "WATER";
    else if (material == CORAL)
        return "CORAL";
    else if (material == SAND)
        return "SAND";
    else if (material == DIRT)
        return "DIRT";
    else if (material == ROCK)
        return "ROCK";
    else if (material == BEDROCK)
        return "BEDROCK";

    // Impossible
    std::cerr << "Unknown material type given." << std::endl;
    return "";
}

ImplicitNaryOperator::ImplicitNaryOperator()
{

}

ImplicitNaryOperator::ImplicitNaryOperator(const ImplicitNaryOperator& other)
    : ImplicitPatch(other)
{
    parent = nullptr;
    setIndex();
    try {
        for (auto* child : other.composables) {
            composables.push_back(child ? child->copy() : nullptr);
            if (composables.back()) composables.back()->parent = this;
        }
    } catch (...) {
        ImplicitNaryOperator::deleteAllChildren();
        throw;
    }
}

ImplicitNaryOperator& ImplicitNaryOperator::operator=(const ImplicitNaryOperator& other)
{
    if (&other == this) return *this;
    // Snapshot before deletion: the source may itself be a descendant.
    ImplicitNaryOperator replacement(other);
    auto* oldParent = parent;
    const int oldIndex = index;
    ImplicitNaryOperator::deleteAllChildren();
    ImplicitPatch::operator=(replacement);
    parent = oldParent;
    index = oldIndex;
    composables.swap(replacement.composables);
    for (auto* child : composables) if (child) child->parent = this;
    return *this;
}

ImplicitNaryOperator::~ImplicitNaryOperator()
{
    ImplicitNaryOperator::deleteAllChildren();
}

float ImplicitNaryOperator::evaluate(const Vector3& pos) const
{
    return getMaterialsAndTotalEvaluation(pos).first;
}

std::map<TerrainTypes, float> ImplicitNaryOperator::getMaterials(const Vector3& pos) const
{
    std::map<TerrainTypes, float> exterior;
    float strongestExterior = 0.f;
    for (int i = int(composables.size()) - 1; i >= 0; --i) {
        if (!composables[i]) continue;
        auto materials = composables[i]->getMaterials(pos);
        float value = 0.f;
        for (const auto& [material, amount] : materials)
            value += (LayerBasedGrid::invisibleLayers.count(material) != 0) ? -amount : amount;
        const float strength = std::abs(value);
        // Keep the existing last-active-child composition, but preserve its
        // material labels and the continuous field just outside the surface.
        if (strength > isovalue) return materials;
        if (strength > strongestExterior) {
            strongestExterior = strength;
            exterior = std::move(materials);
        }
    }
    return exterior;
}

AABBox ImplicitNaryOperator::getSupportBBox() const
{
    Vector3 minPos(false), maxPos(false);
    for (auto& compo : this->composables) {
        if (!compo) continue;
        auto BBox = compo->getSupportBBox();
        minPos = Vector3::min(minPos, BBox.min());
        maxPos = Vector3::max(maxPos, BBox.max());
    }
    if (!minPos.isValid() || !maxPos.isValid()) return AABBox(Vector3(), Vector3());
    return {minPos, maxPos};
}

AABBox ImplicitNaryOperator::getBBox() const
{
    Vector3 minPos(false), maxPos(false);
    for (auto& compo : this->composables) {
        if (!compo) continue;
        auto BBox = compo->getBBox();
        minPos = Vector3::min(minPos, BBox.min());
        maxPos = Vector3::max(maxPos, BBox.max());
    }
    if (!minPos.isValid() || !maxPos.isValid()) return AABBox(Vector3(), Vector3());
    return {minPos, maxPos};
}

void ImplicitNaryOperator::update()
{
    // Nothing to do...
}

std::string ImplicitNaryOperator::toString() const
{
    std::vector<std::string> compoNames;
    for (auto& compo : this->composables)
        if (compo) compoNames.push_back(compo->name);
    return this->name + ": N-operation between " + join(compoNames, ",");
}

nlohmann::json ImplicitNaryOperator::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty() && !snapshotSerialization) {
        content = referenceProperties;
        content["file"] = this->used_json_filename;
    } else {
        content["type"] = "nary";
        content["heightfield"] = dynamic_cast<Implicit2DNary*>(this) != nullptr;
        content["name"] = this->name;
        content["index"] = this->index;
        content["optionalCurve"] = this->optionalCurve;
        content["mirrored"] = this->mirrored;
        std::vector<nlohmann::json> composableJson;
        for (auto& compo : this->composables)
            if (compo) composableJson.push_back(compo->toJson());
        content["composables"] = composableJson;
    }

    content["terrainScaling"] = scaling;
    content["terrainTranslation"] = translation;
    return content;
}

ImplicitPatch *ImplicitNaryOperator::fromJson(nlohmann::json content)
{
    ImplicitNaryOperator* patch = content.value("heightfield", false) ? new Implicit2DNary : new ImplicitNaryOperator;
    std::unique_ptr<ImplicitNaryOperator> owner(patch);
    for (auto& compoJson : content["composables"]) {
        patch->addChild(ImplicitPatch::fromJson(compoJson));
    }
    patch->name = content["name"];
    patch->setIndex(content["index"]); // Advance beyond IDs recorded by older sessions.
    patch->setIndex(); // This imported node is a new live instance, not an alias.
    if (content.contains("optionalCurve"))
        patch->optionalCurve = content["optionalCurve"];
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    patch->update();
    readImplicitTerrainFrame(*patch, content);
    return owner.release();
}

bool ImplicitNaryOperator::contains(PredefinedShapes shape)
{
    for (auto& compo : composables)
        if (compo && compo->contains(shape))
            return true;
    return false;
}

std::vector<ImplicitPatch *> ImplicitNaryOperator::findAll(PredefinedShapes shape)
{
    std::vector<ImplicitPatch*> res;
    for (auto& compo : this->composables)
        if (compo) res = vectorMerge(res, compo->findAll(shape));
    return res;
}

void ImplicitNaryOperator::simplify()
{
    /*
    while (this->composables.size() == 1) {
        auto& previousChild = this->composables[0];
        auto asNary = dynamic_cast<ImplicitNaryOperator*>(previousChild);
        if (asNary == nullptr)
            return;
        this->composables.clear();
        for (auto& subChild : asNary->composables) {
            this->addChild(subChild);
        }
    }
    */
}

void ImplicitNaryOperator::augment()
{
    if (this->contains(PredefinedShapes::Mountain) && this->contains(PredefinedShapes::MountainChain)) {
//        this->composables.push_back(ImplicitPrimitive::createPredefinedShape(Polygon, Vector3(0, 0, 20), 0.f, {}));
    }
}

void ImplicitNaryOperator::updateCache()
{
    ImplicitPatch::updateCache();
    for (auto& compo : this->composables)
        if (compo) compo->updateCache();
}

ImplicitPatch *ImplicitNaryOperator::copy() const
{
    if (const auto* heightfield = dynamic_cast<const Implicit2DNary*>(this)) return new Implicit2DNary(*heightfield);
    return new ImplicitNaryOperator(*this);
}

void ImplicitNaryOperator::addChild(ImplicitPatch *newChild, int index)
{
    if (!newChild) throw std::invalid_argument("Cannot attach a null implicit child");
    if (index < -1 || (index >= 0 && size_t(index) >= composables.size()))
        throw std::out_of_range("Implicit child index outside the existing slots");
    for (ImplicitPatch* ancestor = this; ancestor; ancestor = ancestor->parent)
        if (ancestor == newChild) throw std::invalid_argument("Implicit child would create a cycle");
    const auto existing = std::find(composables.begin(), composables.end(), newChild);
    if (existing != composables.end()) {
        if (index >= 0 && size_t(index) == size_t(existing - composables.begin())) return;
        throw std::invalid_argument("An implicit node cannot own the same child twice");
    }
    // Allocate before detaching so a failed append leaves the source intact.
    if (index == -1) composables.reserve(composables.size() + 1);
    if (auto* previousParent = newChild->parent) {
        auto slot = std::find(previousParent->composables.begin(), previousParent->composables.end(), newChild);
        if (slot == previousParent->composables.end())
            throw std::logic_error("Implicit parent does not contain its child");
        *slot = nullptr; // Preserve fixed-arity slots during wrapping/replacement.
        for (ImplicitPatch* node = previousParent; node; node = node->parent) node->_cached = false;
    }
    if (index == -1) composables.push_back(newChild);
    else {
        if (composables[index] && composables[index]->parent == this) composables[index]->parent = nullptr;
        composables[index] = newChild;
    }
    newChild->parent = this;
    for (ImplicitPatch* node = this; node; node = node->parent) node->_cached = false;
}

void ImplicitNaryOperator::deleteAllChildren()
{
    for (auto& child : composables) {
        if (child)
            delete child;
    }
    this->composables.clear();
}

nlohmann::json ImplicitTranslation::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["translation"] = this->_translation;
    return content;
}

ImplicitPatch *ImplicitTranslation::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitRotation::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["rotation"] = this->_rotation;
    return content;
}

ImplicitPatch *ImplicitRotation::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitScaling::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["scale"] = this->_scale;
    return content;
}

ImplicitPatch *ImplicitScaling::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitWraping::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["distortion"] = this->_distortion;
    return content;
}

ImplicitPatch *ImplicitWraping::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitNoise::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["noise"] = this->_noise;
    return content;
}

ImplicitPatch *ImplicitNoise::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitSpread::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["spreadFactor"] = this->_spreadingFactor;
    return content;
}

ImplicitPatch *ImplicitSpread::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

Implicit2DNary::Implicit2DNary()
    : ImplicitNaryOperator()
{

}

std::map<TerrainTypes, float> Implicit2DNary::getMaterials(const Vector3 &pos) const
{
    return {{DIRT, evaluate(pos)}};
}

float Implicit2DNary::evaluate(const Vector3 &pos) const
{
    const float range = std::max(1.f, getDimensions().z());
    return std::clamp(isovalue + (computeHeight(pos) - pos.z()) * (.5f / range), 0.f, 1.f);
}

GridF Implicit2DNary::getVoxelized(const Vector3& dimensions, const Vector3& scale)
{
    return ImplicitPatch::getVoxelized(dimensions, scale);
}

bool Implicit2DNary::contains(const Vector3 &v)
{
    return computeHeight(v) > v.z();
}

bool Implicit2DNary::checkIsInGround(const Vector3 &position)
{
    return this->contains(position);
}

Vector3 Implicit2DNary::getNormal(const Vector3 &pos) const
{
    const float x0 = computeHeight(pos - Vector3(1, 0));
    const float x1 = computeHeight(pos + Vector3(1, 0));
    const float y0 = computeHeight(pos - Vector3(0, 1));
    const float y1 = computeHeight(pos + Vector3(0, 1));
    return Vector3(x0 - x1, y0 - y1, 2.f).normalized();
}

float Implicit2DNary::getMaxHeight(const Vector3 &pos)
{
    return computeHeight(pos);
}

float Implicit2DNary::getMinHeight(const Vector3 &pos)
{
    return 0.f;
}

float Implicit2DNary::getMinimalHeight(const Vector3 &minBox, const Vector3 &maxBox)
{
    return 0.f;
}

float Implicit2DNary::getMaximalHeight(const Vector3 &minBox, const Vector3 &maxBox)
{
    return ImplicitPatch::getMaximalHeight(minBox, maxBox);
}

void Implicit2DNary::reevaluateAll()
{
    updateCache();
    getVoxelized();
}

namespace {
std::pair<float, float> implicitHeightContributions(const ImplicitPatch* node, const Vector3& point)
{
    if (!node) return {0.f, 0.f};
    if (const auto* primitive = dynamic_cast<const Implicit2DPrimitive*>(node)) {
        const float height = std::abs(primitive->evaluate(point.xy()));
        return LayerBasedGrid::invisibleLayers.count(primitive->material)
            ? std::pair<float, float>{0.f, height} : std::pair<float, float>{height, 0.f};
    }
    if (const auto* unary = dynamic_cast<const ImplicitUnaryOperator*>(node)) {
        if (unary->_noise.x() != 0.f || unary->noiseAt(point.xy()) != 0.f)
            throw std::invalid_argument("Use a 3D implicit root for density-noise operations");
        const Vector3 local = unary->inverseTransform(point.xy());
        const Vector3 zero = unary->applyTransform(local.xy());
        const Vector3 one = unary->applyTransform(local.xy() + Vector3(0, 0, 1));
        if ((one.xy() - zero.xy()).norm() > 1e-5f || one.z() <= zero.z())
            throw std::invalid_argument("Transform does not preserve a single-valued upward heightfield; use a 3D root");
        auto heights = implicitHeightContributions(unary->composableA(), local.xy());
        auto transformed = [&](float h) {
            return h > 0.f ? std::max(0.f, unary->applyTransform(local.xy() + Vector3(0, 0, h)).z()) : 0.f;
        };
        return {transformed(heights.first), transformed(heights.second)};
    }
    if (dynamic_cast<const ImplicitBinaryOperator*>(node))
        throw std::invalid_argument("Use a 3D implicit root for binary composition");
    if (const auto* group = dynamic_cast<const ImplicitNaryOperator*>(node)) {
        std::pair<float, float> result{0.f, 0.f};
        for (const auto* child : group->composables) {
            const auto heights = implicitHeightContributions(child, point);
            result.first = std::max(result.first, heights.first);
            result.second = std::max(result.second, heights.second);
        }
        return result;
    }
    throw std::invalid_argument("A heightfield root requires height-valued primitives");
}
}

float Implicit2DNary::computeHeight(const Vector3 &pos) const
{
    const auto heights = implicitHeightContributions(this, pos.xy());
    return heights.first - heights.second;
}



Implicit2DPrimitive::Implicit2DPrimitive()
    : ImplicitPrimitive()
{

}
float Implicit2DPrimitive::evaluate(const Vector3& pos) const
{
    const auto bounds = getBBox();
    if (!Vector3::isInBox(pos.xy(), bounds.min().xy(), bounds.max().xy())) return 0.f;
    const float localHeight = evalFunction(pos.xy() - position.xy());
    if (!std::isfinite(localHeight) || localHeight <= 0.f) return 0.f;
    const float height = std::max(0.f, position.z() + localHeight);
    return LayerBasedGrid::invisibleLayers.count(material) ? -height : height;
}

std::map<TerrainTypes, float> Implicit2DPrimitive::getMaterials(const Vector3& pos) const
{
    const float height = std::abs(evaluate(pos.xy()));
    if (height <= 0.f) return {{material, 0.f}};
    const float range = std::max(1.f, std::abs(dimensions.z()) + std::abs(position.z()));
    return {{material, std::clamp(isovalue + (height - pos.z()) * (.5f / range), 0.f, 1.f)}};
}
