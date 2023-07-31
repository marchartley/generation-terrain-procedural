#include "ImplicitPatch.h"
#include "Utils/ShapeCurve.h"
//#include "Utils/stb_image.h"
#include "Utils/Utils.h"

float ImplicitPatch::isovalue = .5f;
float ImplicitPatch::zResolution = .1f;
int ImplicitPatch::currentMaxIndex = -1;
std::string ImplicitPatch::json_identifier = "implicitPatches";


ImplicitPatch::ImplicitPatch()
{
    /*
    ImplicitPatch::currentMaxIndex ++;
    this->index = currentMaxIndex;*/
    this->setIndex();

    this->_cachedMaxHeight = Matrix3<float>(1, 1, 1, 0.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;
    this->_cachedMinHeight = Matrix3<float>(1, 1, 1, 0.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

float ImplicitPatch::getMaxHeight(Vector3 pos)
{
    auto AABBox = this->getBBox();
    if (this->_cachedMaxHeight.at((pos - AABBox.min()).xy()) == -1.f) {
        float minHeight = AABBox.min().z;
        float maxHeight = AABBox.max().z;

        this->_cachedMaxHeight.at((pos - AABBox.min()).xy()) = 0;
        for (float z = maxHeight; z > minHeight; z -= ImplicitPatch::zResolution) {
//            float eval = this->evaluate(pos.xy() + Vector3(0, 0, z));
            auto [totalEval, materials] = this->getMaterialsAndTotalEvaluation(pos.xy() + Vector3(0, 0, z));
            float totalPositive = 0.f, totalNegative = 0.f;
            for (auto& [mat, val] : materials) {
                if (isIn(mat, LayerBasedGrid::invisibleLayers)) {
                    totalNegative += val;
                } else {
                    totalPositive += val;
                }
            }
            if (totalEval >= ImplicitPatch::isovalue && totalPositive > totalNegative) {
                this->_cachedMaxHeight.at((pos - AABBox.min()).xy()) = z;
                break;
            }
        }
    }
    return this->_cachedMaxHeight.at((pos - AABBox.min()).xy());
}

float ImplicitPatch::getMinHeight(Vector3 pos)
{
    auto AABBox = this->getBBox();
    if (this->_cachedMinHeight.at((pos - AABBox.min()).xy()) == -1.f) {
        float minHeight = AABBox.min().z;
        float maxHeight = AABBox.max().z;

        this->_cachedMinHeight.at((pos - AABBox.min()).xy()) = 0;
        for (float z = minHeight; z < maxHeight; z += ImplicitPatch::zResolution) {
            float eval = this->evaluate(pos.xy() + Vector3(0, 0, z));
            if (eval >= ImplicitPatch::isovalue) {
                this->_cachedMinHeight.at((pos - AABBox.min()).xy()) = z;
                break;
            }
        }
    }
    return this->_cachedMinHeight.at((pos - AABBox.min()).xy());
}

float ImplicitPatch::getMinimalHeight(AABBox BBox)
{
    return this->getMinimalHeight(BBox.min(), BBox.max());
}

float ImplicitPatch::getMaximalHeight(AABBox BBox)
{
    return this->getMaximalHeight(BBox.min(), BBox.max());
}

float ImplicitPatch::getMinimalHeight(Vector3 minBox, Vector3 maxBox)
{
    auto BBox = this->getBBox();
    if (!minBox.isValid()) minBox = BBox.min();
    if (!maxBox.isValid()) maxBox = BBox.max();

    float minHeight = 10000.f;
    for (int x = minBox.x; x < maxBox.x; x++) {
        for (int y = minBox.y; y < maxBox.y; y++) {
            minHeight = std::min(minHeight, this->getMinHeight(Vector3(x, y)));
        }
    }
    return minHeight;
}

float ImplicitPatch::getMaximalHeight(Vector3 minBox, Vector3 maxBox)
{
    auto BBox = this->getBBox();
    minBox = Vector3::max(minBox, BBox.min());
    maxBox = Vector3::min(maxBox, BBox.max());

    float maxHeight = -10000.f;
    for (int x = minBox.x; x < maxBox.x; x++) {
        for (int y = minBox.y; y < maxBox.y; y++) {
            maxHeight = std::max(maxHeight, this->getMaxHeight(Vector3(x, y)));
        }
    }
    return maxHeight;
}

std::pair<float, std::map<TerrainTypes, float> > ImplicitPatch::getMaterialsAndTotalEvaluation(Vector3 pos)
{
    std::map<TerrainTypes, float> materials = this->getMaterials(pos);
    float totalValue = 0.f;

    for (const auto& [mat, val] : materials)
        totalValue += val;

    return {totalValue, materials};
}

Vector3 ImplicitPatch::getDimensions()
{
    auto AABBox = this->getBBox();
    return AABBox.dimensions();
}

Vector3 ImplicitPatch::getSupportDimensions()
{
    auto AABBox = this->getSupportBBox();
    return AABBox.dimensions();
}

Vector3 ImplicitPatch::getNormal(Vector3 pos)
{
    float delta = 1e-3; // Not too big, not too small...

    float minusX = this->evaluate(pos - Vector3(delta, 0, 0));
    float plusX = this->evaluate(pos + Vector3(delta, 0, 0));
    float minusY = this->evaluate(pos - Vector3(0, delta, 0));
    float plusY = this->evaluate(pos + Vector3(0, delta, 0));
    float minusZ = this->evaluate(pos - Vector3(0, 0, delta));
    float plusZ = this->evaluate(pos + Vector3(0, 0, delta));

    return Vector3(plusX - minusX, plusY - minusY, plusZ - minusZ).normalize();
}

void ImplicitPatch::setIndex(int newIndex)
{
    if (newIndex == -1)
        newIndex = ImplicitPatch::currentMaxIndex + 1;
    this->index = newIndex;
    ImplicitPatch::currentMaxIndex = std::max(this->index, ImplicitPatch::currentMaxIndex);
}

ImplicitPatch *ImplicitPatch::fromJson(nlohmann::json content)
{
    ImplicitPatch* result = nullptr;
    if (content.contains("type")) {
        std::string type = toUpper(content["type"]);
        if (type == "PRIMITIVE")
            result = ImplicitPrimitive::fromJson(content);
        else if (type == "COMPOSE" || type == "BINARY")
            result = ImplicitBinaryOperator::fromJson(content);
        else if (type == "UNARY")
            result = ImplicitUnaryOperator::fromJson(content);
        else if (type == "NARY")
            result = ImplicitNaryOperator::fromJson(content);
    } else if (content.contains("file")) {
        std::string filename = content["file"];
        if (filename.substr(filename.size() - 5) == ".json") {
            nlohmann::json newContent = nlohmann::json::parse(std::ifstream(filename));
            result = ImplicitPatch::fromJson(newContent[ImplicitPatch::json_identifier]);
        } else {
            Vector3 dimensions = json_to_vec3(content["dimensions"]);
            TerrainTypes material = materialFromString(content["material"]);
            result = ImplicitPrimitive::fromHeightmap(filename, dimensions);
            dynamic_cast<ImplicitPrimitive*>(result)->material = material;
        }
        result->used_json_filename = filename;
        return result;
    }
    if (result != nullptr) {
        result->updateCache();
    }
    return result;
}

void ImplicitPatch::updateCache()
{
    this->_cached = false;
    Vector3 dim = Vector3::max(Vector3(1, 1), this->getDimensions().roundedUp());

    this->_cachedMaxHeight = Matrix3<float>(dim.x, dim.y, 1, -1.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;

    this->_cachedMinHeight = Matrix3<float>(dim.x, dim.y, 1, -1.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

bool ImplicitPatch::checkIsInGround(Vector3 position)
{
    auto eval = this->getMaterials(position);
    float groundValue = 0.f;
    float outsideValue = 0.f;
    for (auto& [mat, val] : eval) {
        groundValue += (isIn(mat, LayerBasedGrid::invisibleLayers) ? 0.f : val);
        outsideValue += (isIn(mat, LayerBasedGrid::invisibleLayers) ? val : 0.f);
    }
//    std::cout << groundValue << " " << outsideValue << "\n";
//    if (groundValue > 0)
//        std::cout << groundValue << " " << outsideValue << std::endl;
    return groundValue > ImplicitPatch::isovalue && groundValue >= outsideValue;
}

Mesh ImplicitPatch::getGeometry(Vector3 dimensions)
{
//    Vector3 originalDimensions
    if (!dimensions.isValid())
        dimensions = this->getDimensions();
//    LayerBasedGrid layer(dimensions.x, dimensions.y, 0.f);
//    layer.add(this);
    VoxelGrid voxels;
    voxels.fromImplicit(this);
//    voxels.fromLayerBased(layer);
    return voxels.getGeometry(dimensions);
}

Vector3 ImplicitPatch::getIntersection(Vector3 origin, Vector3 dir, Vector3 minPos, Vector3 maxPos)
{
    auto myAABBox = this->getBBox();
    if (!minPos.isValid())
        minPos = myAABBox.min();
    if (!maxPos.isValid())
        maxPos = myAABBox.max();

    Vector3 currPos = origin;
    float distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
    float distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    // Continue while we are in the grid or we are heading towards the grid
    while((distanceToGrid < 0 || distanceToGridDT < 0) || distanceToGrid > distanceToGridDT)
    {
        if (Vector3::isInBox(currPos, minPos, maxPos)) {
            if (this->checkIsInGround(currPos)) {
                return currPos;
            }
        }
        currPos += dir;
        distanceToGrid = Vector3::signedManhattanDistanceToBoundaries(currPos, minPos, maxPos);
        distanceToGridDT = Vector3::signedManhattanDistanceToBoundaries(currPos + dir, minPos, maxPos);
    }
    return Vector3(false);

}

float ImplicitPatch::getHeight(float x, float y)
{
    return this->getHeight(Vector3(x, y, 0));
}

float ImplicitPatch::getHeight(Vector3 pos)
{
    return this->getMaxHeight(pos);
}

bool ImplicitPatch::contains(Vector3 v)
{
    auto BBox = this->getSupportBBox();
    return BBox.contains(v);
}

bool ImplicitPatch::contains(float x, float y, float z)
{
    return this->contains(Vector3(x, y, z));
}

Vector3 ImplicitPatch::getGlobalPositionOf(Vector3 posInsidePatch)
{
    Vector3 position = posInsidePatch + this->getBBox().min();
    ImplicitPatch* current = this;
    while(current->getParent()) {
        current = current->getParent();
        ImplicitUnaryOperator* asUnary = dynamic_cast<ImplicitUnaryOperator*>(current);
        if (asUnary) {
            position = asUnary->applyTransform(position);
        }
    }
    return position;
}

Vector3 ImplicitPatch::getLocalPositionOf(Vector3 globalPosition)
{
    Vector3 position = globalPosition;
    std::vector<ImplicitPatch*> ancestry;
    ImplicitPatch* current = this;
    while(current->getParent()) {
        current = current->getParent();
        ancestry.push_back(current);
    }
    std::reverse(ancestry.begin(), ancestry.end());
    for (size_t i = 0; i < ancestry.size(); i++) {
        ImplicitUnaryOperator* asUnary = dynamic_cast<ImplicitUnaryOperator*>(ancestry[i]);
        if (asUnary) {
            position = asUnary->inverseTransform(position); // asUnary->unwrapFunction(position);
        }
    }
    return position - this->getBBox().min();
}

Matrix3<float> ImplicitPatch::getVoxelized(Vector3 dimensions, Vector3 scale)
{
    if (_cached) {
        return this->_cachedVoxelized;
    }

    if (!dimensions.isValid())
        dimensions = this->getBBox().max();

    this->_cachedVoxelized = Matrix3<float>(dimensions * scale, -1.f); //LayerBasedGrid::densityFromMaterial(AIR));

    #pragma omp parallel for collapse(3)
    for (int x = 0; x < _cachedVoxelized.sizeX; x++) {
        for (int y = 0; y < _cachedVoxelized.sizeY; y++) {
            for (int z = 0; z < _cachedVoxelized.sizeZ; z++) {
                Vector3 pos = Vector3(x, y, z) * scale;
                std::vector<Vector3> evaluationPoses = {
                    pos + Vector3(0, 0, 0)/*,
                    pos + Vector3(0, 0, 1),
                    pos + Vector3(0, 1, 0),
                    pos + Vector3(0, 1, 1),
                    pos + Vector3(1, 0, 0),
                    pos + Vector3(1, 0, 1),
                    pos + Vector3(1, 1, 0),
                    pos + Vector3(1, 1, 1)*/
                };
                float allPositiveVals = 0.f;
                float allNegativeVals = 0.f;
                std::vector<std::pair<float, std::map<TerrainTypes, float> > > evaluations(evaluationPoses.size());
                for (size_t i = 0; i < evaluations.size(); i++) {
                    auto [totalEval, matVals] = this->getMaterialsAndTotalEvaluation(evaluationPoses[i]);
                    evaluations[i] = {totalEval, matVals};


                    float positiveVal = 0.f;
                    float negativeVal = 0.f;
                    for (const auto& [mat, val] : matVals) {
                        float density = LayerBasedGrid::densityFromMaterial(mat);
                        positiveVal += (density > 0.f ? val : 0.f);
                        negativeVal += (density <= 0.f ? val : 0.f);
                    }
                    totalEval = positiveVal;
                    allPositiveVals += positiveVal;
                    allNegativeVals += negativeVal;
                }

                allPositiveVals = (allPositiveVals - allNegativeVals) / float(evaluations.size());
                if (allPositiveVals >= ImplicitPatch::isovalue) {
                    /*TerrainTypes maxType;
                    float maxVal = 0.f;
                    for (const auto& [mat, val] : matVals) {
                        if (val > maxVal) {
                            maxType = mat;
                            maxVal = val;
                        }
                    }
                    _cachedVoxelized.at(pos) = (isIn(maxType, LayerBasedGrid::invisibleLayers) ? -totalEval : totalEval); // LayerBasedGrid::densityFromMaterial(maxType);
                    */
                    _cachedVoxelized.at(pos) = allPositiveVals;
                }
            }
        }
    }
//    this->_cachedVoxelized.meanSmooth(3, 3, 3, true);
    _cached = true;
//    auto mini = _cachedVoxelized.resize(Vector3(10, 10, 10));
    return _cachedVoxelized;
}

ImplicitPrimitive *ImplicitPatch::createPredefinedShape(PredefinedShapes shape, Vector3 dimensions, float additionalParam, BSpline parametricCurve)
{
    ImplicitPrimitive* primitive = new ImplicitPrimitive();
    primitive->predefinedShape = shape;
    primitive->evalFunction = ImplicitPatch::createPredefinedShapeFunction(shape, dimensions, additionalParam, parametricCurve);
    primitive->optionalCurve = parametricCurve;
    primitive->parametersProvided = {additionalParam};
    primitive->dimensions = dimensions;

    return primitive;
}
ImplicitPatch *ImplicitPatch::createIdentity()
{
    ImplicitPrimitive* primitive = new ImplicitPrimitive();
    primitive->predefinedShape = ImplicitPatch::PredefinedShapes::None;
    primitive->evalFunction = [](Vector3) { return 0.f; };
    primitive->optionalCurve = {};
    primitive->parametersProvided = {0.f};

    return primitive;
}

std::function<float (Vector3)> ImplicitPatch::createPredefinedShapeFunction(PredefinedShapes shape, Vector3 dimensions, float additionalParam, BSpline parametricCurve)
{
    std::function<float(Vector3)> func;
    switch(shape) {
    case Sphere:
        func = ImplicitPatch::createSphereFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Block:
        func = ImplicitPatch::createBlockFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Gaussian:
        func = ImplicitPatch::createGaussianFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Rock:
        func = ImplicitPatch::createRockFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Mountain:
        func = ImplicitPatch::createMountainFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Dune:
        func = ImplicitPatch::createDuneFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Basin:
        func = ImplicitPatch::createBasinFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Cave:
        func = ImplicitPatch::createCaveFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Arch:
        func = ImplicitPatch::createArchFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Noise2D:
        func = ImplicitPatch::createNoise2DFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Cylinder:
        func = ImplicitPatch::createCylinderFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case MountainChain:
        func = ImplicitPatch::createMountainChainFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z, parametricCurve);
        break;
    case Polygon:
        func = ImplicitPatch::createPolygonFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z, parametricCurve);
        break;
    case ImplicitHeightmap:
        // Do it yourself!
        break;
    case ParametricTunnel:
        func = ImplicitPatch::createParametricTunnelFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z, parametricCurve);
        break;
    case None:
        func = ImplicitPatch::createIdentityFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    case Ripple:
        func = ImplicitPatch::createRippleFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    }
    return func;
}




ImplicitPrimitive::ImplicitPrimitive()
{
    this->evalFunction = [](Vector3) -> float { return 0.f; };
    this->parametersProvided = {-1.f}; // Just to have an element at first
}

float ImplicitPrimitive::evaluate(Vector3 pos)
{
    auto [minPos, maxPos] = this->getBBox();
    auto [minSupportPos, maxSupportPos] = this->getSupportBBox();
    if (!Vector3::isInBox(pos, minSupportPos, maxSupportPos))
        return 0.f;
    Vector3 supportWidth = maxSupportPos - minSupportPos;
    Vector3 width = maxPos - minPos;
//    if (this->optionalCurve.points.size() > 2) {
//        ShapeCurve area = this->optionalCurve.points;
//        if (!area.inside(pos)) {
//            return 0.f;
//        }
//    }

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
            evalPos.x = (evalPos.x < width.x ? evalPos.x : 2.f * width.x - evalPos.x);
            evalPos.y = (evalPos.y < width.y ? evalPos.y : 2.f * width.y - evalPos.y);
            evalPos.z = (evalPos.z < width.z ? evalPos.z : 2.f * width.z - evalPos.z);
            evaluation = std::clamp(this->evalFunction(evalPos), 0.f, 1.f) * distanceFalloff * .5f;
        }
    }

    evaluation = std::clamp(evaluation, 0.f, 1.f);

    return evaluation;
}

std::map<TerrainTypes, float> ImplicitPrimitive::getMaterials(Vector3 pos)
{
    return {{this->material, this->evaluate(pos)}};
}

AABBox ImplicitPrimitive::getSupportBBox()
{
    Vector3 margin = (this->supportDimensions - this->dimensions) * .5f;
    if (!this->supportDimensions.isValid() || this->supportDimensions == Vector3())
        margin = this->dimensions * 1.f; // Nice big margin
    return {this->position - margin, this->position + this->dimensions + margin};
}

AABBox ImplicitPrimitive::getBBox()
{
    return {this->position, this->position + this->dimensions};
}

void ImplicitPrimitive::update()
{

    if (this->predefinedShape == ImplicitHeightmap) {
        if (this->cachedHeightmap.size() == 0) {
            this->cachedHeightmap = Heightmap(this->heightmapFilename, this->dimensions.x, this->dimensions.y).heights;
        }
        if (this->cachedHeightmap.height() > 1) {
            this->evalFunction = [=](Vector3 pos) -> float {
                return (this->getBBox().contains(pos) ? clamp(this->cachedHeightmap.interpolate(pos) + ImplicitPatch::isovalue, 0.f, 1.f) : 0.f);
            };
        } else {
            if (this->cachedHeightmap.getDimensions().xy() != this->getDimensions().xy())
                this->cachedHeightmap.resize(this->getDimensions().xy() + Vector3(0, 0, 1.f));
            this->cachedHeightmap = this->cachedHeightmap.normalize() * this->getDimensions().z;
    //        auto cacheCopy = cachedHeightmap;
            this->evalFunction = ImplicitPrimitive::convert2DfunctionTo3Dfunction([=](Vector3 pos) -> float {
    //            auto heightmap = cachedHeightmap;
                auto valueAt = this->cachedHeightmap.interpolate(pos.xy());
                return valueAt;
            });
        }
    } else {
        this->evalFunction = ImplicitPatch::createPredefinedShapeFunction(this->predefinedShape, this->dimensions, this->parametersProvided[0], this->optionalCurve);
    }
}

std::string ImplicitPrimitive::toString()
{
    if (this->name == "") {
        this->name = stringFromPredefinedShape(this->predefinedShape);
    }
    return this->name + ": Primitive " + this->name + " #" + std::to_string(this->index);
}

nlohmann::json ImplicitPrimitive::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        content["file"] = this->used_json_filename;
        content["position"] = vec3_to_json(this->position);
        content["dimensions"] = vec3_to_json(this->dimensions);
        content["material"] = stringFromMaterial(this->material);
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = vec3_to_json(patchOffset);
        content["scale"] = vec3_to_json(patchRescale);*/
    } else {
        content["name"] = this->name;
        content["type"] = "primitive";
        content["position"] = vec3_to_json(this->position);
        content["dimensions"] = vec3_to_json(this->dimensions);
        content["supportDimensions"] = vec3_to_json(this->supportDimensions);
        content["densityValue"] = LayerBasedGrid::densityFromMaterial(this->material);
        content["material"] = stringFromMaterial(this->material);
        content["index"] = this->index;
        content["shape"] = stringFromPredefinedShape(this->predefinedShape);
        content["sigmaValue"] = this->parametersProvided[0];
        content["optionalCurve"] = bspline_to_json(this->optionalCurve);
        content["mirrored"] = this->mirrored;
        content["filename"] = this->heightmapFilename;
    }

    return content;
}

ImplicitPatch *ImplicitPrimitive::fromJson(nlohmann::json content)
{
    ImplicitPrimitive* patch = new ImplicitPrimitive;
    patch->name = content["name"];
    patch->position = json_to_vec3(content["position"]);
    patch->setDimensions(json_to_vec3(content["dimensions"]));
    patch->supportDimensions = json_to_vec3(content["supportDimensions"]);
    patch->material = LayerBasedGrid::materialFromDensity(content["densityValue"]);
    if (content.contains("material"))
        patch->material = materialFromString(content["material"]);
    if (content.contains("optionalCurve"))
        patch->optionalCurve = json_to_bspline(content["optionalCurve"]);
    patch->index = content["index"];
    patch->predefinedShape = predefinedShapeFromString(content["shape"]);
    patch->parametersProvided = {content["sigmaValue"]};
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    if (content.contains("filename"))
        patch->heightmapFilename = content["filename"];
    patch->update();
    return patch;
}

void ImplicitPrimitive::setDimensions(Vector3 newDimensions)
{
    this->dimensions = newDimensions;

    this->updateCache();
}

void ImplicitPrimitive::setSupportDimensions(Vector3 newSupportDimensions)
{
    this->supportDimensions = newSupportDimensions;
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
    ImplicitPrimitive* copy = new ImplicitPrimitive(*this);
    return copy;
}

ImplicitPrimitive *ImplicitPrimitive::fromHeightmap(std::string filename, Vector3 dimensions, ImplicitPrimitive *prim)
{
    std::string ext = toUpper(getExtension(filename));
    if (isIn(ext, {"JPG", "PNG", "TGA", "BMP", "PSD", "GIF", "HDR", "PIC"})) {
        int imgW, imgH, nbChannels;
        unsigned char *data = stbi_load(filename.c_str(), &imgW, &imgH, &nbChannels, STBI_grey); // Load image, force 1 channel
        if (data == NULL)
        {
            std::cerr << "Error : Impossible to load " << filename << "\n";
            std::cerr << "Either file is not found, or type is incorrect. Available file types are : \n";
            std::cerr << "\t- JPG, \n\t- PNG, \n\t- TGA, \n\t- BMP, \n\t- PSD, \n\t- GIF, \n\t- HDR, \n\t- PIC";
            exit (-1);
            return nullptr;
        }
        Matrix3<float> map(imgW, imgH);
        for (int x = 0; x < imgW; x++) {
            for (int y = 0; y < imgH; y++) {
                float value = (float)data[x + y * imgW];
                map.at(x, y) = value;
            }
        }
        stbi_image_free(data);

        if (dimensions.isValid()) {
            map = map.resize(dimensions.x, dimensions.y, 1);
            map = map.normalize() * dimensions.z;
        }
        return ImplicitPrimitive::fromHeightmap(map, filename,  prim);
    } else {
        VoxelGrid voxels = VoxelGrid();
        voxels.retrieveMap(filename);
        Matrix3<float> map = voxels.getVoxelValues();
        return ImplicitPrimitive::fromHeightmap(map, filename, prim);
    }
}

ImplicitPrimitive *ImplicitPrimitive::fromHeightmap(Matrix3<float> heightmap, std::string filename, ImplicitPrimitive *prim)
{
    if (prim == nullptr)
        prim = new ImplicitPrimitive;
    if (heightmap.height() > 1)
        prim->setDimensions(heightmap.getDimensions());
    else
        prim->setDimensions(heightmap.getDimensions().xy() + Vector3(0, 0, heightmap.max()));
    prim->position = Vector3(0, 0, 0);
//    prim->evalFunction = ImplicitPrimitive::convert2DfunctionTo3Dfunction([=](Vector3 pos) -> float {
//        return heightmap.data[heightmap.getIndex(pos.xy())];
//    });

    prim->name = "Heightmap";
    prim->predefinedShape = ImplicitHeightmap;
    prim->cachedHeightmap = heightmap;
    prim->cachedHeightmap.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::DEFAULT_VALUE;
    prim->cachedHeightmap.raiseErrorOnBadCoord = false;
    prim->material = ROCK;
    prim->update();
    return prim;
}






ImplicitBinaryOperator::ImplicitBinaryOperator()
{
    this->composables = std::vector<ImplicitPatch*>(2);
}

float ImplicitBinaryOperator::evaluate(Vector3 pos)
{
    // Get the function, depending on the chosen operator
    float evalA = this->evaluateA(pos);
    float evalB = this->evaluateB(pos);

    return this->evaluateFromAandB(evalA, evalB);
    /*
    if (this->withIntersectionOnB && evalA < ImplicitPatch::isovalue) {
        return evalA; // TODO: Check if we need to use "isovalue" in the condition
    }

    float evaluation = 0.f;
    if (this->composeFunction == STACK) {
        evaluation = std::max(evalA, evalB);
    } else if (this->composeFunction == REPLACE) {
        evaluation = (evalB > ImplicitPatch::isovalue ? evalB : evalA);
    } else if (this->composeFunction == BLEND) {
        evaluation = std::pow(std::pow(evalA, this->blendingFactor) + std::pow(evalB, this->blendingFactor), 1.f/this->blendingFactor);
    } else {
        std::cerr << "Wrong composition function" << std::endl;
    }
    evaluation = std::clamp(evaluation, 0.f, 1.f);
    return evaluation;
    */
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterials(Vector3 pos)
{
    // Get materials on A, get materials on B
    // Take the operator final evaluation
    // Use the ratio "mat on A / final eval" and "mat on B / final eval"
    // to determine final materials
    auto [contributionA, materialsAndEvalsForA] = this->getMaterialsAndTotalEvaluationA(pos);
    auto [contributionB, materialsAndEvalsForB] = this->getMaterialsAndTotalEvaluationB(pos);
    float operationEvaluation = this->evaluateFromAandB(contributionA, contributionB);

    if (this->withIntersectionOnB && contributionA < ImplicitPatch::isovalue) {
        return materialsAndEvalsForA;
    }

    std::map<TerrainTypes, float> result;


    if (operationEvaluation == 0) {
        return {};
    }
/*
    // Contribution A and B are [0, 1]
    // op eval is [0, 1]
    // The contribution of A and B are the gradients of the composite function (?)
    // Or is it "partial derivative" ? Don't know the difference...
    float delta = 1e-3;
    float AminusDelta = std::clamp(contributionA - delta, 0.f, 1.f);
    float AplusDelta = std::clamp(contributionA + delta, 0.f, 1.f);
    float BminusDelta = std::clamp(contributionB - delta, 0.f, 1.f);
    float BplusDelta = std::clamp(contributionB + delta, 0.f, 1.f);
    float deltaA = (AplusDelta - AminusDelta);
    float deltaB = (BplusDelta - BminusDelta);
    float gradX0 = evaluateFromAandB(AplusDelta, contributionB);
    float gradX1 = evaluateFromAandB(AminusDelta, contributionB);
    float gradY0 = evaluateFromAandB(contributionA, BplusDelta);
    float gradY1 = evaluateFromAandB(contributionA, BminusDelta);
    float varyX = (gradX0 - gradX1);
    float varyY = (gradY0 - gradY1);
    float gradX = varyX / deltaA;
    float gradY = varyY / deltaB;


    // No sense to put any material if gradient is null
//    if (gradX <= 0.f && gradY <= 0.f) {
//        return {{TerrainMaterial(), 0.f}};
//    }
    if (gradX <= 0.f && gradY <= 0.f) {
        gradX = 1.f; //return {{TerrainMaterial(), 0.f}};
        gradY = 1.f;
    }

    // Force 0 on null materials
    float ratioA = 0.f;
    float ratioB = 0.f;
    if (contributionA <= 0) {
        ratioA = 0.f;
        ratioB = (operationEvaluation / contributionB);
    } else if (contributionB <= 0) {
        ratioA = (operationEvaluation / contributionA);
        ratioB = 0.f;
    } else {
        ratioA = (gradX / (gradX + gradY)) * (contributionA > 0 ? (operationEvaluation / contributionA) : 0.f);
        ratioB = (gradY / (gradX + gradY)) * (contributionB > 0 ? (operationEvaluation / contributionB) : 0.f);
    }
*/
    float mu_A = contributionA;
    float mu_B = contributionB;
    float rA = mu_A / (mu_A + mu_B);
    float rB = mu_B / (mu_A + mu_B);
    float tempA = rA * mu_A;
    float tempB = rB * mu_B;
    float globalCoef = operationEvaluation / (tempA + tempB);
    float ratioA = rA * globalCoef;
    float ratioB = rB * globalCoef;
/*
    Vector3 gradient = Vector3(gradX, gradY);
//                evaluateFromAandB(contributionA + delta, contributionB) - evaluateFromAandB(contributionA - delta, contributionB),
//                evaluateFromAandB(contributionA, contributionB + delta) - evaluateFromAandB(contributionA, contributionB - delta)
//                );
//    float sumOfContributions = gradient.x + gradient.y; //contributionA + contributionB;
    float ratioA = gradient.x; // / operationEvaluation;
    float ratioB = gradient.y; // / operationEvaluation;
    float sumOfRatios = ratioA + ratioB;
    if (sumOfRatios == 0) {
        return {{TerrainMaterial(), 0.f}};
    }
    // Get ratioA + ratioB = 1
    ratioA = ratioA / sumOfRatios;
    ratioB = ratioB / sumOfRatios;
    // Ratio depends on self evaluation
    if (contributionA > 0)
        ratioA = operationEvaluation * ratioA / contributionA;
    else
        ratioA = 0.f;

    if (contributionB > 0)
        ratioB = operationEvaluation * ratioB / contributionB;
    else
        ratioB = 0.f; */

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

float ImplicitBinaryOperator::evaluateFromAandB(float evalA, float evalB)
{
    if (this->withIntersectionOnB && evalA < ImplicitPatch::isovalue) {
        return evalA; // TODO: Check if we need to use "isovalue" in the condition
    }
    float evaluation = 0.f;
    if (this->composeFunction == STACK) {
        evaluation = std::max(evalA, evalB);
    } else if (this->composeFunction == REPLACE) {
        evaluation = (evalB > ImplicitPatch::isovalue ? evalB : evalA);
    } else if (this->composeFunction == BLEND || this->composeFunction == ONE_SIDE_BLEND) {
        evaluation = std::pow(std::pow(evalA, this->blendingFactor) + std::pow(evalB, this->blendingFactor), 1.f/this->blendingFactor);
    } else {
        std::cerr << "Wrong composition function" << std::endl;
    }
    evaluation = std::clamp(evaluation, 0.f, 1.f);
    return evaluation;
}

float ImplicitBinaryOperator::evaluateA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);

    return this->composableA()->evaluate(evaluationPosA);
}

float ImplicitBinaryOperator::evaluateB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);

    return this->composableB()->evaluate(evaluationPosB);
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterialsA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return this->composableA()->getMaterials(evaluationPosA);
}

std::map<TerrainTypes, float> ImplicitBinaryOperator::getMaterialsB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return this->composableB()->getMaterials(evaluationPosB);
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitBinaryOperator::getMaterialsAndTotalEvaluationA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return this->composableA()->getMaterialsAndTotalEvaluation(evaluationPosA);
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitBinaryOperator::getMaterialsAndTotalEvaluationB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return this->composableB()->getMaterialsAndTotalEvaluation(evaluationPosB);
}

bool ImplicitBinaryOperator::contains(PredefinedShapes shape)
{
    return this->composableA()->contains(shape) || this->composableB()->contains(shape);
}

std::vector<ImplicitPatch *> ImplicitBinaryOperator::findAll(PredefinedShapes shape)
{
    return vectorMerge(this->composableA()->findAll(shape), this->composableB()->findAll(shape));
}

AABBox ImplicitBinaryOperator::getSupportBBox()
{
    auto AABBoxA = this->composableA()->getSupportBBox();
    if (this->withIntersectionOnB) { // No need to go further, we know the limit with the intersection
        return AABBoxA;
    }
    auto AABBoxB = this->composableB()->getSupportBBox();

    if (this->positionalB == PositionalLabel::ABOVE) {
        // If stacked on composable A, composable B can get higher
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_BOTTOM) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_TOP) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::FIXED_POS) {
        // Nothing to do there
    }

    return {Vector3::min(AABBoxA.min(), AABBoxB.min()), Vector3::max(AABBoxA.max(), AABBoxB.max())};
}

AABBox ImplicitBinaryOperator::getBBox()
{
    auto AABBoxA = this->composableA()->getBBox();
    if (this->withIntersectionOnB) { // No need to go further, we know the limit with the intersection
        return AABBoxA;
    }
    auto AABBoxB = this->composableB()->getBBox();

    if (this->positionalB == PositionalLabel::ABOVE) {
        // If stacked on composable A, composable B can get higher
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_BOTTOM) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_TOP) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.dimensions().z;
        AABBoxB.maxi.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::FIXED_POS) {
        // Nothing to do there
    }

    return {Vector3::min(AABBoxA.min(), AABBoxB.min()), Vector3::max(AABBoxA.max(), AABBoxB.max())};
}

void ImplicitBinaryOperator::update()
{
    // Guess we have nothing to do...
}

std::string ImplicitBinaryOperator::toString()
{
    return this->name + ": Operation between " + (composableA() ? "#" + std::to_string(composableA()->index) : "undefined") + " and " + (composableB() ? "#" + std::to_string(composableB()->index) : "undefined");
}

nlohmann::json ImplicitBinaryOperator::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        content["file"] = this->used_json_filename;
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = vec3_to_json(patchOffset);
        content["scale"] = vec3_to_json(patchRescale);*/
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
        content["optionalCurve"] = bspline_to_json(this->optionalCurve);
        content["mirrored"] = this->mirrored;
    }

    return content;
}

ImplicitPatch *ImplicitBinaryOperator::fromJson(nlohmann::json content)
{
    ImplicitBinaryOperator* patch = new ImplicitBinaryOperator;
    patch->addChild(ImplicitPatch::fromJson(content["composableA"]), 0);//patch->composableA() = ImplicitPatch::fromJson(content["composableA"]);
    patch->addChild(ImplicitPatch::fromJson(content["composableB"]), 1);//patch->composableB() = ImplicitPatch::fromJson(content["composableB"]);
    patch->name = content["name"];
    patch->composeFunction = compositionOperationFromString(content["operator"]);
    patch->blendingFactor = content["blendingFactor"];
    patch->positionalB = positionalLabelFromString(content["positionalB"]);
    patch->index = content["index"];
    if (content.contains("useIntersection"))
        patch->withIntersectionOnB = content["useIntersection"];
    if (content.contains("optionalCurve"))
        patch->optionalCurve = json_to_bspline(content["optionalCurve"]);
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    patch->update();
    return patch;
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
    if (this->composableB() != nullptr) {
        std::swap(this->composableA(), this->composableB());
        this->composableA()->updateCache();
        this->composableB()->updateCache();
    }
}

void ImplicitBinaryOperator::deleteAllChildren()
{
    for (auto& child : composables) {
        delete child;
    }
    this->composables.resize(2, nullptr);
}

Vector3 ImplicitBinaryOperator::getEvaluationPositionForComposableA(Vector3 pos)
{
    return pos; // Nothing to do
}

Vector3 ImplicitBinaryOperator::getEvaluationPositionForComposableB(Vector3 pos)
{
    // Get the correct evaluation position for the B composent
    float offsetB = 0.f;
    if (this->positionalB == ABOVE) {
        offsetB = this->composableA()->getMaxHeight(pos);
    } else if (this->positionalB == INSIDE_TOP) {
        offsetB = this->composableA()->getMaxHeight(pos) - this->composableB()->getDimensions().z;
    } else if (this->positionalB == INSIDE_BOTTOM) {
        offsetB = this->composableA()->getMinHeight(pos);
    } else if (this->positionalB == FIXED_POS) {
        offsetB = 0.f;
    } else if (this->positionalB == SMOOTH_ABOVE) {
        float heightB = this->composableB()->getDimensions().z;
        float heightA = heightA = std::min(heightB, this->composableA()->getMaxHeight(pos));
        float minHeightA = this->composableA()->getMinimalHeight(this->composableB()->getBBox());
        float heightInterp = interpolation::wyvill(heightA, minHeightA, (minHeightA + heightB));
//        offsetB = (1.f - heightInterp / heightB) * heightB;
        pos.z *= (heightInterp / heightB);
    } else {
        std::cerr << "Wrong position label" << std::endl;
    }
    pos.z -= offsetB;
    return pos;
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

ImplicitPatch*& ImplicitBinaryOperator::composableA()
{
    return this->composables[0];
}
ImplicitPatch*& ImplicitBinaryOperator::composableB()
{
    return this->composables[1];
}





ImplicitUnaryOperator::ImplicitUnaryOperator()
{
    this->composables = std::vector<ImplicitPatch*>(1);
    this->noiseFunction = [](Vector3) { return 0.f; };
}
/*
    this->wrapFunction = [=](Vector3 pos) {
        for (int i = 0; i < this->transforms.size(); i++)
            pos = this->transforms[i].wrap(pos);
        return pos;
    };
    this->unwrapFunction = [=](Vector3 pos) {
        for (int i = this->transforms.size() - 1; i >= 0; i--)
            pos = this->transforms[i].unwrap(pos);
        return pos;
    };
    this->noiseFunction = [](Vector3) { return 0.f; };
}
*/
float ImplicitUnaryOperator::evaluate(Vector3 pos)
{
//    Vector3 evaluationPos = this->unwrapFunction(pos);
    Vector3 evaluationPos = inverseTransform(pos);
    float evaluation = this->composableA()->evaluate(evaluationPos) + this->noiseFunction(pos);
    evaluation = std::clamp(evaluation, 0.f, 1.f);

    return evaluation;
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
}

void ImplicitUnaryOperator::addRotation(Vector3 rotationAngles) {
    this->_rotation += rotationAngles;
    Matrix rotationMatrix = rotationAngles.toRotationMatrix();
    auto rotationTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + Vector3::fromMatrix(rotationMatrix.product((point - center).toMatrix()));
    };
    auto inverseRotationTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + Vector3::fromMatrix(rotationMatrix.inverse().product((point - center).toMatrix()));
    };
    addTransformation(rotationTransform, inverseRotationTransform);
}

void ImplicitUnaryOperator::addScaling(const Vector3 &scaleFactors) {
    this->_scale *= scaleFactors;
    auto scalingTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + (point - center) * scaleFactors; //Vector3(point.x * scaleFactors.x, point.y * scaleFactors.y, point.z * scaleFactors.z);
    };
    auto inverseScalingTransform = [=](const Vector3& point) {
        Vector3 center = this->composableA()->getBBox().center();
        return center + (point - center) / scaleFactors; // Vector3(point.x / scaleFactors.x, point.y / scaleFactors.y, point.z / scaleFactors.z);
    };
    addTransformation(scalingTransform, inverseScalingTransform);
}

void ImplicitUnaryOperator::addDisplacementField(const Matrix3<Vector3> &displacementField) {
    auto displacementTransform = [&displacementField](const Vector3& point) {
        // Use the displacement field to transform the point
        return point + displacementField.interpolate(point);
    };
    auto inverseDisplacementTransform = [&displacementField](const Vector3& point) {
        // Iteratively "walk" back to the original point
        Vector3 oldPoint;
        Vector3 newPoint = point;
        const int maxIterations = 100;  // Define an appropriate limit for your use case
        int iterations = 0;
        do {
            if (++iterations > maxIterations) {
                throw std::runtime_error("Inverse transformation could not be found");
            }
            oldPoint = newPoint;
            newPoint = oldPoint - displacementField.interpolate(oldPoint);
        } while ((newPoint - oldPoint).length() > 1e-6); // Use a small threshold to stop the iteration
        return newPoint;
    };
    addTransformation(displacementTransform, inverseDisplacementTransform);
}

Vector3 ImplicitUnaryOperator::applyTransform(Vector3 pos) const
{
    for (int i = 0; i < transformations.size(); i++)
        pos = transformations[i].first(pos);
    return pos;
}

Vector3 ImplicitUnaryOperator::inverseTransform(Vector3 pos) const
{
    for (int i = transformations.size() - 1; i >= 0; i--)
        pos = transformations[i].second(pos);
    return pos;
}

std::map<TerrainTypes, float> ImplicitUnaryOperator::getMaterials(Vector3 pos)
{
    Vector3 evaluationPos = this->inverseTransform(pos);
    auto [eval, materials] = this->composableA()->getMaterialsAndTotalEvaluation(evaluationPos);
            if (eval > 0.f) {
        float noiseValue = this->noiseFunction(pos);
        for (auto& [mat, val] : materials) {
            val += (noiseValue * (val / eval));
        }
    }
    return materials;
}
AABBox ImplicitUnaryOperator::getSupportBBox() {
    auto AABBox = this->composableA()->getSupportBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.min(), AABBox.max());

    for (auto& vert : vertices) { // call transform function
        vert = inverseTransform(vert);
    }

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}

AABBox ImplicitUnaryOperator::getBBox()
{
    auto AABBox = this->composableA()->getBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.min(), AABBox.max());

//    std::cout << "Transformations : " << this->transformations.size() << std::endl;

    for (auto& vert : vertices) { // call transform function
        vert = applyTransform(vert);
    }

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}
/*
AABBox ImplicitUnaryOperator::getSupportBBox()
{
    auto AABBox = this->composableA()->getSupportBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.min(), AABBox.max());

    for (auto& vert : vertices) // call transform function
        vert = this->wrapFunction(vert);

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}

AABBox ImplicitUnaryOperator::getBBox()
{
    auto AABBox = this->composableA()->getBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.min(), AABBox.max());

    for (auto& vert : vertices) // call transform function
        vert = this->wrapFunction(vert);

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}*/

std::string ImplicitUnaryOperator::toString()
{
    return this->name + ": Unary operation on " + (composableA() ? "#" + std::to_string(composableA()->index) : "undefined");
}

nlohmann::json ImplicitUnaryOperator::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        content["file"] = this->used_json_filename;
        /*std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = vec3_to_json(patchOffset);
        content["scale"] = vec3_to_json(patchRescale);*/
    } else {
        content["type"] = "unary";
        content["name"] = this->name;
        content["index"] = this->index;
        if (this->composableA())
            content["composableA"] = this->composableA()->toJson();
        content["optionalCurve"] = bspline_to_json(this->optionalCurve);
        content["spreadFactor"] = this->_spreadingFactor;
        content["mirrored"] = this->mirrored;
    }

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
    if (content.contains("translation") && json_to_vec3(content["translation"]) != Vector3()) {
        patch = new ImplicitTranslation();
        translation = json_to_vec3(content["translation"]);
    } else if (content.contains("scale") && json_to_vec3(content["scale"]) != Vector3(1, 1, 1)) {
        patch =  new ImplicitScaling();
        scale = json_to_vec3(content["scale"]);
        if (scale == Vector3())
            scale = Vector3(1.f, 1.f, 1.f);
    } else if (content.contains("rotation") && json_to_vec3(content["rotation"]) != Vector3()) {
        patch =  new ImplicitRotation();
        rotation = json_to_vec3(content["rotation"]);
    } else if (content.contains("noise") && json_to_vec3(content["noise"]) != Vector3()) {
        patch = new ImplicitNoise();
        noise = json_to_vec3(content["noise"]);
    } else if (content.contains("distortion") && json_to_vec3(content["distortion"]) != Vector3()) {
        patch = new ImplicitWraping();
        distortion = json_to_vec3(content["distortion"]);
    } else if (content.contains("spreadFactor") && content["spreadFactor"] != 0) {
        patch = new ImplicitSpread();
        patch->spread(content["spreadFactor"]);
    } else {
        patch = new ImplicitUnaryOperator();
    }
    patch->addChild(ImplicitPatch::fromJson(content["composableA"]));
    patch->name = content["name"];
    patch->index = content["index"];

    if (translation != Vector3()) patch->translate(translation);
    if (scale != Vector3(1.f, 1.f, 1.f)) patch->scale(scale);
    if (rotation != Vector3()) patch->rotate(rotation);
    if (noise != Vector3()) patch->addRandomNoise(noise.x, noise.y, noise.z);
//    if (distortion != Vector3()) patch->addRandomWrap(distortion.x, distortion.y, distortion.z);
    if (distortion != Vector3()) {
        patch->addWrapFunction(Matrix3<Vector3>({
                                                       { Vector3(0.5, 0, 0), Vector3(0.7, -.2, 0), Vector3(1, -1, 0) },
                                                       { Vector3(0.0, 0, 0), Vector3(0.7, 0, 0), Vector3(2, 0, 0) },
                                                       { Vector3(0.5, 0, 0), Vector3(0.7, .2, 0), Vector3(1, 1, 0) }
         }));
    }

    if (content.contains("optionalCurve"))
        patch->optionalCurve = json_to_bspline(content["optionalCurve"]);
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];

    patch->update();
    return patch;
}

bool ImplicitUnaryOperator::contains(PredefinedShapes shape)
{
    return this->composableA()->contains(shape);
}

std::vector<ImplicitPatch *> ImplicitUnaryOperator::findAll(PredefinedShapes shape)
{
    return this->composableA()->findAll(shape);
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
    this->composables.resize(1, nullptr);
}

ImplicitPatch *&ImplicitUnaryOperator::composableA()
{
    return this->composables[0];
}

void ImplicitUnaryOperator::translate(Vector3 translation)
{
    return this->addTranslation(translation);
    /*this->_translation += translation;
    this->transforms.push_back(UnaryOpTranslate(translation));*/
}

void ImplicitUnaryOperator::rotate(Vector3 angles) //float angleX, float angleY, float angleZ)
{
    return this->addRotation(angles); // Vector3(angleX, angleY, angleZ));
    /*
    this->_rotation += Vector3(angleX, angleY, angleZ);

    auto AABBox = this->getBBox();
    Vector3 center = AABBox.center();
    this->transforms.push_back(UnaryOpRotate(Vector3(angleX, angleY, angleZ), center));
    */
}

void ImplicitUnaryOperator::scale(Vector3 scaleFactor)
{
    return this->addScaling(scaleFactor);
    /*
    this->_scale *= scaleFactor;
    auto BBox = this->getBBox();
    Vector3 center = BBox.center();
    this->transforms.push_back(UnaryOpScale(scaleFactor, center));
    */
}

void ImplicitUnaryOperator::addRandomNoise(float amplitude, float period, float offset)
{
    this->_noise = Vector3(amplitude, period, offset);
    this->noiseFunction = [=](Vector3 pos) -> float {
        FastNoiseLite noise;
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseVal = noise.GetNoise(pos.x * period + offset, pos.y * period + offset, pos.z * period + offset) * 1.f;
        return noiseVal * amplitude;
    };
}

void ImplicitUnaryOperator::addRandomWrap(float amplitude, float period, float offset)
{
    this->_distortion = Vector3(amplitude, period, offset);
    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    this->transforms.push_back(UnaryOpWrap(noise, Vector3(amplitude, period, offset)));
}

void ImplicitUnaryOperator::addWrapFunction(Matrix3<Vector3> func)
{
//    AABBox dims = this->getSupportBBox();
    this->transforms.push_back(UnaryOpWrap([=](Vector3 pos) -> Vector3 {
        auto f = func;
        auto supportBBox = this->composableA()->getSupportBBox();
        Vector3 normalizedPos = supportBBox.normalize(pos);
        Vector3 posInMatrix = normalizedPos * (f.getDimensions() - Vector3(1, 1, 1));
        Vector3 distortionValue = f.interpolate(posInMatrix);
//        if (pos.z == 0) //std::abs(normalizedPos.xy().maxComp() - 2.f/3.f) < 0.01f || std::abs(normalizedPos.xy().minComp() - 1.f/3.f) < 0.01f)
//            std::cout << "Pos: " << normalizedPos << " -> " << posInMatrix << ": " << distortionValue << std::endl;
        return distortionValue * supportBBox.dimensions();
    }));
}

void ImplicitUnaryOperator::spread(float factor)
{
    this->_spreadingFactor += factor;
    float newSpreadFactor = this->_spreadingFactor;
    auto BBox = this->getBBox();
    Vector3 center = BBox.center().xy();
    this->transforms.push_back(UnaryOpSpread(BBox, factor));
}

void ImplicitUnaryOperator::addWavelets()
{
    this->transforms.push_back(UnaryOpWrap(
                                   [](Vector3 pos) -> Vector3 {
                                   return Vector3(0, 0, .5f * std::cos(pos.x * 1.f));
                               }));
}



std::function<float (Vector3)> ImplicitPatch::createSphereFunction(float sigma, float width, float depth, float height)
{
    std::function sphere = [width](Vector3 pos) {
        Vector3 center = Vector3(width, width, width) * .5f;
        float sqrDist = (pos - center).norm2();
        float value = std::max(0.f, (sqrDist > width * width ? 0.f : 1 - std::abs(std::sqrt(sqrDist)) / (width)));
        return value;
    };
    return sphere;
}

std::function<float (Vector3)> ImplicitPatch::createBlockFunction(float sigma, float width, float depth, float height)
{
    return [sigma, width, depth, height] (Vector3 pos) {
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

std::function<float (Vector3)> ImplicitPatch::createGaussianFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([width, depth, sigma, height](Vector3 pos) { return normalizedGaussian(Vector3(width, depth, 0), pos.xy(), sigma) * height; });
}

std::function<float (Vector3)> ImplicitPatch::createCylinderFunction(float sigma, float width, float depth, float height)
{
    Vector3 start = Vector3(width * .5f, depth * .5f, height * 0.f);
    Vector3 end = Vector3(width * .5f, depth * .5f, height * 1.f);
//    float radius = sigma;
    return [=] (Vector3 pos) {
        // Line on Z axis
        if (pos.z < start.z || end.z < pos.z)
            return 0.2f;
        BSpline spline({start, end});
        Vector3 closestPoint = spline.estimateClosestPos(pos, 1e-5);
        Vector3 normalizedClosestPoint = ((pos - closestPoint) / (Vector3(width, depth, height)));
        float distance = normalizedClosestPoint.norm();
        float evaluation = std::clamp(1.f - distance, 0.f, 1.f);
        return evaluation;

/*
        // From https://iquilezles.org/articles/distfunctions/
        Vector3 ba = end - start;
        Vector3 pa = pos - start;
        float baba = ba.dot(ba);
        float paba = pa.dot(ba);
        float x = (pa*baba-ba*paba).norm() - radius*baba;
        float y = std::abs(paba-baba*0.5)-baba*0.5;
        float x2 = x*x;
        float y2 = y*y*baba;

        float d = (std::max(x,y)<0.0)?-std::min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));

        float dist = (d > 0 ? 1.f : -1.f) * std::sqrt(std::abs(d))/baba;
        float eval = std::clamp(1 - dist, 0.f, 1.f);
        return eval;*/
    };
}

std::function<float (Vector3)> ImplicitPatch::createRockFunction(float sigma, float width, float depth, float height)
{
    auto sphereFunction = ImplicitPatch::createSphereFunction(sigma, width, depth, height);
    std::function rockFunction = [=] (Vector3 pos) {
        FastNoiseLite noise;
//        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseValue = (noise.GetNoise(pos.x * (100.f / width) + noiseOffset, pos.y * (100.f / width) + noiseOffset, pos.z * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
        return sphereFunction(pos + Vector3(0, 0, 0/*width * .25f*/)) - (noiseValue * .5f);
    };
    return rockFunction;
}

std::function<float (Vector3)> ImplicitPatch::createMountainFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([=](Vector3 pos) {
        Vector3 translatedPos = pos - (Vector3(width, depth, 0) * .5f);
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(width, depth, 1.f) * .5f);
//        return normalizedPos.norm();
        float functionsHeight = std::clamp((1.f - normalizedPos.norm()), 0.f, 1.f) * height;
        return functionsHeight;
    });
}

std::function<float (Vector3)> ImplicitPatch::createDuneFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::createMountainFunction(sigma, width, depth, height);
}

std::function<float (Vector3)> ImplicitPatch::createBasinFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::createBlockFunction(sigma, width, depth, height);
}

std::function<float (Vector3)> ImplicitPatch::createCaveFunction(float sigma, float width, float depth, float height)
{
    std::function caveFunc = [=] (Vector3 pos) {
        Vector3 start = Vector3(width * .5f, depth * 1.f, height * 1.f);
        Vector3 p1 = Vector3(width * .5f, depth * .7f, height * .3f);
        Vector3 p2 = Vector3(width * .5f, depth * .4f, height * .2f);
        Vector3 end = Vector3(width * .5f, depth * .2f, height * .2f);
        BSpline curve = BSpline({start, p1, p2, end});
        float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed
        return 1.f - (curve.estimateDistanceFrom(pos, epsilon) / (width * .5f));
    };
    return caveFunc;
}

std::function<float (Vector3)> ImplicitPatch::createArchFunction(float sigma, float width, float depth, float height)
{
    std::function archFunc = [=] (Vector3 pos) {
        Vector3 start = Vector3(width * .5f, depth * 0.f, height * 0.f);
        Vector3 p1 = Vector3(width * .5f, depth * .3f, height * 1.f);
        Vector3 p2 = Vector3(width * .5f, depth * .7f, height * 1.f);
        Vector3 end = Vector3(width * .5f, depth * 1.f, height * 0.f);
        BSpline curve = BSpline({start, p1, p2, end});
        float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed
        return 1.f - (curve.estimateDistanceFrom(pos, epsilon) / (width * .5f));
    };
    return archFunc;
}

std::function<float (Vector3)> ImplicitPatch::createNoise2DFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([sigma, width, depth, height](Vector3 pos) -> float {
        FastNoiseLite noise;
    //    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseValue = (noise.GetNoise(pos.x * (100.f * sigma / width) + noiseOffset, pos.y * (100.f * sigma / depth) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
        return noiseValue * height;
    });
    //float noiseValue = (noise.GetNoise(pos.x * (100.f / width) + noiseOffset, pos.y * (100.f / width) + noiseOffset, pos.z * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
}

std::function<float (Vector3)> ImplicitPatch::createMountainChainFunction(float sigma, float width, float depth, float height, BSpline _path)
{
//    Vector3 c = Vector3(0.6, 0.8); // std::sin(deg2rad(45.f)), std::cos(deg2rad(45.f)));
    return ImplicitPatch::convert2DfunctionTo3Dfunction([=] (Vector3 pos) -> float {
        BSpline path = _path;
        float closestTime = path.estimateClosestTime(pos);
        Vector3 closestPoint = path.getPoint(closestTime);
        Vector3 vertical(0, 0, 1);
        Vector3 islandDirection = vertical.cross(path.getDirection(closestTime)).normalize();

        Vector3 translatedPos = pos - closestPoint;
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(sigma, sigma, 1.f) * .5f);
        float distance = normalizedPos.norm();
        float functionsHeight = std::clamp((1.f - distance), 0.f, 1.f) * sigma;
        return functionsHeight;
    });
        /*Vector3 p = (pos - closestPoint) / (sigma);
        Vector3 q = Vector3(p.xy().norm(), -p.y);
        float d = (q - c * std::max(q.dot(c), 0.f)).norm();
        float signedDist = d * ((q.x * c.y - q.y * c.x < 0.f) ? -1.f : 1.f);
        return std::clamp(1.f - signedDist, 0.f, 1.f);
        */
        /*
            vec2 q = vec2( length(p.xz), -p.y );
            float d = length(q-c*max(dot(q,c), 0.0));
            return d * ((q.x*c.y-q.y*c.x<0.0)?-1.0:1.0);
        */
        /*
        BSpline path = _path;
        float closestTime = path.estimateClosestTime(pos);
        Vector3 closestPoint = path.getPoint(closestTime);
        Vector3 vertical(0, 0, 1);
        Vector3 islandDirection = vertical.cross(path.getDirection(closestTime)).normalize();
        Vector3 toCurve = (closestPoint - pos);
//        float distance = toCurve.norm();
        float flatDistance = toCurve.xy().norm();
//        toCurve.normalize();
//        float islandFactor = std::pow((toCurve.xy().normalized().dot(islandDirection) + 1.f) * .5f, 1.f); // [-1, 1] -> [0, 1]
        float distanceEval = (flatDistance * toCurve.z); //  * (islandFactor + 1.f); // d * [1, 2]
        float normalizedEval = (distanceEval / sigma);
        float eval = std::clamp(1.f - normalizedEval, 0.f, 1.f);
        return eval;
        */
/*
        Vector3 translatedPos = pos - closestPoint;
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(width, depth, 1.f) * .5f);
//        return normalizedPos.norm();
        float functionsHeight = std::clamp((1.f - normalizedPos.norm()), 0.f, 1.f) * height;
        return functionsHeight;*/
    //    };
}

std::function<float (Vector3)> ImplicitPatch::createPolygonFunction(float sigma, float width, float depth, float height, BSpline path)
{
    ShapeCurve polygon = path;
    for (auto& p : polygon.points)
        p.z = 0.f;
    polygon.points.push_back(polygon.points.front());
    return ImplicitPatch::convert2DfunctionTo3Dfunction([=, _polygon=polygon] (Vector3 pos) -> float {
        ShapeCurve polygon(_polygon);
        return (polygon.contains(pos.xy(), false) ? height : 0.f);
    });
}

std::function<float (Vector3)> ImplicitPatch::createParametricTunnelFunction(float sigma, float width, float depth, float height, BSpline _path)
{
    return [=] (Vector3 pos) -> float {
        BSpline path = _path;
//        float closestTime = path.estimateClosestTime(pos);
//        Vector3 closestPoint = path.getPoint(closestTime);
        float epsilon = 1e-1;
        return 1.f - (path.estimateDistanceFrom(pos, epsilon) / (sigma));
        /*Vector3 newBasisX = path.getFrenetDirection(closestTime); // * sigma;
        Vector3 newBasisY = path.getFrenetBinormal(closestTime); // * sigma;
        Vector3 newBasisZ = path.getFrenetNormal(closestTime); //  * sigma;
        ShapeCurve shape = ShapeCurve({
                                          Vector3(-.5, 0, 0).changeBasis(newBasisX, newBasisY, newBasisZ),
                                          Vector3(-.2, 0, 1).changeBasis(newBasisX, newBasisY, newBasisZ),
                                          Vector3(0.5, 0, 0).changeBasis(newBasisX, newBasisY, newBasisZ),
                                          Vector3(0.2, 0, 1).changeBasis(newBasisX, newBasisY, newBasisZ)
                                      });
        shape.scale(sigma);
        shape.translate(closestPoint);

        float res = ImplicitPatch::isovalue - shape.estimateSignedDistanceFrom(pos);*/
//        float res = /*ImplicitPatch::isovalue - */(closestPoint - pos).norm() / sigma;
//        return std::max(res, 0.f);
        /*
        float distance = (pos - closestPoint).norm() / sigma;
        float functionsHeight = std::clamp((1.f - distance), 0.f, 1.f);

        return functionsHeight;
        */
    };
}

std::function<float (Vector3)> ImplicitPatch::createRippleFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([=] (Vector3 pos) -> float {
        if(!Vector3::isInBox(pos, Vector3(), Vector3(width, depth, height)))
            return 0.f;
        float x = 2 * pos.x / width;
        float y = 2.f * pos.y / depth - 1;
        float maxHeight = (x < 1.f ? x : 1.f/(x*x));
        float damping = std::max(1 - y*y, 0.f);
//        std::cout << x << " " << y << " " << maxHeight << " " << damping << std::endl;
        return height * maxHeight * damping;
    });
}

std::function<float (Vector3)> ImplicitPatch::createIdentityFunction(float sigma, float width, float depth, float height)
{
    return [](Vector3) -> float { return 0.f; };
}

std::function<float (Vector3)> ImplicitPatch::convert2DfunctionTo3Dfunction(std::function<float (Vector3)> func)
{
    std::function _3Dfunction = [func](Vector3 pos) {
        float height = func(pos);
        float inverse_isovalue = .5f - (pos.z / height); //height * .5f - pos.z;
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

/*
ImplicitCSG::ImplicitCSG()
{

}

float ImplicitCSG::evaluate(Vector3 pos)
{
    float valA = this->evaluateA(pos);
    float valB = this->evaluateB(pos);

    float evaluation = 0.f;
    switch (this->operation) {
    case Union:
        evaluation = std::max(valA, valB);
        break;
    case Difference:
        evaluation = valB - valA;
        break;
    case Intersection:
        evaluation = std::min(valA, valB);
        break;
    }
    return evaluation;
}

std::map<TerrainTypes, float> ImplicitCSG::getMaterials(Vector3 pos)
{
    auto matsA = this->getMaterialsA(pos);
    auto matsB = this->getMaterialsB(pos);

    float evaluation = 0.f;
    switch (this->operation) {
    case Union:
        evaluation = std::max(valA, valB);
        break;
    case Difference:
        evaluation = valB - valA;
        break;
    case Intersection:
        evaluation = std::min(valA, valB);
        break;
    }
    return evaluation;
}
*/

UnaryOpTranslate::UnaryOpTranslate(Vector3 translation)
    : UnaryOp()
{
    this->wrap = [=](Vector3 pos) {
        return pos + translation;
    };
    this->unwrap = [=](Vector3 pos) {
        return pos - translation;
    };
}

UnaryOpRotate::UnaryOpRotate(Vector3 rotationAngles, Vector3 center)
    : UnaryOp()
{
    float angleX = rotationAngles.x;
    float angleY = rotationAngles.y;
    float angleZ = rotationAngles.z;
    Matrix Rx (3, 3, std::vector<float>({
        1, 0, 0,
        0, cos(angleX), -sin(angleX),
        0, sin(angleX), cos(angleX)
    }).data());
    Matrix Rz (3, 3, std::vector<float>({
             cos(angleY), 0, -sin(angleY),
             0, 1, 0,
             sin(angleY), 0, cos(angleY)
         }).data());
    Matrix Ry (3, 3, std::vector<float>({
             cos(angleZ), -sin(angleZ), 0,
             sin(angleZ), cos(angleZ), 0,
             0, 0, 1
         }).data());
    Matrix R = Rx.product(Ry).product(Rz);

    angleX *= -1.f;
    angleY *= -1.f;
    angleZ *= -1.f;
    Matrix _Rx (3, 3, std::vector<float>({
            1, 0, 0,
            0, cos(angleX), -sin(angleX),
            0, sin(angleX), cos(angleX)
        }).data());
    Matrix _Rz (3, 3, std::vector<float>({
             cos(angleY), 0, -sin(angleY),
             0, 1, 0,
             sin(angleY), 0, cos(angleY)
         }).data());
    Matrix _Ry (3, 3, std::vector<float>({
             cos(angleZ), -sin(angleZ), 0,
             sin(angleZ), cos(angleZ), 0,
             0, 0, 1
         }).data());
    Matrix _R = _Rx.product(_Ry).product(_Rz);

    //    std::function<Vector3(Vector3)> previousWrapFunction = this->wrapFunction;
    //    std::function<Vector3(Vector3)> previousUnwrapFunction = this->unwrapFunction;
    this->wrap = [=](Vector3 pos) -> Vector3 {
        // Exactly the rotation function, but with the matrix computed only once
        return center + (pos - center).applyTransform(R);
    };
    this->unwrap = [=](Vector3 pos) -> Vector3 {
        // Exactly the rotation function, but with the matrix computed only once
        return center + (pos - center).applyTransform(_R);
    };
}

UnaryOpScale::UnaryOpScale(Vector3 scaling, Vector3 center)
    : UnaryOp()
{
    this->wrap = [=] (Vector3 pos) {
        return center - (pos - center) * scaling;
    };
    this->unwrap = [=] (Vector3 pos) {
        return center + (pos - center) / scaling;
    };
}

UnaryOpWrap::UnaryOpWrap(FastNoiseLite noise, Vector3 strength)
    : UnaryOp()
{
    this->wrap = [=] (Vector3 pos) {
        FastNoiseLite _noise = noise;
        return pos + Vector3(1.f, 1.f, 1.f) * _noise.GetNoise((float)pos.x * strength.y + strength.z, (float)pos.y * strength.y + strength.z, (float)pos.z * strength.y + strength.z) * strength.x;
    };
    this->unwrap = [=] (Vector3 pos) {
        FastNoiseLite _noise = noise;
        return pos - Vector3(1.f, 1.f, 1.f) * _noise.GetNoise((float)pos.x * strength.y + strength.z, (float)pos.y * strength.y + strength.z, (float)pos.z * strength.y + strength.z) * strength.x;
    };
}

UnaryOpWrap::UnaryOpWrap(std::function<Vector3 (Vector3)> func)
{
    this->wrap = [=] (Vector3 pos) {
        Vector3 newPos = pos + func(pos);
        return newPos;
    };
    this->unwrap = [=] (Vector3 pos) {
        return pos - func(pos);
    };
}

UnaryOpWrap::UnaryOpWrap(Matrix3<Vector3> wrapper)
    : wrapper(wrapper)
{
    this->wrap = [=] (Vector3 pos) {
        Vector3 newPos = pos + this->wrapper.at(pos);
        return newPos;
    };
    this->unwrap = [=] (Vector3 pos) {
        return pos - this->wrapper.at(pos);
    };
}

UnaryOpSpread::UnaryOpSpread(AABBox BBox, float spreadFactor)
    : UnaryOp()
{
    Vector3 center = BBox.center();
    Vector3 dimensions = BBox.dimensions();
    this->wrap = [=](Vector3 pos) -> Vector3 {
        Vector3 p = pos;
        float relativeHeight = interpolation::linear(p.z, BBox.min().z, BBox.max().z);
        relativeHeight = interpolation::smooth(relativeHeight);
//        float heightFactor = (1.f - relativeHeight);
        float heightFactor = -(relativeHeight - .5f) * 2.f; // Top = -1, bottom = 1
        float factor = spreadFactor * heightFactor;
        Vector3 toCenter = (p.xy() - center.xy());
        Vector3 displacement = toCenter.normalize() * dimensions * factor;
        return pos + displacement;
    };
    this->unwrap = [=](Vector3 pos) -> Vector3 {
        Vector3 p = pos;
        float relativeHeight = interpolation::linear(p.z, BBox.min().z, BBox.max().z);
        relativeHeight = interpolation::smooth(relativeHeight);
//        float heightFactor = (1.f - relativeHeight);
        float heightFactor = -(relativeHeight - .5f) * 2.f; // Top = -1, bottom = 1
        float factor = spreadFactor * heightFactor;
        Vector3 toCenter = (p.xy() - center.xy());
        Vector3 displacement = toCenter.normalize() * dimensions * .5f * factor;
        return pos - displacement;
    };
}


ImplicitPatch::CompositionFunction compositionOperationFromString(std::string name)
{
    name = toUpper(name);
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

    // Error
    std::cerr << "Wrong operation string given" << std::endl;
    return ImplicitPatch::CompositionFunction::NONE;
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

    // Impossible
    std::cerr << "Unknown operation given." << std::endl;
    return "";
}

ImplicitPatch::PredefinedShapes predefinedShapeFromString(std::string name)
{
    name = toUpper(name);
    if (name == "SPHERE")
        return ImplicitPatch::PredefinedShapes::Sphere;
    else if (name == "BLOCK")
        return ImplicitPatch::PredefinedShapes::Block;
    else if (name == "GAUSSIAN")
        return ImplicitPatch::PredefinedShapes::Gaussian;
    else if (name == "CYLINDER")
        return ImplicitPatch::PredefinedShapes::Cylinder;
    else if (name == "ROCK")
        return ImplicitPatch::PredefinedShapes::Rock;
    else if (name == "MOUNTAIN")
        return ImplicitPatch::PredefinedShapes::Mountain;
    else if (name == "DUNE")
        return ImplicitPatch::PredefinedShapes::Dune;
    else if (name == "BASIN")
        return ImplicitPatch::PredefinedShapes::Basin;
    else if (name == "CAVE")
        return ImplicitPatch::PredefinedShapes::Cave;
    else if (name == "ARCH")
        return ImplicitPatch::PredefinedShapes::Arch;
    else if (name == "NOISE2D")
        return ImplicitPatch::PredefinedShapes::Noise2D;
    else if (name == "MOUNTAINCHAIN")
        return ImplicitPatch::PredefinedShapes::MountainChain;
    else if (name == "POLYGON")
        return ImplicitPatch::PredefinedShapes::Polygon;
    else if (name == "IMPLICIT_HEIGHTMAP")
        return ImplicitPatch::PredefinedShapes::ImplicitHeightmap;
    else if (name == "PARAMETRIC_TUNNEL")
        return ImplicitPatch::PredefinedShapes::ParametricTunnel;
    else if (name == "NONE")
        return ImplicitPatch::PredefinedShapes::None;
    else if (name == "RIPPLE")
        return ImplicitPatch::PredefinedShapes::Ripple;

    // Error
    std::cerr << "Wrong predefined shape string given. (" << name << ")" << std::endl;
    return ImplicitPatch::PredefinedShapes::Sphere;
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

    // Impossible
    std::cerr << "Unknown predefined shape given." << std::endl;
    return "";
}

ImplicitPatch::PositionalLabel positionalLabelFromString(std::string name)
{
    name = toUpper(name);

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
    std::cerr << "Wrong positional label string given. (" << name << ")" << std::endl;
    return ImplicitPatch::PositionalLabel::ABOVE;
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

TerrainTypes materialFromString(std::string name)
{
    name = toUpper(name);

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
    std::cerr << "Unknown material string given. (" << name << ")" << std::endl;
    return LAST;
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

UnaryOp::UnaryOp()
{
    this->wrap = [=] (Vector3 pos) { return pos; };
    this->unwrap = [=] (Vector3 pos) { return pos; };
}

ImplicitNaryOperator::ImplicitNaryOperator()
{

}

float ImplicitNaryOperator::evaluate(Vector3 pos)
{
    float maxVal = 0.f;
    for (auto& compo : this->composables)
        maxVal = maxVal + compo->evaluate(pos); //std::max(maxVal, compo->evaluate(pos));
    return maxVal;
}

std::map<TerrainTypes, float> ImplicitNaryOperator::getMaterials(Vector3 pos)
{
    float maxVal = 0.f;
    std::map<TerrainTypes, float> bestReturn;
    for (auto& compo : this->composables) {
        auto [total, evaluation] = compo->getMaterialsAndTotalEvaluation(pos);
        for (auto& [mat, val] : evaluation) {
            if (bestReturn.count(mat) == 0)
                bestReturn[mat] = 0;
            bestReturn[mat] += val;
        }
    }
    return bestReturn;
    /*
    float maxVal = 0.f;
    std::map<TerrainTypes, float> bestReturn;
    for (auto& compo : this->composables) {
        auto evaluation = compo->getMaterialsAndTotalEvaluation(pos);
        float totalEval = evaluation.first;
        if (totalEval > maxVal) {
            maxVal = std::max(maxVal, totalEval);
            bestReturn = evaluation.second;
        }
    }
    return bestReturn;*/
}

AABBox ImplicitNaryOperator::getSupportBBox()
{
    Vector3 minPos(false), maxPos(false);
    for (auto& compo : this->composables) {
        auto BBox = compo->getSupportBBox();
        minPos = Vector3::min(minPos, BBox.min());
        maxPos = Vector3::max(maxPos, BBox.max());
    }
    return {minPos, maxPos};
}

AABBox ImplicitNaryOperator::getBBox()
{
    Vector3 minPos(false), maxPos(false);
    for (auto& compo : this->composables) {
        auto BBox = compo->getBBox();
        minPos = Vector3::min(minPos, BBox.min());
        maxPos = Vector3::max(maxPos, BBox.max());
    }
    return {minPos, maxPos};
}

void ImplicitNaryOperator::update()
{
    // Nothing to do...
}

std::string ImplicitNaryOperator::toString()
{
    std::string compoNames;
    for (auto& compo : this->composables)
        compoNames += compo->name + ", ";
    return this->name + ": N-operation between " + compoNames;
}

nlohmann::json ImplicitNaryOperator::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        content["file"] = this->used_json_filename;
    } else {
        content["type"] = "nary";
        content["name"] = this->name;
        content["index"] = this->index;
        content["optionalCurve"] = bspline_to_json(this->optionalCurve);
        content["mirrored"] = this->mirrored;
        std::vector<nlohmann::json> composableJson;
        for (auto& compo : this->composables)
            composableJson.push_back(compo->toJson());
        content["composables"] = composableJson;
    }

    return content;
}

ImplicitPatch *ImplicitNaryOperator::fromJson(nlohmann::json content)
{
    ImplicitNaryOperator* patch = new ImplicitNaryOperator;
    for (auto& compoJson : content["composables"]) {
        patch->addChild(ImplicitPatch::fromJson(compoJson));
    }
    patch->name = content["name"];
    patch->index = content["index"];
    if (content.contains("optionalCurve"))
        patch->optionalCurve = json_to_bspline(content["optionalCurve"]);
    if (content.contains("mirrored"))
        patch->mirrored = content["mirrored"];
    patch->update();
    return patch;
}

bool ImplicitNaryOperator::contains(PredefinedShapes shape)
{
    for (auto& compo : composables)
        if (compo->contains(shape))
            return true;
    return false;
}

std::vector<ImplicitPatch *> ImplicitNaryOperator::findAll(PredefinedShapes shape)
{
    std::vector<ImplicitPatch*> res;
    for (auto& compo : this->composables)
        res = vectorMerge(res, compo->findAll(shape));
    return res;
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
        compo->updateCache();
}

ImplicitPatch *ImplicitNaryOperator::copy() const
{
    return new ImplicitNaryOperator(*this);
}

void ImplicitNaryOperator::addChild(ImplicitPatch *newChild, int index)
{
    if (index == -1) {
        this->composables.push_back(newChild);
    } else {
        if (composables[index])
            this->composables[index]->parent = nullptr;
        this->composables[index] = newChild;
    }
    newChild->setParent(this);
//    std::cout << newChild->toString() << " is now child of " << this->toString() << std::endl;
}

void ImplicitNaryOperator::deleteAllChildren()
{
    for (auto& child : composables) {
        delete child;
    }
    this->composables.clear();
}

nlohmann::json ImplicitTranslation::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["translation"] = vec3_to_json(this->_translation);
    return content;
}

ImplicitPatch *ImplicitTranslation::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitRotation::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["rotation"] = vec3_to_json(this->_rotation);
    return content;
}

ImplicitPatch *ImplicitRotation::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitScaling::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["scale"] = vec3_to_json(this->_scale);
    return content;
}

ImplicitPatch *ImplicitScaling::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitWraping::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["distortion"] = vec3_to_json(this->_distortion);
    return content;
}

ImplicitPatch *ImplicitWraping::fromJson(nlohmann::json content)
{
    return ImplicitUnaryOperator::fromJson(content);
}

nlohmann::json ImplicitNoise::toJson()
{
    nlohmann::json content = ImplicitUnaryOperator::toJson();
    content["noise"] = vec3_to_json(this->_noise);
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
