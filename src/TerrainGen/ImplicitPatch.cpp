#include "ImplicitPatch.h"



#ifndef useBefore

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
    if (this->_cachedMaxHeight.at((pos - AABBox.first).xy()) == -1.f) {
        float minHeight = AABBox.first.z;
        float maxHeight = AABBox.second.z;

        this->_cachedMaxHeight.at((pos - AABBox.first).xy()) = 0;
        for (float z = maxHeight; z > minHeight; z -= ImplicitPatch::zResolution) {
            float eval = this->evaluate(pos.xy() + Vector3(0, 0, z));
            if (eval >= ImplicitPatch::isovalue) {
                this->_cachedMaxHeight.at((pos - AABBox.first).xy()) = z;
                break;
            }
        }
    }
    return this->_cachedMaxHeight.at((pos - AABBox.first).xy());
}

float ImplicitPatch::getMinHeight(Vector3 pos)
{
    auto AABBox = this->getBBox();
    if (this->_cachedMinHeight.at((pos - AABBox.first).xy()) == -1.f) {
        float minHeight = AABBox.first.z;
        float maxHeight = AABBox.second.z;

        this->_cachedMinHeight.at((pos - AABBox.first).xy()) = 0;
        for (float z = minHeight; z < maxHeight; z += ImplicitPatch::zResolution) {
            float eval = this->evaluate(pos.xy() + Vector3(0, 0, z));
            if (eval >= ImplicitPatch::isovalue) {
                this->_cachedMinHeight.at((pos - AABBox.first).xy()) = z;
                break;
            }
        }
    }
    return this->_cachedMinHeight.at((pos - AABBox.first).xy());
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
    Vector3 a = AABBox.first, b = AABBox.second;
    return AABBox.second - AABBox.first;
}

Vector3 ImplicitPatch::getSupportDimensions()
{
    auto AABBox = this->getSupportBBox();
    return AABBox.second - AABBox.first;
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
        else if (type == "COMPOSE")
            result = ImplicitOperator::fromJson(content);
        else if (type == "UNARY")
            result = ImplicitUnaryOperator::fromJson(content);
    } else if (content.contains("file")) {
        nlohmann::json newContent = nlohmann::json::parse(std::ifstream(content["file"]));
        result = ImplicitPatch::fromJson(newContent[ImplicitPatch::json_identifier]);
        result->used_json_filename = content["file"];
        return result;
    }
    if (result != nullptr) {
        result->updateCache();
    }
    return result;
}

void ImplicitPatch::updateCache()
{
    Vector3 dim = Vector3::max(Vector3(1, 1), this->getDimensions().roundedUp());

    this->_cachedMaxHeight = Matrix3<float>(dim.x, dim.y, 1, -1.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;

    this->_cachedMinHeight = Matrix3<float>(dim.x, dim.y, 1, -1.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

ImplicitPatch *ImplicitPatch::createPredefinedShape(PredefinedShapes shape, Vector3 dimensions, float additionalParam)
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
    case None:
        func = ImplicitPatch::createIdentityFunction(additionalParam, dimensions.x, dimensions.y, dimensions.z);
        break;
    }
    ImplicitPrimitive* primitive = new ImplicitPrimitive();
    primitive->evalFunction = func;

    return primitive;
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

    float distanceToBBox = Vector3::signedDistanceToBoundaries(pos, minPos, maxPos);
    float distanceFactor = std::max(1.f, (supportWidth - width).maxComp() * .5f); // Don't get a 0 here
    distanceToBBox /= distanceFactor; // Make it depending on the support area
    float distanceFalloff = interpolation::wyvill(std::clamp(distanceToBBox, 0.f, 1.f));
    float evaluation = this->evalFunction(pos - this->position) * distanceFalloff;
    evaluation = std::clamp(evaluation, 0.f, 1.f);

    return evaluation;
}

std::map<TerrainTypes, float> ImplicitPrimitive::getMaterials(Vector3 pos)
{
    return {{this->material, this->evaluate(pos)}};
}

std::pair<Vector3, Vector3> ImplicitPrimitive::getSupportBBox()
{
    Vector3 margin = (this->supportDimensions - this->dimensions) * .5f;
    if (!this->supportDimensions.isValid() || this->supportDimensions == Vector3())
        margin = this->dimensions * .5f; // Nice big margin
    return {this->position - margin, this->position + this->dimensions + margin};
}

std::pair<Vector3, Vector3> ImplicitPrimitive::getBBox()
{
    return {this->position, this->position + this->dimensions};
}

void ImplicitPrimitive::update()
{
    // This is kinda f*cked up...
    ImplicitPrimitive* copy = (ImplicitPrimitive*)ImplicitPatch::createPredefinedShape(this->predefinedShape, this->dimensions, this->parametersProvided[0]);
    copy->setDimensions(this->dimensions);
    copy->setSupportDimensions(this->supportDimensions);
    copy->position = this->position;
    copy->material = this->material;
    copy->parametersProvided = this->parametersProvided;
    copy->index = this->index;
    copy->name = this->name;
    copy->predefinedShape = this->predefinedShape;
    *this = *copy; // Not sure this is legal
    delete copy;
}

std::string ImplicitPrimitive::toString()
{
    return this->name + ": Primitive " + this->name + " #" + std::to_string(this->index);
}

nlohmann::json ImplicitPrimitive::toJson()
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
    patch->index = content["index"];
    patch->predefinedShape = predefinedShapeFromString(content["shape"]);
    patch->parametersProvided = {content["sigmaValue"]};
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






ImplicitOperator::ImplicitOperator()
{

}

float ImplicitOperator::evaluate(Vector3 pos)
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

std::map<TerrainTypes, float> ImplicitOperator::getMaterials(Vector3 pos)
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
        return {{TerrainTypes(), 0.f}};
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

    if (this->composeFunction == REPLACE) {
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
        if (result.count(mat) == 0)
            result[mat] = 0.f;
        result[mat] += val * ratioA;
    }
    for (const auto& [mat, val] : materialsAndEvalsForB) {
        if (result.count(mat) == 0)
            result[mat] = 0.f;
        result[mat] += val * ratioB;
    }

    return result;
}

float ImplicitOperator::evaluateFromAandB(float evalA, float evalB)
{
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
}

float ImplicitOperator::evaluateA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);

    return this->composableA->evaluate(evaluationPosA);
}

float ImplicitOperator::evaluateB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);

    return this->composableB->evaluate(evaluationPosB);
}

std::map<TerrainTypes, float> ImplicitOperator::getMaterialsA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return this->composableA->getMaterials(evaluationPosA);
}

std::map<TerrainTypes, float> ImplicitOperator::getMaterialsB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return this->composableB->getMaterials(evaluationPosB);
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitOperator::getMaterialsAndTotalEvaluationA(Vector3 pos)
{
    Vector3 evaluationPosA = this->getEvaluationPositionForComposableA(pos);
    return this->composableA->getMaterialsAndTotalEvaluation(evaluationPosA);
}

std::pair<float, std::map<TerrainTypes, float>> ImplicitOperator::getMaterialsAndTotalEvaluationB(Vector3 pos)
{
    Vector3 evaluationPosB = this->getEvaluationPositionForComposableB(pos);
    return this->composableB->getMaterialsAndTotalEvaluation(evaluationPosB);
}

std::pair<Vector3, Vector3> ImplicitOperator::getSupportBBox()
{
    auto AABBoxA = this->composableA->getSupportBBox();
    if (this->withIntersectionOnB) { // No need to go further, we know the limit with the intersection
        return AABBoxA;
    }
    auto AABBoxB = this->composableB->getSupportBBox();

    if (this->positionalB == PositionalLabel::ABOVE) {
        // If stacked on composable A, composable B can get higher
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_BOTTOM) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_TOP) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::FIXED_POS) {
        // Nothing to do there
    }

    return {Vector3::min(AABBoxA.first, AABBoxB.first), Vector3::max(AABBoxA.second, AABBoxB.second)};
}

std::pair<Vector3, Vector3> ImplicitOperator::getBBox()
{
    auto AABBoxA = this->composableA->getBBox();
    if (this->withIntersectionOnB) { // No need to go further, we know the limit with the intersection
        return AABBoxA;
    }
    auto AABBoxB = this->composableB->getBBox();

    if (this->positionalB == PositionalLabel::ABOVE) {
        // If stacked on composable A, composable B can get higher
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_BOTTOM) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::INSIDE_TOP) {
        // Same as before "Above"?
        float maxPossibleHeightA = AABBoxA.second.z - AABBoxA.first.z;
        AABBoxB.second.z += maxPossibleHeightA;
    } else if (this->positionalB == PositionalLabel::FIXED_POS) {
        // Nothing to do there
    }

    return {Vector3::min(AABBoxA.first, AABBoxB.first), Vector3::max(AABBoxA.second, AABBoxB.second)};
}

void ImplicitOperator::update()
{
    // Guess we have nothing to do...
}

std::string ImplicitOperator::toString()
{
    return this->name + ": Operation between " + (composableA ? "#" + std::to_string(composableA->index) : "undefined") + " and " + (composableB ? "#" + std::to_string(composableB->index) : "undefined");
}

nlohmann::json ImplicitOperator::toJson()
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
        if (this->composableA)
            content["composableA"] = this->composableA->toJson();
        if (this->composableB)
            content["composableB"] = this->composableB->toJson();
    }

    return content;
}

ImplicitPatch *ImplicitOperator::fromJson(nlohmann::json content)
{
    ImplicitOperator* patch = new ImplicitOperator;
    patch->composableA = ImplicitPatch::fromJson(content["composableA"]);
    patch->composableB = ImplicitPatch::fromJson(content["composableB"]);
    patch->name = content["name"];
    patch->composeFunction = compositionOperationFromString(content["operator"]);
    patch->blendingFactor = content["blendingFactor"];
    patch->positionalB = positionalLabelFromString(content["positionalB"]);
    patch->index = content["index"];
    if (content.contains("useIntersection"))
        patch->withIntersectionOnB = content["useIntersection"];
    patch->update();
    return patch;
}

void ImplicitOperator::updateCache()
{
    ImplicitPatch::updateCache();
    if (this->composableA != nullptr)
        this->composableA->updateCache();
    if (this->composableB != nullptr)
        this->composableB->updateCache();
}

Vector3 ImplicitOperator::getEvaluationPositionForComposableA(Vector3 pos)
{
    return pos; // Nothing to do
}

Vector3 ImplicitOperator::getEvaluationPositionForComposableB(Vector3 pos)
{
    // Get the correct evaluation position for the B composent
    float offsetB = 0.f;
    if (this->positionalB == ABOVE) {
        offsetB = this->composableA->getMaxHeight(pos);
    } else if (this->positionalB == INSIDE_TOP) {
        offsetB = this->composableA->getMaxHeight(pos) - this->composableB->getDimensions().z;
    } else if (this->positionalB == INSIDE_BOTTOM) {
        offsetB = this->composableA->getMinHeight(pos);
    } else if (this->positionalB == FIXED_POS) {
        offsetB = 0.f;
    } else {
        std::cerr << "Wrong position label" << std::endl;
    }
    pos.z -= offsetB;
    return pos;
}





ImplicitUnaryOperator::ImplicitUnaryOperator()
{
    this->wrapFunction = [](Vector3) { return Vector3(0, 0, 0); };
    this->noiseFunction = [](Vector3) { return 0.f; };
}

float ImplicitUnaryOperator::evaluate(Vector3 pos)
{
    Vector3 evaluationPos = pos + this->wrapFunction(pos);
    float evaluation = this->composableA->evaluate(evaluationPos) + this->noiseFunction(pos);
    evaluation = std::clamp(evaluation, 0.f, 1.f);

    return evaluation;
}

std::map<TerrainTypes, float> ImplicitUnaryOperator::getMaterials(Vector3 pos)
{
    Vector3 evaluationPos = pos + this->wrapFunction(pos);
    return this->composableA->getMaterials(evaluationPos);
}

std::pair<Vector3, Vector3> ImplicitUnaryOperator::getSupportBBox()
{
    auto AABBox = this->composableA->getSupportBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.first, AABBox.second);

    for (auto& vert : vertices) // call transform function
        vert -= this->wrapFunction(vert);

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}

std::pair<Vector3, Vector3> ImplicitUnaryOperator::getBBox()
{
    auto AABBox = this->composableA->getBBox();
    auto vertices = Vector3::getAABBoxVertices(AABBox.first, AABBox.second);

    for (auto& vert : vertices) // call transform function
        vert -= this->wrapFunction(vert);

    return {Vector3::min(vertices), Vector3::max(vertices)}; // Get minimal and maximal
}

std::string ImplicitUnaryOperator::toString()
{
    return this->name + ": Unary operation on " + (composableA ? "#" + std::to_string(composableA->index) : "undefined");
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
        content["translation"] = vec3_to_json(this->_translation);
        content["rotation"] = vec3_to_json(this->_rotation);
        if (this->composableA)
            content["composableA"] = this->composableA->toJson();
    }

    return content;
}

ImplicitPatch *ImplicitUnaryOperator::fromJson(nlohmann::json content)
{
    ImplicitUnaryOperator* patch = new ImplicitUnaryOperator;
    patch->composableA = ImplicitPatch::fromJson(content["composableA"]);
    patch->name = content["name"];
    patch->index = content["index"];
    Vector3 rotation = json_to_vec3(content["rotation"]);
    Vector3 translation = json_to_vec3(content["translation"]);
    Vector3 scale;
    if (content.contains("scale"))
        scale = json_to_vec3(content["scale"]);
    if (rotation != Vector3()) patch->rotate(rotation.x, rotation.y, rotation.z);
    if (translation != Vector3()) patch->translate(translation);
    if (scale != Vector3()) patch->scale(scale);
    patch->update();
    return patch;
}

void ImplicitUnaryOperator::translate(Vector3 translation)
{
    this->_translation += translation;
    std::function<Vector3(Vector3)> previousFunction = this->wrapFunction;
    this->wrapFunction = [=](Vector3 pos) -> Vector3 {
        return previousFunction(pos) - translation;
    };
}

void ImplicitUnaryOperator::rotate(float angleX, float angleY, float angleZ)
{
    this->_rotation += Vector3(angleX, angleY, angleZ);

    auto AABBox = this->getBBox();
    Vector3 center = (AABBox.first + AABBox.second) * .5f;
    Matrix Rx (3, 3, std::vector<float>({
                   1, 0, 0,
                   0, cos(-angleX), -sin(-angleX),
                   0, sin(-angleX), cos(-angleX)
               }).data());
    Matrix Rz (3, 3, std::vector<float>({
                    cos(-angleY), 0, -sin(-angleY),
                    0, 1, 0,
                    sin(-angleY), 0, cos(-angleY)
                }).data());
    Matrix Ry (3, 3, std::vector<float>({
                    cos(-angleZ), -sin(-angleZ), 0,
                    sin(-angleZ), cos(-angleZ), 0,
                    0, 0, 1
                }).data());
    Matrix R = Rx.product(Ry).product(Rz);

    std::function<Vector3(Vector3)> previousFunction = this->wrapFunction;
    this->wrapFunction = [=](Vector3 pos) -> Vector3 {
        // Exactly the rotation function, but with the matrix computed only once
        return (((pos + previousFunction(pos)) - center).applyTransform(R) + center) - pos;
    };
    /*
//    std::cout << "Rotation : " << rad2deg(angleX) << "x, " << rad2deg(angleY) << "y, " << rad2deg(angleZ) << "z" << std::endl;
    auto AABBox = this->getBBox();
    Vector3 center = (AABBox.first + AABBox.second) * .5f;
//    Vector3 center = (AABBox.second - AABBox.first) * .5f;
    std::function<Vector3(Vector3)> previousFunction = this->wrapFunction;
    this->wrapFunction = [=](Vector3 pos) -> Vector3 {
//        Vector3 firstTransform = pos + previousFunction(pos);
//        Vector3 fromCenter = (firstTransform - center);
//        Vector3 rotated = fromCenter.rotate(-angleX, -angleY, -angleZ);
//        Vector3 backToCenter = rotated + center;
//        Vector3 backToTranslation = backToCenter - pos;
//        return  backToTranslation;

        return (((pos + previousFunction(pos)) - center).rotate(-angleX, -angleY, -angleZ) + center) - pos;
    };*/
}

void ImplicitUnaryOperator::scale(Vector3 scaleFactor)
{
    this->_scale *= scaleFactor;
    auto BBox = this->getBBox();
    Vector3 center = (BBox.first + BBox.second) * .5f;
    std::function<Vector3(Vector3)> previousFunction = this->wrapFunction;
    this->wrapFunction = [=](Vector3 pos) -> Vector3 {
        return ((previousFunction(pos) - center) * scaleFactor) + center;
    };
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
        Vector3 minPos = Vector3();
        Vector3 maxPos = Vector3(width, depth, height);
        float distanceToBBox = Vector3::signedDistanceToBoundaries(pos, minPos, maxPos);
        float distanceFactor = std::max(1.f, maxPos.maxComp() * .5f); //(supportWidth - width).maxComp() * .5f); // Don't get a 0 here
        distanceToBBox /= distanceFactor; // Make it depending on the support area
        float distanceFalloff = interpolation::wyvill(std::clamp(distanceToBBox, 0.f, 1.f));
        float evaluation = std::clamp(distanceFalloff, 0.f, 1.f); // Maybe...
        return evaluation; // - 0.5f ?
        /*Vector3 boundaries = Vector3(width, depth, height);
        float distToRect = -Vector3::signedDistanceToBoundaries(pos, boundaries * 0.f, boundaries * 1.f);
        float iso = distToRect / boundaries.maxComp();
        return std::clamp(iso + .5f, 0.f, 1.f);*/
    };
//    return ImplicitPatch::convert2DfunctionTo3Dfunction([height](Vector3 pos) { return height; }); // Using the constant 2D function
}

std::function<float (Vector3)> ImplicitPatch::createGaussianFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([width, depth, sigma, height](Vector3 pos) { return normalizedGaussian(Vector3(width, depth, 0), pos.xy(), sigma) * height; });
}

std::function<float (Vector3)> ImplicitPatch::createRockFunction(float sigma, float width, float depth, float height)
{
    auto sphereFunction = ImplicitPatch::createSphereFunction(sigma, width, depth, height);
    std::function rockFunction = [sphereFunction, sigma, width, depth, height] (Vector3 pos) {
        FastNoiseLite noise;
//        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseValue = (noise.GetNoise(pos.x * (100.f / width) + noiseOffset, pos.y * (100.f / width) + noiseOffset, pos.z * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
        return sphereFunction(pos + Vector3(0, 0, width * .25f)) - (noiseValue * .5f);
    };
    return rockFunction;
}

std::function<float (Vector3)> ImplicitPatch::createMountainFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([width, depth, height](Vector3 pos) {
        Vector3 translatedPos = pos - (Vector3(width, depth, 0) * .5f);
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(width, height, 1.f) * .5f);
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
    std::function caveFunc = [width, depth, height] (Vector3 pos) {
        Vector3 start = Vector3(width * .5f, depth * 1.f, height * 1.f);
        Vector3 p1 = Vector3(width * .4f, depth * .7f, height * .7f);
        Vector3 p2 = Vector3(width * .6f, depth * .4f, height * .4f);
        Vector3 end = Vector3(width * .5f, depth * .2f, height * .5f);
        BSpline curve = BSpline({start, p1, p2, end});
        float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed
        return 1.f - (curve.estimateDistanceFrom(pos, epsilon) / (width * .5f));
    };
    return caveFunc;
}

std::function<float (Vector3)> ImplicitPatch::createArchFunction(float sigma, float width, float depth, float height)
{
    std::function archFunc = [width, depth, height] (Vector3 pos) {
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

#else
int ImplicitPatch::currentMaxIndex = -1;
ImplicitPatch::ImplicitPatch()
    : ImplicitPatch(Vector3(), Vector3(), Vector3(), [](Vector3 _pos) {return 0.f; })
{

}

ImplicitPatch::ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, evalFunction, [&](Vector3 pos) { return this->evaluate(pos); }, densityValue)
{
    this->compositionFunction = [&](Vector3 pos) {
        float eval = this->evaluate(pos);
        return eval;
    };
}

ImplicitPatch::ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, evalFunction, compositionFunction, [&](Vector3 _pos) { return 1.f; }, densityValue)
{
    this->alphaDistanceFunction = [&](Vector3 _pos) {
        float distPower = 2.f;
        float minDistance = 0.8f;
        float distance = -Vector3::signedDistanceToBoundaries(_pos, this->position/*boundMin*/, this->dimension);
        float maxDistance = ((this->dimension - this->position/*boundMin*/) * .5f).maxComp();
        float ratio = 1.f - std::clamp(distance / maxDistance, 0.f, 1.f);
        return 1.f - (ratio <= minDistance ? 0.f : std::pow((ratio - minDistance) / (1.f - minDistance), distPower));
    };
}

ImplicitPatch::ImplicitPatch(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : position(pos)/*, boundMin(boundMin)*/, dimension(boundMax), evalFunction(evalFunction), densityValue(densityValue), compositionFunction(compositionFunction), alphaDistanceFunction(alphaDistanceFunction)
{
/*    this->_cachedMaxHeight = Matrix3<float>(boundMax - boundMin, -1000.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;
*/
    ImplicitPatch::currentMaxIndex ++;
    this->index = ImplicitPatch::currentMaxIndex;
    this->initCachedAttributes();
}

ImplicitPatch::~ImplicitPatch()
{
}

void ImplicitPatch::deleteComposables()
{
    if (this->composableA != nullptr) {
        this->composableA->deleteComposables();
        delete this->composableA;
    } if (this->composableB != nullptr) {
        this->composableB->deleteComposables();
        delete this->composableB;
    }
}

ImplicitPatch ImplicitPatch::clone()
{
    ImplicitPatch copy = *this;
    ImplicitPatch::currentMaxIndex ++;
    copy.index = ImplicitPatch::currentMaxIndex;
    return copy;
}

void ImplicitPatch::cleanCache()
{
    this->_cachedMaxHeight.reset(-1000.f);
    this->_cachedMinHeight.reset(-1000.f);
    if (this->composableA)
        this->composableA->cleanCache();
    if (this->composableB)
        this->composableB->cleanCache();
}

float ImplicitPatch::getMaxHeight(Vector3 position)
{
    if (this->_cachedMaxHeight.at(position.xy()) < -1.f) { // Unknown height yet
        if (this->compositionOperation == NONE) {
            float resolution = 1.f;
            float maxHeight = this->dimension.z;
            // TODO: check if the next condition is sufficint / optimal
            while (maxHeight > 0.f/*(maxHeight >= this->position.z || maxHeight >= position.z)*/ && evaluate(Vector3(position.x, position.y, maxHeight)) < .5f)
                maxHeight -= resolution;
            this->_cachedMaxHeight.at(position.xy()) = maxHeight /*- std::min(this->position.z, position.z)*/;
        }
        else if (this->compositionOperation == STACK) {
            this->_cachedMaxHeight.at(position.xy()) = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)) + this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position));
        }
        else if (this->compositionOperation == BLEND) {
            // TODO : CORRECT THIS (?)
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }
        else if (this->compositionOperation == REPLACE) {
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }/*
        else if (this->compositionOperation == NEG_STACKING) {
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }
        else if (this->compositionOperation == STACK_IN_WATER) {
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }*/
    }
    return this->_cachedMaxHeight.at(position.xy());
}

float ImplicitPatch::getMinHeight(Vector3 position)
{
    if (this->_cachedMinHeight.at(position.xy()) < -1.f) { // Unknown height yet
        if (this->compositionOperation == NONE) {
            float resolution = 1.f;
            float minHeight = this->position.z;
            // TODO: check if the next condition is sufficint / optimal
            while (minHeight < this->getDimensions().z /*(maxHeight >= this->position.z || maxHeight >= position.z)*/ && evaluate(Vector3(position.x, position.y, minHeight)) < .5f)
                minHeight += resolution;
            if (minHeight < this->getDimensions().z)
                this->_cachedMinHeight.at(position.xy()) = 100000.f;
            else
                this->_cachedMinHeight.at(position.xy()) = minHeight /*- std::min(this->position.z, position.z)*/;
        } else {
            float minHeightA = this->composableA->getMinHeight(position);
            if (minHeightA < this->composableA->getDimensions().z)
                this->_cachedMinHeight.at(position.xy()) = minHeightA;
            else
                this->_cachedMinHeight.at(position.xy()) = this->composableA->position.z + this->composableA->getDimensions().z + this->composableB->getMinHeight(position);
        }
        /*
        else if (this->compositionOperation == STACK) {
            float minHeightA = this->composableA->getMinHeight(position);
            if (minHeightA < this->composableA->getDimensions().z)
                this->_cachedMinHeight.at(position.xy()) = minHeightA;
            else
                this->_cachedMinHeight.at(position.xy()) = this->composableB->getMinHeight(position);
        }
        else if (this->compositionOperation == BLEND) {
            // TODO : CORRECT THIS (?)
            this->_cachedMinHeight.at(position.xy()) = std::max(this->composableA->getMinHeight(position), this->composableB->getMinHeight(position));
        }
        else if (this->compositionOperation == REPLACE) {
            this->_cachedMinHeight.at(position.xy()) = std::max(this->composableA->getMinHeight(position), this->composableB->getMinHeight(position));
        }
        else if (this->compositionOperation == NEG_STACKING) {
            this->_cachedMinHeight.at(position.xy()) = std::max(this->composableA->getMinHeight(position), this->composableB->getMinHeight(position));
        }
        else if (this->compositionOperation == STACK_IN_WATER) {
            this->_cachedMinHeight.at(position.xy()) = std::max(this->composableA->getMinHeight(position), this->composableB->getMinHeight(position));
        }*/
    }
    return this->_cachedMinHeight.at(position.xy());
}

float ImplicitPatch::evaluate(Vector3 position)
{
//    position -= this->position;
    if (!Vector3::isInBox(position, Vector3() - dimension, 2.f * dimension))
        return 0.f;

    float evaluationValue = 0.f;
    if (this->compositionOperation == NONE) {
        evaluationValue =  this->evalFunction(position);
        Vector3 normalizedPosition = position / this->getDimensions();
        Vector3 centeredPosition = normalizedPosition - Vector3(.5f, .5f, .5f);
        Vector3 absPos = centeredPosition.abs();
        float maxComp = absPos.maxComp();
        float distance = std::clamp(maxComp - .5f, 0.f, 1.f);
        float distanceFalloff = interpolation::wyvill(distance);
        float isoFalloff = distanceFalloff;
        evaluationValue = evaluationValue * isoFalloff;
    } else {
        Vector3 posA = this->getPositionToEvaluateComposableA(position);
        Vector3 posB = this->getPositionToEvaluateComposableB(position);
        float evalA = this->composableA->evaluate(posA);
        float evalB = this->composableB->evaluate(posB);
        evaluationValue = this->evaluateFromValues(evalA, evalB);
    }
    return std::clamp(evaluationValue, 0.f, 1.f);
}

float ImplicitPatch::evaluateFromValues(float evaluationA, float evaluationB)
{
    float evaluationValue = 0.f;
    switch(this->compositionOperation) {
        case CompositionFunction::STACK: {
            evaluationValue = std::max(evaluationA, evaluationB);
            break;
        }
        case CompositionFunction::BLEND: {
            evaluationValue = std::pow(std::pow(evaluationA, this->blendingFactor) + std::pow(evaluationB, this->blendingFactor), 1.f/this->blendingFactor);
            break;
        }
        case CompositionFunction::REPLACE: {
            evaluationValue = (evaluationB >= 0.5f ? evaluationB : evaluationA);
            break;
        }/*
        case CompositionFunction::NEG_STACKING: {
            evaluationValue = (evaluationB >= 0.5f ? evaluationB : evaluationA);
            break;
        }
        case CompositionFunction::STACK_IN_WATER: {
            evaluationValue = (evaluationB >= 0.5f ? evaluationB : evaluationA);
            break;
        }*/
        case CompositionFunction::NONE: {
            break; // Cannot pass here
        }
    }
    return std::clamp(evaluationValue, 0.f, 1.f);
}

std::map<float, float> ImplicitPatch::getDensityAtPosition(Vector3 position, Vector3 worldPosition)
{
    if (!worldPosition.isValid())
        worldPosition = position;

    Vector3 minPos = this->getMinPosition();
    Vector3 maxPos = this->getMaxPosition();
    Vector3 margin = maxPos - minPos;
    if (!Vector3::isInBox(worldPosition, minPos - margin, maxPos + margin)) {
        return {};
    }

    std::map<float, float> densityAndEvalA;
    std::map<float, float> densityAndEvalB;
    std::map<float, float> finalEvalutation;

    if (this->compositionOperation == NONE) { // Not a composition => one unique material to return
        finalEvalutation = {{this->densityValue, this->evaluate(position)}};
    } else { // If it's a composite, need to merge the result of the two branches
        float totalA = 0.f;
        float totalB = 0.f;
//        float densityA = 0.f; // Should not be necessary
//        float densityB = 0.f; // Should not be necessary
        if (this->composableA && !this->composableA->isIdentity()) {
            densityAndEvalA = this->composableA->getDensityAtPosition(this->getPositionToEvaluateComposableA(position), worldPosition);
            for (const auto& [densA, valueA] : densityAndEvalA) {
                totalA += valueA;
//                densityA = densA; // Should not be necessary
            }
        }
        if (this->composableB && !this->composableB->isIdentity()) {
            densityAndEvalB = this->composableB->getDensityAtPosition(this->getPositionToEvaluateComposableB(position), worldPosition);
            for (const auto& [densB, valueB] : densityAndEvalB) {
                totalB += valueB;
//                densityB = densB; // Should not be necessary
            }
        }
        float totalQuantity = totalA + totalB;
        if (totalQuantity <= 0.f) {
            finalEvalutation = {{-1.f, 0.f}};
        } else {
            // Multiply the densities in patches A and B to match the composition value after evaluation
            float myEval = this->evaluateFromValues(totalA, totalB);
            float ratio = myEval / totalQuantity;
            for (const auto& [densA, valueA] : densityAndEvalA) {
                if (finalEvalutation.count(densA) == 0)
                    finalEvalutation[densA] = 0.f;
                finalEvalutation[densA] += valueA * ratio;
            }
            for (const auto& [densB, valueB] : densityAndEvalB) {
                if (finalEvalutation.count(densB) == 0)
                    finalEvalutation[densB] = 0.f;
                finalEvalutation[densB] += valueB * ratio;
            }
    //        finalEvalutation = {{(totalA > totalB ? densityA : densityB), myEval}};
                    /*
            if (myEval >= 0.5f) {
                finalEvalutation = {{(totalA > totalB ? densityA : densityB), myEval}};
            } else {
                finalEvalutation = {{-1.f, myEval}};
            }*/
            /*
            for (const auto& [densA, valueA] : densityAndEvalA) {
                if (finalEvalutation.count(densA) == 0)
                    finalEvalutation[densA] = 0;
                finalEvalutation[densA] += valueA;
            }
            for (const auto& [densB, valueB] : densityAndEvalB) {
                if (finalEvalutation.count(densB) == 0)
                    finalEvalutation[densB] = 0;
                finalEvalutation[densB] += valueB;
            }
            // TODO : CORRECT THIS FUNCTION, PLEASE!!
            float ourEvaluation = 1.f; //this->evaluateFromValues(totalA, totalB);
            for (auto& [dens, val] : finalEvalutation) {
                val *= ourEvaluation;
            }*/
        }
    }
    return finalEvalutation;
}

Vector3 ImplicitPatch::getPositionToEvaluateComposableA(Vector3 pos)
{
    return pos - (this->composableA->position + this->position);
//    return pos; // No modification needed
}

Vector3 ImplicitPatch::getPositionToEvaluateComposableB(Vector3 pos)
{
    if (this->compositionOperation == NONE) {
        // Doesn't really make sense to get here
    } else {
        if (this->positioning == ABOVE) {
            float heightA = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(pos));
            pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
            pos.z -= heightA;
        } else if (this->positioning == INSIDE_TOP) {
            float heightA = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(pos));
            pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
            pos.z -= heightA;
//            pos.z -= this->composableB->getMaxHeight(pos);
        } else if (this->positioning == INSIDE_BOTTOM) {
            float heightA = this->composableA->getMinHeight(this->getPositionToEvaluateComposableA(pos));
            pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
            pos.z -= heightA;
        } else if (this->positioning == FIXED) {
            // Nothing to do
        }
    }
    /*
    if (this->compositionOperation == NONE) {
        // Doesn't really make sense to get here
    } else if (this->compositionOperation == STACK) {
        float heightA = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(pos));
        pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
        pos.z -= heightA;
    } else if (this->compositionOperation == BLEND) {
        // Nothing to do
    } else if (this->compositionOperation == REPLACE) {
        // Nothing to do
    }
    */
    /*else if (this->compositionOperation == NEG_STACKING) {
        float heightA = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(pos));
        pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
        pos.z -= (heightA + this->composableB->getDimensions().z);
    } else if (this->compositionOperation == STACK_IN_WATER) {
        float heightA = this->composableA->getMinHeight(this->getPositionToEvaluateComposableA(pos));
        pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
        pos.z -= heightA;
    }*/
    return pos - (this->composableB->position + this->position);
}

Vector3 ImplicitPatch::getMinPosition()
{
    if (!this->isOperation()) {
        if (this->isIdentity())
            return Vector3(false);
        return this->position;
    } else {
        Vector3 minA = this->composableA->getMinPosition();
        Vector3 minB = this->composableB->getMinPosition();
        if (minA.isValid() && minB.isValid()) {
            return Vector3::min(minA, minB);
        } else {
            if (!minA.isValid() && minB.isValid()) {
                return minB;
            } else if (minA.isValid() && !minB.isValid()) {
                return minA;
            } else {
                return Vector3(false);
            }
        }
    }

}

Vector3 ImplicitPatch::getMaxPosition()
{
    if (!this->isOperation()) {
        if (this->isIdentity())
            return Vector3(false);
        return this->position + this->dimension;
    } else {
        Vector3 maxA = this->composableA->getMaxPosition();
        Vector3 maxB = this->composableB->getMaxPosition();
        if (maxA.isValid() && maxB.isValid()) {
            return Vector3::max(maxA, maxB);
        } else {
            if (!maxA.isValid() && maxB.isValid()) {
                return maxB;
            } else if (maxA.isValid() && !maxB.isValid()) {
                return maxA;
            } else {
                return Vector3(false);
            }
        }
    }

}

void ImplicitPatch::defineFunctionsBasedOnPredefinedShape()
{
    if (this->shape == Sphere)
        this->evalFunction = ImplicitPatch::createSphereFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Block)
        this->evalFunction = ImplicitPatch::createBlockFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Gaussian)
        this->evalFunction = ImplicitPatch::createGaussianFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Rock)
        this->evalFunction = ImplicitPatch::createRockFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Mountain)
        this->evalFunction = ImplicitPatch::createMountainFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Dune)
        this->evalFunction = ImplicitPatch::createDuneFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Basin)
        this->evalFunction = ImplicitPatch::createBasinFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Cave)
        this->evalFunction = ImplicitPatch::createCaveFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
    if (this->shape == Arch)
        this->evalFunction = ImplicitPatch::createArchFunction(this->sigmaValue, this->dimension.x, this->dimension.y, this->dimension.z);
}

nlohmann::json ImplicitPatch::toJson()
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content[ImplicitPatch::json_identifier]);

        Vector3 patchOffset = this->position - oldValues->position;
        Vector3 patchRescale = this->getDimensions() / oldValues->getDimensions();
        content["path"] = this->used_json_filename;
        content["offset"] = vec3_to_json(patchOffset);
        content["scale"] = vec3_to_json(patchRescale);
    } else {
        content["name"] = this->name;
        content["operator"] = stringFromCompositionOperation(this->compositionOperation);
        content["blendingFactor"] = this->blendingFactor;
        content["position"] = vec3_to_json(this->position);
        content["dimension"] = vec3_to_json(this->dimension);
        content["densityValue"] = this->densityValue;
        content["index"] = this->index;
        content["shape"] = stringFromPredefinedShape(this->shape);
        content["sigmaValue"] = this->sigmaValue;

        if (this->composableA)
            content["composableA"] = this->composableA->toJson();
        if (this->composableB)
            content["composableB"] = this->composableB->toJson();
    }

    return content;
}

ImplicitPatch *ImplicitPatch::fromJson(nlohmann::json json_content)
{
    if (json_content.contains("path")) {
        std::string filename = json_content["path"];
        std::ifstream file(filename);
        nlohmann::json sub_json_content = nlohmann::json::parse(file);
        if (sub_json_content.contains(ImplicitPatch::json_identifier)) {
            Vector3 patchOffset = Vector3(0, 0, 0);
            Vector3 patchRescale = Vector3(1, 1, 1);
            if (sub_json_content["implicitPatch"].contains("offset")) {
                patchOffset = json_to_vec3(sub_json_content["implicitPatch"]["offset"]);
            }
            if (sub_json_content["implicitPatch"].contains("scale")) {
                patchRescale = json_to_vec3(sub_json_content["implicitPatch"]["scale"]);
            }
            ImplicitPatch* result = ImplicitPatch::fromJson(sub_json_content[ImplicitPatch::json_identifier]);
            result->used_json_filename = filename; // Save the file in this object
            result->dimension *= patchRescale;
            result->position += patchOffset;
            return result;
        } else {
            std::cerr << "No patches defined in file " << filename << std::endl;
            return nullptr;
        }
    } else {
        ImplicitPatch* patch = new ImplicitPatch;
        patch->name = json_content["name"];
        patch->compositionOperation = compositionOperationFromString(json_content["operator"]);
        patch->blendingFactor = json_content["blendingFactor"];
        patch->position = json_to_vec3(json_content["position"]);
        patch->dimension = json_to_vec3(json_content["dimension"]);
        patch->densityValue = json_content["densityValue"];
        patch->index = json_content["index"];
        patch->shape = predefinedShapeFromString(json_content["shape"]);
        patch->sigmaValue = json_content["sigmaValue"];

        if (json_content.contains("composableA"))
            patch->composableA = ImplicitPatch::fromJson(json_content["composableA"]);
        if (json_content.contains("composableB"))
            patch->composableB = ImplicitPatch::fromJson(json_content["composableB"]);
        patch->defineFunctionsBasedOnPredefinedShape();

        return patch;
    }

}

std::string ImplicitPatch::toString()
{
    if (this->used_json_filename.empty()) {
        return this->name + " #" + std::to_string(this->index);
    } else {
        // Gives the filename if it's coming from a JSON file
        return QString::fromStdString(this->used_json_filename).split("\\").last().split("/").last().split(".json").front().toStdString();
    }
}

ImplicitPatch ImplicitPatch::createStack(ImplicitPatch *P1, ImplicitPatch *P2, PositionalLabel pos)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::STACK;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(P1->boundMin - relativePos1, P2->boundMin - relativePos2);
    newPatch.dimension = Vector3::max(P1->dimension + relativePos1, P2->dimension + relativePos2);
    newPatch.dimension.z = P1->dimension.z + P2->dimension.z;
    newPatch.initCachedAttributes(P1, P2, pos);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createReplacement(ImplicitPatch *P1, ImplicitPatch *P2, PositionalLabel pos)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::REPLACE;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2, pos);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createBlending(ImplicitPatch *P1, ImplicitPatch *P2, PositionalLabel pos)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::BLEND;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2, pos);
    return newPatch;
}

/*ImplicitPatch ImplicitPatch::createNegStacking(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::NEG_STACKING;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createStackInWater(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::STACK_IN_WATER;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2);
    return newPatch;
}*/

void ImplicitPatch::initCachedAttributes(ImplicitPatch *composableA, ImplicitPatch *composableB, PositionalLabel pos)
{
    this->positioning = pos;
    this->composableA = composableA;
    this->composableB = composableB;
    this->_cachedMaxHeight = Matrix3<float>(this->getDimensions().x, this->getDimensions().y, 1, -1000.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;
    this->_cachedMinHeight = Matrix3<float>(this->getDimensions().x, this->getDimensions().y, 1, -1000.f);
    this->_cachedMinHeight.raiseErrorOnBadCoord = false;
    this->_cachedMinHeight.defaultValueOnBadCoord = 0.f;
}

std::function<float (Vector3)> ImplicitPatch::createSphereFunction(float sigma, float width, float depth, float height)
{
    std::function sphere = [width](Vector3 pos) {
        Vector3 center = Vector3(width, width, width) * .5f;
        float sqrDist = (pos - center).norm2();
        float value = std::max(0.f, (sqrDist > width * width ? 0.f : 1 - std::abs(std::sqrt(sqrDist)) / (width)));
        if ((pos.xy() - center.xy()).norm2() < 5.f) {
            int a = 0;
        }
        return value;
    };
    return sphere;
}

std::function<float (Vector3)> ImplicitPatch::createBlockFunction(float sigma, float width, float depth, float height)
{
    return [sigma, width, depth, height] (Vector3 pos) {
        Vector3 boundaries = Vector3(width, depth, height);
        float distToRect = -Vector3::signedDistanceToBoundaries(pos, boundaries * 0.f, boundaries * 1.f);
        float iso = distToRect / boundaries.maxComp();
        return std::clamp(iso + .5f, 0.f, 1.f);
    };
//    return ImplicitPatch::convert2DfunctionTo3Dfunction([height](Vector3 pos) { return height; }); // Using the constant 2D function
}

std::function<float (Vector3)> ImplicitPatch::createGaussianFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([width, depth, sigma, height](Vector3 pos) { return normalizedGaussian(Vector3(width, depth, 0), pos.xy(), sigma) * height; });
}

std::function<float (Vector3)> ImplicitPatch::createRockFunction(float sigma, float width, float depth, float height)
{
    auto sphereFunction = ImplicitPatch::createSphereFunction(sigma, width, depth, height);
    std::function rockFunction = [sphereFunction, sigma, width, depth, height] (Vector3 pos) {
        FastNoiseLite noise;
//        noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise.SetFractalType(FastNoiseLite::FractalType_FBm);
        float noiseOffset = sigma * 1000.f; // noise.GetNoise(sigma * 1000.f, sigma * 1000.f);
        float noiseValue = (noise.GetNoise(pos.x * (100.f / width) + noiseOffset, pos.y * (100.f / width) + noiseOffset, pos.z * (100.f / width) + noiseOffset) + 1.f) * .5f; // Between 0 and 1
        return sphereFunction(pos + Vector3(0, 0, width * .5f)) - (noiseValue * .5f);
    };
    return rockFunction;
}

std::function<float (Vector3)> ImplicitPatch::createMountainFunction(float sigma, float width, float depth, float height)
{
    return ImplicitPatch::convert2DfunctionTo3Dfunction([width, depth, height](Vector3 pos) {
        Vector3 translatedPos = pos - (Vector3(width, depth, 0) * .5f);
        Vector3 normalizedPos = translatedPos.xy() / (Vector3(width, height, 1.f) * .5f);
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
    std::function caveFunc = [width, depth, height] (Vector3 pos) {
        Vector3 start = Vector3(width * .5f, depth * 1.f, height * 1.f);
        Vector3 p1 = Vector3(width * .4f, depth * .7f, height * .7f);
        Vector3 p2 = Vector3(width * .6f, depth * .4f, height * .4f);
        Vector3 end = Vector3(width * .5f, depth * .2f, height * .5f);
        BSpline curve = BSpline({start, p1, p2, end});
        float epsilon = 1e-1; // Keep a coarse epsilon just to speed up the process, no real precision needed
        return 1.f - (curve.estimateDistanceFrom(pos, epsilon) / (width * .5f));
    };
    return caveFunc;
}

std::function<float (Vector3)> ImplicitPatch::createArchFunction(float sigma, float width, float depth, float height)
{
    std::function archFunc = [width, depth, height] (Vector3 pos) {
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

#endif


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
    else if (name == "NONE")
        return ImplicitPatch::PredefinedShapes::None;

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
    else if (shape == ImplicitPatch::PredefinedShapes::None)
        return "None";

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

    // Impossible
    std::cerr << "Unknown positional label given." << std::endl;
    return "";
}

TerrainTypes materialFromString(std::string name)
{
    name = toUpper(name);

    if (name == "AIR")
        return AIR;
    else if (name == "STRONG_WATER")
        return STRONG_WATER;
    else if (name == "LIGHT_WATER")
        return LIGHT_WATER;
    else if (name == "TURBULENT_WATER")
        return TURBULENT_WATER;
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
    else if (material == STRONG_WATER)
        return "STRONG_WATER";
    else if (material == LIGHT_WATER)
        return "LIGHT_WATER";
    else if (material == TURBULENT_WATER)
        return "TURBULENT_WATER";
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
