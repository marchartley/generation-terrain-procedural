#include "ImplicitPatch.h"

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
    this->_cachedMaxHeight = Matrix3<float>(boundMax - boundMin, -1000.f);
    this->_cachedMaxHeight.raiseErrorOnBadCoord = false;
    this->_cachedMaxHeight.defaultValueOnBadCoord = 0.f;

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
        }
        else if (this->compositionOperation == NEG_STACKING) {
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }
        else if (this->compositionOperation == STACK_IN_WATER) {
            this->_cachedMaxHeight.at(position.xy()) = std::max(this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(position)), this->composableB->getMaxHeight(this->getPositionToEvaluateComposableB(position)));
        }
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
        if ((position.xy() - Vector3(15, 15)).norm2() < 3.f) {
            int a =0;
        }
        if ((position.xy() - Vector3(25, 25)).norm2() < 3.f) {
            int a =0;
        }
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
        }
        case CompositionFunction::NEG_STACKING: {
            evaluationValue = (evaluationB >= 0.5f ? evaluationB : evaluationA);
            break;
        }
        case CompositionFunction::STACK_IN_WATER: {
            evaluationValue = (evaluationB >= 0.5f ? evaluationB : evaluationA);
            break;
        }
        case CompositionFunction::NONE: {
            break; // Cannot pass here
        }
    }
    return std::clamp(evaluationValue, 0.f, 1.f);
}

std::map<float, float> ImplicitPatch::getDensityAtPosition(Vector3 position)
{
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
            densityAndEvalA = this->composableA->getDensityAtPosition(this->getPositionToEvaluateComposableA(position));
            for (const auto& [densA, valueA] : densityAndEvalA) {
                totalA += valueA;
//                densityA = densA; // Should not be necessary
            }
        }
        if (this->composableB && !this->composableB->isIdentity()) {
            densityAndEvalB = this->composableB->getDensityAtPosition(this->getPositionToEvaluateComposableB(position));
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
    return pos - (this->composableA->position - this->position);
//    return pos; // No modification needed
}

Vector3 ImplicitPatch::getPositionToEvaluateComposableB(Vector3 pos)
{
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
    } else if (this->compositionOperation == NEG_STACKING) {
        float heightA = this->composableA->getMaxHeight(this->getPositionToEvaluateComposableA(pos));
        pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
        pos.z -= (heightA + this->composableB->getDimensions().z);
    } else if (this->compositionOperation == STACK_IN_WATER) {
        float heightA = this->composableA->getMinHeight(this->getPositionToEvaluateComposableA(pos));
        pos.z +=  this->composableB->position.z; // Stick the patch to the ground.
        pos.z -= heightA;
    }
    return pos - (this->composableB->position - this->position);
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

nlohmann::json ImplicitPatch::toJson() const
{
    nlohmann::json content;
    if (!this->used_json_filename.empty()) {
        std::ifstream file(this->used_json_filename);
        nlohmann::json old_json_content = nlohmann::json::parse(file);
        ImplicitPatch* oldValues = ImplicitPatch::fromJson(old_json_content["implicitPatches"]);

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
        if (sub_json_content.contains("implicitPatches")) {
            Vector3 patchOffset = Vector3(0, 0, 0);
            Vector3 patchRescale = Vector3(1, 1, 1);
            if (sub_json_content["implicitPatch"].contains("offset")) {
                patchOffset = json_to_vec3(sub_json_content["implicitPatch"]["offset"]);
            }
            if (sub_json_content["implicitPatch"].contains("scale")) {
                patchRescale = json_to_vec3(sub_json_content["implicitPatch"]["scale"]);
            }
            ImplicitPatch* result = ImplicitPatch::fromJson(sub_json_content["implicitPatches"]);
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

ImplicitPatch ImplicitPatch::createStack(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::STACK;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(P1->boundMin - relativePos1, P2->boundMin - relativePos2);
    newPatch.dimension = Vector3::max(P1->dimension + relativePos1, P2->dimension + relativePos2);
    newPatch.dimension.z = P1->dimension.z + P2->dimension.z;
    newPatch.initCachedAttributes(P1, P2);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createReplacement(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::REPLACE;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createBlending(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.compositionOperation = CompositionFunction::BLEND;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch.initCachedAttributes(P1, P2);
    return newPatch;
}

ImplicitPatch ImplicitPatch::createNegStacking(ImplicitPatch *P1, ImplicitPatch *P2)
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
}

void ImplicitPatch::initCachedAttributes(ImplicitPatch *composableA, ImplicitPatch *composableB)
{
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



ImplicitPatch::CompositionFunction compositionOperationFromString(std::string name)
{
    name = toUpper(name);
    if (name == "NONE")
        return ImplicitPatch::CompositionFunction::NONE;
    if (name == "STACK")
        return ImplicitPatch::CompositionFunction::STACK;
    if (name == "BLEND")
        return ImplicitPatch::CompositionFunction::BLEND;
    if (name == "REPLACE")
        return ImplicitPatch::CompositionFunction::REPLACE;
    if (name == "NEG_STACKING")
        return ImplicitPatch::CompositionFunction::NEG_STACKING;
    if (name == "STACK_IN_WATER")
        return ImplicitPatch::CompositionFunction::STACK_IN_WATER;

    // Error
    std::cerr << "Wrong operation string given" << std::endl;
    return ImplicitPatch::CompositionFunction::NONE;
}

std::string stringFromCompositionOperation(ImplicitPatch::CompositionFunction operation)
{
    if (operation == ImplicitPatch::CompositionFunction::NONE)
        return "None";
    if (operation == ImplicitPatch::CompositionFunction::STACK)
        return "Stack";
    if (operation == ImplicitPatch::CompositionFunction::BLEND)
        return "Blend";
    if (operation == ImplicitPatch::CompositionFunction::REPLACE)
        return "Replace";
    if (operation == ImplicitPatch::CompositionFunction::NEG_STACKING)
        return "Neg_stacking";
    if (operation == ImplicitPatch::CompositionFunction::STACK_IN_WATER)
        return "Stack_in_water";

    // Impossible
    std::cerr << "Unknown operation given." << std::endl;
    return "";
}

ImplicitPatch::PredefinedShapes predefinedShapeFromString(std::string name)
{
    name = toUpper(name);
    if (name == "NONE")
        return ImplicitPatch::PredefinedShapes::None;
    if (name == "SPHERE")
        return ImplicitPatch::PredefinedShapes::Sphere;
    if (name == "BLOCK")
        return ImplicitPatch::PredefinedShapes::Block;
    if (name == "GAUSSIAN")
        return ImplicitPatch::PredefinedShapes::Gaussian;
    if (name == "ROCK")
        return ImplicitPatch::PredefinedShapes::Rock;
    if (name == "MOUNTAIN")
        return ImplicitPatch::PredefinedShapes::Mountain;
    if (name == "DUNE")
        return ImplicitPatch::PredefinedShapes::Dune;
    if (name == "BASIN")
        return ImplicitPatch::PredefinedShapes::Basin;
    if (name == "CAVE")
        return ImplicitPatch::PredefinedShapes::Cave;
    if (name == "ARCH")
        return ImplicitPatch::PredefinedShapes::Arch;

    // Error
    std::cerr << "Wrong predefined shape string given" << std::endl;
    return ImplicitPatch::PredefinedShapes::None;
}

std::string stringFromPredefinedShape(ImplicitPatch::PredefinedShapes shape)
{
    if (shape == ImplicitPatch::PredefinedShapes::None)
        return "None";
    if (shape == ImplicitPatch::PredefinedShapes::Sphere)
        return "Sphere";
    if (shape == ImplicitPatch::PredefinedShapes::Block)
        return "Block";
    if (shape == ImplicitPatch::PredefinedShapes::Gaussian)
        return "Gaussian";
    if (shape == ImplicitPatch::PredefinedShapes::Rock)
        return "Rock";
    if (shape == ImplicitPatch::PredefinedShapes::Mountain)
        return "Mountain";
    if (shape == ImplicitPatch::PredefinedShapes::Dune)
        return "Dune";
    if (shape == ImplicitPatch::PredefinedShapes::Basin)
        return "Basin";
    if (shape == ImplicitPatch::PredefinedShapes::Cave)
        return "Cave";
    if (shape == ImplicitPatch::PredefinedShapes::Arch)
        return "Arch";

    // Impossible
    std::cerr << "Unknown predefined shape given." << std::endl;
    return "";
}
