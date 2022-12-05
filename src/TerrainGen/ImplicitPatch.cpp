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
    this->_cachedHeights = Matrix3<float>(boundMax - boundMin, -1000.f);
    this->_cachedHeights.raiseErrorOnBadCoord = false;
    this->_cachedHeights.defaultValueOnBadCoord = 0.f;

    ImplicitPatch::currentMaxIndex ++;
    this->index = ImplicitPatch::currentMaxIndex;
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

float ImplicitPatch::getMaxHeight(Vector3 position)
{
    if (this->_cachedHeights.at(position.xy()) < -1.f || true) { // Unknown height yet
        if (this->compositionOperation == NONE) {
            float resolution = 1.f;
            float maxHeight = this->dimension.z;
            // TODO: check if the next condition is sufficint / optimal
            while (maxHeight > 0.f/*(maxHeight >= this->position.z || maxHeight >= position.z)*/ && evaluate(Vector3(position.x, position.y, maxHeight)) < .5f)
                maxHeight -= resolution;
            this->_cachedHeights.at(position.xy()) = maxHeight /*- std::min(this->position.z, position.z)*/;
        }
        else if (this->compositionOperation == STACK) {
            this->_cachedHeights.at(position.xy()) = this->composableA->getMaxHeight(position) + this->composableB->getMaxHeight(position);
        }
        else if (this->compositionOperation == BLEND) {
            // TODO : CORRECT THIS (?)
            this->_cachedHeights.at(position.xy()) = std::max(this->composableA->getMaxHeight(position), this->composableB->getMaxHeight(position));
        }
        else if (this->compositionOperation == REPLACE) {
            this->_cachedHeights.at(position.xy()) = std::max(this->composableA->getMaxHeight(position), this->composableB->getMaxHeight(position));
        }
    }
    return this->_cachedHeights.at(position.xy());
}

float ImplicitPatch::evaluate(Vector3 position)
{
    position -= this->position;
    if (!Vector3::isInBox(position, Vector3() - dimension, 2.f * dimension))
        return 0.f;

    float evaluationValue = 0.f;
    if (this->compositionOperation == NONE) {
        evaluationValue =  this->evalFunction(position);
        float distance = std::clamp(((position / this->getDimensions()) - Vector3(.5, .5, .5)).abs().maxComp() - .5f, 0.f, 1.f);
        float distanceFalloff = interpolation::wyvill(/*1.f - */distance);
        evaluationValue = evaluationValue * distanceFalloff;
    } else {
        float evalA = this->composableA->evaluate(this->getPositionToEvaluateComposableA(position));
        float evalB = this->composableB->evaluate(this->getPositionToEvaluateComposableB(position));
        evaluationValue = this->evaluateFromValues(evalA, evalB);
        /*
        switch(this->compositionOperation) {
            case CompositionFunction::STACK: {
                evaluationValue = evalA + evalB;
                break;
            }
            case CompositionFunction::BLEND: {
                float blendPower = 5.f;
                evaluationValue = std::pow(std::pow(evalA, blendPower) + std::pow(evalB, blendPower), 1.f/blendPower);
                break;
            }
            case CompositionFunction::REPLACE: {
                evaluationValue = (evalB >= 0.5f ? evalB : evalA);
                break;
            }
            case CompositionFunction::NONE: {
                break; // Cannot pass here
            }
        }*/
    }
    return std::clamp(evaluationValue, 0.f, 1.f);
//    return interpolation::wyvill(1.f - std::clamp(evaluationValue, 0.f, 1.f));
    /*else if (this->compositionOperation == STACK) {
        Vector3 posA =  position; // - (this->composableA->position - this->position);
        Vector3 posB =  position; // - (this->composableB->position - this->position);
        posB.z +=  this->composableB->position.z; // Stick the patch to the ground.
        float evalA = this->composableA->evaluate(posA);
        float heightA = this->composableA->getMaxHeight(posA);
        posB.z -= heightA;
        float evalB = this->composableB->evaluate(posB);
        evaluationValue =  evalA + evalB;
    } else if (this->compositionOperation == BLEND) {
        float blendPower = 2.f;
        Vector3 posA =  position; // - (this->composableA->position - this->position);
        Vector3 posB =  position; // - (this->composableB->position - this->position);
        float evalA = this->composableA->evaluate(posA);
        float evalB = this->composableB->evaluate(posB);
        evaluationValue = std::pow(std::pow(evalA, blendPower) + std::pow(evalB, blendPower), 1.f/blendPower);
//        evaluationValue =  this->composableA->evaluate(position) + this->composableB->evaluate(position - Vector3(0, 0, this->composableA->getMaxHeight(position)));
    } else if (this->compositionOperation == REPLACE) {
        Vector3 posA =  position; // - (this->composableA->position - this->position);
        Vector3 posB =  position; // - (this->composableB->position - this->position);
        float evalB = this->composableB->evaluate(posB);
        evaluationValue =  (evalB >= 0.5f ? evalB : this->composableA->evaluate(posA));
    }
    return interpolation::wyvill(1.f - std::clamp(evaluationValue, 0.f, 1.f));
    */
}

float ImplicitPatch::evaluateFromValues(float evaluationA, float evaluationB)
{
    float evaluationValue = 0.f;
    switch(this->compositionOperation) {
        case CompositionFunction::STACK: {
            evaluationValue = evaluationA + evaluationB;
            break;
        }
        case CompositionFunction::BLEND: {
            float blendPower = 5.f;
            evaluationValue = std::pow(std::pow(evaluationA, blendPower) + std::pow(evaluationB, blendPower), 1.f/blendPower);
            break;
        }
        case CompositionFunction::REPLACE: {
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
        float densityA = 0.f; // Should not be necessary
        float densityB = 0.f; // Should not be necessary
        if (this->composableA && !this->composableA->isIdentity()) {
            densityAndEvalA = this->composableA->getDensityAtPosition(this->getPositionToEvaluateComposableA(position));
            for (const auto& [densA, valueA] : densityAndEvalA) {
                totalA += valueA;
                densityA = densA; // Should not be necessary
            }
        }
        if (this->composableB && !this->composableB->isIdentity()) {
            densityAndEvalB = this->composableB->getDensityAtPosition(this->getPositionToEvaluateComposableB(position));
            for (const auto& [densB, valueB] : densityAndEvalB) {
                totalB += valueB;
                densityB = densB; // Should not be necessary
            }
        }
        float myEval = this->evaluateFromValues(totalA, totalB);
        if (myEval >= 0.5f) {
            finalEvalutation = {{(totalA > totalB ? densityA : densityB), myEval}};
        } else {
            finalEvalutation = {{-1.f, myEval}};
        }
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
    return finalEvalutation;
}

Vector3 ImplicitPatch::getPositionToEvaluateComposableA(Vector3 pos)
{
    return pos; // No modification needed
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
    }
    return pos;
}

ImplicitPatch ImplicitPatch::createStack(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.composableA = P1;
    newPatch.composableB = P2;
    newPatch.compositionOperation = CompositionFunction::STACK;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(P1->boundMin - relativePos1, P2->boundMin - relativePos2);
    newPatch.dimension = Vector3::max(P1->dimension + relativePos1, P2->dimension + relativePos2);
    newPatch.dimension.z = P1->dimension.z + P2->dimension.z;
    newPatch._cachedHeights = Matrix3<float>(newPatch.getDimensions().x, newPatch.getDimensions().y, 1, -1000.f);
    newPatch._cachedHeights.raiseErrorOnBadCoord = false;
    newPatch._cachedHeights.defaultValueOnBadCoord = 0.f;
    return newPatch;
}

ImplicitPatch ImplicitPatch::createReplacement(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.composableA = P1;
    newPatch.composableB = P2;
    newPatch.compositionOperation = CompositionFunction::REPLACE;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch._cachedHeights = Matrix3<float>(newPatch.getDimensions().x, newPatch.getDimensions().y, 1, -1000.f);
    newPatch._cachedHeights.raiseErrorOnBadCoord = false;
    newPatch._cachedHeights.defaultValueOnBadCoord = 0.f;
    return newPatch;
}

ImplicitPatch ImplicitPatch::createBlending(ImplicitPatch *P1, ImplicitPatch *P2)
{
    ImplicitPatch newPatch;
    newPatch.composableA = P1;
    newPatch.composableB = P2;
    newPatch.compositionOperation = CompositionFunction::BLEND;

    newPatch.position = Vector3(); // Vector3::min(P1->position, P2->position);
    Vector3 relativePos1 = (P1->position - newPatch.position);
    Vector3 relativePos2 = (P2->position - newPatch.position);
//    newPatch.boundMin = Vector3::min(relativePos1 + P1->boundMin, relativePos2 + P2->boundMin);
    newPatch.dimension = Vector3::max(relativePos1 + P1->dimension, relativePos2 + P2->dimension);
    newPatch._cachedHeights = Matrix3<float>(newPatch.getDimensions().x, newPatch.getDimensions().y, 1, -1000.f);
    newPatch._cachedHeights.raiseErrorOnBadCoord = false;
    newPatch._cachedHeights.defaultValueOnBadCoord = 0.f;
    return newPatch;
}

std::function<float (Vector3)> ImplicitPatch::createSphereFunction(float radius)
{
    std::function sphere = [radius](Vector3 pos) {
        Vector3 center(radius, radius, radius);
        float sqrDist = (pos - center).norm2();
        float value = std::max(0.f, (sqrDist > 4 * radius * radius ? 0.f : 1 - std::abs(std::sqrt(sqrDist)) / (2.f * radius)));
        return value;
    };
    return sphere;
}

std::function<float (Vector3)> ImplicitPatch::createBlockFunction(float width, float depth, float height)
{
    /*std::function levelFunction = [height](Vector3 pos) {
        return std::max(0.f, 1.f - std::abs(height * .5f - pos.z));
    };
    return levelFunction;*/
    return convert2DfunctionTo3Dfunction([height](Vector3 pos) { return height; }); // Using the constant 2D function
}

std::function<float (Vector3)> ImplicitPatch::createGaussianFunction(float sigma, float width, float depth, float height)
{
//    std::function gauss = [sigma, height, this](Vector3 pos) {
//        return (normalizedGaussian(this->functionSize.xy(), pos.xy(), sigma) * height) - pos.z;
//    };
    return convert2DfunctionTo3Dfunction([width, depth, sigma, height](Vector3 pos) { return normalizedGaussian(Vector3(width, depth, 0), pos.xy(), sigma) * height; });
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

//Patch::~Patch()
//{

//}

/*
Patch2D::Patch2D()
    : Patch2D(Vector3(0, 0, 0), Vector3(-1, -1, 0), Vector3(1, 1, 0), [](Vector3 _pos) {return 0.f; })
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, heightFunction, densityValue)
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, heightFunction, compositionFunction, densityValue)
{

}

Patch2D::Patch2D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> heightFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, heightFunction, compositionFunction, alphaDistanceFunction, densityValue)
{

}

float Patch2D::getMaxHeight(Vector3 position)
{
    Vector3 p = position; // - this->position;
    if (boundMin.x > p.x || p.x > dimension.x || boundMin.y > p.y || p.y > dimension.y)
        return 0.f;
//    if (this->shape.grow(.1f).translate(Vector3(.001f, 0.f)).inside(p))
    return this->evalFunction(p.xy());
//    return 0.f;
}

float Patch2D::evaluate(Vector3 position)
{
    Vector3 p = position; // - this->position;
//    if (!shape.inside(p.xy())) return 0.f;

    if (boundMin.x > p.x || p.x > dimension.x || boundMin.y > p.y || p.y > dimension.y)
        return 0.f;

    float val = this->getMaxHeight(p.xy());
    return (val > p.z ? 1.f : 0.f);
}

Patch3D::Patch3D()
    : Patch3D(Vector3(0, 0, 0), Vector3(-1, -1, -1), Vector3(1, 1, 1), [](Vector3 _pos) {return 0.f; })
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, evalFunction, densityValue)
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, evalFunction, compositionFunction, densityValue)
{

}

Patch3D::Patch3D(Vector3 pos, Vector3 boundMin, Vector3 boundMax, std::function<float (Vector3)> evalFunction, std::function<float (Vector3)> compositionFunction, std::function<float (Vector3)> alphaDistanceFunction, float densityValue)
    : ImplicitPatch(pos, boundMin, boundMax, evalFunction, compositionFunction, alphaDistanceFunction, densityValue)
{
}

Patch3D::Patch3D(const Patch3D &copy)
{
    this->position = copy.position;
    this->boundMin = copy.boundMin;
    this->dimension = copy.dimension;
    this->evalFunction = copy.evalFunction;
    this->compositionFunction = copy.compositionFunction;
    this->alphaDistanceFunction = copy.alphaDistanceFunction;
    this->densityValue = copy.densityValue;
    float test = compositionFunction(Vector3(10,10,10));
    std::cout << test << std::endl;
    this->_cachedHeights = Matrix3<float>(dimension - boundMin, -1000.f);
    this->_cachedHeights.raiseErrorOnBadCoord = false;
    this->_cachedHeights.defaultValueOnBadCoord = 0.f;
}

float Patch3D::getMaxHeight(Vector3 position)
{
    if (this->_cachedHeights.at(position.xy()) < -1.f) {
        float resolution = 0.1;
        float maxHeight = this->dimension.z;
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
    if (boundMin.x > p.x || p.x > dimension.x || boundMin.y > p.y || p.y > dimension.y || boundMin.z > p.z || p.z > dimension.z)
        return 0.f;
    return (this->evalFunction(p) > 0.f ? 1.f : 0.f);
}

Patch3D Patch3D::stack(Patch3D *P1, Patch3D *P2)
{
    Vector3 P1center = P1->position;
    Vector3 P2center = P2->position;

    Vector3 boundMin = Vector3::min(P1->boundMin + P1center, P2->boundMin + P2center);
    Vector3 boundMax = Vector3::max(P1->dimension + P1center, P2->dimension + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->dimension.z - P1->boundMin.z) + (P2->dimension.z - P2->boundMin.z);

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
    Vector3 boundMax = Vector3::max(P1->dimension + P1center, P2->dimension + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->dimension.z - P1->boundMin.z) + (P2->dimension.z - P2->boundMin.z);

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
    Vector3 boundMax = Vector3::max(P1->dimension + P1center, P2->dimension + P2center);
    // Stacking means that the height can double (at most)
    boundMax.z = boundMin.z + (P1->dimension.z - P1->boundMin.z) + (P2->dimension.z - P2->boundMin.z);

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
*/
