#include "EnvCurve.h"
#include "EnvObject/EnvironmentalScene.h"

#include "serialization/Serializer.h"

EnvCurve::EnvCurve()
    : EnvObject()
{}

float EnvCurve::getSqrDistance(const Vector3 &position)
{
    return (position - this->curve.estimateClosestPos(position)).norm2();
}

std::map<std::string, Vector3> EnvCurve::getAllProperties(const Vector3 &position) const
{
    float closestTime = this->curve.estimateClosestTime(position);
    Vector3 closestPos = this->curve.getPoint(closestTime);
    return {
        {"default", closestPos},
        {"center", this->curve.center()},
        {"start", this->curve.points.front()},
        {"end", this->curve.points.back()},
        {"inside", ((position - closestPos).norm2() < width * width ? Vector3(true) : Vector3(false))},
        {"normal", this->curve.getNormal(closestTime)},
        {"dir", this->curve.getDirection(closestTime)},
        {"curvature", Vector3(this->curve.getCurvature(closestTime), 0, 0)}
    };
}

EnvCurve *EnvCurve::clone()
{
    EnvCurve* self = new EnvCurve;
    *self = *this;
    return self;
}

bool EnvCurve::placeInTerrain(const Vector3 &seedPosition)
{
    BSpline initialCurve;
    if (this->curveFollow == EnvCurve::SKELETON) {
        initialCurve = ContinuousCurveOptimizer::getSkeletonCurve(seedPosition, this->fitnessFunction, this->length);
    } else if (this->curveFollow == EnvCurve::ISOVALUE) {
        initialCurve = ContinuousCurveOptimizer::followIsolevel(seedPosition, this->fitnessFunction, this->length);
    } else if (this->curveFollow == EnvCurve::GRADIENTS) {
        initialCurve = ContinuousCurveOptimizer::getExactLengthCurveFollowingGradients(seedPosition, this->fitnessFunction, this->length);
    }
    this->snake.position = seedPosition;
    return this->placeInTerrain(initialCurve);
}

bool EnvCurve::placeInTerrain(const BSpline &seedCurve)
{
    if (seedCurve.empty()) {
        return false;
    }
    this->curve = seedCurve;
    this->curve.resamplePoints(10);
    Vector3 position = this->curve[this->curve.size() / 2];
    this->curve.translate(-position);
    this->translate(position.xy());
    this->recomputeEvaluationPoints();
    this->fitnessScoreAtCreation = this->evaluate();
    if (this->fitnessScoreAtCreation < this->minScore)
        return false;
    this->spawnTime = this->scene->currentTime;
    return true;
}

void EnvCurve::improvePositionning(float steps)
{
    this->snake.contour = this->curve;
    this->snake.position = this->curve.getPoint(.5f);
    this->updateCurve(this->snake.runSegmentation(int(steps)));
}

void EnvCurve::recomputeEvaluationPoints()
{
    this->evaluationPositions = curve.points;
}

void EnvCurve::applyDeposition(EnvMaterial& material)
{
    if (this->materialDepositionRate.count(material.name) == 0) return;
    auto depositionProperties = this->materialDepositionRate[material.name];
    if (depositionProperties.rate == 0 || depositionProperties.radius == 0) return;

    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = this->curve; //.getPath(100);
    translatedCurve.translate(Vector3(depositionProperties.radius, depositionProperties.radius, 0) - box.min());

    GridF deposition = GridF(box.dimensions().x() + depositionProperties.radius * 2.f, box.dimensions().y() + depositionProperties.radius * 2.f);

    deposition.iterateParallel([&] (const Vector3& pos) {
        float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
        // float amount = normalizedGaussian(depositionProperties.radius * .25f, distToCurve);
        float amount = (distToCurve < depositionProperties.radius * depositionProperties.radius ? 1.f : 0.f);

        deposition.at(pos) = amount;
    });
    material.currentState.add(deposition * depositionProperties.rate, box.min().xy() - Vector3(depositionProperties.radius, depositionProperties.radius));

    /*if (_cachedAbsorptionDepositionField.empty()) {
        _cachedAbsorptionDepositionField = GridF(box.dimensions().x() + depositionProperties.radius * 2.f, box.dimensions().y() + depositionProperties.radius * 2.f);

        _cachedAbsorptionDepositionField.iterateParallel([&] (const Vector3& pos) {
            float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
            // float amount = normalizedGaussian(depositionProperties.radius * .25f, distToCurve);
            float amount = (distToCurve < depositionProperties.radius ? 1.f : 0.f);

            _cachedAbsorptionDepositionField.at(pos) = amount;
        });
    }
    material.currentState.add(_cachedAbsorptionDepositionField * this->materialDepositionRate[material.name].rate, box.min().xy() - Vector3(depositionProperties.radius, depositionProperties.radius));
    */
}

void EnvCurve::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate.count(material.name) == 0) return;
    auto absorptionProperties = this->materialAbsorptionRate[material.name];
    if (absorptionProperties.rate == 0 || absorptionProperties.radius == 0) return;

    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = this->curve; //.getPath(100);
    translatedCurve.translate(Vector3(absorptionProperties.radius, absorptionProperties.radius, 0) - box.min());

    GridF absorption = GridF(box.dimensions().x() + absorptionProperties.radius * 2.f, box.dimensions().y() + absorptionProperties.radius * 2.f);

    absorption.iterateParallel([&] (const Vector3& pos) {
        float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
        // float amount = normalizedGaussian(absorptionProperties.radius * .25f, distToCurve);
        float amount = (distToCurve < absorptionProperties.radius * absorptionProperties.radius ? 1.f : 0.f);

        absorption.at(pos) = amount;
    });
    material.currentState.add(absorption * absorptionProperties.rate * -1.f, box.min().xy() - Vector3(absorptionProperties.radius, absorptionProperties.radius));
    material.currentState.iterateParallel([&] (size_t i) { material.currentState[i] = std::max(material.currentState[i], 0.f); });
    /*
    if (this->materialAbsorptionRate[material.name].rate == 0) return;
    displayProcessTime("Absorption... ", [&]() {
        AABBox box = AABBox(this->curve.points);
        BSpline translatedCurve = BSpline(this->curve.getPath(100));
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();

        if (_cachedAbsorptionDepositionField.empty()) {
            _cachedAbsorptionDepositionField = GridF(box.dimensions().x() + width * 2.f, box.dimensions().y() + width * 2.f);

            _cachedAbsorptionDepositionField.iterateParallel([&] (const Vector3& pos) {
                float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
                float amount = normalizedGaussian(width * .25f, distToCurve);

                _cachedAbsorptionDepositionField.at(pos) = amount;
            });
        }
        material.currentState.add(_cachedAbsorptionDepositionField * this->materialAbsorptionRate[material.name].rate, box.min().xy() - Vector3(width, width));

        material.currentState.iterateParallel([&] (size_t i) {
            material.currentState[i] = std::max(material.currentState[i], 0.f);
        });
    }, false);
    */
}

void EnvCurve::applyDepositionOnDeath()
{
    for (auto& [materialName, depos] : materialDepositionOnDeath) {
        auto& material = this->scene->materials[materialName];
        if (depos.rate == 0) return;

        AABBox box = AABBox(this->curve.points);
        BSpline translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();
        GridF sand = GridF(box.dimensions().x() + width * 2.f, box.dimensions().y() + width * 2.f);

        sand.iterateParallel([&] (const Vector3& pos) {
            sand.at(pos) = normalizedGaussian(width * .25f, translatedCurve.estimateSqrDistanceFrom(pos)) * depos.rate;
        });
        material.currentState.add(sand, box.min().xy() - Vector3(width, width));
    }
}

GridV3& EnvCurve::computeFlowModification(GridV3& waterFlow)
{
    std::vector<RelativeKelvinlet> relativeFlowsStarting;
    std::vector<RelativeKelvinlet> relativeFlowsEnding;
    std::vector<RelativeKelvinlet> relativeCurveFlow;
    for (size_t i = 0; i < startingPointKelvinlets.size(); i++) {
        if (startingPointKelvinlets[i]->valid())
            relativeFlowsStarting.push_back(RelativeKelvinlet(startingPointKelvinlets[i], this->curve.points.front()));
    }
    for (size_t i = 0; i < endingPointKelvinlets.size(); i++) {
        if (endingPointKelvinlets[i]->valid())
            relativeFlowsEnding.push_back(RelativeKelvinlet(endingPointKelvinlets[i], this->curve.points.back()));
    }
    for (size_t i = 0; i < curveKelvinlets.size(); i++) {
        if (curveKelvinlets[i]->valid())
            relativeCurveFlow.push_back(RelativeKelvinlet(curveKelvinlets[i], Vector3()));
    }

    const Vector3 initialFlowStarting = waterFlow.interpolate(this->curve.points.front());
    float flowAngleStarting = curve.getDerivative(0).getSignedAngleWith(Vector3(1, 0, 0)); //initialFlowStarting.getSignedAngleWith(Vector3(1, 0, 0));
    float flowStrengthStarting = initialFlowStarting.length();

    const Vector3 initialFlowEnding = waterFlow.interpolate(this->curve.points.back());
    float flowAngleEnding = curve.getDerivative(1).getSignedAngleWith(Vector3(1, 0, 0)); // initialFlowEnding.getSignedAngleWith(Vector3(1, 0, 0));
    float flowStrengthEnding = initialFlowEnding.length();

    waterFlow.iterateParallel([&](const Vector3& p) {
        for (const auto& relativeK : relativeFlowsStarting) {
            waterFlow[p] += relativeK.evaluate(p, flowAngleStarting, flowStrengthStarting, true);
        }
        for (const auto& relativeK : relativeFlowsEnding) {
            waterFlow[p] += relativeK.evaluate(p, flowAngleEnding, flowStrengthEnding, true);
        }
        for (const auto& k : relativeCurveFlow) {
            waterFlow[p] += k.evaluate(p, 0.f, waterFlow[p].norm());
        }
    });
    /*
    if (this->flowEffect == Vector3()) return {waterFlow, GridF()};

    float growingState = computeGrowingState2();

    if (this->_cachedFlowModif.empty()) {
        Vector3 objectWidth = Vector3(width, width, 0);
        Vector3 halfWidth = objectWidth * .5f;
        BSpline translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p.z() = 0;
        AABBox box = AABBox(translatedCurve.points);
        box.expand({box.min() - halfWidth, box.max() + halfWidth});

        GrabKelvinletCurve k;
        k.radialScale = width * .05f;
        k.force = 10.f;

        k.curve = translatedCurve;

        GridV3 flow = waterFlow;
        flow.iterateParallel([&](const Vector3i& p) {
            Vector3 displacement = k.evaluate(p);
            flow(p) += displacement.xy();
        });
        this->_cachedFlowModif = flow;
    }
    return {waterFlow.add(_cachedFlowModif * growingState, Vector3()), GridF()}; //{flow, GridF(flow.getDimensions(), 1.f)};
    */
    return waterFlow;
}


/*
ImplicitPatch* EnvCurve::createImplicitPatch(const GridF& _heights, ImplicitPrimitive *previousPrimitive)
{
    if (this->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    AABBox box(this->curve.points);
    float growingState = this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    float height = this->height * growingState;
    Vector3 offset(this->width, this->width, height * .5f);

    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        // BSpline translatedCurve = previousPrimitive->optionalCurve;
        // box = AABBox(translatedCurve.points);
        // box.expand({box.min() - offset, box.max() + offset});
        // std::cout << "Not nullptr box: " << box << std::endl;
        // *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, false);
    } else {
        BSpline translatedCurve = this->curve;
        GridF heights = _heights;
        heights.raiseErrorOnBadCoord = false;
        heights.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        for (Vector3& p : translatedCurve) {
            p.z() = heights(p.xy());
        }
        if (height == 0){
            box = AABBox({box.center()});
            offset = Vector3();
        } else {
            box = AABBox(translatedCurve.points);
            box.expand({box.min() - offset, box.max() + offset});
        }
        std::cout << "Nullptr box: " << box << std::endl;
        // translatedCurve.translate(-(box.min() - offset * .5f));
        // translatedCurve.translate(-(box.min()).xy());
        // translatedCurve.translate(Vector3(0, 0, -offset.z() * 0.5f));
        // box = AABBox(translatedCurve.points);
        // box.expand({box.min() - offset, box.max() + offset});
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, false);
    }
    patch->position = (box.min() - offset.xy() * .5f).xy();
    patch->supportDimensions = box.dimensions() + offset;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}*/
ImplicitPatch* EnvCurve::createImplicitPatch(const GridF& _heights, ImplicitPrimitive *previousPrimitive)
{
    if (!geometryNeedsUpdate) return this->_patch;
    if (this->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    AABBox box(this->curve.points);
    float growingState = 1.f; // this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    // float height = this->height;
    float radius = this->width * growingState;
    Vector3 offset(this->width, this->width, radius * .5f);

    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        /*BSpline translatedCurve = previousPrimitive->optionalCurve;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, false);*/

        BSpline translatedCurve = this->curve;
        GridF heights = _heights;
        heights.raiseErrorOnBadCoord = false;
        heights.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        float maxHeight = 0;
        for (Vector3& p : translatedCurve) {
            p.z() = heights(p.xy());
            maxHeight = std::max(maxHeight, p.z());
        }
        if (radius == 0){
            box = AABBox({box.center()});
            offset = Vector3();
        }
        box = AABBox(translatedCurve.points);
        box.expand({box.min(), box.max() + offset * 1.f + Vector3(0, 0, maxHeight + 10)});
        translatedCurve.translate(-(box.min() - offset * .5f));
        *patch = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, radius, translatedCurve, false);
        patch->position = (box.min() - offset.xy() * .5f).xy();
        patch->supportDimensions = box.dimensions() + offset;

    } else {
        BSpline translatedCurve = this->curve;
        GridF heights = _heights;
        heights.raiseErrorOnBadCoord = false;
        heights.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::MIRROR_VALUE;
        float maxHeight = 0;
        for (Vector3& p : translatedCurve) {
            p.z() = heights(p.xy());
            maxHeight = std::max(maxHeight, p.z());
        }
        if (radius == 0){
            box = AABBox({box.center()});
            offset = Vector3();
        }
        box = AABBox(translatedCurve.points);
        box.expand({box.min(), box.max() + offset * 1.f + Vector3(0, 0, maxHeight + 10)});
        translatedCurve.translate(-(box.min() - offset * .5f));
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, radius, translatedCurve, false);
        patch->position = (box.min() - offset.xy() * .5f).xy();
        patch->supportDimensions = box.dimensions() + offset;
    }
    patch->material = this->material;
    patch->name = this->name;
    this->_patch = patch;
    this->geometryNeedsUpdate = false;
    return patch;
}

/*GridF EnvCurve::createHeightfield()
{
    return GridF();
}*/

EnvCurve &EnvCurve::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    for (auto& evaluationPosition : evaluationPositions)
        evaluationPosition.translate(translation);
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
    this->geometryNeedsUpdate = true;
    return *this;
}


void EnvCurve::updateCurve(const BSpline &newCurve)
{
    /*
    float evaluationPointClosestTime = this->curve.estimateClosestTime(this->evaluationPosition);
    Vector3 relativeDisplacementFromEvaluationToCurve = (this->evaluationPosition - this->curve.getPoint(evaluationPointClosestTime));
    this->evaluationPosition = newCurve.getPoint(evaluationPointClosestTime) + relativeDisplacementFromEvaluationToCurve;*/
    this->evaluationPositions = newCurve.points;
    for (auto& k : this->curveKelvinlets) {
        if (auto asKelvinletCurve = dynamic_cast<KelvinletCurve*>(k)) {
            asKelvinletCurve->curve = newCurve;
        }
    }
    this->curve = newCurve;
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
    this->_cachedAbsorptionDepositionField.clear();
}

/*
nlohmann::json EnvCurve::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["curve"] = this->curve;
    return json;
}
*/
