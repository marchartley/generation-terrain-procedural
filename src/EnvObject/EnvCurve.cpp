#include "EnvCurve.h"

EnvCurve::EnvCurve()
    : EnvObject()
{

}

EnvCurve *EnvCurve::fromJSON(nlohmann::json content)
{
    return dynamic_cast<EnvCurve*>(EnvObject::fromJSON(content));
}
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

EnvCurve *EnvCurve::instantiate(std::string objectName)
{
    return dynamic_cast<EnvCurve*>(EnvObject::instantiate(objectName));
}

void EnvCurve::recomputeEvaluationPoints()
{
    this->evaluationPositions = curve.points;
}

void EnvCurve::applyDeposition(EnvMaterial& material)
{
    if (this->materialDepositionRate[material.name] == 0) return;

    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = BSpline(this->curve.getPath(100));
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();

    if (_cachedAbsorptionDepositionField.empty()) {
        _cachedAbsorptionDepositionField = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

        _cachedAbsorptionDepositionField.iterateParallel([&] (const Vector3& pos) {
            float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
            float amount = normalizedGaussian(width * .25f, distToCurve);

            _cachedAbsorptionDepositionField.at(pos) = amount;
        });
    }
    material.currentState.add(_cachedAbsorptionDepositionField * this->materialDepositionRate[material.name], box.min().xy() - Vector3(width, width));
}

void EnvCurve::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate[material.name] == 0) return;
    displayProcessTime("Absorption... ", [&]() {
        AABBox box = AABBox(this->curve.points);
        BSpline translatedCurve = BSpline(this->curve.getPath(100));
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();

        if (_cachedAbsorptionDepositionField.empty()) {
            _cachedAbsorptionDepositionField = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

            _cachedAbsorptionDepositionField.iterateParallel([&] (const Vector3& pos) {
                float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
                float amount = normalizedGaussian(width * .25f, distToCurve);

                _cachedAbsorptionDepositionField.at(pos) = amount;
            });
        }
        material.currentState.add(_cachedAbsorptionDepositionField * this->materialAbsorptionRate[material.name], box.min().xy() - Vector3(width, width));

        material.currentState.iterateParallel([&] (size_t i) {
            material.currentState[i] = std::max(material.currentState[i], 0.f);
        });
    }, false);
}

void EnvCurve::applyDepositionOnDeath()
{
    for (auto& [materialName, amount] : materialDepositionOnDeath) {
        auto& material = EnvObject::materials[materialName];
        if (amount == 0) return;

        AABBox box = AABBox(this->curve.points);
        BSpline translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();
        GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

        sand.iterateParallel([&] (const Vector3& pos) {
            sand.at(pos) = normalizedGaussian(width * .25f, translatedCurve.estimateSqrDistanceFrom(pos)) * amount;
        });
        material.currentState.add(sand, box.min().xy() - Vector3(width, width));
    }
}

std::pair<GridV3, GridF> EnvCurve::computeFlowModification()
{
    if (this->flowEffect == Vector3()) return {EnvObject::flowfield, GridF()};

    float growingState = computeGrowingState2();

    if (this->_cachedFlowModif.empty()) {
        Vector3 objectWidth = Vector3(width, width, 0);
        Vector3 halfWidth = objectWidth * .5f;
        BSpline translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p.z = 0;
        AABBox box = AABBox(translatedCurve.points);
        box.expand({box.min() - halfWidth, box.max() + halfWidth});


        /*
        PinchKelvinletCurve k;
        k.radialScale = width * .05f;
        k.force = 10.f;
        */
        /*
        ScaleKelvinletCurve k2;
        k2.radialScale = width * .1f;
        k2.force = -10.f;
        k2.mu = .9f;
        k2.v = 0;
        k2.curve = translatedCurve;
        */
        TranslateKelvinletCurve k;
        k.radialScale = width * .05f;
        k.force = 10.f;

        /*
        TwistKelvinletCurve k;
        k.radialScale = width * .05f;
        k.force = 10.f;
        */

        k.curve = translatedCurve;

        GridV3 flow = EnvObject::flowfield;
        flow.iterateParallel([&](const Vector3& p) {
            Vector3 displacement = k.evaluate(p);
            flow(p) += displacement.xy();
        });
        this->_cachedFlowModif = flow;
    }
    return {EnvObject::flowfield.add(_cachedFlowModif * growingState, Vector3()), GridF()}; //{flow, GridF(flow.getDimensions(), 1.f)};
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
            p.z = heights(p.xy());
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
        // translatedCurve.translate(Vector3(0, 0, -offset.z * 0.5f));
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
    float height = this->height;
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
            p.z = heights(p.xy());
            maxHeight = std::max(maxHeight, p.z);
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
            p.z = heights(p.xy());
            maxHeight = std::max(maxHeight, p.z);
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
    this->curve = newCurve;
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
}

nlohmann::json EnvCurve::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["curve"] = bspline_to_json(this->curve);
    return json;
}
