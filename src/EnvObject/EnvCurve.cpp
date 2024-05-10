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

void EnvCurve::applyDeposition(EnvMaterial& material)
{
    if (this->materialDepositionRate[material.name] == 0) return;

    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = BSpline(this->curve.getPath(100));
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();
    GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
        float amount = normalizedGaussian(width * .25f, distToCurve) * this->materialDepositionRate[material.name];

        sand.at(pos) = amount;
    });
    material.currentState.add(sand, box.min().xy() - Vector3(width, width));
}

void EnvCurve::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate[material.name] == 0) return;
    displayProcessTime("Absorption... ", [&]() {
        AABBox box;
        BSpline translatedCurve;
        GridF sand;

        box = AABBox(this->curve.points);
        translatedCurve = BSpline(this->curve.getPath(100));
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();
        sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

        sand.iterateParallel([&] (const Vector3& pos) {
            float distToCurve = translatedCurve.estimateSqrDistanceFrom(pos, true);
            float amount = normalizedGaussian(width * .25f, distToCurve) * this->materialAbsorptionRate[material.name];
            sand.at(pos) = amount;
        });

        material.currentState.add(-sand, box.min().xy() - Vector3(width, width));

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
//    return {GridV3(EnvObject::flowfield.getDimensions()), GridF(EnvObject::flowfield.getDimensions())};
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
    return {flow, GridF(flow.getDimensions(), 1.f)};

    /*
    GridV3 flow(EnvObject::flowfield.getDimensions());
    GridF occupancy(flow.getDimensions());

    flow.iterateParallel([&] (const Vector3& pos) {
        auto initialFlow = EnvObject::flowfield(pos).normalized();
        float closestTime = translatedCurve.estimateClosestTime(pos);
        Vector3 closestPos = translatedCurve.getPoint(closestTime);
        float sqrDist = (closestPos - pos).norm2();
        if (sqrDist > (width * .5f) * (width * .5f)) return;
        float gauss = 1.f; // normalizedGaussian(width * .1f, sqrDist);
        Vector3 impact = this->flowEffect * gauss;
        auto [direction, normal, binormal] = translatedCurve.getFrenetFrame(closestTime);
        direction *= direction.dot(initialFlow);
        normal *= -normal.dot(initialFlow);
        binormal *= binormal.dot(initialFlow);
        flow(pos) = impact.changedBasis(direction, normal, binormal);
        occupancy(pos) = 1.f;
    });
    return {flow, occupancy};*/
}

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
        BSpline translatedCurve = previousPrimitive->optionalCurve;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, false);
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
        }
        translatedCurve.translate(-(box.min() - offset * .5f));
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, false);
    }
    patch->position = (box.min() - offset.xy() * .5f).xy();
    patch->supportDimensions = box.dimensions() + offset;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

GridF EnvCurve::createHeightfield() const
{
    return GridF();
}

EnvCurve &EnvCurve::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    this->evaluationPosition.translate(translation);
    return *this;
}


void EnvCurve::updateCurve(const BSpline &newCurve)
{
    float evaluationPointClosestTime = this->curve.estimateClosestTime(this->evaluationPosition);
    Vector3 relativeDisplacementFromEvaluationToCurve = (this->evaluationPosition - this->curve.getPoint(evaluationPointClosestTime));
    this->curve = newCurve;
    this->evaluationPosition = this->curve.getPoint(evaluationPointClosestTime) + relativeDisplacementFromEvaluationToCurve;
}

nlohmann::json EnvCurve::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["curve"] = bspline_to_json(this->curve);
    return json;
}
