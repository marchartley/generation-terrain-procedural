#include "EnvArea.h"

EnvArea::EnvArea()
    : EnvObject()
{

}

EnvObject *EnvArea::fromJSON(nlohmann::json content)
{
    return EnvObject::fromJSON(content);
}

float EnvArea::getSqrDistance(const Vector3 &position, std::string complement)
{
    if (complement == "center")
        return (position - this->curve.center()).norm2();
    else if (complement == "border") // Just keep the absolute value
        return this->curve.estimateSqrDistanceFrom(position);
    else if (complement == "start") // Yeah, that doesn't make sense...
        return 0.f;
    else if (complement == "end") // Yeah, that doesn't make sense...
        return 0.f;

    float dist = this->curve.estimateSignedDistanceFrom(position);
    return (dist * dist) * (dist > 0 ? 1.f : -1.f); // Keep the sign
}

Vector3 EnvArea::getVector(const Vector3 &position, std::string complement)
{
    if (complement == "direction")
        return this->curve.getDirection(curve.estimateClosestTime(position));
    else if (complement == "normal")
        return this->curve.getNormal(curve.estimateClosestTime(position));
    else if (complement == "binormal")
        return this->curve.getBinormal(curve.estimateClosestTime(position));
    return Vector3();
}

Vector3 EnvArea::getNormal(const Vector3 &position)
{
    return this->curve.getNormal(this->curve.estimateClosestTime(position));
}

Vector3 EnvArea::getDirection(const Vector3 &position)
{
    return Vector3::invalid();
}

Vector3 EnvArea::getProperty(const Vector3& position, std::string prop) const
{
    if (prop == "center") {
        return this->curve.center();
    } else if (prop == "start") {
        return Vector3::invalid();
    } else if (prop == "end") {
        return Vector3::invalid();
    } else if (prop == "inside") {
        return (this->curve.containsXY(position, false) ? Vector3(true) : Vector3(false));
    }
    return this->curve.estimateClosestPos(position); // Default
}

std::map<std::string, Vector3> EnvArea::getAllProperties(const Vector3 &position) const
{
    float closestTime = this->curve.estimateClosestTime(position);
    Vector3 closestPos = this->curve.getPoint(closestTime);
    return {
        {"default", closestPos},
        {"center", this->curve.center()},
        {"start", Vector3::invalid()},
        {"end", Vector3::invalid()},
        {"inside", (this->curve.containsXY(position, false) ? Vector3(true) : Vector3(false))},
        {"normal", this->curve.getNormal(closestTime)},
        {"dir", Vector3::invalid()},
        {"curvature", Vector3(this->curve.getCurvature(closestTime), 0, 0)}
    };
}

EnvArea *EnvArea::clone()
{
    EnvArea* self = new EnvArea;
    *self = *this;
    return self;
}

void EnvArea::applyDeposition(EnvMaterial& material)
{
    if (this->materialDepositionRate[material.name] == 0) return;
    AABBox box = AABBox(this->curve.points);
    ShapeCurve translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();
    GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        bool inside = translatedCurve.contains(pos);
        sand(pos) = (inside ? 1.f : 0.f) * this->materialDepositionRate[material.name]; //gaussian(width, translatedCurve.estimateSqrDistanceFrom(Vector3(x, y, 0)));
    });
    material.currentState.add(sand.meanSmooth(width, width, 1), box.min() - Vector3(width, width));
}

void EnvArea::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate[material.name] == 0) return;
    return;
}

std::pair<GridV3, GridF> EnvArea::computeFlowModification()
{
    Vector3 objectWidth = Vector3(width, width, 0);
    Vector3 halfWidth = objectWidth * .5f;
    ShapeCurve translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p.z = 0;
    AABBox box = AABBox(translatedCurve.points);
    box.expand({box.min() - halfWidth, box.max() + halfWidth});


    GridV3 flow = EnvObject::flowfield; //(EnvObject::flowfield.getDimensions());
    GridF occupancy(EnvObject::flowfield.getDimensions());

    GridF dist(flow.getDimensions());
    flow.iterateParallel([&] (const Vector3& pos) {
        dist(pos) = (box.contains(pos) && translatedCurve.contains(pos, true) ? 1.f : 0.f);
    });
    dist = dist.toDistanceMap(true, false);
    GridV3 grad = dist.gradient() * -1.f;
    for (auto& v : grad)
        v.normalize();

    float timePrepare = 0.f;
    float timeApply = 0.f;

    flow.iterateParallel([&] (const Vector3& pos) {
        if (!box.contains(pos) || !translatedCurve.contains(pos, true))
            return;

        Vector3 impact, direction, normal, binormal;

        timePrepare += timeIt([&]() {
            float closestTime = translatedCurve.estimateClosestTime(pos);
            Vector3 closestPos = translatedCurve.getPoint(closestTime);

            float distanceToBorder = (pos - closestPos).norm();
            float distFactor = clamp(distanceToBorder / (width * .5f), 0.f, 1.f); // On border = 1, at w/2 = 0, more inside = 0
            Vector3 previousFlow = EnvObject::flowfield(pos);
            // Change the order of the Frenet Frame to get the direction in the direction of the "outside" and the normal is along the borders
    //            auto [normal, direction, binormal] = translatedCurve.getFrenetFrame(closestTime);
            // We will use the distance map to get the direction, then we know (0, 0, 1) is the binormal (2D shape), so normal is cross product.
            direction = grad(pos);
            binormal = Vector3(0, 0, 1);
            normal = direction.cross(binormal); // Direction and binormal are normalized.
            if (normal.dot(previousFlow) > 0) {
                normal *= -1.f;
            }
            impact = this->flowEffect;// * distFactor;
        });
        timeApply += timeIt([&]() {
            flow.at(pos) = impact.changedBasis(direction, normal, binormal);// + impact * previousFlow;
            occupancy.at(pos) = 1.f;
        });
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvArea::createImplicitPatch(const GridF &heights, ImplicitPrimitive* previousPrimitive)
{
    if (this->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    BSpline translatedCurve = this->curve;
    for (Vector3& p : translatedCurve) {
        p.z = heights(p.xy());
    }
    AABBox box(translatedCurve.points);
    Vector3 offset(this->width, this->width, this->height * this->computeGrowingState());
    translatedCurve.translate(-(box.min() - offset * .5f));
    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        translatedCurve = previousPrimitive->optionalCurve;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->height * this->computeGrowingState(), translatedCurve, false);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->height * this->computeGrowingState(), translatedCurve, false);
    }
    patch->position = (box.min() - offset.xy() * .5f).xy();
    patch->position.z += heights(patch->position.xy());
    patch->supportDimensions = box.dimensions() + offset;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

GridF EnvArea::createHeightfield() const
{
    return GridF();
}

EnvObject &EnvArea::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    this->evaluationPosition.translate(translation);
    return *this;
}

void EnvArea::updateCurve(const BSpline &newCurve)
{
    float evaluationPointClosestTime = this->curve.estimateClosestTime(this->evaluationPosition);
    Vector3 relativeDisplacementFromEvaluationToCurve = (this->evaluationPosition - this->curve.getPoint(evaluationPointClosestTime));
    this->curve = newCurve;
    this->evaluationPosition = this->curve.getPoint(evaluationPointClosestTime) + relativeDisplacementFromEvaluationToCurve;
}
