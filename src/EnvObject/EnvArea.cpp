#include "EnvArea.h"

EnvArea::EnvArea()
    : EnvObject()
{

}

EnvArea *EnvArea::fromJSON(nlohmann::json content)
{
    return dynamic_cast<EnvArea*>(EnvObject::fromJSON(content));
}
float EnvArea::getSqrDistance(const Vector3 &position)
{
    return (position - this->curve.estimateClosestPos(position)).norm2() * (this->curve.containsXY(position, false) ? -1.f : 1.f);
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

EnvArea *EnvArea::instantiate(std::string objectName)
{
    return dynamic_cast<EnvArea*>(EnvObject::instantiate(objectName));
}

void EnvArea::recomputeEvaluationPoints()
{
    if (evaluateInside) {
        this->evaluationPositions = curve.randomPointsInside(curve.size());
    } else {
        this->evaluationPositions = curve.points;
    }
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

void EnvArea::applyDepositionOnDeath()
{
    for (auto& [materialName, amount] : materialDepositionOnDeath) {
        auto& material = EnvObject::materials[materialName];
        if (amount == 0) return;
        AABBox box = AABBox(this->curve.points);
        ShapeCurve translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p = p + Vector3(width, width, 0) - box.min();
        GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

        sand.iterateParallel([&] (const Vector3& pos) {
            bool inside = translatedCurve.contains(pos);
            sand(pos) = (inside ? 1.f : 0.f) * amount;
        });
        material.currentState.add(sand.meanSmooth(width, width, 1), box.min() - Vector3(width, width));
    }
}

std::pair<GridV3, GridF> EnvArea::computeFlowModification()
{
    if (flowEffect == Vector3()) return {EnvObject::flowfield, GridF()};

    float growingState = computeGrowingState2();

    if (_cachedFlowModif.empty()) {
        Vector3 objectWidth = Vector3(width, width, 0);
        Vector3 halfWidth = objectWidth * .5f;
        ShapeCurve translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p.z = 0;
        AABBox box = AABBox(translatedCurve.points);
        box.expand({box.min() - halfWidth, box.max() + halfWidth});


        GridV3 flow = GridV3(EnvObject::flowfield.getDimensions());

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
                flow.at(pos) += impact.changedBasis(direction, normal, binormal);// + impact * previousFlow;
                // occupancy.at(pos) = 1.f;
            });
        });
        // return {flow, occupancy};
        this->_cachedFlowModif = flow;
    }
    return {EnvObject::flowfield.add(_cachedFlowModif * growingState, Vector3()), GridF()};
}

ImplicitPatch* EnvArea::createImplicitPatch(const GridF &heights, ImplicitPrimitive* previousPrimitive)
{
    if (!geometryNeedsUpdate) return this->_patch;

    if (this->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    BSpline translatedCurve = this->curve;
    for (Vector3& p : translatedCurve) {
        p.z = heights(p.xy());
    }
    AABBox box(translatedCurve.points);
    float growingState = 1.f; //this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    Vector3 offset(this->width, this->width, this->height * growingState);
    translatedCurve.translate(-(box.min() - offset * .5f));
    box.expand(box.max() + Vector3(0, 0, this->height * growingState));
    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        translatedCurve = previousPrimitive->optionalCurve;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->height * growingState, translatedCurve, true);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->height * growingState, translatedCurve, true);
    }
    patch->position = (box.min() - offset.xy() * .5f).xy();
    // patch->position.z += heights(patch->position.xy());
    patch->supportDimensions = box.dimensions() + offset;
    patch->material = this->material;
    patch->name = this->name;
    this->_patch = patch;
    this->geometryNeedsUpdate = false;
    return patch;
}
/*
GridF EnvArea::createHeightfield()
{
    if (_cachedHeightfield.empty()) {
        auto patch = dynamic_cast<ImplicitPrimitive*>(_patch);
        GridF heights(patch->getDimensions().x, patch->getDimensions().y, 1);
        heights.iterateParallel([&] (const Vector3& p) {
            for (int i = 1; i < height; i++) {
                heights(p) += (patch->evaluate(p.xy() + patch->position.xy() + Vector3(0, 0, i)) > 0 ? 1.f : 0.f);
            }
        });
        _cachedHeightfield = heights;
    }
    return _cachedHeightfield;
}
*/
EnvArea &EnvArea::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    for (auto& p : evaluationPositions)
        p.translate(translation);
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
    this->geometryNeedsUpdate = true;
    return *this;
}

void EnvArea::updateCurve(const BSpline &newCurve)
{
    if (evaluateInside) {
        std::cerr << "Need to implement 'updateCurve' for inside...!" << std::endl;
        for (auto& evaluationPosition : this->evaluationPositions) {
            float evaluationPointClosestTime = this->curve.estimateClosestTime(evaluationPosition);
            Vector3 relativeDisplacementFromEvaluationToCurve = (evaluationPosition - this->curve.getPoint(evaluationPointClosestTime));
            evaluationPosition = newCurve.getPoint(evaluationPointClosestTime) + relativeDisplacementFromEvaluationToCurve;
        }
    } else {
        this->evaluationPositions = newCurve.points;
    }
    this->curve = newCurve;
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
}

nlohmann::json EnvArea::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["curve"] = bspline_to_json(this->curve);
    return json;
}
