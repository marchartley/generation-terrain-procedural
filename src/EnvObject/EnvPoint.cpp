#include "EnvPoint.h"

EnvPoint::EnvPoint()
    : EnvObject()
{

}

EnvPoint *EnvPoint::fromJSON(nlohmann::json content)
{
    return dynamic_cast<EnvPoint*>(EnvObject::fromJSON(content));
}

float EnvPoint::getSqrDistance(const Vector3 &position)
{
    return (position - this->position).norm2();
}

std::map<std::string, Vector3> EnvPoint::getAllProperties(const Vector3 &position) const
{
    Vector3 diff = (this->position - position);
    return {
        {"default", this->position},
        {"center", this->position},
        {"start", this->position},
        {"end", this->position},
        {"inside", (diff.norm2() < this->radius * this->radius ? Vector3(true) : Vector3(false))},
        {"normal", diff.normalized()},
        {"dir", Vector3::invalid()},
        {"curvature", Vector3(this->radius, 0, 0)}
    };
}

EnvPoint *EnvPoint::clone()
{
    EnvPoint* self = new EnvPoint;
    *self = *this;
    return self;
}

EnvPoint *EnvPoint::instantiate(std::string objectName)
{
    return dynamic_cast<EnvPoint*>(EnvObject::instantiate(objectName));
}

void EnvPoint::recomputeEvaluationPoints()
{
    this->evaluationPositions = {position};
}

void EnvPoint::applyDeposition(EnvMaterial& material)
{
    // if (toLower(this->name) == "coralpolyp") std::cout << "Deposing " << material.name << "? " << (this->materialDepositionRate[material.name] == 0 ? "No" : "Yes") << std::endl;
    if (this->materialDepositionRate[material.name] == 0) return;
    float growingState = this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    GridF sand = GridF::normalizedGaussian(radius, radius, 1, radius * .25f) * this->materialDepositionRate[material.name] * growingState /** (EnvObject::flowfield(this->position).norm() * 10.f)*/;
    material.currentState.add(sand, this->position - sand.getDimensions() * .5f);
}

void EnvPoint::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate[material.name] == 0) return;
    if (this->needsForGrowth.count(material.name) == 0) return;

    float growingState = this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    float sandMissing = (needsForGrowth[material.name] - currentSatisfaction[material.name]) * (growingState + .01f) / (1.01f);

    GridF sandAbsorbed = GridF::normalizedGaussian(radius, radius, 1, radius * .1f) * sandMissing;
    Vector3 subsetStart = (this->position - sandAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 0);
    Vector3 subsetEnd = (this->position + sandAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 1);
    GridF currentSand = material.currentState.subset(subsetStart, subsetEnd);
    sandAbsorbed = sandAbsorbed.min(currentSand, 0, 0, 0);
    this->currentSatisfaction[material.name] += sandAbsorbed.sum();
    material.currentState.add(-sandAbsorbed, this->position - sandAbsorbed.getDimensions() * .5f);
}

void EnvPoint::applyDepositionOnDeath()
{
    for (auto& [materialName, amount] : materialDepositionOnDeath) {
        auto& material = EnvObject::materials[materialName];
        if (amount == 0) return;
        GridF sand = GridF::normalizedGaussian(radius, radius, 1, radius * .25f) * amount;
        material.currentState.add(sand, this->position - sand.getDimensions() * .5f);
    }
}

std::pair<GridV3, GridF> EnvPoint::computeFlowModification()
{
    if (flowEffect == Vector3()) return {EnvObject::flowfield, GridF()};
    float growingState = this->computeGrowingState2();
    if (_cachedFlowModif.size() == 0) {
        // float growingState = this->computeGrowingState();

        ScaleKelvinlet k;
        k.pos = this->position;
        k.radialScale = this->radius * .2f;
        k.scale = 10.f * /*growingState * */this->flowEffect.x;
        k.mu = .9f;
        k.v = 0.f;

        // GridV3 flow = EnvObject::flowfield;
        GridV3 flow(EnvObject::flowfield.getDimensions());
        flow.iterateParallel([&](const Vector3& p) {
            Vector3 displacement = k.evaluate(p);
            flow(p) += displacement;
        });
        _cachedFlowModif = flow;
        // return {flow, GridF(flow.getDimensions(), 1.f)};
        /*
        GridV3 gradients(2.f*radius, 2.f*radius, 1);
        Vector3 center = Vector3(radius, radius);
        gradients.iterateParallel([&](const Vector3& pos) {
            Vector3 dir = pos - center;
            float mag = std::max(0.f, 1.f - (dir.magnitude() / radius));
            gradients(pos) = dir.setMag(mag * currentGrowth);
        });

        GridV3 flow = GridV3(EnvObject::flowfield.getDimensions()).paste(gradients * 1.f, this->position - Vector3(radius, radius));
    //    GridV3 flow = GridV3(EnvObject::flowfield.getDimensions()).paste(GridV3(gauss.getDimensions(), Vector3(1, 0, 0) * this->flowEffect) * gauss, this->position - Vector3(radius, radius));
        GridF occupancy(flow.getDimensions());
        occupancy.iterateParallel([&] (const Vector3& pos) {
            occupancy(pos) = ((pos - this->position).norm2() < radius * radius ? 1.f : 0.f);
        });
        return {flow, occupancy};*/
    }
    return {EnvObject::flowfield.add(_cachedFlowModif * growingState, Vector3()), GridF(EnvObject::flowfield.getDimensions())}; // , this->position - _cachedFlowModif.getDimensions().xy() * .5f);
}

ImplicitPatch* EnvPoint::createImplicitPatch(const GridF &heights, ImplicitPrimitive *previousPrimitive)
{
    if (!geometryNeedsUpdate) return this->_patch;
    if (this->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    ImplicitPrimitive* patch;
    float growingState = 1.f; // this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    Vector3 dimensions = Vector3(radius * growingState, radius * growingState, radius * growingState* this->height);
    Vector3 patchPosition = this->position.xy() - Vector3(radius, radius) * .5f;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, dimensions, 0, {}, false);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, dimensions, radius * .25f * growingState, {}, false);
        patchPosition.z = heights(this->position.xy());
    }

    patch->position = patchPosition;
    patch->supportDimensions = dimensions;
    patch->material = this->material;
    patch->name = this->name;
    this->_patch = patch;
    this->geometryNeedsUpdate = false;
    return patch;
}

/*GridF EnvPoint::createHeightfield()
{
    return GridF();
}*/

EnvPoint &EnvPoint::translate(const Vector3 &translation)
{
    this->position.translate(translation);
    // this->evaluationPosition.translate(translation);
    for (auto& p : evaluationPositions)
        p.translate(translation);
    this->_cachedFlowModif.clear();
    this->geometryNeedsUpdate = true;
    return *this;
}

nlohmann::json EnvPoint::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["position"] = vec3_to_json(this->position);
    return json;
}


/*float EnvPoint::estimateShadowing(const GridV3& flow, const Vector3& pos) {
    Vector3 currents = flow.interpolate(this->position);
    Vector3 toPos = (pos - this->position);

    return (std::max(0.f, std::abs(std::pow(currents.normalized().dot(toPos.normalized()), 5.f)) * sign(currents.dot(toPos)) - (toPos.norm() / 100.f))); // > 0.5f ? 1.f : 0.f);
}*/
