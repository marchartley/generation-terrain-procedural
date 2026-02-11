#include "EnvPoint.h"
#include "EnvObject/EnvironmentalScene.h"

#include "serialization/Serializer.h"

EnvPoint::EnvPoint()
    : EnvObject()
{

}

/*
EnvPoint *EnvPoint::fromJSON(nlohmann::json content)
{
    return dynamic_cast<EnvPoint*>(EnvObject::fromJSON(content));
}
*/
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

/*EnvPoint *EnvPoint::instantiate(std::string objectName)
{
    return dynamic_cast<EnvPoint*>(this->scene->instantiate(objectName));
}*/

bool EnvPoint::placeInTerrain(const Vector3 &seedPosition)
{
    if (!seedPosition.isValid() || seedPosition == Vector3()) { // Not sure if second test is really needed...
        // std::cout << "WTF pos = " << seedPosition << std::endl;
        return false;
    }
    this->position = seedPosition;
    // this->translate(seedPosition.xy());
    this->recomputeEvaluationPoints();
    this->fitnessScoreAtCreation = this->evaluate();
    if (this->fitnessScoreAtCreation < this->minScore) {
        // std::cout << "Evaluation of " << name << " for " << this->position << " : " << fitnessScoreAtCreation << " / " << this->minScore << std::endl;
        return false;
    }
    this->spawnTime = this->scene->currentTime;
    return true;
}

bool EnvPoint::placeInTerrain(const BSpline &seedCurve)
{
    if (seedCurve.empty())
        return false;
    return this->placeInTerrain(seedCurve.points.back());
}

void EnvPoint::improvePositionning(float maxDistance)
{
    this->position += gradientFromFieldFunction(this->fitnessFunction)(this->position).normalize() * maxDistance;
}

void EnvPoint::recomputeEvaluationPoints()
{
    this->evaluationPositions = {position};
}

void EnvPoint::applyDeposition(EnvMaterial& material)
{
    if (this->materialDepositionRate.count(material.name) == 0) return;
    auto depositionProperties = this->materialDepositionRate[material.name];
    if (depositionProperties.rate == 0 || depositionProperties.radius == 0) return;

    AABBox box = AABBox(this->position - Vector3(1, 1, 0) * depositionProperties.radius, this->position + Vector3(1, 1, 0) * depositionProperties.radius);
    GridF deposition = GridF(box.dimensions().x(), box.dimensions().y());
    Vector3i center = deposition.getDimensions() / 2;
    deposition.iterateParallel([&] (const Vector3& pos) {
        float distToCenter = (pos - center).norm2();
        // float amount = normalizedGaussian(depositionProperties.radius * .25f, distToCurve);
        float amount = (distToCenter < depositionProperties.radius * depositionProperties.radius ? 1.f : 0.f);

        deposition.at(pos) = amount;
    });
    material.currentState.add(deposition * depositionProperties.rate, box.min().xy());

    // if (toLower(this->name) == "coralpolyp") std::cout << "Deposing " << material.name << "? " << (this->materialDepositionRate[material.name] == 0 ? "No" : "Yes") << std::endl;
    /*if (this->materialDepositionRate[material.name].rate == 0) return;
    float growingState = this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    GridF sand = GridF::normalizedGaussian(radius, radius, 1, radius * .25f) * this->materialDepositionRate[material.name].rate * growingState;
    material.currentState.add(sand, this->position - sand.getDimensions() * .5f);
    */
}

void EnvPoint::applyAbsorption(EnvMaterial& material)
{
    if (this->materialAbsorptionRate.count(material.name) == 0) return;
    auto absorptionProperties = this->materialAbsorptionRate[material.name];
    if (absorptionProperties.rate == 0 || absorptionProperties.radius == 0) return;

    AABBox box = AABBox(this->position - Vector3(1, 1, 0) * absorptionProperties.radius, this->position + Vector3(1, 1, 0) * absorptionProperties.radius);
    GridF absorption = GridF(box.dimensions().x(), box.dimensions().y());
    Vector3i center = absorption.getDimensions() / 2;
    absorption.iterateParallel([&] (const Vector3& pos) {
        float distToCenter = (pos - center).norm2();
        // float amount = normalizedGaussian(absorptionProperties.radius * .25f, distToCurve);
        float amount = (distToCenter < absorptionProperties.radius * absorptionProperties.radius ? 1.f : 0.f);

        absorption.at(pos) = amount;
    });
    material.currentState.add(absorption * absorptionProperties.rate * -1.f, box.min().xy());
    material.currentState.iterateParallel([&] (size_t i) { material.currentState[i] = std::max(material.currentState[i], 0.f); });
    /*
    if (this->materialAbsorptionRate[material.name].rate == 0) return;
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
    */
}

void EnvPoint::applyDepositionOnDeath()
{
    for (auto& [materialName, depos] : materialDepositionOnDeath) {
        auto& material = this->scene->materials[materialName];
        if (depos.rate == 0) return;
        GridF sand = GridF::normalizedGaussian(radius, radius, 1, radius * .25f) * depos.rate;
        material.currentState.add(sand, this->position - sand.getDimensions() * .5f);
    }
}

GridV3 EnvPoint::computeFlowModification()
{
    /*
    if (flowEffect == Vector3()) return {this->scene->flowfield, GridF()};
    float growingState = this->computeGrowingState2();
    if (_cachedFlowModif.size() == 0) {
        // float growingState = this->computeGrowingState();

        ScaleKelvinlet k;
        k.pos = this->position;
        k.radialScale = this->radius * .2f;
        k.force = 10.f * this->flowEffect.x();
        k.mu = .9f;
        k.v = 0.f;

        // GridV3 flow = this->scene->flowfield;
        GridV3 flow(this->scene->flowfield.getDimensions());
        flow.iterateParallel([&](const Vector3i& p) {
            Vector3 displacement = k.evaluate(p);
            flow(p) += displacement;
        });
        _cachedFlowModif = flow;
    }
    return {this->scene->flowfield.add(_cachedFlowModif * growingState, Vector3()), GridF(this->scene->flowfield.getDimensions())}; // , this->position - _cachedFlowModif.getDimensions().xy() * .5f);
    */
    return this->scene->flowfield;
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
        patchPosition.z() = heights(this->position.xy());
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
/*
nlohmann::json EnvPoint::toJSON() const
{
    auto json = EnvObject::toJSON();
    json["position"] = this->position;
    return json;
}
*/

/*float EnvPoint::estimateShadowing(const GridV3& flow, const Vector3& pos) {
    Vector3 currents = flow.interpolate(this->position);
    Vector3 toPos = (pos - this->position);

    return (std::max(0.f, std::abs(std::pow(currents.normalized().dot(toPos.normalized()), 5.f)) * sign(currents.dot(toPos)) - (toPos.norm() / 100.f))); // > 0.5f ? 1.f : 0.f);
}*/
