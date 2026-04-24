#include "EnvPoint.h"
#include "EnvObject/EnvironmentalScene.h"

#include "serialization/Serializer.h"

EnvPoint::EnvPoint()
    : EnvObject()
{
}

EnvPoint *EnvPoint::clone() const
{
    auto newDefinition = new EnvPoint();
    *newDefinition = *this;
    return newDefinition;
}

EnvObjectInstance* EnvPoint::instantiate()
{
    auto newObject = new EnvPointInstance(this);
    return newObject;
}

void EnvPoint::clearKelvinlets()
{
    for (auto& k : mainKelvinlets) {
        delete k;
    }
    mainKelvinlets.clear();
}

EnvPointInstance::EnvPointInstance()
    : EnvPointInstance(nullptr)
{

}
EnvPointInstance::EnvPointInstance(EnvPoint *definition)
    : EnvObjectInstance(definition)
{

}

float EnvPointInstance::getSqrDistance(const Vector3 &position)
{
    return (position - this->position).norm2();
}

std::map<std::string, Vector3> EnvPointInstance::getAllProperties(const Vector3 &position) const
{
    Vector3 diff = (this->position - position);
    return {
        {"default", this->position},
        {"center", this->position},
        {"start", this->position},
        {"end", this->position},
        {"inside", (diff.norm2() < this->getDefinition()->radius * this->getDefinition()->radius ? Vector3(true) : Vector3(false))},
        {"normal", diff.normalized()},
        {"dir", Vector3::invalid},
        {"curvature", Vector3(this->getDefinition()->radius, 0, 0)}
    };
}

EnvPointInstance *EnvPointInstance::clone()
{
    EnvPointInstance* self = new EnvPointInstance;
    *self = *this;
    return self;
}

bool EnvPointInstance::placeInTerrain(const Vector3 &seedPosition)
{
    if (!seedPosition.isValid() || seedPosition == Vector3()) { // Not sure if second test is really needed...
        // std::cout << "WTF pos = " << seedPosition << std::endl;
        return false;
    }
    this->position = seedPosition;
    // this->translate(seedPosition.xy());
    this->recomputeEvaluationPoints();
    this->fitnessScoreAtCreation = this->evaluate();
    if (this->fitnessScoreAtCreation < this->getDefinition()->minScore) {
        // std::cout << "Evaluation of " << name << " for " << this->position << " : " << fitnessScoreAtCreation << " / " << this->minScore << std::endl;
        return false;
    }
    this->spawnTime = this->scene->currentTime;
    return true;
}

bool EnvPointInstance::placeInTerrain(const BSpline &seedCurve)
{
    if (seedCurve.empty())
        return false;
    return this->placeInTerrain(seedCurve.back());
}

void EnvPointInstance::improvePositionning(float maxDistance)
{
    Vector3 newPosition = this->position + gradientFromFieldFunction(this->getDefinition()->fitnessFunction)(this->position).normalize() * maxDistance;
    if (evaluate(newPosition) >= evaluate(this->position)) {
        this->position = newPosition;
    }
}

void EnvPointInstance::recomputeEvaluationPoints()
{
    this->evaluationPositions = {position};
}

void EnvPointInstance::applyDeposition(EnvMaterial& material)
{
    if (this->getDefinition()->materialDepositionRate.count(material.name) == 0) return;
    auto depositionProperties = this->getDefinition()->materialDepositionRate[material.name];
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
}

void EnvPointInstance::applyAbsorption(EnvMaterial& material)
{
    if (this->getDefinition()->materialAbsorptionRate.count(material.name) == 0) return;
    auto absorptionProperties = this->getDefinition()->materialAbsorptionRate[material.name];
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
    // material.currentState.iterateParallel([&] (size_t i) { material.currentState[i] = std::max(material.currentState[i], 0.f); });
}

void EnvPointInstance::applyDepositionOnDeath()
{
    for (auto& [materialName, depos] : this->getDefinition()->materialDepositionOnDeath) {
        auto& material = this->scene->materials[materialName];
        if (depos.rate == 0) return;
        GridF sand = GridF::normalizedGaussian(depos.radius, depos.radius, 1, depos.radius * .25f) * depos.rate;
        material.currentState.add(sand, this->position - sand.getDimensions() * .5f);
    }
}

GridV3& EnvPointInstance::computeFlowModification(GridV3& waterFlow, float scale)
{
    std::vector<RelativeKelvinlet> relativeFlows;
    for (size_t i = 0; i < this->getDefinition()->mainKelvinlets.size(); i++) {
        if (this->getDefinition()->mainKelvinlets[i]->valid()) {
            relativeFlows.push_back(RelativeKelvinlet(this->getDefinition()->mainKelvinlets[i], this->position, scale));
        }
    }

    const Vector3 initialFlow = waterFlow.interpolate(this->position);
    float flowAngle = initialFlow.getSignedAngleWith(Vector3(1, 0, 0));
    float flowStrength = initialFlow.length();
    waterFlow.iterateParallel([&](const Vector3& p) {
        for (const auto& relativeK : relativeFlows) {
            waterFlow[p] += relativeK.evaluate(p, flowAngle, flowStrength, true);
        }
    });
    return waterFlow;
}

ImplicitPatch* EnvPointInstance::createImplicitPatch(const GridF &heights, ImplicitPrimitive *previousPrimitive)
{
    if (!geometryNeedsUpdate) return this->_patch;
    if (this->getDefinition()->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    ImplicitPrimitive* patch;
    float growingState = 1.f; // this->computeGrowingState2();
    // float growingState = this->computeGrowingState();
    Vector3 dimensions = Vector3(this->getDefinition()->radius * growingState, this->getDefinition()->radius * growingState, this->getDefinition()->radius * growingState* this->getDefinition()->height);
    Vector3 patchPosition = this->position.xy() - Vector3(this->getDefinition()->radius, this->getDefinition()->radius) * .5f;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->getDefinition()->implicitShape, dimensions, 0, {}, false);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->getDefinition()->implicitShape, dimensions, this->getDefinition()->radius * .25f * growingState, {}, false);
        patchPosition.z() = heights(this->position.xy());
    }

    patch->position = patchPosition;
    patch->supportDimensions = dimensions;
    patch->material = this->getDefinition()->material;
    patch->name = this->getDefinition()->name;
    this->_patch = patch;
    this->geometryNeedsUpdate = false;
    return patch;
}

EnvPointInstance &EnvPointInstance::translate(const Vector3 &translation)
{
    this->position.translate(translation);
    // this->evaluationPosition.translate(translation);
    for (auto& p : evaluationPositions)
        p.translate(translation);
    this->_cachedFlowModif.clear();
    this->geometryNeedsUpdate = true;
    return *this;
}
