#include "EnvArea.h"
#include "EnvObject/EnvironmentalScene.h"

#include "serialization/Serializer.h"

EnvArea::EnvArea()
    : EnvObject()
{

}

EnvObjectInstance* EnvArea::instantiate()
{
    auto newObject = new EnvAreaInstance(this);
    return newObject;
}

EnvArea *EnvArea::clone() const
{
    auto newDefinition = new EnvArea();
    *newDefinition = *this;
    return newDefinition;
}

void EnvArea::clearKelvinlets()
{
    for (auto& k : curveKelvinlets) {
        delete k;
    }
    curveKelvinlets.clear();
}
EnvAreaInstance::EnvAreaInstance()
    : EnvAreaInstance(nullptr)
{

}
EnvAreaInstance::EnvAreaInstance(EnvArea* definition)
    : EnvObjectInstance(definition)
{

}

float EnvAreaInstance::getSqrDistance(const Vector3 &position)
{
    return (position - this->curve.estimateClosestPos(position)).norm2() * (this->curve.containsXY(position, false) ? -1.f : 1.f);
}

std::map<std::string, Vector3> EnvAreaInstance::getAllProperties(const Vector3 &position) const
{
    float closestTime = this->curve.estimateClosestTime(position);
    Vector3 closestPos = this->curve.getPoint(closestTime);
    return {
        {"default", closestPos},
        {"center", this->curve.center()},
        {"start", Vector3::invalid},
        {"end", Vector3::invalid},
        {"inside", (this->curve.containsXY(position, false) ? Vector3(true) : Vector3(false))},
        {"normal", this->curve.getNormal(closestTime)},
        {"dir", Vector3::invalid},
        {"curvature", Vector3(this->curve.getCurvature(closestTime), 0, 0)}
    };
}

EnvAreaInstance *EnvAreaInstance::clone()
{
    EnvAreaInstance* self = new EnvAreaInstance;
    *self = *this;
    return self;
}

bool EnvAreaInstance::placeInTerrain(const Vector3 &seedPosition)
{
    ShapeCurve initialCurve = ContinuousAreaOptimizer::getAreaOptimizedShape(seedPosition, this->getDefinition()->fitnessFunction, this->getDefinition()->length * this->getDefinition()->width);
    this->snake.position = seedPosition;
    return this->placeInTerrain(initialCurve);
}

bool EnvAreaInstance::placeInTerrain(const BSpline &seedCurve)
{
    ShapeCurve initialCurve = ShapeCurve(seedCurve);
    initialCurve.close().resamplePoints();
    if (seedCurve.empty()) {
        return false;
    }
    Vector3 position = initialCurve.centroid(); // The optimisation process might have moved the evaluation position greatly
    this->curve = initialCurve;
    this->curve.translate(-position);
    this->curve.resamplePoints(10);
    this->translate(position.xy());
    this->recomputeEvaluationPoints();
    this->fitnessScoreAtCreation = this->evaluate();
    if (this->fitnessScoreAtCreation < this->getDefinition()->minScore)
        return false;
    this->spawnTime = this->scene->currentTime;
    return true;
}

void EnvAreaInstance::improvePositionning(float steps)
{
    this->snake.contour = this->curve;
    this->snake.position = this->curve.center();
    this->updateCurve(this->snake.runSegmentation(steps));
}

void EnvAreaInstance::recomputeEvaluationPoints()
{
    if (this->getDefinition()->evaluateInside) {
        this->evaluationPositions = curve.randomPointsInside(curve.size());
    } else {
        this->evaluationPositions = curve.getPoints();
    }
}

void EnvAreaInstance::applyDeposition(EnvMaterial& material)
{
    if (this->getDefinition()->materialDepositionRate.count(material.name) == 0) return;
    auto depositionProperties = this->getDefinition()->materialDepositionRate[material.name];
    if (depositionProperties.rate == 0 || depositionProperties.radius == 0) return;
    ShapeCurve translatedCurve = this->curve;
    translatedCurve = translatedCurve.grow(depositionProperties.radius);
    AABBox box = AABBox(translatedCurve.AABBox());
    translatedCurve.translate(Vector3(depositionProperties.radius, depositionProperties.radius, 0) - box.min());
    // for (auto& p : translatedCurve)
        // p = p + Vector3(depositionProperties.radius, depositionProperties.radius, 0) - box.min();
    GridF sand = GridF(box.dimensions().x() + depositionProperties.radius * 2.f, box.dimensions().y() + depositionProperties.radius * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        bool inside = translatedCurve.contains(pos);
        sand(pos) = (inside ? 1.f : 0.f) * depositionProperties.rate; //gaussian(width, translatedCurve.estimateSqrDistanceFrom(Vector3(x, y, 0)));
    });
    material.currentState.add(sand, box.min() - Vector3(depositionProperties.radius, depositionProperties.radius));
    // material.currentState.add(sand.meanSmooth(depositionProperties.radius, depositionProperties.radius, 1), box.min() - Vector3(depositionProperties.radius, depositionProperties.radius));
}

void EnvAreaInstance::applyAbsorption(EnvMaterial& material)
{
    if (this->getDefinition()->materialAbsorptionRate.count(material.name) == 0) return;
    auto absorptionProperties = this->getDefinition()->materialAbsorptionRate[material.name];
    if (absorptionProperties.rate == 0 || absorptionProperties.radius == 0) return;
    ShapeCurve translatedCurve = this->curve;
    translatedCurve = translatedCurve.grow(absorptionProperties.radius);
    AABBox box = AABBox(translatedCurve.AABBox());
    translatedCurve.translate(Vector3(absorptionProperties.radius, absorptionProperties.radius, 0) - box.min());
    // for (auto& p : translatedCurve)
    // p = p + Vector3(absorptionProperties.radius, absorptionProperties.radius, 0) - box.min();
    GridF sand = GridF(box.dimensions().x() + absorptionProperties.radius * 2.f, box.dimensions().y() + absorptionProperties.radius * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        bool inside = translatedCurve.contains(pos);
        sand(pos) = -1.f * (inside ? 1.f : 0.f) * absorptionProperties.rate; //gaussian(width, translatedCurve.estimateSqrDistanceFrom(Vector3(x, y, 0)));
    });
    material.currentState.add(sand, box.min() - Vector3(absorptionProperties.radius, absorptionProperties.radius));
    // material.currentState.add(sand.meanSmooth(absorptionProperties.radius, absorptionProperties.radius, 1), box.min() - Vector3(absorptionProperties.radius, absorptionProperties.radius));
    material.currentState.iterateParallel([&](size_t i) { material.currentState[i] = std::max(material.currentState[i], 0.f); });
}

void EnvAreaInstance::applyDepositionOnDeath()
{
    for (auto& [materialName, depos] : this->getDefinition()->materialDepositionOnDeath) {
        auto& material = this->scene->materials[materialName];
        if (depos.rate == 0) return;
        AABBox box = AABBox(this->curve.AABBox());
        ShapeCurve translatedCurve = this->curve;
        for (auto& p : translatedCurve)
            p = p + Vector3(depos.radius, depos.radius, 0) - box.min();
        GridF sand = GridF(box.dimensions().x() + depos.radius * 2.f, box.dimensions().y() + depos.radius * 2.f);

        sand.iterateParallel([&] (const Vector3& pos) {
            bool inside = translatedCurve.contains(pos);
            sand(pos) = (inside ? 1.f : 0.f) * depos.rate;
        });
        material.currentState.add(sand.meanSmooth(depos.radius, depos.radius, 1), box.min() - Vector3(depos.radius, depos.radius));
    }
}

GridV3& EnvAreaInstance::computeFlowModification(GridV3& waterFlow, float scale)
{
    std::vector<KelvinletCurve*> evaluatedCurveKelvinlets;
    for (size_t i = 0; i < this->getDefinition()->curveKelvinlets.size(); i++) {
        auto& k = this->getDefinition()->curveKelvinlets[i];
        k->setCurve(&this->curve);
        if (this->getDefinition()->curveKelvinlets[i]->valid()) {
            evaluatedCurveKelvinlets.push_back(k);
            // if (auto asCurveKelvinlet = dynamic_cast<KelvinletCurve*>(this->getDefinition()->curveKelvinlets[i])) {
                // evaluatedCurveKelvinlets.push_back(asCurveKelvinlet);
            // }
        }
    }
    waterFlow.iterateParallel([&](const Vector3& p) {
        for (const auto& k : evaluatedCurveKelvinlets) {
            waterFlow[p] += k->evaluate(p) * waterFlow[p].norm();
        }
    });
    return waterFlow;
}

ImplicitPatch* EnvAreaInstance::createImplicitPatch(const GridF &heights, ImplicitPrimitive* previousPrimitive)
{
    if (!geometryNeedsUpdate) return this->_patch;

    if (this->getDefinition()->implicitShape == ImplicitPatch::PredefinedShapes::None) {
        previousPrimitive = nullptr;
        return nullptr;
    }
    BSpline translatedCurve = this->curve;
    for (Vector3& p : translatedCurve) {
        p.z() = heights(p.xy());
    }
    AABBox box(translatedCurve.AABBox());
    float growingState = 1.f;
    Vector3 offset(this->getDefinition()->width, this->getDefinition()->width, this->getDefinition()->height * growingState);
    translatedCurve.translate(-(box.min() - offset * .5f));
    box.expand(box.max() + Vector3(0, 0, this->getDefinition()->height * growingState));
    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        translatedCurve = previousPrimitive->optionalCurve;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->getDefinition()->implicitShape, box.dimensions() + offset, this->getDefinition()->height * growingState, translatedCurve, true);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->getDefinition()->implicitShape, box.dimensions() + offset, this->getDefinition()->height * growingState, translatedCurve, true);
    }
    patch->position = (box.min() - offset.xy() * .5f).xy();
    // patch->position.z() += heights(patch->position.xy());
    patch->supportDimensions = box.dimensions() + offset;
    patch->material = this->getDefinition()->material;
    patch->name = this->getDefinition()->name;
    this->_patch = patch;
    this->geometryNeedsUpdate = false;
    return patch;
}
/*
GridF EnvAreaInstance::createHeightfield()
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
EnvAreaInstance &EnvAreaInstance::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    for (auto& p : evaluationPositions)
        p.translate(translation);
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
    this->geometryNeedsUpdate = true;
    return *this;
}

void EnvAreaInstance::updateCurve(const BSpline &newCurve)
{
    if (this->getDefinition()->evaluateInside) {
        std::cerr << "Need to implement 'updateCurve' for inside...!" << std::endl;
        for (auto& evaluationPosition : this->evaluationPositions) {
            float evaluationPointClosestTime = this->curve.estimateClosestTime(evaluationPosition);
            Vector3 relativeDisplacementFromEvaluationToCurve = (evaluationPosition - this->curve.getPoint(evaluationPointClosestTime));
            evaluationPosition = newCurve.getPoint(evaluationPointClosestTime) + relativeDisplacementFromEvaluationToCurve;
        }
    } else {
        this->evaluationPositions = newCurve.getPoints();
    }
    this->curve = newCurve;
    this->_cachedFlowModif.clear();
    this->_cachedHeightfield.clear();
}
