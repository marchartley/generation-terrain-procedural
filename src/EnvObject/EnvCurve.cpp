#include "EnvCurve.h"

EnvCurve::EnvCurve()
    : EnvObject()
{

}

EnvObject *EnvCurve::fromJSON(nlohmann::json content)
{
    return EnvObject::fromJSON(content);
}

float EnvCurve::getSqrDistance(const Vector3 &position, std::string complement)
{
    if (complement == "front")
        return (position - this->curve.points.front()).norm2();
    else if (complement == "end")
        return (position - this->curve.points.back()).norm2();
    else if (complement == "center")
        return (position - this->curve.center()).norm2();
    else
        return (position - this->curve.estimateClosestPos(position)).norm2();
}

Vector3 EnvCurve::getVector(const Vector3 &position, std::string complement)
{
    if (complement == "direction")
        return this->curve.getDirection(curve.estimateClosestTime(position));
    else if (complement == "normal")
        return this->curve.getNormal(curve.estimateClosestTime(position));
    else if (complement == "binormal")
        return this->curve.getBinormal(curve.estimateClosestTime(position));
    return Vector3();
}

Vector3 EnvCurve::getNormal(const Vector3 &position)
{
    return this->curve.getNormal(this->curve.estimateClosestTime(position));
}

Vector3 EnvCurve::getDirection(const Vector3 &position)
{
    return this->curve.getDirection(this->curve.estimateClosestTime(position));
}

Vector3 EnvCurve::getProperty(const Vector3& position, std::string prop) const
{
    if (prop == "center") {
        return this->curve.center();
    } else if (prop == "start") {
        return this->curve.points.front();
    } else if (prop == "end") {
        return this->curve.points.back();
    } else if (prop == "inside") {
        return ((position - this->curve.estimateClosestPos(position)).norm2() < this->width * this->width ? Vector3(true) : Vector3(false));
    }
    return this->curve.estimateClosestPos(position); // Default
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

void EnvCurve::applySandDeposit()
{
    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();
    GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        sand.at(pos) = normalizedGaussian(width * .5f, translatedCurve.estimateSqrDistanceFrom(pos)) * this->sandEffect;
    });
//    std::cout << "Deposition of " << sand.sum() << " (" << sand.min() << " -- " << sand.max() << ") at " << box.min() << " for " << sand << std::endl;
//    sand *= this->sandEffect;
    EnvObject::sandDeposit.add(sand, box.min().xy() - Vector3(width, width));
}

void EnvCurve::applySandAbsorption()
{
    return;
}

void EnvCurve::applyPolypDeposit()
{
    if (this->polypEffect == 0) return;
    AABBox box = AABBox(this->curve.points);
    BSpline translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();
    GridF polyp = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

    polyp.iterateParallel([&] (const Vector3& pos) {
        polyp.at(pos) = normalizedGaussian(width * .5f, translatedCurve.estimateSqrDistanceFrom(pos)) * this->polypEffect;
    });
//    std::cout << "Deposition of " << polyp.sum() << " (" << polyp.min() << " -- " << polyp.max() << ") at " << box.min() << " for " << polyp << std::endl;
//    polyp *= this->polypEffect;
    EnvObject::polypDeposit.add(polyp, box.min().xy() - Vector3(width, width));
}

void EnvCurve::applyPolypAbsorption()
{
    return;
}

std::pair<GridV3, GridF> EnvCurve::computeFlowModification()
{
    Vector3 objectWidth = Vector3(width, width, 0);
    Vector3 halfWidth = objectWidth * .5f;
    BSpline translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p.z = 0;
    AABBox box = AABBox(translatedCurve.points);
    box.expand({box.min() - halfWidth, box.max() + halfWidth});


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
        /*if (impact.y != 0) { // Add an impact on the normal direction: need to change the normal to be in the right side of the curve
            if ((closestPos - pos).dot(normal) > 0)
                normal *= -1.f;
        }*/
        direction *= direction.dot(initialFlow);
        normal *= -normal.dot(initialFlow);
        binormal *= binormal.dot(initialFlow);
        flow(pos) = impact.changedBasis(direction, normal, binormal);
        occupancy(pos) = 1.f;
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvCurve::createImplicitPatch(ImplicitPrimitive *previousPrimitive)
{
    BSpline translatedCurve = this->curve;
    AABBox box(this->curve.points);
    float height = this->height * this->computeGrowingState();
    Vector3 offset(this->width, this->width, height * .5f);
    if (height == 0){
        box = AABBox({box.center()});
        offset = Vector3();
    }
    translatedCurve.translate(-(box.min() - offset * .5f));

    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, true);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, height, translatedCurve, true);
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

EnvObject &EnvCurve::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    return *this;
}
