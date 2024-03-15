#include "EnvPoint.h"

EnvPoint::EnvPoint()
    : EnvObject()
{

}

EnvObject *EnvPoint::fromJSON(nlohmann::json content)
{
    return EnvObject::fromJSON(content);
}

float EnvPoint::getSqrDistance(const Vector3 &position, std::string complement)
{
    return (position - this->position).norm2();
}

Vector3 EnvPoint::getVector(const Vector3 &position, std::string complement)
{
    return Vector3();
}

Vector3 EnvPoint::getNormal(const Vector3 &position)
{
    return (position - this->position).normalized();
}

Vector3 EnvPoint::getDirection(const Vector3 &position)
{
    return Vector3(false);
}

Vector3 EnvPoint::getProperty(const Vector3& position, std::string prop) const
{
    if (prop == "center") {
        return this->position;
    } else if (prop == "start") {
        return this->position;
    } else if (prop == "end") {
        return this->position;
    } else if (prop == "inside") {
        return ((position - this->position).norm2() < this->radius * this->radius ? Vector3(true) : Vector3(false));
    }
    return this->position; // Default
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

void EnvPoint::applySandDeposit()
{
    GridF sand = GridF::normalizedGaussian(radius, radius, 1, radius * 1.f) * this->sandEffect * this->computeGrowingState();
    EnvObject::sandDeposit.add(sand, this->position - sand.getDimensions() * .5f);

}

void EnvPoint::applySandAbsorption()
{
    if (this->needsForGrowth.count("sand") == 0) return;

    float sandMissing = (needsForGrowth["sand"] - currentSatisfaction["sand"]) * (this->computeGrowingState() + .01f) / (1.01f);

    GridF sandAbsorbed = GridF::normalizedGaussian(radius, radius, 1, radius * .1f) * sandMissing;
    Vector3 subsetStart = (this->position - sandAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 0);
    Vector3 subsetEnd = (this->position + sandAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 1);
    GridF currentSand = EnvObject::sandDeposit.subset(subsetStart, subsetEnd);
    sandAbsorbed = sandAbsorbed.min(currentSand, 0, 0, 0);
    this->currentSatisfaction["sand"] += sandAbsorbed.sum();
    EnvObject::sandDeposit.add(-sandAbsorbed, this->position - sandAbsorbed.getDimensions() * .5f);
}

void EnvPoint::applyPolypDeposit()
{
    if (this->polypEffect == 0) return;
    GridF polyp = GridF::normalizedGaussian(radius, radius, 1, radius * 1.f) * this->polypEffect * this->computeGrowingState();
    EnvObject::polypDeposit.add(polyp, this->position - polyp.getDimensions() * .5f);
}

void EnvPoint::applyPolypAbsorption()
{
    if (this->needsForGrowth.count("polyp") == 0) return;

    float polypMissing = (needsForGrowth["polyp"] - currentSatisfaction["polyp"]) * (this->computeGrowingState() + .01f) / (1.01f);

    GridF polypAbsorbed = GridF::normalizedGaussian(radius, radius, 1, radius * .1f) * polypMissing;
    Vector3 subsetStart = (this->position - polypAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 0);
    Vector3 subsetEnd = (this->position + polypAbsorbed.getDimensions() * .5f).xy() + Vector3(0, 0, 1);
    GridF currentPolyp = EnvObject::polypDeposit.subset(subsetStart, subsetEnd);
    polypAbsorbed = polypAbsorbed.min(currentPolyp, 0, 0, 0);
    this->currentSatisfaction["polyp"] += polypAbsorbed.sum();
    EnvObject::polypDeposit.add(-polypAbsorbed, this->position - polypAbsorbed.getDimensions() * .5f);
}

std::pair<GridV3, GridF> EnvPoint::computeFlowModification()
{
//    GridV3 flowSubset = EnvObject::flowfield.subset(Vector3(this->position.x - radius, this->position.y - radius), Vector3(this->position.x + radius, this->position.y + radius));
    //GridF gauss = GridF::gaussian(2.f*radius, 2.f*radius, 1, radius/* * .5f*/).normalize();

    float currentGrowth = 1.f; //computeGrowingState();

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
    return {flow, occupancy};
}

ImplicitPatch* EnvPoint::createImplicitPatch(ImplicitPrimitive *previousPrimitive)
{

    ImplicitPrimitive* patch;
    float currentGrowth = this->computeGrowingState();
    Vector3 dimensions = Vector3(radius * currentGrowth, radius * currentGrowth, radius * currentGrowth * this->height);
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, dimensions, 0, {}, true);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, dimensions, radius * .25f * currentGrowth, {}, true);
    }

    patch->supportDimensions = dimensions;
    patch->position = this->position.xy() - Vector3(radius, radius) * .5f;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

GridF EnvPoint::createHeightfield() const
{
    return GridF();
}

EnvObject &EnvPoint::translate(const Vector3 &translation)
{
    this->position += translation;
    return *this;
}


float EnvPoint::estimateShadowing(const GridV3& flow, const Vector3& pos) {
    Vector3 currents = flow.interpolate(this->position);
    Vector3 toPos = (pos - this->position);

    return (std::max(0.f, std::abs(std::pow(currents.normalized().dot(toPos.normalized()), 5.f)) * sign(currents.dot(toPos)) - (toPos.norm() / 100.f))); // > 0.5f ? 1.f : 0.f);
}
