#include "EnvObject.h"

#include <fstream>

GridV3 initFlow();

GridV3 EnvObject::flowfield = initFlow();
GridV3 EnvObject::terrainNormals;
GridF EnvObject::sandDeposit;
std::map<std::string, EnvObject*> EnvObject::availableObjects;
std::vector<EnvObject*> EnvObject::instantiatedObjects;
float EnvObject::flowImpactFactor = .99f;

GridV3 initFlow() {
    if (EnvObject::flowfield.empty()) {
        EnvObject::flowfield = GridV3(100, 100, 1, Vector3(0, 0, 0));
        EnvObject::flowfield.raiseErrorOnBadCoord = false;
        EnvObject::flowfield.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
    }
    return EnvObject::flowfield;
}

EnvObject::EnvObject()
{
    initFlow();
}

EnvObject::~EnvObject()
{

}

void EnvObject::readFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    auto json = nlohmann::json::parse(toLower(content));

    for (auto& obj : json) {
        if (!obj.contains("type")) {
            throw std::domain_error("No type given for Environmental Object defined as " + nlohmann::to_string(obj));
        }
        std::string objName = obj["name"];

        if (obj["type"] == "point")
            EnvObject::availableObjects[objName] = EnvPoint::fromJSON(obj);
        else if (obj["type"] == "curve")
            EnvObject::availableObjects[objName] = EnvCurve::fromJSON(obj);
        else if (obj["type"] == "area")
            EnvObject::availableObjects[objName] = EnvArea::fromJSON(obj);
        else {
            throw std::domain_error("Unrecognized type for Environmental Object defined as " + nlohmann::to_string(obj));
        }
    }
}

EnvObject *EnvObject::fromJSON(nlohmann::json content)
{
    EnvObject* obj = nullptr;
    std::string objName = content["name"];
    std::string objType = content["type"];
    float sandEffect = content["sand"];
    TerrainTypes material = materialFromString(content["material"]);
    ImplicitPatch::PredefinedShapes shape = predefinedShapeFromString(content["geometry"]);
    Vector3 dimensions = json_to_vec3(content["dimensions"]);

    Vector3 flowEffect;
    if (objType == "point") {
        auto asPoint = new EnvPoint;
        asPoint->radius = dimensions.x;
        obj = asPoint;
        flowEffect = Vector3(content["flow"], content["flow"], content["flow"]);
    } else if (objType == "curve") {
        auto asCurve = new EnvCurve;
        asCurve->width = dimensions.x;
        asCurve->length = dimensions.y;
        obj = asCurve;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "area") {
        auto asArea = new EnvArea;
        asArea->width = dimensions.x;
        obj = asArea;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    }

    obj->name = toLower(objName);
    obj->flowEffect = flowEffect;
    obj->sandEffect = sandEffect;
    obj->fittingFunction = EnvObject::parseFittingFunction(content["rule"]);
    obj->material = material;
    obj->implicitShape = shape;
    return obj;
}

std::function<float (Vector3)> EnvObject::parseFittingFunction(std::string formula)
{
    std::vector<std::string> tokens = split(toLower(formula), " ");

    std::vector<std::pair<std::string, std::function<float (Vector3)>>> operations;

    std::string currentOperation = "+";

    for (auto token : tokens) {
//        float tokenScore = 0.f;
        // Extract the distance function
        if (startsWith(token, "d2(")) {
            std::string object = token.substr(std::string("d2(").size());
            object = object.substr(0, object.size() - std::string(")").size());
            operations.push_back({currentOperation, [=](const Vector3& pos) {
                return EnvObject::getSqrDistanceTo(object, pos).first;
            }});
        }
        else if (startsWith(token, "<")) {
            std::string objects = token.substr(std::string("<").size());
            objects = objects.substr(0, objects.size() - std::string(">").size());
            std::string object1 = objects.substr(0, objects.find(","));
            std::string object2 = objects.substr(objects.find(",") + 1);
            operations.push_back({currentOperation, [=](const Vector3& pos) {
                auto [vec1, obj1] = EnvObject::getVectorOf(object1, pos);
                auto [vec2, obj2] = EnvObject::getVectorOf(object2, pos);
                return vec1.dot(vec2);
            }});
        } else if (std::atof(token.c_str())) { // Is a number
            operations.push_back({currentOperation, [=](const Vector3& pos) {
                return std::atof(token.c_str());
            }});
        } else {
            currentOperation = token;
            continue;
        }


    }

    std::function<float (Vector3)> func = [=](const Vector3& pos) -> float {
        float score = 0.f;
        for (auto& [sign, function] : operations) {
            float funcRes = function(pos);
            if (sign == "+")
                score += funcRes;
            else if (sign == "-")
                score -= funcRes;
            else if (sign == "/")
                score /= funcRes;
            else if (sign == "*")
                score *= funcRes;
        }
        return score;
    };
    return func;
}

std::pair<std::string, std::string> EnvObject::extractNameAndComplement(std::string variable)
{
    auto splitted = split(variable, ".");
    if (splitted.size() == 0) splitted = {"", ""};
    else if (splitted.size() == 1) splitted = {splitted[0], ""};
    return {splitted[0], splitted[1]};
}

std::pair<float, EnvObject *> EnvObject::getSqrDistanceTo(std::string objectName, const Vector3 &position)
{
    auto [name, complement] = EnvObject::extractNameAndComplement(objectName);
    float minDist = std::numeric_limits<float>::max();
    EnvObject* bestElem = nullptr;
    for (auto& instance : EnvObject::instantiatedObjects) {
        if (instance->name != name) continue;
        float distance = instance->getSqrDistance(position, complement);
        if (distance < minDist) {
            minDist = distance;
            bestElem = instance;
        }
    }
    return {minDist, bestElem};
}

std::pair<Vector3, EnvObject *> EnvObject::getVectorOf(std::string objectName, const Vector3 &position)
{
    auto [name, complement] = EnvObject::extractNameAndComplement(objectName);
    if (name == "current") {
        return {EnvObject::flowfield.at(position), nullptr};
    }
    auto object = EnvObject::getSqrDistanceTo(objectName, position).second;
    return {object->getVector(position, complement), object};
}

EnvObject *EnvObject::instantiate(std::string objectName)
{
    auto object = EnvObject::availableObjects[objectName]->clone();
    EnvObject::instantiatedObjects.push_back(object);
    return object;
}

void EnvObject::removeAllObjects()
{
    for (auto& object : EnvObject::instantiatedObjects) {
        delete object;
    }
    EnvObject::instantiatedObjects.clear();
}

void EnvObject::applyEffects()
{
    GridV3 totalNewFlow(EnvObject::flowfield.getDimensions());
    GridF totalOccupancy(EnvObject::flowfield.getDimensions());
    for (auto& object : EnvObject::instantiatedObjects) {
        const auto [flow, occupancy] = object->computeFlowModification();
        totalNewFlow += flow;
        totalOccupancy += occupancy;
    }
    totalNewFlow.iterateParallel([&] (size_t i) {
        totalNewFlow[i] = (totalOccupancy[i] != 0.f ? totalNewFlow[i] / totalOccupancy[i] : totalNewFlow[i]);
    });
    EnvObject::flowfield = EnvObject::flowfield * (1.f - EnvObject::flowImpactFactor) + totalNewFlow * EnvObject::flowImpactFactor;

    for (auto& object : EnvObject::instantiatedObjects) {
        object->applySandDeposit();
    }
    EnvObject::sandDeposit = EnvObject::sandDeposit.wrapWith(EnvObject::flowfield * 3.f);
}

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

EnvPoint *EnvPoint::clone()
{
    EnvPoint* self = new EnvPoint;
    *self = *this;
    return self;
}

void EnvPoint::applySandDeposit()
{
    GridF sand = GridF::gaussian(radius, radius, 1, radius * .5f) * this->sandEffect;
    EnvObject::sandDeposit.add(sand, sand.getDimensions() * .5f);

}

std::pair<GridV3, GridF> EnvPoint::computeFlowModification()
{
//    GridV3 flowSubset = EnvObject::flowfield.subset(Vector3(this->position.x - radius, this->position.y - radius), Vector3(this->position.x + radius, this->position.y + radius));
    GridF gauss = GridF::gaussian(radius, radius, 1, radius * .5f).normalize();
    GridV3 flow = GridV3(EnvObject::flowfield.getDimensions()).paste(GridV3(gauss.getDimensions(), Vector3(1, 0, 0) * this->flowEffect) * gauss, this->position - Vector3(radius, radius));
    GridF occupancy(flow.getDimensions());
    occupancy.iterateParallel([&] (const Vector3& pos) {
        occupancy(pos) = ((pos - this->position).norm2() < radius * radius ? 1.f : 0.f);
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvPoint::createImplicitPatch()
{
    ImplicitPrimitive* patch = ImplicitPatch::createPredefinedShape(this->implicitShape, Vector3(radius, radius, radius), 0);
    patch->position = this->position;
    patch->material = this->material;
    return patch;
}

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
        sand.at(pos) = gaussian(width * .5f, translatedCurve.estimateSqrDistanceFrom(pos)) * this->sandEffect;
    });
//    sand *= this->sandEffect;
    EnvObject::sandDeposit.add(sand, box.min() /*- sand.getDimensions() * .5f*/);
}

std::pair<GridV3, GridF> EnvCurve::computeFlowModification()
{
    Vector3 objectWidth = Vector3(width, width, 0);
    Vector3 halfWidth = objectWidth * .5f;
    ShapeCurve translatedCurve = this->curve;
    for (auto& p : translatedCurve)
        p.z = 0;
    AABBox box = AABBox(translatedCurve.points);
    box.expand({box.min() - halfWidth, box.max() + halfWidth});


    GridV3 flow(EnvObject::flowfield.getDimensions());
    GridF occupancy(flow.getDimensions());

    flow.iterateParallel([&] (const Vector3& pos) {
        float closestTime = translatedCurve.estimateClosestTime(pos);
        Vector3 closestPos = translatedCurve.getPoint(closestTime);
        float sqrDist = (closestPos - pos).norm2();
        if (sqrDist > (width * .5f) * (width * .5f)) return;
        float gauss = normalizedGaussian(width * .5f, sqrDist);
        Vector3 impact = this->flowEffect * gauss;
        auto [direction, normal, binormal] = translatedCurve.getFrenetFrame(closestTime);
        if (impact.y != 0) { // Add an impact on the normal direction: need to change the normal to be in the right side of the curve
            if ((translatedCurve.getPoint(closestTime) - pos).dot(normal) > 0)
                normal *= -1.f;
        }
        flow(pos) = impact.changedBasis(direction, normal, binormal);
        occupancy(pos) = 1.f;
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvCurve::createImplicitPatch()
{
    BSpline translatedCurve = this->curve;
    AABBox box(this->curve.points);
    translatedCurve.translate(box.min());
    ImplicitPrimitive* patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + Vector3(0, 0, this->width*2.f), this->width, translatedCurve);
    patch->position = box.min();
    patch->material = this->material;
    return patch;
}

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
        return (position - this->area.center()).norm2();
    else if (complement == "border") // Just keep the absolute value
        return this->area.estimateSqrDistanceFrom(position);

    float dist = this->area.estimateSignedDistanceFrom(position);
    return (dist * dist) * (dist > 0 ? 1.f : -1.f); // Keep the sign
}

Vector3 EnvArea::getVector(const Vector3 &position, std::string complement)
{
    if (complement == "direction")
        return this->area.getDirection(area.estimateClosestTime(position));
    else if (complement == "normal")
        return this->area.getNormal(area.estimateClosestTime(position));
    else if (complement == "binormal")
        return this->area.getBinormal(area.estimateClosestTime(position));
    return Vector3();
}

EnvArea *EnvArea::clone()
{
    EnvArea* self = new EnvArea;
    *self = *this;
    return self;
}

void EnvArea::applySandDeposit()
{
    AABBox box = AABBox(this->area.points);
    ShapeCurve translatedCurve = this->area;
    for (auto& p : translatedCurve)
        p = p + Vector3(width, width, 0) - box.min();
    GridF sand = GridF(box.dimensions().x + width * 2.f, box.dimensions().y + width * 2.f);

    sand.iterateParallel([&] (const Vector3& pos) {
        bool inside = translatedCurve.contains(pos);
        sand(pos) = (inside ? 1.f : 0.f) * sandEffect; //gaussian(width, translatedCurve.estimateSqrDistanceFrom(Vector3(x, y, 0)));
    });
//    sand *= this->sandEffect;
    //EnvObject::sandDeposit.add(sand, box.min());
}

std::pair<GridV3, GridF> EnvArea::computeFlowModification()
{
    Vector3 objectWidth = Vector3(width, width, 0);
    Vector3 halfWidth = objectWidth * .5f;
    ShapeCurve translatedCurve = this->area;
    for (auto& p : translatedCurve)
        p.z = 0;
    AABBox box = AABBox(translatedCurve.points);
    box.expand({box.min() - halfWidth, box.max() + halfWidth});


    GridV3 flow(EnvObject::flowfield.getDimensions());
    GridF occupancy(EnvObject::flowfield.getDimensions());

    GridF dist(flow.getDimensions());
    flow.iterateParallel([&] (const Vector3& pos) {
        dist(pos) = (box.contains(pos) && translatedCurve.contains(pos, false) ? 1.f : 0.f);
    });
    dist = dist.toDistanceMap(true, false);
    GridV3 grad = dist.meanSmooth().gradient() * -1.f;
    for (auto& v : grad)
        v.normalize();

    flow.iterateParallel([&] (const Vector3& pos) {
        if (!box.contains(pos) || !translatedCurve.contains(pos))
            return;

        float closestTime = translatedCurve.estimateClosestTime(pos);
        Vector3 closestPos = translatedCurve.getPoint(closestTime);

        float distanceToBorder = (pos - closestPos).norm();
        float distFactor = clamp(distanceToBorder / (width * .5f), 0.f, 1.f); // On border = 1, at w/2 = 0, more inside = 0
        Vector3 previousFlow = EnvObject::flowfield(pos);
        // Change the order of the Frenet Frame to get the direction in the direction of the "outside" and the normal is along the borders
//            auto [normal, direction, binormal] = translatedCurve.getFrenetFrame(closestTime);
        // We will use the distance map to get the direction, then we know (0, 0, 1) is the binormal (2D shape), so normal is cross product.
        Vector3 direction = grad(pos);
        Vector3 binormal = Vector3(0, 0, 1);
        Vector3 normal = direction.cross(binormal); // Direction and binormal are normalized.
        Vector3 impact = this->flowEffect * distFactor;
        // Ignore the normal for now, I can't find any logic with it...
//            if (impact.y != 0) { // Add an impact on the normal direction: need to change the normal to be in the right side of the curve
//                if ((translatedCurve.getPoint(closestTime) - Vector3(x, y, 0)).dot(direction) > 0)
//                    normal *= -1.f;
//            }
        flow.at(pos) = impact.changedBasis(direction, normal, binormal);// + impact * previousFlow;
        occupancy.at(pos) = 1.f;
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvArea::createImplicitPatch()
{
    BSpline translatedCurve = this->area;
    AABBox box(this->area.points);
    translatedCurve.translate(box.min());
    ImplicitPrimitive* patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + Vector3(0, 0, this->width * 2.f), this->width, translatedCurve);
    patch->position = box.min();
    patch->material = this->material;
    return patch;
}
