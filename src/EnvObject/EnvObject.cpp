#include "EnvObject.h"
#include "ExpressionParser.h"

#include <fstream>

GridV3 initFlow();

GridV3 EnvObject::flowfield = initFlow();
GridV3 EnvObject::terrainNormals;
GridF EnvObject::sandDeposit;
std::map<std::string, EnvObject*> EnvObject::availableObjects;
std::vector<EnvObject*> EnvObject::instantiatedObjects;
float EnvObject::flowImpactFactor = .8f;
int EnvObject::currentMaxID = -1;

std::map<std::string, GridV3> EnvObject::allVectorProperties;
std::map<std::string, GridF> EnvObject::allScalarProperties;

GridV3 initFlow() {
    if (EnvObject::flowfield.empty()) {
        EnvObject::flowfield = GridV3(100, 100, 1, Vector3(0, 0, 0));
        EnvObject::flowfield.raiseErrorOnBadCoord = false;
        EnvObject::flowfield.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        EnvObject::sandDeposit = GridF(EnvObject::flowfield.getDimensions());
        EnvObject::sandDeposit.raiseErrorOnBadCoord = false;
        EnvObject::sandDeposit.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
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

    for (auto& [name, obj] : EnvObject::availableObjects) {
        obj->fittingFunction = EnvObject::parseFittingFunction(obj->s_FittingFunction, obj->name);
    }
//    precomputeTerrainProperties(Heightmap());
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
        asArea->height = dimensions.z;
        obj = asArea;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    }

    obj->name = toLower(objName);
    obj->flowEffect = flowEffect;
    obj->sandEffect = sandEffect;
    obj->s_FittingFunction = content["rule"];
//    obj->fittingFunction = EnvObject::parseFittingFunction(content["rule"]);
    obj->material = material;
    obj->implicitShape = shape;
    obj->inputDimensions = dimensions;
    return obj;
}

std::function<float (Vector3)> EnvObject::parseFittingFunction(std::string formula, std::string currentObject)
{
    if (formula == "")
        formula = "0.0";

    std::map<std::string, Variable> variables;
    for (auto& [name, obj] : EnvObject::availableObjects) {
        variables[name] = Vector3();
        variables[name + ".center"] = Vector3();
        variables[name + ".start"] = Vector3();
        variables[name + ".end"] = Vector3();
        variables[name + ".normal"] = Vector3();
        variables[name + ".dir"] = Vector3();
    }

    variables["current.center"] = Vector3();
    variables["current.start"] = Vector3();
    variables["current.end"] = Vector3();
    variables["current.normal"] = Vector3();
    variables["current.dir"] = Vector3();
    variables["current.vel"] = float();

    ExpressionParser parser;
    variables["pos"] = Vector3();
    parser.validate(formula, variables);
    std::set<std::string> neededVariables = parser.extractAllVariables(formula);
    auto _func = parser.parse(formula, variables);
    return [&, _func, neededVariables, currentObject](Vector3 pos) -> float {
        std::map<std::string, Variable> vars;
        for (auto& [prop, map] : EnvObject::allVectorProperties) {
            if (neededVariables.count(prop) && !map(pos).isValid()) { // A variable is needed but there are no value attributed (eg : object not instantiated yet)
                if (prop == currentObject || startsWith(prop, currentObject + ".")) {
                    std::cout << prop << " -> " << currentObject << std::endl;
                    vars[prop] = pos;
                } else {
                    return std::numeric_limits<float>::max(); // Return max value
                }
            } else {
                vars[prop] = map(pos);
            }
        }
        for (auto& [prop, map] : EnvObject::allScalarProperties) {
            vars[prop] = map(pos);
        }
        vars["pos"] = pos;
        float score = _func(vars);
        return score;
    };

/*
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
            std::string object1 = objects.substr(0, objects.find("dot"));
            std::string object2 = objects.substr(objects.find("dot") + 1);
            operations.push_back({currentOperation, [=](const Vector3& pos) {
                auto [vec1, obj1] = EnvObject::getVectorOf(object1, pos);
                auto [vec2, obj2] = EnvObject::getVectorOf(object2, pos);
                if (obj1 == nullptr || obj2 == nullptr)
                    return 10000000.f;
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
        return (score >= 0 ? score : 100000.f);
    };
    return func;
    */
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
    if (object == nullptr) return {Vector3::invalid(), nullptr};
    return {object->getVector(position, complement), object};
}

EnvObject *EnvObject::instantiate(std::string objectName)
{
    EnvObject::currentMaxID++;
    auto object = EnvObject::availableObjects[objectName]->clone();
    object->ID = EnvObject::currentMaxID;
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
    EnvObject::sandDeposit = EnvObject::sandDeposit.wrapWith(EnvObject::flowfield.meanSmooth(3, 3, 1) * 10.f);
}

void EnvObject::beImpactedByEvents()
{
    // TODO!!!
    for (auto& obj : EnvObject::instantiatedObjects) {
        obj->growingState = std::min(obj->growingState + .2f, 1.f);
    }
}

float EnvObject::evaluate(const Vector3 &position)
{
    return this->fittingFunction(position);
}

void EnvObject::precomputeTerrainProperties(const Heightmap &heightmap)
{
    Vector3 terrainDimensions = Vector3(heightmap.getDimensions().x, heightmap.getDimensions().y, 1);
    GridF initialScalarPropertyMap(terrainDimensions, 0.f);
    GridV3 initialVectorPropertyMap(terrainDimensions, Vector3::invalid());

    // Initialize the maps
    for (auto& [name, obj] : EnvObject::availableObjects) {
        EnvObject::allVectorProperties[name] = initialVectorPropertyMap;
        EnvObject::allVectorProperties[name + ".center"] = initialVectorPropertyMap;
        EnvObject::allVectorProperties[name + ".start"] = initialVectorPropertyMap;
        EnvObject::allVectorProperties[name + ".end"] = initialVectorPropertyMap;
        EnvObject::allVectorProperties[name + ".normal"] = initialVectorPropertyMap;
        EnvObject::allVectorProperties[name + ".dir"] = initialVectorPropertyMap;
    }
    EnvObject::allVectorProperties["current.dir"] = initialVectorPropertyMap;
    EnvObject::allScalarProperties["current.vel"] = initialScalarPropertyMap;


    // Evaluate at each point
    initialScalarPropertyMap.iterateParallel([&](const Vector3& pos) {
        for (auto& [name, obj] : EnvObject::availableObjects) {
            auto [distance, object] = EnvObject::getSqrDistanceTo(name, pos);
            if (object == nullptr) continue;

            EnvObject::allVectorProperties[name](pos) = object->getProperty(pos, "default");
            EnvObject::allVectorProperties[name + ".center"](pos) = object->getProperty(pos, "center");
            EnvObject::allVectorProperties[name + ".start"](pos) = object->getProperty(pos, "start");
            EnvObject::allVectorProperties[name + ".end"](pos) = object->getProperty(pos, "end");
            EnvObject::allVectorProperties[name + ".normal"](pos) = object->getNormal(pos);
            EnvObject::allVectorProperties[name + ".dir"](pos) = object->getDirection(pos);
        }
        Vector3 waterFlow = EnvObject::flowfield(pos);
        EnvObject::allVectorProperties["current.dir"](pos) = waterFlow.normalized();
        EnvObject::allScalarProperties["current.vel"](pos) = waterFlow.length();
    });
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
    }
    return this->position; // Default
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

ImplicitPatch* EnvPoint::createImplicitPatch(ImplicitPrimitive *previousPrimitive)
{

    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, Vector3(radius, radius, radius * growingState), 0);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, Vector3(radius, radius, radius * growingState), 0);
    }

    patch->position = this->position;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

EnvObject &EnvPoint::translate(const Vector3 &translation)
{
    this->position += translation;
    return *this;
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
    }
    return this->curve.estimateClosestPos(position); // Default
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

ImplicitPatch* EnvCurve::createImplicitPatch(ImplicitPrimitive *previousPrimitive)
{
    BSpline translatedCurve = this->curve;
    AABBox box(this->curve.points);
    Vector3 offset(this->width, this->width, this->width);
    translatedCurve.translate(-(box.min() - offset * .5f));

    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->width * growingState, translatedCurve);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->width * growingState, translatedCurve);
    }
    patch->position = box.min() - offset.xy() * .5f;
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

EnvObject &EnvCurve::translate(const Vector3 &translation)
{
    this->curve.translate(translation);
    return *this;
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
    else if (complement == "start") // Yeah, that doesn't make sense...
        return 0.f;
    else if (complement == "end") // Yeah, that doesn't make sense...
        return 0.f;

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

Vector3 EnvArea::getNormal(const Vector3 &position)
{
    return this->area.getNormal(this->area.estimateClosestTime(position));
}

Vector3 EnvArea::getDirection(const Vector3 &position)
{
    return Vector3::invalid();
}

Vector3 EnvArea::getProperty(const Vector3& position, std::string prop) const
{
    if (prop == "center") {
        return this->area.center();
    } else if (prop == "start") {
        return Vector3::invalid();
    } else if (prop == "end") {
        return Vector3::invalid();
    }
    return this->area.estimateClosestPos(position); // Default
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
    EnvObject::sandDeposit.add(sand.meanSmooth(width, width, 1), box.min() - Vector3(width, width));
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

    float timePrepare = 0.f;
    float timeApply = 0.f;

    flow.iterate([&] (const Vector3& pos) {
        if (!box.contains(pos) || !translatedCurve.contains(pos))
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
            impact = this->flowEffect * distFactor;
            // Ignore the normal for now, I can't find any logic with it...
    //            if (impact.y != 0) { // Add an impact on the normal direction: need to change the normal to be in the right side of the curve
    //                if ((translatedCurve.getPoint(closestTime) - Vector3(x, y, 0)).dot(direction) > 0)
    //                    normal *= -1.f;
    //            }
        });
        timeApply += timeIt([&]() {
            flow.at(pos) = impact.changedBasis(direction, normal, binormal);// + impact * previousFlow;
            occupancy.at(pos) = 1.f;
        });
//        std::cout << "Prepare: " <<showTime(timePrepare) << " Apply: " << showTime(timeApply) << std::endl;
    });
    return {flow, occupancy};
}

ImplicitPatch* EnvArea::createImplicitPatch(ImplicitPrimitive* previousPrimitive)
{
    BSpline translatedCurve = this->area;
    AABBox box(this->area.points);
    Vector3 offset(this->width, this->width, this->height * growingState);
    translatedCurve.translate(-(box.min() - offset * .5f));
    ImplicitPrimitive* patch;
    if (previousPrimitive != nullptr) {
        patch = previousPrimitive;
        *previousPrimitive = *ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->width * growingState, translatedCurve);
    } else {
        patch = ImplicitPatch::createPredefinedShape(this->implicitShape, box.dimensions() + offset, this->width * growingState, translatedCurve);
    }
    patch->position = box.min() - offset.xy() * .5f;
//    patch->position = box.center();
    patch->material = this->material;
    patch->name = this->name;
    return patch;
}

EnvObject &EnvArea::translate(const Vector3 &translation)
{
    this->area.translate(translation);
    return *this;
}
