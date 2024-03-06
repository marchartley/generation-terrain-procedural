#include "EnvObject.h"
#include "ExpressionParser.h"

#include "Graphics/DisplayGraphics.h"
#include <fstream>

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

GridV3 initFlow();

GridV3 EnvObject::flowfield = initFlow();
GridV3 EnvObject::initialFlowfield;
GridV3 EnvObject::terrainNormals;
GridF EnvObject::sandDeposit;
std::map<std::string, EnvObject*> EnvObject::availableObjects;
std::vector<EnvObject*> EnvObject::instantiatedObjects;
float EnvObject::flowImpactFactor = .5f;
int EnvObject::currentMaxID = -1;

std::map<std::string, GridV3> EnvObject::allVectorProperties;
std::map<std::string, GridF> EnvObject::allScalarProperties;

GridV3 initFlow() {
    if (EnvObject::flowfield.empty()) {
        EnvObject::initialFlowfield = GridV3(100, 100, 1, Vector3(1, 0, 0));
        EnvObject::initialFlowfield.raiseErrorOnBadCoord = false;
        EnvObject::initialFlowfield.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        EnvObject::flowfield = EnvObject::initialFlowfield;
        EnvObject::sandDeposit = GridF(EnvObject::flowfield.getDimensions(), 0.f);
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

    for (auto& obj : EnvObject::instantiatedObjects) {
        auto name = obj->name;
        obj->flowEffect = EnvObject::availableObjects[name]->flowEffect;
        obj->sandEffect = EnvObject::availableObjects[name]->sandEffect;
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
        asCurve->height = dimensions.z;
        obj = asCurve;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "area") {
        auto asArea = new EnvArea;
        asArea->width = dimensions.x;
        asArea->height = dimensions.z;
        obj = asArea;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    }


    if (content.contains("needs")) {
        obj->needsForGrowth = content["needs"].get<std::map<std::string, float>>();
        for (auto [key, val] : obj->needsForGrowth)
            obj->currentSatisfaction[key] = 0;
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

float EnvObject::computeGrowingState()
{
    bool verbose = false;
    std::ostringstream oss;

    if (this->needsForGrowth.count("age")) {
        currentSatisfaction["age"] = this->age;
    }

    float totalScore = 1.f; // Start expecting to be "adult"
    for (auto [key, value] : needsForGrowth) {
        if (value == 0) continue;
        float score = std::clamp(currentSatisfaction[key] / needsForGrowth[key], 0.f, 1.f);
        totalScore = std::min(totalScore, score);

        oss << key << " (" << currentSatisfaction[key] << "/" << needsForGrowth[key] << "), ";
    }
    oss << " => total: " << totalScore * 100 << "%";
    if (verbose)
        std::cout << this->name << " state : " << oss.str() << std::endl;
    return totalScore;
}

std::function<float (Vector3)> EnvObject::parseFittingFunction(std::string formula, std::string currentObject)
{
    formula = toLower(formula);
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
        variables[name + ".inside"] = float();
        variables[name + ".curvature"] = float();
    }

    variables["current"] = Vector3();
    variables["current.center"] = Vector3();
    variables["current.start"] = Vector3();
    variables["current.end"] = Vector3();
    variables["current.normal"] = Vector3();
    variables["current.dir"] = Vector3();
    variables["current.vel"] = float();
    variables["sand"] = float();
    variables["depth"] = float();

    ExpressionParser parser;
    variables["pos"] = Vector3();
    if (!parser.validate(formula, variables, false)) {
        throw std::runtime_error("The formula " + formula + " is not valid!!");
    }
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
                    std::cout << "Missing data for " << prop << std::endl;
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
    EnvObject::updateFlowfield();
    EnvObject::updateSedimentation();
}

void EnvObject::updateSedimentation()
{
    EnvObject::sandDeposit = EnvObject::sandDeposit.meanSmooth(5, 5, 1, false); // Diffuse sand
    for (auto& object : EnvObject::instantiatedObjects) {
        object->applySandAbsorption();
        object->applySandDeposit();
    }
    EnvObject::sandDeposit = EnvObject::sandDeposit.warpWith(EnvObject::flowfield.meanSmooth(3, 3, 1, true) * 2.f);
    EnvObject::sandDeposit.iterateParallel([&] (const Vector3& pos) {
        if (pos.x > 5) return;
        EnvObject::sandDeposit(pos) = std::max(.1f, EnvObject::sandDeposit(pos));
    });
}

void EnvObject::updateFlowfield()
{
    EnvObject::flowfield = EnvObject::initialFlowfield;
    GridV3 totalNewFlow(EnvObject::flowfield.getDimensions());
    GridF totalOccupancy(EnvObject::flowfield.getDimensions());
    for (auto& object : EnvObject::instantiatedObjects) {
        auto [flow, occupancy] = object->computeFlowModification();
        totalNewFlow += flow;
        totalOccupancy += occupancy;

        flow = flow * occupancy + EnvObject::flowfield * (1.f - occupancy);
//        flow.iterateParallel([&] (size_t i) {
//            flow[i] = (occupancy[i] != 0.f ? flow[i] / occupancy[i] : EnvObject::flowfield[i]);
//        });
        EnvObject::flowfield = EnvObject::flowfield * (1.f - EnvObject::flowImpactFactor) + flow * EnvObject::flowImpactFactor;
    }
    EnvObject::flowfield = EnvObject::flowfield.meanSmooth(3, 3, 1, true);
    /*totalNewFlow.iterateParallel([&] (size_t i) {
        totalNewFlow[i] = (totalOccupancy[i] != 0.f ? totalNewFlow[i] / totalOccupancy[i] : totalNewFlow[i]);
    });
    EnvObject::flowfield = EnvObject::flowfield * (1.f - EnvObject::flowImpactFactor) + totalNewFlow * EnvObject::flowImpactFactor;
    */
}

void EnvObject::beImpactedByEvents()
{
    // TODO!!!
    for (auto& obj : EnvObject::instantiatedObjects) {
        obj->age += 1.f;
//        obj->growingState = std::min(obj->growingState + .2f, 1.f);
    }
}

float EnvObject::evaluate(const Vector3 &position)
{
    return this->fittingFunction(position);
}

void EnvObject::precomputeTerrainProperties(const Heightmap &heightmap)
{

    displayProcessTime("Computing terrain properties... ", [&]() {
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
            EnvObject::allScalarProperties[name + ".inside"] = initialScalarPropertyMap;
            EnvObject::allScalarProperties[name + ".curvature"] = initialScalarPropertyMap;
        }
        EnvObject::allVectorProperties["current"] = initialVectorPropertyMap;
        EnvObject::allVectorProperties["current.dir"] = initialVectorPropertyMap;
        EnvObject::allScalarProperties["current.vel"] = initialScalarPropertyMap;
        EnvObject::allScalarProperties["sand"] = initialScalarPropertyMap;
        EnvObject::allScalarProperties["depth"] = initialScalarPropertyMap;


        // Evaluate at each point
        for (auto& [name, obj] : EnvObject::availableObjects) {
            displayProcessTime("Computing properties for " + name + "... ", [&]() {
                EnvObject::recomputeTerrainPropertiesForObject(heightmap, name);
            }, false);
        }
        EnvObject::recomputeFlowAndSandProperties(heightmap);
    });
}

void EnvObject::recomputeTerrainPropertiesForObject(const Heightmap &heightmap, std::string objectName)
{
    auto name = objectName;
    EnvObject::flowfield.iterateParallel([&](const Vector3& pos) {
        auto [distance, object] = EnvObject::getSqrDistanceTo(name, pos);
        if (object == nullptr) return;
        auto allProperties = object->getAllProperties(pos);
        EnvObject::allVectorProperties[name](pos) = allProperties["default"]; // object->getProperty(pos, "default");
        EnvObject::allVectorProperties[name + ".center"](pos) = allProperties["center"]; // object->getProperty(pos, "center");
        EnvObject::allVectorProperties[name + ".start"](pos) = allProperties["start"]; // object->getProperty(pos, "start");
        EnvObject::allVectorProperties[name + ".end"](pos) = allProperties["end"]; // object->getProperty(pos, "end");
        EnvObject::allScalarProperties[name + ".inside"](pos) = (allProperties["inside"].isValid() ? 1.f : 0.f); // (object->getProperty(pos, "inside").isValid() ? 1.f : 0.f);
        EnvObject::allVectorProperties[name + ".normal"](pos) = allProperties["normal"]; // object->getNormal(pos);
        EnvObject::allVectorProperties[name + ".dir"](pos) = allProperties["dir"]; // object->getDirection(pos);
        EnvObject::allScalarProperties[name + ".curvature"](pos) = (allProperties["curvature"].x < 1e5 ? allProperties["curvature"].x : -1.f);
    });
}

void EnvObject::recomputeFlowAndSandProperties(const Heightmap &heightmap)
{
    EnvObject::flowfield.iterateParallel([&](const Vector3& pos) {
        Vector3 waterFlow = EnvObject::flowfield(pos);
        EnvObject::allVectorProperties["current"](pos) = waterFlow;
        EnvObject::allVectorProperties["current.dir"](pos) = waterFlow.normalized();
        EnvObject::allScalarProperties["current.vel"](pos) = waterFlow.length();
    });
    EnvObject::allScalarProperties["sand"] = EnvObject::sandDeposit;
    EnvObject::allScalarProperties["depth"] = (heightmap.properties->waterLevel * heightmap.getSizeZ()) - heightmap.heights;
}

void EnvObject::reset()
{
    for (auto& obj : EnvObject::instantiatedObjects)
        delete obj;
    EnvObject::instantiatedObjects.clear();
    for (auto& [name, obj] : EnvObject::availableObjects)
        delete obj;
    EnvObject::availableObjects.clear();
    initFlow();

}
