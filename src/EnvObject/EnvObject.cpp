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
//GridF EnvObject::sandDeposit;
//GridF EnvObject::polypDeposit;
//std::map<std::string, GridF> EnvObject::materialDeposit;
std::map<std::string, EnvMaterial> EnvObject::materials;
std::map<std::string, EnvObject*> EnvObject::availableObjects;
std::vector<EnvObject*> EnvObject::instantiatedObjects;
float EnvObject::flowImpactFactor = .9f;
int EnvObject::currentMaxID = -1;
std::vector<MaterialsTransformation> EnvObject::transformationRules;

std::map<std::string, GridV3> EnvObject::allVectorProperties;
std::map<std::string, GridF> EnvObject::allScalarProperties;

GridV3 initFlow() {
    if (EnvObject::flowfield.empty()) {
        EnvObject::initialFlowfield = GridV3(100, 100, 1, Vector3(0, 0, 0));
        EnvObject::initialFlowfield.raiseErrorOnBadCoord = false;
        EnvObject::initialFlowfield.returned_value_on_outside = RETURN_VALUE_ON_OUTSIDE::REPEAT_VALUE;
        EnvObject::flowfield = EnvObject::initialFlowfield;
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

void EnvObject::readEnvObjectsFile(std::string filename)
{

    /*EnvMaterial sand("sand", 1.f, 2.f, 2.f, EnvObject::flowfield.getDimensions());
    EnvMaterial polyp("polyp", 1.f, .1f, 0.f, EnvObject::flowfield.getDimensions());
    EnvMaterial pebbles("pebbles", .5f, 0.5f, 10.f, EnvObject::flowfield.getDimensions());
    EnvObject::materials["sand"] = sand;
    EnvObject::materials["polyp"] = polyp;
    EnvObject::materials["pebbles"] = pebbles;
    */
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EnvObject::readEnvObjectsFileContent(content);
}

void EnvObject::readEnvObjectsFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));
    for (auto& obj : json) {
        std::string objName = obj["name"];
        if (startsWith(objName, "--")) continue; // Ignore some objects if the name starts with "--"
        if (!obj.contains("type")) {
            throw std::domain_error("No type given for Environmental Object defined as " + nlohmann::to_string(obj));
        }

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
//        obj->sandEffect = EnvObject::availableObjects[name]->sandEffect;
        obj->materialAbsorptionRate = EnvObject::availableObjects[name]->materialAbsorptionRate;
        obj->materialDepositionRate = EnvObject::availableObjects[name]->materialDepositionRate;
    }
    //    precomputeTerrainProperties(Heightmap());
}

void EnvObject::readEnvMaterialsFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EnvObject::readEnvMaterialsFileContent(content);
}

void EnvObject::readEnvMaterialsFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));
    for (auto& mat : json) {
        std::string matName = mat["name"];
        if (startsWith(matName, "--")) continue; // Ignore some materials if the name starts with "--"

        float diffusionSpeed = mat["diffusionspeed"];
        float waterTransport = mat["watertransport"];
        float mass = mat["mass"];
        float decay = mat["decay"];

        EnvMaterial material;
        material.name = matName;
        material.diffusionSpeed = diffusionSpeed;
        material.waterTransport = waterTransport;
        material.mass = mass;
        material.decay = decay;

        if (EnvObject::materials.count(matName) != 0) {
            material.currentState = EnvObject::materials[matName].currentState;
        } else {
            material.currentState = GridF(EnvObject::flowfield.getDimensions(), 0.f);
        }

        EnvObject::materials[matName] = material;
    }
}

void EnvObject::readEnvMaterialsTransformationsFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EnvObject::readEnvMaterialsTransformationsFileContent(content);
}

void EnvObject::readEnvMaterialsTransformationsFileContent(std::string content)
{
    std::vector<MaterialsTransformation> rules;
//    std::string sline;
    auto lines = split(content, "\n");
//    while (std::getline(content, sline)) {
    for (std::string sline : lines) {
        if (sline.empty() || sline[0] == '#') continue; // Comments with "#"
        std::istringstream line(sline);
        std::map<std::string, float> inputs, outputs;
        std::string value, word, operation;
        bool transformationValid = true;
        // Get inputs
        while (operation != "=") {
            line >> value;
            line >> word;
            line >> operation;
            inputs[word] = std::stof(value);
            transformationValid &= (EnvObject::materials.count(word) != 0);
        }
        // Get outputs
        while (true) {
            line >> value;
            line >> word;
            outputs[word] = std::stof(value);
            transformationValid &= (EnvObject::materials.count(word) != 0);
            if (!(line >> operation)) break;
        }
        if (transformationValid)
            rules.push_back({inputs, outputs});
        else {
            std::cerr << "Transformation not valid : " << sline << std::endl;
        }
    }
    EnvObject::transformationRules = rules;
}

EnvObject *EnvObject::fromJSON(nlohmann::json content)
{
    EnvObject* obj = nullptr;
    std::string objName = content["name"];
    std::string objType = content["type"];
    std::map<std::string, float> materialDepositionRate;
    std::map<std::string, float> materialAbsorptionRate;
    TerrainTypes material = materialFromString(content["material"]);
    ImplicitPatch::PredefinedShapes shape = predefinedShapeFromString(content["geometry"]);
    Vector3 dimensions = json_to_vec3(content["dimensions"]);

    Vector3 flowEffect;
    if (objType == "point") {
        auto asPoint = new EnvPoint;
        asPoint->radius = dimensions.x;
        asPoint->height = dimensions.z;
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

    if (content.contains("depositionrate")) {
        auto deposits = content["depositionrate"].get<std::map<std::string, float>>();
        for (auto& [mat, val] : deposits) {
            materialDepositionRate[mat] = val;
            materialAbsorptionRate[mat] = val;
        }
    }
    if (content.contains("absorptionrate")) {
        auto absorbs = content["absorptionrate"].get<std::map<std::string, float>>();
        for (auto& [mat, val] : absorbs) {
            materialAbsorptionRate[mat] = val;
        }
    }
    if(content.contains("sand")) {
        materialDepositionRate["sand"] = float(content["sand"]);
        materialAbsorptionRate["sand"] = float(content["sand"]);
    }
    if(content.contains("polyp")) {
        materialDepositionRate["polyp"] = float(content["polyp"]);
        materialAbsorptionRate["polyp"] = float(content["polyp"]);
    }
    if(content.contains("pebbles")) {
        materialDepositionRate["pebbles"] = float(content["pebbles"]);
        materialAbsorptionRate["pebbles"] = float(content["pebbles"]);
    }


    if (content.contains("needs")) {
        obj->needsForGrowth = content["needs"].get<std::map<std::string, float>>();
        for (auto [key, val] : obj->needsForGrowth)
            obj->currentSatisfaction[key] = 0;
    }

    obj->name = toLower(objName);
    obj->flowEffect = flowEffect;
    obj->materialAbsorptionRate = materialAbsorptionRate;
    obj->materialDepositionRate = materialDepositionRate;
//    obj->materialAbsorptionRate["sand"] = sandEffect;
//    obj->materialDepositionRate["sand"] = sandEffect;
//    obj->materialAbsorptionRate["polyp"] = polypEffect;
//    obj->materialDepositionRate["polyp"] = polypEffect;
//    obj->sandEffect = sandEffect;
//    obj->polypEffect = polypEffect;
    obj->s_FittingFunction = content["rule"];
    obj->material = material;
    obj->implicitShape = shape;
    obj->inputDimensions = dimensions;
    if (dimensions.z == 0) obj->inputDimensions = Vector3();
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
//    variables["sand"] = float();
//    variables["polyp"] = float();
    variables["depth"] = float();
    for (auto& [matName, material] : EnvObject::materials)
        variables[matName] = float();

    ExpressionParser parser;
    variables["pos"] = Vector3();
    if (!parser.validate(formula, variables, false)) {
        throw std::runtime_error("The formula " + formula + " is not valid for object '" + currentObject + "'");
    }
    std::set<std::string> neededVariables = parser.extractAllVariables(formula);
    auto _func = parser.parse(formula, variables);
    return [&, _func, neededVariables, currentObject](Vector3 pos) -> float {
        std::map<std::string, Variable> vars;
        for (auto& [prop, map] : EnvObject::allVectorProperties) {
            /*if (neededVariables.count(prop) && !map(pos).isValid()) { // A variable is needed but there are no value attributed (eg : object not instantiated yet)
                if (prop == currentObject || startsWith(prop, currentObject + ".")) {
                    std::cout << prop << " -> " << currentObject << std::endl;
                    vars[prop] = pos;
                } else {
                    std::cout << "Missing data for " << prop << std::endl;
                    return std::numeric_limits<float>::max(); // Return max value
                }
            } else {
                vars[prop] = map(pos);
            }*/
            vars[prop] = map(pos);
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

bool EnvObject::applyEffects(const GridF& heights)
{
    EnvObject::updateFlowfield();
    return EnvObject::updateSedimentation(heights);
//    EnvObject::applyMaterialsTransformations();
}

bool EnvObject::updateSedimentation(const GridF& heights)
{
    bool bigChangesInAtLeastOneMaterialDistribution = false;
    GridV3 heightsGradients = heights.gradient();
    auto smoothFluids = EnvObject::flowfield.meanSmooth(3, 3, 1, true);
    for (auto& [name, material] : EnvObject::materials) {
        float startingAmount = material.currentState.sum();
        std::vector<std::pair<float, float>> depoAbso(EnvObject::instantiatedObjects.size());
        if (material.diffusionSpeed < 1.f) {
            if (random_gen::generate() < material.diffusionSpeed) {
                material.currentState = material.currentState.meanSmooth(3, 3, 1, false); // Diffuse a little bit
            }
        } else {
            material.currentState = material.currentState.meanSmooth(material.diffusionSpeed * 2 + 1, material.diffusionSpeed * 2 + 1, 1, false); // Diffuse
        }
        for (size_t i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
            float start = material.currentState.sum();
            auto& object = EnvObject::instantiatedObjects[i];
            object->applyAbsorption(material);
            depoAbso[i].second = material.currentState.sum() - start;
        }
        material.currentState = material.currentState.warpWith((smoothFluids * material.waterTransport) - (heightsGradients * material.mass));

        for (size_t i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
            auto& object = EnvObject::instantiatedObjects[i];
            float start = material.currentState.sum();
            object->applyDeposition(material);
            depoAbso[i].first = material.currentState.sum() - start;
        }

        material.currentState.iterateParallel([&](size_t i) {
//            material.currentState[i] = std::clamp(material.currentState[i], 0.f, 1.f); // Limited between 0 and 1 ?
        });
        material.currentState *= material.decay;

        float endingAmount = material.currentState.sum();

        if (std::abs(endingAmount - startingAmount) > 1e-3) bigChangesInAtLeastOneMaterialDistribution = true;
    }
    return bigChangesInAtLeastOneMaterialDistribution;
}

void EnvObject::applyMaterialsTransformations()
{
    displayProcessTime("Filling compact materials... ", [&]() {
        std::set<std::string> neededMaterials;
        for (size_t iRule = 0; iRule < transformationRules.size(); iRule++) {
            auto [input, output] = transformationRules[iRule];
            for (auto [inMaterial, inDose] : input) {
                neededMaterials.insert(inMaterial);
            }
            for (auto [outMaterial, outDose] : output) {
                neededMaterials.insert(outMaterial);
            }
        }
        std::map<std::string, float> initialState; // Loop the map creation only once
        for (const auto& matName : neededMaterials)
            initialState.insert({matName, 0.f});
        Matrix3<std::map<std::string, float>> allMaterials(EnvObject::flowfield.getDimensions(), initialState);
        allMaterials.iterateParallel([&] (size_t i) {
            for (const auto& [matName, amount] : allMaterials[i]) {
                allMaterials[i][matName] = EnvObject::materials[matName].currentState[i];
            }

            for (size_t iRule = 0; iRule < transformationRules.size(); iRule++) {
                const auto& [input, output] = transformationRules[iRule];
                float maxTransform = 10000.f;
                for (const auto& [inMaterial, inDose] : input) {
                    float inAmount = allMaterials[i][inMaterial];
                    float transformVal = inAmount / inDose;
                    maxTransform = std::min(maxTransform, transformVal);
                }
                if (maxTransform > 1e-3) {
                    for (const auto& [inMaterial, inDose] : input) {
                        allMaterials[i][inMaterial] -= inDose * maxTransform;
                    }
                    for (const auto& [outMaterial, outDose] : output) {
                        allMaterials[i][outMaterial] += outDose * maxTransform;
                    }
                }
            }
        });

        for (auto& matName : neededMaterials) {
            auto& mat = EnvObject::materials[matName];
            mat.currentState.iterateParallel([&](size_t i) {
                mat.currentState[i] = allMaterials[i][matName];
            });
        }
    }, false);
}

void EnvObject::updateFlowfield()
{
    EnvObject::flowfield = EnvObject::initialFlowfield;
//    for (int i = int(EnvObject::instantiatedObjects.size()) - 1; i >= 0; i--) {
    for (int i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
        auto& object = EnvObject::instantiatedObjects[i];
        auto [flow, occupancy] = object->computeFlowModification();
//        float currentGrowth = object->computeGrowingState();
//        occupancy *= currentGrowth;

        EnvObject::flowfield = flow;
//        EnvObject::flowfield += flow * EnvObject::flowImpactFactor;
//        flow = flow * occupancy + EnvObject::flowfield * (1.f - occupancy);
//        EnvObject::flowfield = EnvObject::flowfield * (1.f - EnvObject::flowImpactFactor) + flow * EnvObject::flowImpactFactor;
    }
    EnvObject::flowfield = EnvObject::flowfield.meanSmooth(3, 3, 1, true);
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
    return this->fittingFunction(position.xy());
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
//        EnvObject::allScalarProperties["sand"] = initialScalarPropertyMap;
//        EnvObject::allScalarProperties["polyp"] = initialScalarPropertyMap;
        EnvObject::allScalarProperties["depth"] = initialScalarPropertyMap;
        for (const auto& [matName, material] : EnvObject::materials) {
            EnvObject::allScalarProperties[matName] = initialScalarPropertyMap;
        }


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
        EnvObject::allVectorProperties[name](pos) = allProperties["default"];
        EnvObject::allVectorProperties[name + ".center"](pos) = allProperties["center"];
        EnvObject::allVectorProperties[name + ".start"](pos) = allProperties["start"];
        EnvObject::allVectorProperties[name + ".end"](pos) = allProperties["end"];
        EnvObject::allScalarProperties[name + ".inside"](pos) = (allProperties["inside"].isValid() ? 1.f : 0.f);
        EnvObject::allVectorProperties[name + ".normal"](pos) = allProperties["normal"];
        EnvObject::allVectorProperties[name + ".dir"](pos) = allProperties["dir"];
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
    for (auto& [matName, material] : EnvObject::materials) {
        EnvObject::allScalarProperties[matName] = material.currentState;
    }
    EnvObject::allScalarProperties["depth"] = (heightmap.properties->waterLevel * heightmap.getSizeZ()) - heightmap.heights;
}

void EnvObject::reset()
{
    for (auto& obj : EnvObject::instantiatedObjects)
        delete obj;
    EnvObject::instantiatedObjects.clear();
//    for (auto& [name, obj] : EnvObject::availableObjects)
//        delete obj;
//    EnvObject::availableObjects.clear();
    initFlow();
    for (auto& [matName, mat] : materials) {
        mat.currentState.reset();
    }

}

#include "Utils/Delaunay.h"
GraphObj EnvObject::sceneToGraph()
{
    GraphObj graph;
    std::vector<GraphNodeObj*> nodes(EnvObject::instantiatedObjects.size());
    std::vector<Vector3> positions(nodes.size());

    for (int i = 0; i < nodes.size(); i++) {
        auto& obj = EnvObject::instantiatedObjects[i];
        positions[i] = dynamic_cast<EnvPoint*>(obj)->position;
//        nodes[i] = graph.addNode(new GraphNodeObj(obj, positions[i], i));
    }

    graph = Delaunay().fromVoronoi(Voronoi(positions)).graph.cast<EnvObject*>();

    return graph;
}
