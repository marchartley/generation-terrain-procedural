#include "EnvObject.h"
#include "ExpressionParser.h"

#include "Graphics/Plotter.h"
#include <fstream>

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

GridV3 initFlow(bool force = false);

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
int EnvObject::currentTime = 0;
Scenario EnvObject::scenario;

std::map<std::string, GridV3> EnvObject::allVectorProperties;
std::map<std::string, GridF> EnvObject::allScalarProperties;

GridV3 initFlow(bool force) {
    if (force || EnvObject::flowfield.empty()) {
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
        obj->fitnessFunction = EnvObject::parseFittingFunction(obj->s_FitnessFunction, obj->name);
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
        float virtualHeight = mat["virtualheight"];

        EnvMaterial material;
        material.name = matName;
        material.diffusionSpeed = diffusionSpeed;
        material.waterTransport = waterTransport;
        material.mass = mass;
        material.decay = decay;
        material.virtualHeight = virtualHeight;

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

void EnvObject::readScenarioFile(std::string filename)
{
    std::ifstream file(filename);
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EnvObject::readScenarioFileContent(content);
}

void EnvObject::readScenarioFileContent(std::string content)
{
    auto json = nlohmann::json::parse(toLower(content));

    scenario = Scenario();

    auto objects = json["objects"].get<std::map<std::string, nlohmann::json>>();
    for (auto [name, obj] : objects) {
        float proba = obj["proba"];
        int amount = (obj.contains("amount") ? obj["amount"].get<int>() : -1);

        scenario.addObject(name, proba, amount);
    }

    auto events = json["events"];
    for (auto& event : events) {
        std::string type = toLower(event["type"]);
        float startTime = event["start"];
        float endTime = event["end"];

        if (type == "waterlevel") {
            float amount = event["amount"];
            scenario.waterLevelEvents.push_back(WaterLevelEvent(amount, startTime, endTime));
        } else if (type == "storm") {
            Vector3 position = json_to_vec3(event["position"]);
            Vector3 direction = json_to_vec3(event["direction"]);
            float sigma = event["sigma"];
            scenario.stormEvents.push_back(StormEvent(position, direction, sigma, startTime, endTime));
        } else if (type == "subsidence") {
            Vector3 position = json_to_vec3(event["position"]);
            float amount = event["amount"];
            float sigma = event["sigma"];
            scenario.subsidenceEvents.push_back(SubsidenceEvent(position, amount, sigma, startTime, endTime));
        } else if (type == "tectonic") {
            Vector3 direction = json_to_vec3(event["direction"]);
            float sigma = event["sigma"];
            scenario.tectonicEvents.push_back(TectonicEvent(direction, sigma, startTime, endTime));
        } else {
            std::cerr << "The event " << type << " is not recognized..." << std::endl;
        }
    }

    auto parameters = json["simulation"];
    float duration = parameters["end"];
    float dt = parameters["dt"];
    float waterLevel = parameters["waterlevel"];

    scenario.duration = duration;
    scenario.dt = dt;
    scenario.waterLevel = waterLevel;

}

EnvObject *EnvObject::fromJSON(nlohmann::json content)
{
    EnvObject* obj = nullptr;
    std::string objName = content["name"];
    std::string objType = content["type"];
    std::map<std::string, float> materialDepositionRate;
    std::map<std::string, float> materialAbsorptionRate;
    std::map<std::string, float> materialDepositionOnDeath;
    TerrainTypes material = materialFromString(content["material"]);
    ImplicitPatch::PredefinedShapes shape = predefinedShapeFromString(content["geometry"]);
    Vector3 dimensions = json_to_vec3(content["dimensions"]);
    HeightmapFrom heightFrom = (!content.contains("heightfrom") || content["heightfrom"] == "surface" ? SURFACE : (content["heightfrom"] == "water" ? WATER : GROUND));
    float minScore = (content.contains("minscore") ? content["minscore"].get<float>() : 0.f);

    Vector3 flowEffect;
    if (objType == "point") {
        auto asPoint = new EnvPoint;
        asPoint->radius = dimensions.x;
        asPoint->height = dimensions.z;
        obj = asPoint;
        if (content["flow"].is_number())
            flowEffect = Vector3(content["flow"], content["flow"], content["flow"]);
        else
            flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "curve") {
        auto asCurve = new EnvCurve;
        asCurve->width = dimensions.x;
        asCurve->length = dimensions.y;
        asCurve->height = dimensions.z;
        if (content.contains("follows")) {
            if (content["follows"] == "isovalue") asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::ISOVALUE;
            else if (content["follows"] == "gradients") asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::GRADIENTS;
            else if (content["follows"] == "skeleton") asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::SKELETON;
            else std::cerr << "Value for 'follow' in object " << objName << " not recognized. Should be 'isovalue', 'gradients' or 'skeleton'. Got " << content["follows"] << std::endl;
        }
        obj = asCurve;
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "area") {
        auto asArea = new EnvArea;
        asArea->width = dimensions.x;
        asArea->length = dimensions.y;
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
    if (content.contains("ondeath")) {
        auto depos = content["ondeath"].get<std::map<std::string, float>>();
        for (auto& [mat, val] : depos) {
            materialDepositionOnDeath[mat] = val;
        }
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
    obj->materialDepositionOnDeath = materialDepositionOnDeath;
    obj->s_FittingFunction = content["rule"];
    if (content.contains("fitness")) {
        obj->s_FitnessFunction = content["fitness"];
    } else {
        obj->s_FitnessFunction = obj->s_FittingFunction;
    }
    obj->material = material;
    obj->implicitShape = shape;
    obj->inputDimensions = dimensions;
    obj->recomputeEvaluationPoints();
    obj->heightFrom = heightFrom;
    obj->minScore = minScore;
    if (dimensions.z == 0) obj->inputDimensions = Vector3();
    return obj;
}

nlohmann::json EnvObject::toJSON() const
{
    nlohmann::json json;

    json["name"] = this->name;
    json["ID"] = this->ID;
    json["age"] = this->age;
    json["needs"] = this->currentSatisfaction;
    // json["evaluationPosition"] = vec3_to_json(this->evaluationPosition);
    /*std::vector<nlohmann::json> positions;
    for (auto& p : this->evaluationPositions)
        positions.push_back(vec3_to_json(p));
    json["evaluationPositions"] = positions;*/
    json["fitnessScoreAtCreation"] = this->fitnessScoreAtCreation;

    return json;
}

float EnvObject::computeGrowingState()
{
    if (this->createdManually) return 1.f;
    return 1.f; // Yep, let's say that it is always mature....

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

float EnvObject::computeGrowingState2()
{
    if (this->createdManually) return 1.f;
    return (this->evaluate() > 0.f ? 1.f : 0.f); // Let's say that it is either dead or alive...

    // float newFitnessEvaluation = this->evaluate(this->evaluationPosition);
    float newFitnessEvaluation = this->evaluate();
    // std::cout << this->fitnessScoreAtCreation << " / " << newFitnessEvaluation << std::endl;
    // if (newFitnessEvaluation <= 0) return 0;
    if (this->fitnessScoreAtCreation <= 1e-6) return 0; // This is a problem...
    return std::clamp(newFitnessEvaluation / this->fitnessScoreAtCreation, 0.f, 1.f);
}

GridF EnvObject::createHeightfield()
{
    if (_patch) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(_patch)) {
            if (_cachedHeightfield.empty() && patch && patch->predefinedShape != ImplicitPatch::NONE) {
                auto previousMaterial = patch->material;
                patch->material = SAND; // Temporarly be solid to get some height (which is then depth...)
                GridF heights(patch->getDimensions().x, patch->getDimensions().y, 1);
                heights.iterateParallel([&] (const Vector3& p) {
                    float resolution = .5f;
                    for (float z = 1; z < patch->getDimensions().z * 1.f; z += resolution) {
                        heights(p) += (patch->evaluate(p.xy() + patch->position.xy() + Vector3(0, 0, z)) > 0 ? resolution : 0.f);
                    }
                });
                patch->material = previousMaterial;
                _cachedHeightfield = heights;
            }
        }
    }
    return _cachedHeightfield;
}

std::function<float (Vector3)> EnvObject::parseFittingFunction(std::string formula, std::string currentObject, bool removeSelfInstances, EnvObject *myObject)
{
    formula = toLower(formula);
    if (formula == "")
        return [](Vector3 _) { return 0.f; };

    std::map<std::string, Variable> variables;
    for (auto& [name, obj] : EnvObject::availableObjects) {
        variables[name] = Vector3::invalid();
        variables[name + ".center"] = Vector3::invalid();
        variables[name + ".start"] = Vector3::invalid();
        variables[name + ".end"] = Vector3::invalid();
        variables[name + ".normal"] = Vector3::invalid();
        variables[name + ".dir"] = Vector3::invalid();
        variables[name + ".inside"] = float();
        variables[name + ".curvature"] = float();
    }

    variables["current"] = Vector3::invalid();
    variables["current.center"] = Vector3::invalid();
    variables["current.start"] = Vector3::invalid();
    variables["current.end"] = Vector3::invalid();
    variables["current.normal"] = Vector3::invalid();
    variables["current.dir"] = Vector3::invalid();
    variables["current.vel"] = float();
    variables["current.gradient"] = Vector3::invalid();
    variables["depth"] = float();
    variables["depth.gradient"] = Vector3::invalid();
    variables["fracture"] = float();
    variables["fracture.gradient"] = Vector3::invalid();
    for (auto& [matName, material] : EnvObject::materials) {
        variables[matName] = float();
        variables[matName + ".gradient"] = Vector3::invalid();
    }

    ExpressionParser parser;
    variables["currenttime"] = float();
    variables["spawntime"] = float();
    variables["pos"] = Vector3::invalid();
    if (!parser.validate(formula, variables, false)) {
        throw std::runtime_error("The formula " + formula + " is not valid for object '" + currentObject + "'");
    }
    std::set<std::string> neededVariables = parser.extractAllVariables(formula);
    auto _func = parser.parse(formula, variables);
    return [&, formula, _func, neededVariables, currentObject, removeSelfInstances, myObject](Vector3 pos) -> float {
        // ExpressionParser parser;
        std::map<std::string, Variable> vars;
        for (auto& [prop, map] : EnvObject::allVectorProperties) {
            if (removeSelfInstances && (startsWith(prop, currentObject + ".") || startsWith(prop, currentObject))) {
                vars[prop] = Vector3::invalid();
            } else {
                vars[prop] = map(pos);
            }
        }
        for (auto& [prop, map] : EnvObject::allScalarProperties) {
            if (removeSelfInstances && (startsWith(prop, currentObject + ".") || startsWith(prop, currentObject))) {
                vars[prop] = float();
            } else {
                vars[prop] = map(pos);
            }
        }
        if (myObject && myObject->_patch) {
            // if (auto patch = dynamic_cast<ImplicitPrimitive*>(myObject->_patch)) {
            vars["depth"] = std::get<float>(vars["depth"]) + (myObject->createHeightfield().at(pos)); //.at(pos - patch->position.xy()));
            // } else {
                // std::cout << "NO" << std::endl;
            // }
        }
        vars["pos"] = pos;
        vars["currenttime"] = float(EnvObject::currentTime);
        if (myObject != nullptr)
            vars["spawntime"] = float(std::min(EnvObject::currentTime, myObject->spawnTime));
        else
            vars["spawntime"] = float(EnvObject::currentTime);
        /*
        for (std::string neededVar : neededVariables) {
            if (std::holds_alternative<float>(vars[neededVar])) {
                std::cout << neededVar << ": " << std::get<float>(vars[neededVar]) << std::endl;
            } else if (std::holds_alternative<Vector3>(vars[neededVar])) {
                std::cout << neededVar << ": " << std::get<Vector3>(vars[neededVar]) << std::endl;
            } else {
                std::cout << neededVar << ": unknown" << std::endl;
            }
        }*/
        float score = _func(vars);

        /*
        std::cout << "Values used for " << currentObject << ":\n";
        for (const std::string& usedVar : neededVariables) {
            if (std::holds_alternative<float>(vars[usedVar]))
                std::cout << usedVar << " = " << std::get<float>(vars[usedVar]) << std::endl;
            else
                std::cout << usedVar << " = " << std::get<Vector3>(vars[usedVar]) << std::endl;
        }*/
        return score;
    };
}

EnvObject *EnvObject::findClosest(std::string objectName, const Vector3 &pos)
{
    float minDist = std::numeric_limits<float>::max();
    EnvObject* bestElem = nullptr;
    for (auto& instance : EnvObject::instantiatedObjects) {
        if (instance->name != objectName) continue;
        float distance = instance->getSqrDistance(pos);
        if (distance < minDist) {
            minDist = distance;
            bestElem = instance;
        }
    }
    return bestElem;
}
/*
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
*/
EnvObject *EnvObject::instantiate(std::string objectName)
{
    if (EnvObject::availableObjects.count(objectName) == 0) {
        return nullptr;
    }
    EnvObject::currentMaxID++;
    auto object = EnvObject::availableObjects[objectName]->clone();
    object->ID = EnvObject::currentMaxID;
    EnvObject::instantiatedObjects.push_back(object);
    return object;
}

void EnvObject::removeObject(EnvObject *obj)
{
    if (obj) {
        auto& list = EnvObject::instantiatedObjects;
        list.erase(std::find(list.begin(), list.end(), obj));
    }
}

void EnvObject::removeAllObjects()
{
    for (auto& object : EnvObject::instantiatedObjects) {
        delete object;
    }
    EnvObject::instantiatedObjects.clear();
}

bool EnvObject::applyEffects(const GridF& heights, const GridV3& userFlow)
{
    EnvObject::updateFlowfield(userFlow);
    return EnvObject::updateSedimentation(heights);
//    EnvObject::applyMaterialsTransformations();
    // return false;
}

bool EnvObject::updateSedimentation(const GridF& heights)
{
    bool bigChangesInAtLeastOneMaterialDistribution = false;
    GridV3 heightsGradients = heights.gradient();
    auto smoothFluids = EnvObject::flowfield.meanSmooth(3, 3, 1, true);

    std::vector<std::string> names;
    for (auto& [name, material] : EnvObject::materials) {
        names.push_back(name);
    }
    #pragma omp parallel for
    for(int i = 0; i < names.size(); i++) {
        auto& material = materials[names[i]];
        if (material.isStable) {
            // std::cout << material.name << " is stable" << std::endl;
            // continue;
        } else {
            // std::cout << material.name << " NOT stable" << std::endl;
        }

        float startingAmount = material.currentState.sum();
        for (size_t i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
            auto& object = EnvObject::instantiatedObjects[i];
            // if (object->materialAbsorptionRate[material.name] != 0) {
                #pragma omp critical
                object->applyAbsorption(material);
            // }
        }

        for (size_t i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
            auto& object = EnvObject::instantiatedObjects[i];
            // if (object->materialDepositionRate[material.name] != 0) {
                #pragma omp critical
                object->applyDeposition(material);
            // }
        }

        material.update(smoothFluids, heightsGradients, EnvObject::scenario.dt);
        // material.currentState *= material.decay;

        float endingAmount = material.currentState.sum();

        if (std::abs(endingAmount - startingAmount) > 1e-3) {
            #pragma omp critical
            bigChangesInAtLeastOneMaterialDistribution = true;
        } else {
            // std::cout << material.name << " diff : " << std::abs(endingAmount - startingAmount) << std::endl;
            material.isStable = true;
        }
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

void EnvObject::updateFlowfield(const GridV3 &userFlow)
{
    EnvObject::flowfield = EnvObject::initialFlowfield;
    if (!userFlow.empty())
        EnvObject::flowfield += userFlow;
    for (int i = 0; i < EnvObject::instantiatedObjects.size(); i++) {
        auto& object = EnvObject::instantiatedObjects[i];
        auto [flow, occupancy] = object->computeFlowModification();
        EnvObject::flowfield = flow;
    }
    EnvObject::flowfield = EnvObject::flowfield.meanSmooth(3, 3, 1, true);
}

void EnvObject::beImpactedByEvents()
{
    for (auto& obj : EnvObject::instantiatedObjects) {
        obj->age += 1.f;
    }
}

float EnvObject::evaluate(const Vector3 &position)
{
    return this->fitnessFunction(position.xy());
}

float EnvObject::evaluate()
{
    // Should only be at one point...

    if (evaluationPositions.empty()) {
        std::cerr << "Object " << name << " has no evaluation point..." << std::endl;
        return 0;
    }
    float totalScore = 0;
    for (auto& p : evaluationPositions) {
        float score = evaluate(p);
        if (score != score) {
            std::cerr << "NaN found while evaluating " << this->name << " at " << p << std::endl;
        } else {
            totalScore += score;
        }
    }
    return totalScore / float(evaluationPositions.size());
}

void EnvObject::die()
{
    this->applyDepositionOnDeath();
}

void EnvObject::precomputeTerrainProperties(const GridF& heightmap, float waterLevel, float maxHeight)
{

    displayProcessTime("Computing terrain properties... ", [&]() {
        Vector3 terrainDimensions = heightmap.getDimensions();
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
        EnvObject::allVectorProperties["current.gradient"] = initialVectorPropertyMap;

        EnvObject::allScalarProperties["depth"] = initialScalarPropertyMap;
        EnvObject::allVectorProperties["depth.gradient"] = initialVectorPropertyMap;

        EnvObject::allScalarProperties["fracture"] = initialScalarPropertyMap;
        EnvObject::allVectorProperties["fracture.gradient"] = initialVectorPropertyMap;

        for (const auto& [matName, material] : EnvObject::materials) {
            EnvObject::allScalarProperties[matName] = initialScalarPropertyMap;
            EnvObject::allVectorProperties[matName + ".gradient"] = initialVectorPropertyMap;
        }


        // Evaluate at each point
        for (auto& [name, obj] : EnvObject::availableObjects) {
            displayProcessTime("Computing properties for " + name + "... ", [&]() {
                EnvObject::recomputeTerrainPropertiesForObject(name);
            }, false);
        }
        EnvObject::recomputeFlowAndSandProperties(heightmap, waterLevel, maxHeight);
    });
}

void EnvObject::recomputeTerrainPropertiesForObject(std::string objectName)
{
    auto name = objectName;
    EnvObject::flowfield.iterateParallel([&](const Vector3& pos) {
//        auto [distance, object] = EnvObject::getSqrDistanceTo(name, pos);
        EnvObject* object = EnvObject::findClosest(objectName, pos);
        if (object == nullptr) {
            EnvObject::allVectorProperties[name](pos) = Vector3::invalid();
            EnvObject::allVectorProperties[name + ".center"](pos) = Vector3::invalid();
            EnvObject::allVectorProperties[name + ".start"](pos) = Vector3::invalid();
            EnvObject::allVectorProperties[name + ".end"](pos) = Vector3::invalid();
            EnvObject::allScalarProperties[name + ".inside"](pos) = 0.f;
            EnvObject::allVectorProperties[name + ".normal"](pos) = Vector3::invalid();
            EnvObject::allVectorProperties[name + ".dir"](pos) = Vector3::invalid();
            EnvObject::allScalarProperties[name + ".curvature"](pos) = 0.f;
        } else {
            auto allProperties = object->getAllProperties(pos);
            EnvObject::allVectorProperties[name](pos) = allProperties["default"];
            EnvObject::allVectorProperties[name + ".center"](pos) = allProperties["center"];
            EnvObject::allVectorProperties[name + ".start"](pos) = allProperties["start"];
            EnvObject::allVectorProperties[name + ".end"](pos) = allProperties["end"];
            EnvObject::allScalarProperties[name + ".inside"](pos) = (allProperties["inside"].isValid() ? 1.f : 0.f);
            EnvObject::allVectorProperties[name + ".normal"](pos) = allProperties["normal"];
            EnvObject::allVectorProperties[name + ".dir"](pos) = allProperties["dir"];
            EnvObject::allScalarProperties[name + ".curvature"](pos) = (allProperties["curvature"].x < 1e5 ? allProperties["curvature"].x : -1.f);
        }
    });
}

void EnvObject::recomputeFlowAndSandProperties(const GridF& heightmap, float waterLevel, float maxHeight)
{
    EnvObject::flowfield.iterateParallel([&](const Vector3& pos) {
        Vector3 waterFlow = EnvObject::flowfield(pos);
        EnvObject::allVectorProperties["current"](pos) = waterFlow;
        EnvObject::allVectorProperties["current.dir"](pos) = waterFlow.normalized();
        EnvObject::allScalarProperties["current.vel"](pos) = waterFlow.length();
    });
    EnvObject::allVectorProperties["current.gradient"] = EnvObject::allScalarProperties["current.vel"].gradient();
    for (auto& [matName, material] : EnvObject::materials) {
        EnvObject::allScalarProperties[matName] = material.currentState;
        EnvObject::allVectorProperties[matName + ".gradient"] = material.currentState.gradient();
    }
    EnvObject::allScalarProperties["depth"] = ((waterLevel * maxHeight) - heightmap.gaussianSmooth(1.f, true));
    EnvObject::allVectorProperties["depth.gradient"] = EnvObject::allScalarProperties["depth"].gradient();

    EnvObject::allScalarProperties["fracture"] = EnvObject::scenario.computeTectonic(EnvObject::allScalarProperties["fracture"].getDimensions());
    EnvObject::allVectorProperties["fracture.gradient"] = EnvObject::allScalarProperties["fracture"].gradient();
}

void EnvObject::reset()
{
    std::set<std::string> destroyedObjects;
    for (auto& obj : EnvObject::instantiatedObjects) {
        destroyedObjects.insert(obj->name);
        delete obj;
    }
    EnvObject::instantiatedObjects.clear();
    for (auto name : destroyedObjects) {
        EnvObject::recomputeTerrainPropertiesForObject(name);
    }
    initFlow(true);
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
