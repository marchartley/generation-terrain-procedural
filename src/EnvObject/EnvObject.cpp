#include "EnvObject.h"
#include "ExpressionParser.h"

#include "GUIElements/ImageViewer.h"
#include <fstream>

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

#include "EnvObject/EnvironmentalScene.h"

EnvObject::EnvObject()
{

}

EnvObject::~EnvObject()
{

}


EnvObject* EnvObject::fromJSON(nlohmann::json content)
{
    EnvObject* obj = nullptr;
    std::string objName = content["name"];
    std::string objType = content["type"];
    std::map<std::string, float> materialDepositionRate;
    std::map<std::string, float> materialAbsorptionRate;
    std::map<std::string, float> materialDepositionOnDeath;
    TerrainTypes material = materialFromString(content["material"]);
    ImplicitPatch::PredefinedShapes shape = predefinedShapeFromString(content["geometry"]);
    Vector3 dimensions = json_to_vec3<float>(content["dimensions"]);
    HeightmapFrom heightFrom = (!content.contains("heightfrom") || content["heightfrom"] == "surface" ? SURFACE : (content["heightfrom"] == "water" ? WATER : GROUND));
    float minScore = (content.contains("minscore") ? content["minscore"].get<float>() : 0.f);

    // bool snakeDefined = false;
    SnakeSegmentationImplicit snakeParameters;

    Vector3 flowEffect;
    if (objType == "point") {
        auto asPoint = new EnvPoint;
        asPoint->radius = dimensions.x();
        asPoint->height = dimensions.z();
        obj = asPoint;
        if (content["flow"].is_number())
            flowEffect = Vector3(content["flow"], content["flow"], content["flow"]);
        else
            flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "curve") {
        auto asCurve = new EnvCurve;
        asCurve->width = dimensions.x();
        asCurve->length = dimensions.y();
        asCurve->height = dimensions.z();
        obj = asCurve;
        if (content.contains("follows")) {
            if (content["follows"] == "isovalue") {
                asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::ISOVALUE;
                snakeParameters = ContinuousCurveOptimizer::getSnakeForMinLengthCurveFollowingIsolevel(obj->fittingFunction, asCurve->length);
            }
            else if (content["follows"] == "gradients") {
                asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::GRADIENTS;
                snakeParameters = ContinuousCurveOptimizer::getSnakeForExactLengthCurveFollowingGradients(obj->fittingFunction, asCurve->length);
            }
            else if (content["follows"] == "skeleton") {
                asCurve->curveFollow = EnvCurve::CURVE_FOLLOW::SKELETON;
                snakeParameters = ContinuousCurveOptimizer::getSnakeForSkeletonCurve(obj->fittingFunction, asCurve->length);
            }
            else {
                std::cerr << "Value for 'follow' in object " << objName << " not recognized. Should be 'isovalue', 'gradients' or 'skeleton'. Got " << content["follows"] << std::endl;
            }
        }
        flowEffect = Vector3(content["flow"]["direction"], content["flow"]["normal"], content["flow"]["binormal"]);
    } else if (objType == "area") {
        auto asArea = new EnvArea;
        asArea->width = dimensions.x();
        asArea->length = dimensions.y();
        asArea->height = dimensions.z();
        obj = asArea;
        snakeParameters = ContinuousAreaOptimizer::getSnakeForAreaOptimizedShape(obj->fittingFunction, asArea->width * asArea->length);
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

    if (content.contains("snake")) {
        auto snakeParamsJSON = content["snake"];
        snakeParameters.connectivityCost = snakeParamsJSON["connectivity"];
        snakeParameters.curvatureCost = snakeParamsJSON["curvature"];
        snakeParameters.lengthCost = snakeParamsJSON["length"];
        snakeParameters.areaCost = snakeParamsJSON["area"];
        snakeParameters.imageCost = snakeParamsJSON["image"];
        snakeParameters.imageInsideCoef = snakeParamsJSON["imageinside"];
        snakeParameters.imageBordersCoef = snakeParamsJSON["imageborders"];
        snakeParameters.nbCatapillars = snakeParamsJSON["catapillars"];
        snakeParameters.targetLength = (snakeParamsJSON.contains("targetlength") ? snakeParamsJSON["targetlength"].get<float>() : snakeParameters.targetLength);
        snakeParameters.targetArea = (snakeParamsJSON.contains("targetarea") ? snakeParamsJSON["targetarea"].get<float>() : snakeParameters.targetArea);
        // snakeDefined = true;
    }
    snakeParameters.positionCost = 1.f;

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
    obj->snake = snakeParameters;
    // obj->snakeDefined = snakeDefined;
    obj->material = material;
    obj->implicitShape = shape;
    obj->inputDimensions = dimensions;
    obj->recomputeEvaluationPoints();
    obj->heightFrom = heightFrom;
    obj->minScore = minScore;
    if (dimensions.z() == 0) obj->inputDimensions = Vector3();
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
    /*
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
    return totalScore;*/
}

float EnvObject::computeGrowingState2()
{
    if (this->createdManually) return 1.f;
    // return this->evaluate();
    return (this->evaluate() > this->minScore ? 1.f : 0.f); // Let's say that it is either dead or alive...

    /*
    // float newFitnessEvaluation = this->evaluate(this->evaluationPosition);
    float newFitnessEvaluation = this->evaluate();
    // std::cout << this->fitnessScoreAtCreation << " / " << newFitnessEvaluation << std::endl;
    // if (newFitnessEvaluation <= 0) return 0;
    if (this->fitnessScoreAtCreation <= 1e-6) return 0; // This is a problem...
    return std::clamp(newFitnessEvaluation / this->fitnessScoreAtCreation, 0.f, 1.f);
    */
}

GridF EnvObject::createHeightfield()
{
    if (_patch) {
        if (auto patch = dynamic_cast<ImplicitPrimitive*>(_patch)) {
            if (_cachedHeightfield.empty() && patch && patch->predefinedShape != ImplicitPatch::NONE) {
                auto previousMaterial = patch->material;
                patch->material = SAND; // Temporarly be solid to get some height (which is then depth...)
                GridF heights(patch->getDimensions().x(), patch->getDimensions().y(), 1);
                heights.iterateParallel([&] (const Vector3& p) {
                    float resolution = .5f;
                    for (float z = 1; z < patch->getDimensions().z() * 1.f; z += resolution) {
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

std::function<float (const Vector3&)> EnvObject::parseFittingFunction(std::string formula, std::string currentObject, EnvironmentalScene* scene, bool removeSelfInstances, EnvObject *myObject)
{
    formula = toLower(formula);
    if (formula == "")
        return [](const Vector3&) { return 0.f; };

    std::map<std::string, Variable> variables;
    for (auto& [name, obj] : scene->availableObjects) {
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
    for (auto& [matName, material] : scene->materials) {
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
    const std::string depthStr = "depth";
    const std::string currenttimeStr = "currenttime";
    const std::string spawntimeStr = "spawntime";
    const std::string posStr = "pos";
    auto _func = parser.parse(formula, variables);
    return [&, formula, _func, neededVariables, currentObject, removeSelfInstances, myObject, depthStr, currenttimeStr, spawntimeStr, posStr, scene](Vector3 pos) -> float {
        // ExpressionParser parser;
        std::map<std::string, Variable> vars;
        // displayProcessTime("Variables: ", [&]() {
        for (auto& [prop, map] : scene->allVectorProperties) {
                if (!isIn(prop, neededVariables)) continue;
            if (removeSelfInstances && (startsWith(prop, currentObject + ".") || startsWith(prop, currentObject))) {
                vars[prop] = Vector3::invalid();
            } else {
                vars[prop] = map(pos);
            }
        }
        for (auto& [prop, map] : scene->allScalarProperties) {
            if (!isIn(prop, neededVariables)) continue;
            if (removeSelfInstances && (startsWith(prop, currentObject + ".") || startsWith(prop, currentObject))) {
                vars[prop] = float();
            } else {
                vars[prop] = map(pos);
            }
        }
        if (myObject && myObject->_patch && isIn(depthStr, neededVariables)) {
            // if (auto patch = dynamic_cast<ImplicitPrimitive*>(myObject->_patch)) {
            vars["depth"] = std::get<float>(vars["depth"]) + (myObject->createHeightfield().at(pos)); //.at(pos - patch->position.xy()));
            // } else {
                // std::cout << "NO" << std::endl;
            // }
        }
        if (isIn(posStr, neededVariables))
            vars["pos"] = pos;
        if (isIn(currenttimeStr, neededVariables))
            vars["currenttime"] = float(scene->currentTime);
        if (isIn(spawntimeStr, neededVariables)) {
            if (myObject != nullptr)
                vars["spawntime"] = float(std::min(scene->currentTime, myObject->spawnTime));
            else
                vars["spawntime"] = float(scene->currentTime);
        }
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
        // });
        float score;
        // displayProcessTime("Score: ", [&](){
        score = _func(vars);
        // });

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

float EnvObject::evaluate(const Vector3 &position)
{
    return this->fitnessFunction(position.xy());
}

float EnvObject::evaluate()
{
    if (this->premature) return 1.f;
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
