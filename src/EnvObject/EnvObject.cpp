#include "EnvObject.h"
#include "ExpressionParser.h"

#include "GUIElements/ImageViewer.h"
#include <fstream>

#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"

#include "EnvObject/EnvironmentalScene.h"

#include "serialization/Serializer.h"

EnvObject::EnvObject()
{
    this->snakeParameters = new SnakeSegmentationParameters;
    this->snakeField = new SnakeImageFieldImplicit;
}

EnvObject::~EnvObject()
{

}

EnvObjectInstance::EnvObjectInstance()
    : EnvObjectInstance(nullptr)
{

}
EnvObjectInstance::EnvObjectInstance(EnvObject* definition)
    : definition(definition), scene(definition->scene)
{
    this->snake = SnakeSegmentation(definition->snakeParameters, definition->snakeField);
}

float EnvObjectInstance::computeGrowingState()
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

float EnvObjectInstance::computeGrowingState2()
{
    if (this->createdManually) return 1.f;
    // return this->evaluate();
    return (this->evaluate() > this->getDefinition()->minScore ? 1.f : 0.f); // Let's say that it is either dead or alive...

    /*
    // float newFitnessEvaluation = this->evaluate(this->evaluationPosition);
    float newFitnessEvaluation = this->evaluate();
    // std::cout << this->fitnessScoreAtCreation << " / " << newFitnessEvaluation << std::endl;
    // if (newFitnessEvaluation <= 0) return 0;
    if (this->fitnessScoreAtCreation <= 1e-6) return 0; // This is a problem...
    return std::clamp(newFitnessEvaluation / this->fitnessScoreAtCreation, 0.f, 1.f);
    */
}

GridF EnvObjectInstance::createHeightfield()
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

std::function<float (const Vector3&)> EnvObject::parseFittingFunction(std::string formula, std::string currentObject, EnvironmentalScene* scene, bool removeSelfInstances, EnvObjectInstance *myObject)
{
    formula = toLower(formula);
    if (formula == "")
        return [](const Vector3&) { return 0.f; };

    std::map<std::string, Variable> variables;
    for (auto& [name, obj] : scene->availableObjects) {
        variables[name] = Vector3::invalid;
        variables[name + ".center"] = Vector3::invalid;
        variables[name + ".start"] = Vector3::invalid;
        variables[name + ".end"] = Vector3::invalid;
        variables[name + ".normal"] = Vector3::invalid;
        variables[name + ".dir"] = Vector3::invalid;
        variables[name + ".inside"] = float();
        variables[name + ".curvature"] = float();
    }

    variables["current"] = Vector3::invalid;
    variables["current.center"] = Vector3::invalid;
    variables["current.start"] = Vector3::invalid;
    variables["current.end"] = Vector3::invalid;
    variables["current.normal"] = Vector3::invalid;
    variables["current.dir"] = Vector3::invalid;
    variables["current.vel"] = float();
    variables["current.gradient"] = Vector3::invalid;
    variables["depth"] = float();
    variables["depth.gradient"] = Vector3::invalid;
    variables["fracture"] = float();
    variables["fracture.gradient"] = Vector3::invalid;
    for (auto& [matName, material] : scene->materials) {
        variables[matName] = float();
        variables[matName + ".gradient"] = Vector3::invalid;
    }

    ExpressionParser parser;
    variables["currenttime"] = float();
    variables["spawntime"] = float();
    variables["pos"] = Vector3::invalid;
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
                vars[prop] = Vector3::invalid;
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
        for (const std::string& neededVar : neededVariables) {
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

void EnvObject::updateFittingFunction()
{
    if (this->s_FittingFunction != "") {
        this->snakeField->imageField = this->fittingFunction;
        this->snakeField->gradientField = gradientFromFieldFunction(this->fittingFunction);
    } else {
        std::cout << "Warning: '" << name << "' has no fitting function defined, using the fitness." << std::endl;

        this->snakeField->imageField = this->fitnessFunction;
        this->snakeField->gradientField = gradientFromFieldFunction(this->fitnessFunction);
    }
}

float EnvObjectInstance::evaluate(const Vector3 &position)
{
    return this->getDefinition()->fitnessFunction(position.xy());
}

float EnvObjectInstance::evaluate()
{
    if (this->premature) return 1.f;
    // Should only be at one point...

    if (evaluationPositions.empty()) {
        std::cerr << "Object " << this->getDefinition()->name << " has no evaluation point..." << std::endl;
        return 0;
    }
    float totalScore = 0;
    for (auto& p : evaluationPositions) {
        float score = evaluate(p);
        if (score != score) {
            std::cerr << "NaN found while evaluating " << this->getDefinition()->name << " at " << p << std::endl;
        } else {
            totalScore += score;
        }
    }
    return totalScore / float(evaluationPositions.size());
}

void EnvObjectInstance::die()
{
    this->applyDepositionOnDeath();
}
