#include "EnvObjSerializer.h"

#include "Serializer.h"



void to_json(nlohmann::json& json, const EnvObject& obj)
{
    json["name"] = obj.name;
    json["ID"] = obj.ID;
    json["age"] = obj.age;
    json["fitnessScoreAtCreation"] = obj.fitnessScoreAtCreation;
}

void to_json(nlohmann::json& json, const EnvPoint& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["position"] = obj.position;
}

void to_json(nlohmann::json& json, const EnvCurve& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["curve"] = obj.curve;
}

void to_json(nlohmann::json& json, const EnvArea& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["curve"] = obj.curve;
}

void to_json(nlohmann::json& json, const EnvMaterial& material)
{
    json["name"] = material.name;
    json["data"] = material.currentState;
}


void to_json(nlohmann::json &json, const DepositionRate &depos)
{
    json["radius"] = depos.radius;
    json["rate"] = depos.rate;
}

void to_json(nlohmann::json &json, const AbsorptionRate &absorb)
{
    json["radius"] = absorb.radius;
    json["rate"] = absorb.rate;
}








void from_json(const nlohmann::json& json, EnvObject& obj)
{
    std::map<std::string, DepositionRate> materialDepositionRate;
    std::map<std::string, AbsorptionRate> materialAbsorptionRate;
    std::map<std::string, DepositionRate> materialDepositionOnDeath;

    if (json.contains("deposition-rate")) {
        auto deposits = json["deposition-rate"].get<std::map<std::string, DepositionRate>>();
        for (auto& [mat, val] : deposits) {
            materialDepositionRate[mat] = val;
        }
    }
    if (json.contains("absorption-rate")) {
        auto absorbs = json["absorption-rate"].get<std::map<std::string, AbsorptionRate>>();
        for (auto& [mat, val] : absorbs) {
            materialAbsorptionRate[mat] = val;
        }
    }
    if (json.contains("on-death-deposition")) {
        auto depos = json["on-death-deposition"].get<std::map<std::string, DepositionRate>>();
        for (auto& [mat, val] : depos) {
            materialDepositionOnDeath[mat] = val;
        }
    }


    obj.name = toLower(json["name"]);
    obj.materialAbsorptionRate = materialAbsorptionRate;
    obj.materialDepositionRate = materialDepositionRate;
    obj.materialDepositionOnDeath = materialDepositionOnDeath;
    obj.s_FittingFunction = json["rule"];
    if (json.contains("fitness")) {
        obj.s_FitnessFunction = json["fitness"];
    } else {
        obj.s_FitnessFunction = obj.s_FittingFunction;
    }
    if (json.contains("snake"))
        obj.snake = json["snake"];
    obj.material = materialFromString(json["material"]);
    obj.implicitShape = predefinedShapeFromString(json["geometry"]);
    obj.heightFrom = (!json.contains("heightfrom") || json["heightfrom"] == "surface" ? EnvObject::SURFACE : (json["heightfrom"] == "water" ? EnvObject::WATER : EnvObject::GROUND));
    obj.minScore = (json.contains("minscore") ? json["minscore"].get<float>() : 0.f);
}

void from_json(const nlohmann::json& json, EnvPoint& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));
    obj.radius = json["radius"];
    obj.height = json["height"];
    // obj.flowEffect = json["flow"];
    obj.recomputeEvaluationPoints();
}

void from_json(const nlohmann::json& json, EnvCurve& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));

    obj.width = json["width"];
    obj.length = json["length"];
    obj.height = json["height"];
    if (json.contains("follows")) {
        if (json["follows"] == "isovalue") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::ISOVALUE;
            obj.snake = ContinuousCurveOptimizer::getSnakeForMinLengthCurveFollowingIsolevel(obj.fittingFunction, obj.length);
        }
        else if (json["follows"] == "gradients") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::GRADIENTS;
            obj.snake = ContinuousCurveOptimizer::getSnakeForExactLengthCurveFollowingGradients(obj.fittingFunction, obj.length);
        }
        else if (json["follows"] == "skeleton") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::SKELETON;
            obj.snake = ContinuousCurveOptimizer::getSnakeForSkeletonCurve(obj.fittingFunction, obj.length);
        }
        else {
            std::cerr << "Value for 'follow' in object " << obj.name << " not recognized. Should be 'isovalue', 'gradients' or 'skeleton'. Got " << json["follows"] << std::endl;
        }
    }
    // obj.flowEffect = json["flow"];
    obj.recomputeEvaluationPoints();
}

void from_json(const nlohmann::json& json, EnvArea& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));

    obj.width = json["width"];
    obj.length = json["length"];
    obj.height = json["height"];
    obj.snake = ContinuousAreaOptimizer::getSnakeForAreaOptimizedShape(obj.fittingFunction, obj.width * obj.length);
    // obj.flowEffect = json["flow"];
    obj.recomputeEvaluationPoints();
}

void from_json(const nlohmann::json& json, EnvMaterial& material)
{
    material.name = json["name"];
    material.currentState = json["data"];
}

void from_json(const nlohmann::json &json, DepositionRate &depos)
{
    depos.radius = json["radius"];
    depos.rate = json["rate"];
}

void from_json(const nlohmann::json &json, AbsorptionRate &absorb)
{
    absorb.radius = json["radius"];
    absorb.rate = json["rate"];
}
