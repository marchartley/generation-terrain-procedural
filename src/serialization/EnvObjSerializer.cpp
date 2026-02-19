#include "EnvObjSerializer.h"

#include "Serializer.h"




void to_json(nlohmann::json& json, const EnvObject& obj)
{
    json["deposition-rate"] = obj.materialDepositionRate;
    json["absorption-rate"] = obj.materialAbsorptionRate;
    json["on-death-deposition"] = obj.materialDepositionOnDeath;

    json["name"] = obj.name;
    json["fitness"] = obj.s_FitnessFunction;
    json["fitting"] = obj.s_FittingFunction;
    json["snake"] = obj.snake;
    json["material"] = obj.material;
    json["geometry"] = obj.implicitShape;
    json["height-from"] = obj.heightFrom;
    json["min-score"] = obj.minScore;
}

void to_json(nlohmann::json& json, const EnvPoint& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = obj.mainKelvinlets;
    json["radius"] = obj.radius;
    json["type"] = "Point";
}

void to_json(nlohmann::json& json, const EnvCurve& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = nlohmann::json({
      {"starting-effect", obj.startingPointKelvinlets},
      {"ending-effect", obj.endingPointKelvinlets},
      {"curve-effect", obj.curveKelvinlets}
    });
    json["length"] = obj.length;
    json["height"] = obj.height;
    json["width"] = obj.width;
    json["type"] = "Curve";
}

void to_json(nlohmann::json& json, const EnvArea& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = nlohmann::json({
        {"curve-effect", obj.curveKelvinlets}
    });
    json["length"] = obj.length;
    json["height"] = obj.height;
    json["width"] = obj.width;
    json["type"] = "Area";
}

void to_json(nlohmann::json& json, const EnvMaterial& material)
{
    json["name"] = material.name;
    // json["data"] = material.currentState;
    json["decay"] = material.decay;
    json["diffusion-speed"] = material.diffusionSpeed;
    json["mass"] = material.mass;
    json["water-transport"] = material.waterTransport;
    json["virtual-height"] = material.virtualHeight;
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

void to_json(nlohmann::json &json, const EnvObject::HeightmapFrom &heightfrom)
{
    if (heightfrom == EnvObject::SURFACE) { json = "Surface"; return; }
    if (heightfrom == EnvObject::WATER) { json = "Water"; return; }
    if (heightfrom == EnvObject::GROUND) { json = "Ground"; return; }
}


void to_json(nlohmann::json &json, const ImplicitPatch::PredefinedShapes& predefinedShape)
{
    json = stringFromPredefinedShape(predefinedShape);
}


void to_json(nlohmann::json &json, const TerrainTypes& material)
{
    json = stringFromMaterial(material);
}

/*
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
*/







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
    obj.heightFrom = json.value("height-from", EnvObject::SURFACE);
    obj.minScore = (json.contains("min-score") ? json["min-score"].get<float>() : 0.f);
}

void from_json(const nlohmann::json& json, EnvPoint& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));
    obj.radius = json["radius"];
    obj.height = json["height"];

    if (json.contains("flow-effect")) {
        auto depos = json["flow-effect"];
    }
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

void from_json(const nlohmann::json &json, EnvObject::HeightmapFrom &heightfrom)
{
    const std::string identifier = toLower(json.get<std::string>());
    if (identifier == "Surface") { heightfrom = EnvObject::SURFACE; return; }
    if (identifier == "Water") { heightfrom = EnvObject::WATER; return; }
    if (identifier == "Ground") { heightfrom = EnvObject::GROUND; return; }
}


void from_json(const nlohmann::json &json, ImplicitPatch::PredefinedShapes& predefinedShape)
{
    predefinedShape = predefinedShapeFromString(json);
}


void from_json(const nlohmann::json &json, TerrainTypes& material)
{
    material = materialFromString(json);
}















void to_json(nlohmann::json& json, const EnvObject* envObj)
{
    if (auto obj = dynamic_cast<const EnvPoint*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvCurve*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvArea*>(envObj))
        to_json(json, *obj);
}



EnvObject* make_envobj_from_json(const nlohmann::json& j)
{
    const std::string type = toLower(j["type"]);

    if (type == "point")  return new EnvPoint();
    if (type == "curve") return new EnvCurve();
    if (type == "area") return new EnvArea();

    throw std::runtime_error("Unknown Environmental Object type: " + type);
}


void from_json(const nlohmann::json& json, EnvObject*& envObject)
{
    envObject = make_envobj_from_json(json);
    from_json(json, *envObject);
}














