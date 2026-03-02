#ifndef ENVOBJSERIALIZER_H
#define ENVOBJSERIALIZER_H

#include "EnvObject/EnvObject.h"
#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"
#include "EnvObject/EnvMaterial.h"

#include "Utils/json.h"

template <class Json>
void to_json(Json& json, const EnvObject& obj);
template <class Json>
void to_json(Json& json, const EnvPoint& obj);
template <class Json>
void to_json(Json& json, const EnvCurve& obj);
template <class Json>
void to_json(Json& json, const EnvArea& obj);
template <class Json>
void to_json(Json& json, const EnvObjectInstance& obj);
template <class Json>
void to_json(Json& json, const EnvPointInstance& obj);
template <class Json>
void to_json(Json& json, const EnvCurveInstance& obj);
template <class Json>
void to_json(Json& json, const EnvAreaInstance& obj);
template <class Json>
void to_json(Json& json, const EnvMaterial& material);

template <class Json>
void to_json(Json& json, const DepositionRate& depos);
template <class Json>
void to_json(Json& json, const AbsorptionRate& absorb);

template <class Json>
void to_json(Json& json, const EnvObject::HeightmapFrom& heightfrom);

template <class Json>
void to_json(Json& json, const ImplicitPatch::PredefinedShapes& predefinedShape);

template <class Json>
void to_json(Json& json, const TerrainTypes& material);

template <class Json>
void from_json(const Json& json, EnvObject& obj);
template <class Json>
void from_json(const Json& json, EnvPoint& obj);
template <class Json>
void from_json(const Json& json, EnvCurve& obj);
template <class Json>
void from_json(const Json& json, EnvArea& obj);
template <class Json>
void from_json(const Json& json, EnvObjectInstance& obj);
template <class Json>
void from_json(const Json& json, EnvPointInstance& obj);
template <class Json>
void from_json(const Json& json, EnvCurveInstance& obj);
template <class Json>
void from_json(const Json& json, EnvAreaInstance& obj);
template <class Json>
void from_json(const Json& json, EnvMaterial& material);

template <class Json>
void from_json(const Json& json, DepositionRate& depos);
template <class Json>
void from_json(const Json& json, AbsorptionRate& absorb);

template <class Json>
void from_json(const Json& json, EnvObject::HeightmapFrom& heightfrom);

template <class Json>
void from_json(const Json& json, ImplicitPatch::PredefinedShapes& predefinedShape);

template <class Json>
void from_json(const Json& json, TerrainTypes& material);





template <class Json>
void to_json(Json& json, const EnvObject* envObject);

template <class Json>
void from_json(const Json& json, EnvObject*& envObject);


template <class Json>
void to_json(Json& json, const EnvObjectInstance* envObject);

template <class Json>
void from_json(const Json& json, EnvObjectInstance*& envObject);
























#include "Serializer.h"





template <class Json>
void to_json(Json& json, const EnvObject& obj)
{
    json["name"] = obj.name;
    json["fitness"] = obj.s_FitnessFunction;
    json["fitting"] = obj.s_FittingFunction;
    json["snake"] = obj.snakeParameters;
    json["material"] = obj.material;
    json["geometry"] = obj.implicitShape;
    json["height-from"] = obj.heightFrom;
    json["min-score"] = obj.minScore;
    json["deposition-rate"] = obj.materialDepositionRate;
    json["absorption-rate"] = obj.materialAbsorptionRate;
    json["on-death-deposition"] = obj.materialDepositionOnDeath;
    json["height"] = obj.height;
}


template <class Json>
void to_json(Json& json, const EnvPoint& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = obj.mainKelvinlets;
    json["radius"] = obj.radius;
    json["type"] = "Point";
}


template <class Json>
void to_json(Json& json, const EnvCurve& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = nlohmann::json({
        {"starting-effect", obj.startingPointKelvinlets},
        {"ending-effect", obj.endingPointKelvinlets},
        {"curve-effect", obj.curveKelvinlets}
    });
    json["length"] = obj.length;
    json["width"] = obj.width;
    json["type"] = "Curve";
}


template <class Json>
void to_json(Json& json, const EnvArea& obj)
{
    to_json(json, static_cast<const EnvObject&>(obj));
    json["flow-effect"] = nlohmann::json({
        {"curve-effect", obj.curveKelvinlets},
        {"flow-attenuation", obj.flowAttenuation}
    });
    json["length"] = obj.length;
    json["width"] = obj.width;
    json["type"] = "Area";
}


template <class Json>
void to_json(Json& json, const EnvObjectInstance& obj)
{
    json["name"] = obj.getDefinition()->name;
    json["age"] = obj.age;
    json["manual"] = obj.createdManually;
    json["fitness-at-creation"] = obj.fitnessScoreAtCreation;
    json["orientation"] = obj.storedOrientation;
    json["id"] = obj.ID;
}


template <class Json>
void to_json(Json& json, const EnvPointInstance& obj)
{
    to_json(json, static_cast<const EnvObjectInstance&>(obj));
    json["position"] = obj.position;
}


template <class Json>
void to_json(Json& json, const EnvCurveInstance& obj)
{
    to_json(json, static_cast<const EnvObjectInstance&>(obj));
    json["curve"] = obj.curve;
}


template <class Json>
void to_json(Json& json, const EnvAreaInstance& obj)
{
    to_json(json, static_cast<const EnvObjectInstance&>(obj));
    json["curve"] = obj.curve;
}

template <class Json>
void to_json(Json& json, const EnvMaterial& material)
{
    json["name"] = material.name;
    // json["data"] = material.currentState;
    json["decay"] = material.decay;
    json["diffusion-speed"] = material.diffusionSpeed;
    json["mass"] = material.mass;
    json["water-transport"] = material.waterTransport;
    json["virtual-height"] = material.virtualHeight;
}



template <class Json>
void to_json(Json&json, const DepositionRate &depos)
{
    json["radius"] = depos.radius;
    json["rate"] = depos.rate;
}


template <class Json>
void to_json(Json&json, const AbsorptionRate &absorb)
{
    json["radius"] = absorb.radius;
    json["rate"] = absorb.rate;
}


template <class Json>
void to_json(Json&json, const EnvObject::HeightmapFrom &heightfrom)
{
    if (heightfrom == EnvObject::SURFACE) { json = "Surface"; return; }
    if (heightfrom == EnvObject::WATER) { json = "Water"; return; }
    if (heightfrom == EnvObject::GROUND) { json = "Ground"; return; }
}



template <class Json>
void to_json(Json&json, const ImplicitPatch::PredefinedShapes& predefinedShape)
{
    json = stringFromPredefinedShape(predefinedShape);
}



template <class Json>
void to_json(Json&json, const TerrainTypes& material)
{
    json = stringFromMaterial(material);
}






template <class Json>
void from_json(const Json& json, EnvObject& obj)
{
    std::map<std::string, DepositionRate> materialDepositionRate;
    std::map<std::string, AbsorptionRate> materialAbsorptionRate;
    std::map<std::string, DepositionRate> materialDepositionOnDeath;

    if (json.contains("deposition-rate")) {
        auto deposits = json.at("deposition-rate").template get<std::map<std::string, DepositionRate>>();
        for (auto& [mat, val] : deposits) {
            materialDepositionRate[mat] = val;
        }
    }
    if (json.contains("absorption-rate")) {
        auto absorbs = json.at("absorption-rate").template get<std::map<std::string, AbsorptionRate>>();
        for (auto& [mat, val] : absorbs) {
            materialAbsorptionRate[mat] = val;
        }
    }
    if (json.contains("on-death-deposition")) {
        auto depos = json.at("on-death-deposition").template get<std::map<std::string, DepositionRate>>();
        for (auto& [mat, val] : depos) {
            materialDepositionOnDeath[mat] = val;
        }
    }


    obj.name = toLower(json.at("name"));
    obj.materialAbsorptionRate = materialAbsorptionRate;
    obj.materialDepositionRate = materialDepositionRate;
    obj.materialDepositionOnDeath = materialDepositionOnDeath;
    if (json.contains("rule"))
        obj.s_FitnessFunction = json.at("rule");
    else
        obj.s_FitnessFunction = json.at("fitness");

    if (json.contains("fitting")) {
        obj.s_FittingFunction = json.at("fitting");
    } else {
        obj.s_FittingFunction = "";
    }
    // if (json.contains("snake"))
        // obj.snakeParameters = json.at("snake");
    obj.snakeParameters = new SnakeSegmentationParameters;
    from_json(json.at("snake"), *obj.snakeParameters);

    obj.material = json.at("material");
    obj.implicitShape = json.at("geometry");
    obj.heightFrom = json.value("height-from", EnvObject::SURFACE);
    obj.minScore = json.at("min-score");
    obj.height = json.at("height");
}

template <class Json>
void from_json(const Json& json, EnvPoint& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));
    obj.radius = json.at("radius");
    if (json.contains("flow-effect")) {
        obj.mainKelvinlets = json.at("flow-effect");
    }
}

template <class Json>
void from_json(const Json& json, EnvCurve& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));

    obj.width = json.at("width");
    obj.length = json.at("length");
    /*if (json.contains("follows")) {
        if (json.at("follows") == "isovalue") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::ISOVALUE;
            obj.snake = ContinuousCurveOptimizer::getSnakeForMinLengthCurveFollowingIsolevel(obj.fittingFunction, obj.length);
        }
        else if (json.at("follows") == "gradients") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::GRADIENTS;
            obj.snake = ContinuousCurveOptimizer::getSnakeForExactLengthCurveFollowingGradients(obj.fittingFunction, obj.length);
        }
        else if (json.at("follows") == "skeleton") {
            obj.curveFollow = EnvCurve::CURVE_FOLLOW::SKELETON;
            obj.snake = ContinuousCurveOptimizer::getSnakeForSkeletonCurve(obj.fittingFunction, obj.length);
        }
        else {
            std::cerr << "Value for 'follow' in object " << obj.name << " not recognized. Should be 'isovalue', 'gradients' or 'skeleton'. Got " << json.at("follows") << std::endl;
        }
    }*/
    if (json.contains("flow-effect")) {
        auto flow = json.at("flow-effect");
        obj.startingPointKelvinlets = flow.at("starting-effect");
        obj.endingPointKelvinlets = flow.at("ending-effect");
        obj.curveKelvinlets = flow.at("curve-effect");
    }

    // obj.recomputeEvaluationPoints();
}

template <class Json>
void from_json(const Json& json, EnvArea& obj)
{
    from_json(json, static_cast<EnvObject&>(obj));

    obj.width = json.at("width");
    obj.length = json.at("length");
    // obj.snake = ContinuousAreaOptimizer::getSnakeForAreaOptimizedShape(obj.fittingFunction, obj.width * obj.length);

    if (json.contains("flow-effect")) {
        auto flow = json.at("flow-effect");
        obj.curveKelvinlets = flow.at("curve-effect");
        obj.flowAttenuation = flow.at("flow-attenuation");
    }

    // obj.recomputeEvaluationPoints();
}





template <class Json>
void from_json(const Json& json, EnvObjectInstance& obj, EnvironmentalScene*& scene)
{
    obj.ID = json.at("id");
    obj.createdManually = json.at("manual");
    obj.fitnessScoreAtCreation = json.at("fitness-at-creation");
    obj.scene = scene;
    obj.definition = obj.scene->availableObjects[json.at("name")];
}

template <class Json>
void from_json(const Json& json, EnvPointInstance& obj, EnvironmentalScene*& scene)
{
    from_json(json, static_cast<EnvObjectInstance&>(obj), scene);
    obj.position = json.at("position");
}

template <class Json>
void from_json(const Json& json, EnvCurveInstance& obj, EnvironmentalScene*& scene)
{
    from_json(json, static_cast<EnvObjectInstance&>(obj), scene);
    obj.curve = json.at("curve");
}

template <class Json>
void from_json(const Json& json, EnvAreaInstance& obj, EnvironmentalScene*& scene)
{
    from_json(json, static_cast<EnvObjectInstance&>(obj), scene);
    obj.curve = json.at("curve");
}




template <class Json>
void from_json(const Json& json, EnvMaterial& material)
{
    material.name = json.at("name");
    // material.currentState = json.at("data");
    material.decay = json.at("decay");
    material.diffusionSpeed = json.at("diffusion-speed");
    material.mass = json.at("mass");
    material.waterTransport = json.at("water-transport");
    material.virtualHeight = json.at("virtual-height");
}

template <class Json>
void from_json(const Json &json, DepositionRate &depos)
{
    depos.radius = json.at("radius");
    depos.rate = json.at("rate");
}

template <class Json>
void from_json(const Json &json, AbsorptionRate &absorb)
{
    absorb.radius = json.at("radius");
    absorb.rate = json.at("rate");
}

template <class Json>
void from_json(const Json &json, EnvObject::HeightmapFrom &heightfrom)
{
    const std::string identifier = toLower(json);
    if (identifier == "Surface") { heightfrom = EnvObject::SURFACE; return; }
    if (identifier == "Water") { heightfrom = EnvObject::WATER; return; }
    if (identifier == "Ground") { heightfrom = EnvObject::GROUND; return; }
}


template <class Json>
void from_json(const Json &json, ImplicitPatch::PredefinedShapes& predefinedShape)
{
    predefinedShape = predefinedShapeFromString(json);
}


template <class Json>
void from_json(const Json &json, TerrainTypes& material)
{
    material = materialFromString(json);
}
















template <class Json>
void to_json(Json& json, const EnvObject* envObj)
{
    if (auto obj = dynamic_cast<const EnvPoint*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvCurve*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvArea*>(envObj))
        to_json(json, *obj);
}


template <class Json>
EnvObject* make_envobj_from_json(const Json& json)
{
    const std::string type = toLower(json["type"]);

    if (type == "point")  return new EnvPoint();
    if (type == "curve") return new EnvCurve();
    if (type == "area") return new EnvArea();

    throw std::runtime_error("Unknown Environmental Object type: " + type);
}


template <class Json>
void from_json(const Json& json, EnvObject*& envObject)
{
    envObject = make_envobj_from_json(json);
    from_json(json, *envObject);
}








template <class Json>
void to_json(Json& json, const EnvObjectInstance* envObj)
{
    if (auto obj = dynamic_cast<const EnvPointInstance*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvCurveInstance*>(envObj))
        to_json(json, *obj);
    if (auto obj = dynamic_cast<const EnvAreaInstance*>(envObj))
        to_json(json, *obj);
}


template <class Json>
EnvObjectInstance* make_envobj_instance_from_json(const Json& json, EnvironmentalScene*& scene)
{
    return scene->instantiate(json.at("name"));
}


template <class Json>
void from_json(const Json& json, EnvObjectInstance*& envObject, EnvironmentalScene*& scene)
{
    envObject = make_envobj_instance_from_json(json, scene);
    from_json(json, *envObject, scene);
}


#endif // ENVOBJSERIALIZER_H
