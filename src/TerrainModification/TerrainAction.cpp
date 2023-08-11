#include "TerrainAction.h"

TerrainAction::TerrainAction(std::string actionType, std::map<std::string, std::string> parameters, GridF exactEffect)
    : actionType(actionType), parameters(parameters), exactEffect(exactEffect)
{

}

nlohmann::json TerrainAction::serialize()
{
    nlohmann::json json;

    json["type"] = actionType;
    json["parameters"] = parameters;
    /*std::string exactEffectStr = "";
    for (const auto& val : exactEffect) exactEffectStr += std::to_string(val) + " ";
    json["exactEffect"] = nlohmann::json::object({});
    json["exactEffect"]["size"] = nlohmann::json::object({
                                                             {"x", exactEffect.sizeX},
                                                             {"y", exactEffect.sizeY},
                                                             {"z", exactEffect.sizeZ}
                                                         });
    json["exactEffect"]["data"] = exactEffectStr;*/
    return json;
}

TerrainAction TerrainAction::fromSerialized(nlohmann::json& actionSerialized)
{
    std::string actionType;
    std::map<std::string, std::string> parameters;
    GridF exactEffect;

    actionType = actionSerialized.at("type");
    parameters = actionSerialized.at("parameters");
    /*
    auto matrix = actionSerialized.at("exactEffect");
    exactEffect = GridF(matrix.at("size").at("x").get<int>(),
                                 matrix.at("size").at("y").get<int>(),
                                 matrix.at("size").at("z").get<int>());
    std::string matrixDataStr = matrix.at("data").get<std::string>();
    std::vector<float> matrixData(exactEffect.size());

    std::istringstream iss(matrixDataStr);
    std::copy(std::istream_iterator<float>(iss),
            std::istream_iterator<float>(),
            std::back_inserter(exactEffect.data));
    */
    return TerrainAction(actionType, parameters, exactEffect);
}
