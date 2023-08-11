#ifndef TERRAINACTION_H
#define TERRAINACTION_H

#include <map>
#include "DataStructure/Matrix3.h"
#include "Utils/json.h"


class TerrainAction
{
public:
    TerrainAction(std::string actionType, std::map<std::string, std::string> parameters, GridF exactEffect);

    nlohmann::json serialize();
    static TerrainAction fromSerialized(nlohmann::json& actionSerialized);

    std::string actionType;
    std::map<std::string, std::string> parameters;
    GridF exactEffect;
};

#endif // TERRAINACTION_H
