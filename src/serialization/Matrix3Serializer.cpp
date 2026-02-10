#include "Matrix3Serializer.h"

void to_json(nlohmann::json &json, const GridF &grid)
{
    json = stringifyGridF(grid, false);
}

void to_json(nlohmann::json &json, const GridV3 &grid)
{
    json = stringifyGridV3(grid, false);
}

void from_json(const nlohmann::json &json, GridF &grid)
{
    grid = loadGridF(json, false);
}

void from_json(const nlohmann::json &json, GridV3 &grid)
{
    grid = loadGridV3(json, false);
}
