#include "SnakeSerializer.h"

#include "Serializer.h"

void to_json(nlohmann::json &json, const SnakeSegmentation& snake)
{
    json["connectivity-cost"] = snake.connectivityCost;
    json["curvature-cost"] = snake.curvatureCost;
    json["length-cost"] = snake.lengthCost;
    json["area-cost"] = snake.areaCost;
    json["image-cost"] = snake.imageCost;
    json["image-inside-cost"] = snake.imageInsideCoef;
    json["image-borders-cost"] = snake.imageBordersCoef;
    json["catapillars"] = snake.nbCatapillars;
    json["target-length"] = snake.targetLength;
    json["target-area"] = snake.targetArea;
    json["position-cost"] = snake.positionCost;
    json["closed"] = snake.collapseFirstAndLastPoint;
    json["position"] = snake.position;
    json["slope-cost"] = snake.slopeCost;
}

void from_json(const nlohmann::json &json, SnakeSegmentation& snake)
{
    if (json.contains("connectivity-cost"))
        snake.connectivityCost = json["connectivity-cost"];
    if (json.contains("curvature-cost"))
        snake.curvatureCost = json["curvature-cost"];
    if (json.contains("length-cost"))
        snake.lengthCost = json["length-cost"];
    if (json.contains("area-cost"))
        snake.areaCost = json["area-cost"];
    if (json.contains("image-cost"))
        snake.imageCost = json["image-cost"];
    if (json.contains("image-inside-cost"))
        snake.imageInsideCoef = json["image-inside-cost"];
    if (json.contains("image-borders-cost"))
        snake.imageBordersCoef = json["image-borders-cost"];
    if (json.contains("catapillars"))
        snake.nbCatapillars = json["catapillars"];
    if (json.contains("target-length"))
        snake.targetLength = json["target-length"];
    if (json.contains("target-area"))
        snake.targetArea = json["target-area"];
    if (json.contains("position-cost"))
        snake.positionCost = json["position-cost"];
    if (json.contains("closed"))
        snake.collapseFirstAndLastPoint = json["closed"];
    if (json.contains("position"))
        snake.position = json["position"];
    if (json.contains("slope-cost"))
        snake.slopeCost = json["slope-cost"];
}
