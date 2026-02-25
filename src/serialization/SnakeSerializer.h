#ifndef SNAKESERIALIZER_H
#define SNAKESERIALIZER_H

#include "EnvObject/SnakeSegmentation.h"

#include "Utils/json.h"

template <class Json>
void to_json(Json& json, const SnakeSegmentation& snake);

template <class Json>
void from_json(const Json& json, SnakeSegmentation& snake);





#include "Serializer.h"

template <class Json>
void to_json(Json &json, const SnakeSegmentation& snake)
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

template <class Json>
void from_json(const Json &json, SnakeSegmentation& snake)
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

#endif // SNAKESERIALIZER_H
