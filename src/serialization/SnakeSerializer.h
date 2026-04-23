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
void to_json(Json &json, const SnakeSegmentationParameters& snake)
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
    json["slope-cost"] = snake.slopeCost;
}

template <class Json>
void from_json(const Json &json, SnakeSegmentationParameters& snake)
{
    if (json.contains("connectivity-cost"))
        snake.connectivityCost = json.at("connectivity-cost");
    if (json.contains("curvature-cost"))
        snake.curvatureCost = json.at("curvature-cost");
    if (json.contains("length-cost"))
        snake.lengthCost = json.at("length-cost");
    if (json.contains("area-cost"))
        snake.areaCost = json.at("area-cost");
    if (json.contains("image-cost"))
        snake.imageCost = json.at("image-cost");
    if (json.contains("image-inside-cost"))
        snake.imageInsideCoef = json.at("image-inside-cost");
    if (json.contains("image-borders-cost"))
        snake.imageBordersCoef = json.at("image-borders-cost");
    if (json.contains("catapillars"))
        snake.nbCatapillars = json.at("catapillars");
    if (json.contains("target-length"))
        snake.targetLength = json.at("target-length");
    if (json.contains("target-area"))
        snake.targetArea = json.at("target-area");
    if (json.contains("position-cost"))
        snake.positionCost = json.at("position-cost");
    if (json.contains("closed"))
        snake.collapseFirstAndLastPoint = json.at("closed");
    if (json.contains("slope-cost"))
        snake.slopeCost = json.at("slope-cost");
}


template <class Json>
void to_json(Json& json, const SnakeSegmentation& snake) {
    json["parameters"] = snake.params;
    if (dynamic_cast<SnakeImageFieldImplicit*>(snake.field)) {
        json["type"] = "implicit";
    } else if (dynamic_cast<SnakeImageFieldExplicit*>(snake.field)) {
        json["type"] = "explicit";
    } else {
        json["type"] = "undefined";
    }
}

template <class Json>
void from_json(const Json& json, SnakeSegmentation& snake) {
    if (json.contains("type")) {
        std::string type = json.at("type");
        if (type == "implicit") {
            snake.field = new SnakeImageFieldImplicit;
        } else if (type == "explicit") {
            snake.field = new SnakeImageFieldExplicit;
        } else {
            throw std::invalid_argument("Snake type unknown: '" + type + "'");
        }
    } else {
        snake.field = new SnakeImageFieldImplicit; // Default, most versatile
    }
    snake.params = new SnakeSegmentationParameters;
    from_json(json, *snake.params);
}

#endif // SNAKESERIALIZER_H
