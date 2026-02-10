#ifndef SNAKESERIALIZER_H
#define SNAKESERIALIZER_H

#include "EnvObject/SnakeSegmentation.h"

#include "Utils/json.h"

void to_json(nlohmann::json& json, const SnakeSegmentation& snake);

void from_json(const nlohmann::json& json, SnakeSegmentation& snake);

#endif // SNAKESERIALIZER_H
