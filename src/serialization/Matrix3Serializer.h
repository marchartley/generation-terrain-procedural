#ifndef MATRIX3SERIALIZER_H
#define MATRIX3SERIALIZER_H

#include "DataStructure/Matrix3.h"
#include "Utils/json.h"

void to_json(nlohmann::json& json, const GridF& grid);
void to_json(nlohmann::json& json, const GridV3& grid);

void from_json(const nlohmann::json& json, GridF& grid);
void from_json(const nlohmann::json& json, GridV3& grid);

#endif // MATRIX3SERIALIZER_H
