#ifndef BSPLINESERIALIZER_H
#define BSPLINESERIALIZER_H

#include "Utils/BSpline.h"

#include "Utils/json.h"


void to_json(nlohmann::json& json, const BSpline& spline);
void from_json(const nlohmann::json& json, BSpline& spline);

#endif // BSPLINESERIALIZER_H
