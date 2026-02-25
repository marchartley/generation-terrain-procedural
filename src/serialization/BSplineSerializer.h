#ifndef BSPLINESERIALIZER_H
#define BSPLINESERIALIZER_H

#include "Utils/BSpline.h"

#include "Utils/json.h"


template <class Json>
void to_json(Json& json, const BSpline& spline);
template <class Json>
void from_json(const Json& json, BSpline& spline);


#include "Serializer.h"

template <class Json>
void to_json(Json &json, const BSpline &spline)
{
    std::vector<Json> points(spline.points.size());
    for (size_t i = 0; i < spline.points.size(); i++) {
        points[i] = spline.points[i];
    }
    json["points"] = points;
    json["closed"] = spline.closed;
}


template <class Json>
void from_json(const Json& json, BSpline& spline) {
    auto points = json.at("points");
    spline.points.resize(points.size());
    for (size_t i = 0; i < points.size(); i++)
        spline.points[i] = points[i];
    if (json.at("closed").template get<bool>())
        spline.close();
}
#endif // BSPLINESERIALIZER_H
