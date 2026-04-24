#ifndef BSPLINESERIALIZER_H
#define BSPLINESERIALIZER_H

#include "Curves/BSpline.h"

#include "Utils/json.h"


template <class Json>
void to_json(Json& json, const BSpline& spline);
template <class Json>
void from_json(const Json& json, BSpline& spline);


#include "Serializer.h"

template <class Json>
void to_json(Json &json, const BSpline &spline)
{
    auto p = spline.getPoints();
    std::vector<Json> points(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        points[i] = p[i];
    }
    json["points"] = points;
    json["closed"] = spline.isClosed();
}


template <class Json>
void from_json(const Json& json, BSpline& spline) {
    auto points = json.at("points");
    // spline.points.resize(points.size());
    for (size_t i = 0; i < points.size(); i++)
        spline.addPoint(points[i]); // spline.points[i] = points[i];
    if (json.at("closed").template get<bool>())
        spline.close();
}
#endif // BSPLINESERIALIZER_H
