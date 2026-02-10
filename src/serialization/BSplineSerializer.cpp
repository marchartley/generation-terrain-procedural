#include "BSplineSerializer.h"

#include "Serializer.h"

void to_json(nlohmann::json &json, const BSpline &spline)
{
    std::vector<nlohmann::json> points(spline.points.size());
    for (size_t i = 0; i < spline.points.size(); i++) {
        points[i] = spline.points[i];
    }
    json["points"] = points;
    json["closed"] = spline.closed;
}


void from_json(const nlohmann::json& json, BSpline& spline) {
    auto points = json.at("points");
    spline.points.resize(points.size());
    for (size_t i = 0; i < points.size(); i++)
        spline.points[i] = points[i].get<Vector3>();
    if (json.at("closed").get<bool>())
        spline.close();
}
