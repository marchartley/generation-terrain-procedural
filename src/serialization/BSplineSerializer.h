#ifndef BSPLINESERIALIZER_H
#define BSPLINESERIALIZER_H

#include "Curves/Curves.h"

#include "Utils/json.h"


template <class Json>
void to_json(Json& json, const BSpline& spline);
template <class Json>
void from_json(const Json& json, BSpline& spline);


#include "Serializer.h"

template <class Json>
void to_json(Json& json, const Curve& curve) {
    json["closed"] = curve.isClosed();
}

template <class Json>
void to_json(Json &json, const BSpline& curve)
{
    to_json(json, static_cast<const Curve&>(curve));
    auto p = curve.getPath();
    std::vector<Json> points(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        points[i] = p[i];
    }
    json["points"] = points;
    json["type"] = "catmull";
}

template <class Json>
void to_json(Json &json, const Polyline& curve)
{
    to_json(json, static_cast<const Curve&>(curve));
    auto p = curve.getPath();
    std::vector<Json> points(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        points[i] = p[i];
    }
    json["points"] = points;
    json["type"] = "polyline";
}

template <class Json>
void to_json(Json &json, const BezierCurve& curve)
{
    to_json(json, static_cast<const Curve&>(curve));
    auto p = curve.getPoints();
    std::vector<Json> points(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        points[i] = p[i];
    }
    auto h = curve.getHandles();
    std::vector<Json> handles(p.size());
    for (size_t i = 0; i < h.size(); i++) {
        handles[i] = h[i];
    }
    json["points"] = points;
    json["handles"] = handles;
    json["type"] = "bezier";
}





template <class Json>
void from_json(const Json& json, Curve& curve)
{
    if (json.at("closed").template get<bool>())
        curve.close();
}

template <class Json>
void from_json(const Json& json, BSpline& curve) {
    auto points = json.at("points").template get<std::vector<Vector3>>();
    curve = BSpline(points);
    from_json(json, static_cast<Curve&>(curve));
}

template <class Json>
void from_json(const Json& json, Polyline& curve) {
    auto points = json.at("points").template get<std::vector<Vector3>>();
    curve = Polyline(points);
    from_json(json, static_cast<Curve&>(curve));
}

template <class Json>
void from_json(const Json& json, BezierCurve& curve) {
    auto points = json.at("points").template get<std::vector<Vector3>>();
    auto handles = json.at("handles").template get<std::vector<Vector3>>();
    curve = BezierCurve(points, handles);
    from_json(json, static_cast<Curve&>(curve));
}




template <class Json>
Curve* make_curve_from_json(const Json& json) {
    Curve* c;
    bool foundType = false;
    if (json.contains("type")) {
        const std::string type = toLower(json.at("type"));
        if (type == "catmull") {
            c = new BSpline;
            from_json(json, *(dynamic_cast<BSpline*>(c)));
            foundType = true;
        } else if (type == "polyline") {
            c = new Polyline;
            from_json(json, *(dynamic_cast<Polyline*>(c)));
            foundType = true;
        } else if (type == "bezier") {
            c = new BezierCurve;
            from_json(json, *(dynamic_cast<BezierCurve*>(c)));
            foundType = true;
        }
    }
    if (!foundType) {
        c = new BezierCurve;
        from_json(json, *(dynamic_cast<BezierCurve*>(c)));
    }
    return c;
}

template <class Json>
void from_json(const Json& json, Curve*& curve) {
    curve = make_curve_from_json(json);
}
#endif // BSPLINESERIALIZER_H
