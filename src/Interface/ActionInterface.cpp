#include "ActionInterface.h"

ActionInterface::ActionInterface(std::string actionTypeName, QWidget *parent)
    : CustomInteractiveObject(parent), actionType(actionTypeName)
{

}

ActionInterface::~ActionInterface()
{
}

nlohmann::json vec3_to_json(const Vector3& vec) {
    return nlohmann::json({{"x", vec.x}, {"y", vec.y}, {"z", vec.z}});
}

Vector3 json_to_vec3(nlohmann::json json)
{
    return Vector3(json.at("x").get<float>(), json.at("y").get<float>(), json.at("z").get<float>());
}

nlohmann::json bspline_to_json(const BSpline& spline) {
    std::vector<nlohmann::json> points;
    for (const auto& p : spline.points) {
        points.push_back(vec3_to_json(p));
    }
    return nlohmann::json({
                              {"points", points},
                              {"closed", spline.closed}
                          });
}
BSpline json_to_bspline(nlohmann::json json) {
    BSpline spline;
    for (auto& point : json.at("points"))
        spline.points.push_back(json_to_vec3(point));
    if (json.at("closed").get<bool>())
        spline.close();
    return spline;
}
nlohmann::json karst_profile_to_json(KarstHoleProfile profile) {
    return nlohmann::json({
                              {"outlines", bspline_to_json(profile.vertices) },
                              {"scaling", vec3_to_json(profile.scaling) }
                          });
}
KarstHoleProfile json_to_karst_profile(nlohmann::json json) {
    KarstHoleProfile profile(json_to_bspline(json.at("outlines")));
    profile.scaling = json_to_vec3(json.at("scaling"));
    return profile;
}
