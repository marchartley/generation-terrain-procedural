#ifndef VEC3SERIALIZER_H
#define VEC3SERIALIZER_H

#include "DataStructure/Vector3.h"
#include "Utils/json.h"

template <class T>
void to_json(nlohmann::json& json, const Vec3<T>& vec) {
    json = nlohmann::json({{"x", vec.x()}, {"y", vec.y()}, {"z", vec.z()}});
}
template <class T>
void from_json(const nlohmann::json& json, Vec3<T>& vec) {
    if (json.contains("x")) {
        vec = Vec3<T>(json.at("x").get<T>(), json.at("y").get<T>(), (json.contains("z") ? json.at("z").get<T>() : 0));
    } else if (json.contains("r")) {
        vec = Vec3<T>(json.at("r").get<T>(), json.at("g").get<T>(), json.at("b").get<T>());
    }
}

#endif // VEC3SERIALIZER_H
