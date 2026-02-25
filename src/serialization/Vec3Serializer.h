#ifndef VEC3SERIALIZER_H
#define VEC3SERIALIZER_H

#include "DataStructure/Vector3.h"
#include "Utils/json.h"

template <class Json, class T>
void to_json(Json& json, const Vec3<T>& vec) {
    json = Json({{"x", vec.x()}, {"y", vec.y()}, {"z", vec.z()}});
}
template <class Json, class T>
void from_json(const Json& json, Vec3<T>& vec) {
    if (json.contains("x")) {
        vec = Vec3<T>(json.at("x").template get<T>(), json.at("y").template get<T>(), (json.contains("z") ? json.at("z").template get<T>() : 0));
    } else if (json.contains("r")) {
        vec = Vec3<T>(json.at("r").template get<T>(), json.at("g").template get<T>(), json.at("b").template get<T>());
    }
}

#endif // VEC3SERIALIZER_H
