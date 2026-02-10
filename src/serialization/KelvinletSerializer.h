#ifndef KELVINLETSERIALIZER_H
#define KELVINLETSERIALIZER_H

#include "DataStructure/Kelvinlet.h"
#include "Utils/json.h"


void to_json(nlohmann::json& json, const Kelvinlet& kelvinlet);
void to_json(nlohmann::json& json, const KelvinletPoint& kelvinlet);
void to_json(nlohmann::json& json, const KelvinletCurve& kelvinlet);
void to_json(nlohmann::json& json, const GrabKelvinlet& kelvinlet);
void to_json(nlohmann::json& json, const ScaleKelvinlet& kelvinlet);
void to_json(nlohmann::json& json, const TwistKelvinlet& kelvinlet);
void to_json(nlohmann::json& json, const PinchKelvinlet& kelvinlet);
void to_json(nlohmann::json& json, const GrabKelvinletCurve& kelvinlet);
void to_json(nlohmann::json& json, const ScaleKelvinletCurve& kelvinlet);
void to_json(nlohmann::json& json, const TwistKelvinletCurve& kelvinlet);
void to_json(nlohmann::json& json, const PinchKelvinletCurve& kelvinlet);

void from_json(const nlohmann::json& json, Kelvinlet& kelvinlet);
void from_json(const nlohmann::json& json, KelvinletPoint& kelvinlet);
void from_json(const nlohmann::json& json, KelvinletCurve& kelvinlet);
void from_json(const nlohmann::json& json, GrabKelvinlet& kelvinlet);
void from_json(const nlohmann::json& json, ScaleKelvinlet& kelvinlet);
void from_json(const nlohmann::json& json, TwistKelvinlet& kelvinlet);
void from_json(const nlohmann::json& json, PinchKelvinlet& kelvinlet);
void from_json(const nlohmann::json& json, GrabKelvinletCurve& kelvinlet);
void from_json(const nlohmann::json& json, ScaleKelvinletCurve& kelvinlet);
void from_json(const nlohmann::json& json, TwistKelvinletCurve& kelvinlet);
void from_json(const nlohmann::json& json, PinchKelvinletCurve& kelvinlet);

#endif // KELVINLETSERIALIZER_H
