#include "KelvinletSerializer.h"

#include "Serializer.h"


void to_json(nlohmann::json& json, const Kelvinlet& kelvinlet)
{
    json["type"] = kelvinlet.getShortName();
    json["epsLimit"] = kelvinlet.epsLimit;
    json["mu"] = kelvinlet.mu;
    json["radialScale"] = kelvinlet.radialScale;
    json["v"] = kelvinlet.v;
}
void to_json(nlohmann::json& json, const KelvinletPoint& kelvinlet)
{
    to_json(json, static_cast<const Kelvinlet&>(kelvinlet));
    json["pos"] = kelvinlet.pos;
}

void to_json(nlohmann::json& json, const KelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const Kelvinlet&>(kelvinlet));
    json["curve"] = kelvinlet.curve;
}

void to_json(nlohmann::json& json, const GrabKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const ScaleKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const TwistKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const PinchKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const GrabKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const ScaleKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const TwistKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

void to_json(nlohmann::json& json, const PinchKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}




void from_json(const nlohmann::json& json, Kelvinlet& kelvinlet)
{
    kelvinlet.epsLimit = json["epsLimit"];
    kelvinlet.mu = json["mu"];
    kelvinlet.radialScale = json["radialScale"];
    kelvinlet.v = json["v"];
}

void from_json(const nlohmann::json& json, KelvinletPoint& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.pos = json["pos"];
}

void from_json(const nlohmann::json& json, KelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.curve = json["curve"];
}

void from_json(const nlohmann::json& json, GrabKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, ScaleKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, TwistKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, PinchKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, GrabKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, ScaleKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, TwistKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

void from_json(const nlohmann::json& json, PinchKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}
