#ifndef KELVINLETSERIALIZER_H
#define KELVINLETSERIALIZER_H

#include "DataStructure/Kelvinlet.h"
#include "Utils/json.h"


template <class Json>
void to_json(Json& json, const Kelvinlet& kelvinlet);
template <class Json>
void to_json(Json& json, const KelvinletPoint& kelvinlet);
template <class Json>
void to_json(Json& json, const KelvinletCurve& kelvinlet);
template <class Json>
void to_json(Json& json, const GrabKelvinlet& kelvinlet);
template <class Json>
void to_json(Json& json, const ScaleKelvinlet& kelvinlet);
template <class Json>
void to_json(Json& json, const TwistKelvinlet& kelvinlet);
template <class Json>
void to_json(Json& json, const PinchKelvinlet& kelvinlet);
template <class Json>
void to_json(Json& json, const GrabKelvinletCurve& kelvinlet);
template <class Json>
void to_json(Json& json, const ScaleKelvinletCurve& kelvinlet);
template <class Json>
void to_json(Json& json, const TwistKelvinletCurve& kelvinlet);
template <class Json>
void to_json(Json& json, const PinchKelvinletCurve& kelvinlet);

template <class Json>
void to_json(Json &json, const RelativeKelvinlet& relativeKelvinlet);


template <class Json>
void from_json(const Json& json, Kelvinlet& kelvinlet);
template <class Json>
void from_json(const Json& json, KelvinletPoint& kelvinlet);
template <class Json>
void from_json(const Json& json, KelvinletCurve& kelvinlet);
template <class Json>
void from_json(const Json& json, GrabKelvinlet& kelvinlet);
template <class Json>
void from_json(const Json& json, ScaleKelvinlet& kelvinlet);
template <class Json>
void from_json(const Json& json, TwistKelvinlet& kelvinlet);
template <class Json>
void from_json(const Json& json, PinchKelvinlet& kelvinlet);
template <class Json>
void from_json(const Json& json, GrabKelvinletCurve& kelvinlet);
template <class Json>
void from_json(const Json& json, ScaleKelvinletCurve& kelvinlet);
template <class Json>
void from_json(const Json& json, TwistKelvinletCurve& kelvinlet);
template <class Json>
void from_json(const Json& json, PinchKelvinletCurve& kelvinlet);

template <class Json>
void from_json(const Json &json, RelativeKelvinlet& relativeKelvinlet);


template <class Json>
void to_json(Json& json, const Kelvinlet* kelvinlet);

template <class Json>
void from_json(const Json& json, Kelvinlet*& kelvinlet);












#include "Serializer.h"

template <class Json>
void to_json(Json& json, const Kelvinlet& kelvinlet)
{
    json["type"] = kelvinlet.getShortName();
    json["epsLimit"] = kelvinlet.epsLimit;
    json["mu"] = kelvinlet.mu;
    json["radialScale"] = kelvinlet.radialScale;
    json["v"] = kelvinlet.v;
}
template <class Json>
void to_json(Json& json, const KelvinletPoint& kelvinlet)
{
    to_json(json, static_cast<const Kelvinlet&>(kelvinlet));
    json["pos"] = kelvinlet.pos;
}

template <class Json>
void to_json(Json& json, const KelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const Kelvinlet&>(kelvinlet));
    json["curve"] = kelvinlet.curve;
}

template <class Json>
void to_json(Json& json, const GrabKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const ScaleKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const TwistKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const PinchKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const GrabKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const ScaleKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const TwistKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const PinchKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force;
}

template <class Json>
void to_json(Json& json, const RelativeKelvinlet& relativeKelvinlet)
{
    json["anchor"] = relativeKelvinlet.anchorPoint;
    json["kelvinlet"] = relativeKelvinlet.kelvinlet;
}




template <class Json>
void from_json(const Json& json, Kelvinlet& kelvinlet)
{
    kelvinlet.epsLimit = json["epsLimit"];
    kelvinlet.mu = json["mu"];
    kelvinlet.radialScale = json["radialScale"];
    kelvinlet.v = json["v"];
}

template <class Json>
void from_json(const Json& json, KelvinletPoint& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.pos = json["pos"];
}

template <class Json>
void from_json(const Json& json, KelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.curve = json["curve"];
}

template <class Json>
void from_json(const Json& json, GrabKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, ScaleKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, TwistKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, PinchKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, GrabKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, ScaleKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, TwistKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, PinchKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.force = json["force"];
}

template <class Json>
void from_json(const Json& json, RelativeKelvinlet& relativeKelvinlet)
{
    relativeKelvinlet.anchorPoint = json["anchor"];
    relativeKelvinlet.kelvinlet = json["kelvinlet"].template get<Kelvinlet*>();
}
















template <class Json>
void to_json(Json& json, const KelvinletPoint* kelvinlet);
template <class Json>
void to_json(Json& json, const KelvinletCurve* kelvinlet);


template <class Json>
void to_json(Json& json, const Kelvinlet* kelvinlet)
{
    if (auto k = dynamic_cast<const KelvinletPoint*>(kelvinlet))
        to_json(json, k);
    if (auto k = dynamic_cast<const KelvinletCurve*>(kelvinlet))
        to_json(json, k);
}


template <class Json>
void to_json(Json& json, const KelvinletPoint* kelvinlet)
{
    if (auto k = dynamic_cast<const GrabKelvinlet*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const TwistKelvinlet*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const ScaleKelvinlet*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const PinchKelvinlet*>(kelvinlet))
        to_json(json, *k);
}

template <class Json>
void to_json(Json& json, const KelvinletCurve* kelvinlet)
{
    if (auto k = dynamic_cast<const GrabKelvinletCurve*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const TwistKelvinletCurve*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const ScaleKelvinletCurve*>(kelvinlet))
        to_json(json, *k);
    if (auto k = dynamic_cast<const PinchKelvinletCurve*>(kelvinlet))
        to_json(json, *k);
}

template <class Json>
Kelvinlet* make_kelvinlet_from_json(const Json& j)
{
    const std::string type = toLower(j["type"]);

    const bool isCurve = j.contains("curve");
    const bool isPoint = j.contains("pos");

    if (isPoint) {
        if (type == toLower(GrabKelvinlet().getShortName()))  return new GrabKelvinlet();
        if (type == toLower(ScaleKelvinlet().getShortName())) return new ScaleKelvinlet();
        if (type == toLower(TwistKelvinlet().getShortName())) return new TwistKelvinlet();
        if (type == toLower(PinchKelvinlet().getShortName())) return new PinchKelvinlet();
    }
    else if (isCurve) {
        if (type == toLower(GrabKelvinletCurve().getShortName()))  return new GrabKelvinletCurve();
        if (type == toLower(ScaleKelvinletCurve().getShortName())) return new ScaleKelvinletCurve();
        if (type == toLower(TwistKelvinletCurve().getShortName())) return new TwistKelvinletCurve();
        if (type == toLower(PinchKelvinletCurve().getShortName())) return new PinchKelvinletCurve();
    }

    throw std::runtime_error("Unknown Kelvinlet type: " + type);
}


template <class Json>
void from_json(const Json& json, Kelvinlet*& kelvinlet)
{
    kelvinlet = make_kelvinlet_from_json(json);
    from_json(json, *kelvinlet);
}


#endif // KELVINLETSERIALIZER_H
