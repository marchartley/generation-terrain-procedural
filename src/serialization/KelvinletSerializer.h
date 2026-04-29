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
    json["eps-limit"] = kelvinlet.epsLimit;
    json["mu"] = kelvinlet.mu;
    json["radial-scale"] = kelvinlet.radialScale;
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
    json["curve"] = *(kelvinlet.curve);
}

template <class Json>
void to_json(Json& json, const GrabKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const ScaleKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const TwistKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const PinchKelvinlet& kelvinlet)
{
    to_json(json, static_cast<const KelvinletPoint&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const GrabKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const ScaleKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const TwistKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force();
}

template <class Json>
void to_json(Json& json, const PinchKelvinletCurve& kelvinlet)
{
    to_json(json, static_cast<const KelvinletCurve&>(kelvinlet));
    json["force"] = kelvinlet.force();
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
    kelvinlet.epsLimit = json.at("eps-limit");
    kelvinlet.mu = json.at("mu");
    kelvinlet.radialScale = json.at("radial-scale");
    kelvinlet.v = json.at("v");
}

template <class Json>
void from_json(const Json& json, KelvinletPoint& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.pos = json.at("pos");
}

template <class Json>
void from_json(const Json& json, KelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<Kelvinlet&>(kelvinlet));
    kelvinlet.curve = std::shared_ptr<Curve>(make_curve_from_json(json.at("curve")));
}

template <class Json>
void from_json(const Json& json, GrabKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, ScaleKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, TwistKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, PinchKelvinlet& kelvinlet)
{
    from_json(json, static_cast<KelvinletPoint&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, GrabKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, ScaleKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, TwistKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, PinchKelvinletCurve& kelvinlet)
{
    from_json(json, static_cast<KelvinletCurve&>(kelvinlet));
    kelvinlet.setForce(json.at("force"));
}

template <class Json>
void from_json(const Json& json, RelativeKelvinlet& relativeKelvinlet)
{
    relativeKelvinlet.anchorPoint = json.at("anchor");
    relativeKelvinlet.kelvinlet = json.at("kelvinlet").template get<Kelvinlet*>();
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
Kelvinlet* make_kelvinlet_from_json(const Json& json)
{
    Kelvinlet* kelvinlet;
    const std::string type = toLower(json["type"]);

    const bool isCurve = json.contains("curve");
    const bool isPoint = json.contains("pos");

    if (isPoint) {
        if (type == toLower(GrabKelvinlet().getShortName()))  {
            kelvinlet = new GrabKelvinlet();
            from_json(json, *(dynamic_cast<GrabKelvinlet*>(kelvinlet)));
        }
        if (type == toLower(ScaleKelvinlet().getShortName())) {
            kelvinlet = new ScaleKelvinlet();
            from_json(json, *(dynamic_cast<ScaleKelvinlet*>(kelvinlet)));
        }
        if (type == toLower(TwistKelvinlet().getShortName())) {
            kelvinlet = new TwistKelvinlet();
            from_json(json, *(dynamic_cast<TwistKelvinlet*>(kelvinlet)));
        }
        if (type == toLower(PinchKelvinlet().getShortName())) {
            kelvinlet = new PinchKelvinlet();
            from_json(json, *(dynamic_cast<PinchKelvinlet*>(kelvinlet)));
        }
    }
    else if (isCurve) {
        if (type == toLower(GrabKelvinletCurve().getShortName()))  {
            kelvinlet = new GrabKelvinletCurve();
            from_json(json, *(dynamic_cast<GrabKelvinletCurve*>(kelvinlet)));
        }
        if (type == toLower(ScaleKelvinletCurve().getShortName())) {
            kelvinlet = new ScaleKelvinletCurve();
            from_json(json, *(dynamic_cast<ScaleKelvinletCurve*>(kelvinlet)));
        }
        if (type == toLower(TwistKelvinletCurve().getShortName())) {
            kelvinlet = new TwistKelvinletCurve();
            from_json(json, *(dynamic_cast<TwistKelvinletCurve*>(kelvinlet)));
        }
        if (type == toLower(PinchKelvinletCurve().getShortName())) {
            kelvinlet = new PinchKelvinletCurve();
            from_json(json, *(dynamic_cast<PinchKelvinletCurve*>(kelvinlet)));
        }
    }
    return kelvinlet;
    // throw std::runtime_error("Unknown Kelvinlet type: " + type);
}


template <class Json>
void from_json(const Json& json, Kelvinlet*& kelvinlet)
{
    kelvinlet = make_kelvinlet_from_json(json);
    // from_json(json, *kelvinlet);
}


#endif // KELVINLETSERIALIZER_H
