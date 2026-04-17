#ifndef KELVINLET_H
#define KELVINLET_H

#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"

class Kelvinlet
{
public:
    Kelvinlet();
    virtual ~Kelvinlet() = default;

    virtual Vector3 evaluate(const Vector3& p) const = 0;
    virtual bool valid() const = 0;

    virtual Kelvinlet* clone() const = 0;
    virtual Kelvinlet& reset() = 0;
    virtual std::string getShortName() const = 0;

    virtual Kelvinlet& translate(const Vector3& translation) = 0;
    virtual Kelvinlet& scale(float scaling) { this->radialScale *= scaling; return *this; }

//    inline float strength(const Vector3& p) const { return densityFunction((p - pos).norm()); }
    inline float regularDistance(float r) const { return std::sqrt(r*r + radialScale*radialScale); }
    inline float densityFunction(float r) const { return (15 * std::pow(radialScale, 4)) / (8 * M_PI * std::pow(regularDistance(r), 7)); }

    float radialScale = 1; // epsilon
    float mu = .5f; // 1.f; // Elastic shear modulus
    float v = 0.5f; //.0f; // Poisson ratio (v = 1/2 is special case of Stokeslets [incompressible])

    float epsLimit = 1000.f; // 100.f;
};

class KelvinletPoint : public Kelvinlet {
public:
    KelvinletPoint();
    virtual ~KelvinletPoint() = default;

    virtual Vector3 evaluate(const Vector3& p) const = 0;
    virtual bool valid() const { return pos.isValid(); }

    virtual KelvinletPoint* clone() const = 0;
    virtual KelvinletPoint& reset() { pos = Vector3::invalid; return *this; }
    virtual KelvinletPoint& translate(const Vector3& translation) { this->pos += translation; return *this; }
    virtual KelvinletPoint& scale(float scaling) { Kelvinlet::scale(scaling); this->pos *= scaling; return *this; }

    Vector3 pos = Vector3::invalid;
};

class GrabKelvinlet : public KelvinletPoint {
public:
    GrabKelvinlet();
    virtual ~GrabKelvinlet() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual GrabKelvinlet& scale(float scaling) { KelvinletPoint::scale(scaling); this->force *= scaling; return *this; }

    virtual GrabKelvinlet* clone() const { return new GrabKelvinlet(*this); }
    virtual std::string getShortName() const { return "Grab"; }
    virtual bool valid() const { return KelvinletPoint::valid() && force.isValid() && force.norm2() > 0; }
    virtual GrabKelvinlet& reset() { KelvinletPoint::reset(); force = Vector3::invalid; return *this; }

    Vector3 force = Vector3::invalid;
};

class TwistKelvinlet : public KelvinletPoint
{
public:
    TwistKelvinlet();
    virtual ~TwistKelvinlet() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual TwistKelvinlet& scale(float scaling) { KelvinletPoint::scale(scaling); this->force *= scaling; return *this; }

    virtual TwistKelvinlet* clone() const { return new TwistKelvinlet(*this); }
    virtual std::string getShortName() const { return "Twist"; }
    virtual bool valid() const { return KelvinletPoint::valid() && force.isValid() && force.norm2() > 0; }
    virtual TwistKelvinlet& reset() { KelvinletPoint::reset(); force = Vector3::invalid; return *this; }

    Vector3 force = Vector3::invalid;
};

class ScaleKelvinlet : public KelvinletPoint
{
public:
    ScaleKelvinlet();
    virtual ~ScaleKelvinlet() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual ScaleKelvinlet& scale(float scaling) { KelvinletPoint::scale(scaling); this->force *= scaling; return *this; }

    virtual ScaleKelvinlet* clone() const { return new ScaleKelvinlet(*this); }
    virtual std::string getShortName() const { return "Scale"; }
    virtual ScaleKelvinlet& reset() { KelvinletPoint::reset(); force = 0; return *this; }
    virtual bool valid() const { return KelvinletPoint::valid() && force != 0; }

    float force = 0;
};

class PinchKelvinlet : public KelvinletPoint
{
public:
    PinchKelvinlet();
    virtual ~PinchKelvinlet() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual PinchKelvinlet& scale(float scaling) { KelvinletPoint::scale(scaling); this->force *= scaling; return *this; }

    virtual PinchKelvinlet* clone() const { return new PinchKelvinlet(*this); }
    virtual std::string getShortName() const { return "Pinch"; }
    virtual bool valid() const { return KelvinletPoint::valid() && force.isValid() && force.norm2() > 0; }
    virtual PinchKelvinlet& reset() { KelvinletPoint::reset(); force = Vector3::invalid; return *this; }

    Vector3 force = Vector3::invalid;
};

class KelvinletCurve : public Kelvinlet {
public:
    KelvinletCurve();
    virtual ~KelvinletCurve() = default;

    virtual Vector3 evaluate(const Vector3& p) const = 0;
    virtual bool valid() const { return !curve.empty(); }

    virtual KelvinletCurve* clone() const = 0;
    virtual KelvinletCurve& reset() { curve = BSpline(); return *this; }
    virtual KelvinletCurve& translate(const Vector3& translation) { this->curve.translate(translation); return *this; }
    virtual KelvinletCurve& scale(float scaling) { Kelvinlet::scale(scaling); this->curve.scale(scaling); return *this; }

    BSpline curve;
};


class GrabKelvinletCurve : public KelvinletCurve {
public:
    GrabKelvinletCurve();
    virtual ~GrabKelvinletCurve() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual GrabKelvinletCurve& scale(float scaling) { KelvinletCurve::scale(scaling); this->force *= scaling; return *this; }

    virtual GrabKelvinletCurve* clone() const { return new GrabKelvinletCurve(*this); }
    virtual std::string getShortName() const { return "Grab Curve"; }
    virtual GrabKelvinletCurve& reset() { KelvinletCurve::reset(); force = 0; return *this; }
    virtual bool valid() const { return KelvinletCurve::valid() && force != 0; }

    float force = 0;
};

class TwistKelvinletCurve : public KelvinletCurve
{
public:
    TwistKelvinletCurve();
    virtual ~TwistKelvinletCurve() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual TwistKelvinletCurve& scale(float scaling) { KelvinletCurve::scale(scaling); this->force *= scaling; return *this; }

    virtual TwistKelvinletCurve* clone() const { return new TwistKelvinletCurve(*this); }
    virtual std::string getShortName() const { return "Twist Curve"; }
    virtual TwistKelvinletCurve& reset() { KelvinletCurve::reset(); force = 0; return *this; }
    virtual bool valid() const { return KelvinletCurve::valid() && force != 0; }

    float force = 0;
};

class ScaleKelvinletCurve : public KelvinletCurve
{
public:
    ScaleKelvinletCurve();
    virtual ~ScaleKelvinletCurve() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual ScaleKelvinletCurve& scale(float scaling) { KelvinletCurve::scale(scaling); this->force *= scaling; return *this; }

    virtual ScaleKelvinletCurve* clone() const { return new ScaleKelvinletCurve(*this); }
    virtual std::string getShortName() const { return "Scale Curve"; }
    virtual ScaleKelvinletCurve& reset() { KelvinletCurve::reset(); force = 0; return *this; }
    virtual bool valid() const { return KelvinletCurve::valid() && force != 0; }

    float force = 0;
};

class PinchKelvinletCurve : public KelvinletCurve
{
public:
    PinchKelvinletCurve();
    virtual ~PinchKelvinletCurve() = default;

    virtual Vector3 evaluate(const Vector3& p) const;

    virtual PinchKelvinletCurve& scale(float scaling) { KelvinletCurve::scale(scaling); this->force *= scaling; return *this; }

    virtual PinchKelvinletCurve* clone() const { return new PinchKelvinletCurve(*this); }
    virtual std::string getShortName() const { return "Pinch Curve"; }
    virtual PinchKelvinletCurve& reset() { KelvinletCurve::reset(); force = 0; return *this; }
    virtual bool valid() const { return KelvinletCurve::valid() && force != 0; }

    float force = 0;
};


struct RelativeKelvinlet {
    RelativeKelvinlet() {}
    RelativeKelvinlet(Kelvinlet* kelvinlet, const Vector3& anchorPoint, float zoom = 1.f);

    Kelvinlet* kelvinlet = nullptr;
    Vector3 anchorPoint = Vector3::invalid;
    float zoom = 1.f;


    Vector3 evaluate(const Vector3& p, float angle = 0, float scaleFactor = 1.f, bool KelvinletPositionIsRelativeToAnchor = false) const;
};

#endif // KELVINLET_H
