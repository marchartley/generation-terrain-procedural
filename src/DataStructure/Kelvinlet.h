#ifndef KELVINLET_H
#define KELVINLET_H

#include "DataStructure/Vector3.h"
#include "Utils/BSpline.h"

class Kelvinlet
{
public:
    Kelvinlet();

    virtual Vector3 evaluate(const Vector3& p) const = 0;

//    inline float strength(const Vector3& p) const { return densityFunction((p - pos).norm()); }
    inline float regularDistance(float r) const { return std::sqrt(r*r + radialScale*radialScale); }
    inline float densityFunction(float r) const { return (15 * std::pow(radialScale, 4)) / (8 * M_PI * std::pow(regularDistance(r), 7)); }

    float radialScale; // epsilon
    float mu = .5f; // 1.f; // Elastic shear modulus
    float v = 0.5f; //.0f; // Poisson ratio (v = 1/2 is special case of Stokeslets [incompressible])

    float epsLimit = 10000.f; // 100.f;
};

class KelvinletPoint : public Kelvinlet {
public:
    KelvinletPoint();

    virtual Vector3 evaluate(const Vector3& p) const = 0;

    Vector3 pos;
};

class TranslateKelvinlet : public KelvinletPoint {
public:
    TranslateKelvinlet();

    virtual Vector3 evaluate(const Vector3& p) const;

    Vector3 force;
};

class TwistKelvinlet : public KelvinletPoint
{
public:
    TwistKelvinlet();

    virtual Vector3 evaluate(const Vector3& p) const;

    Vector3 force;
};

class ScaleKelvinlet : public KelvinletPoint
{
public:
    ScaleKelvinlet();

    virtual Vector3 evaluate(const Vector3& p) const;

    float scale;
};

class PinchKelvinlet : public KelvinletPoint
{
public:
    PinchKelvinlet();

    virtual Vector3 evaluate(const Vector3& p) const;

    Vector3 force;
};

class KelvinletCurve : public Kelvinlet {
public:
    KelvinletCurve();

    virtual Vector3 evaluate(const Vector3& p) const = 0;

    BSpline curve;
};


class TranslateKelvinletCurve : public KelvinletCurve {
public:
    TranslateKelvinletCurve();

    virtual Vector3 evaluate(const Vector3& p) const;

    float force;
};

class TwistKelvinletCurve : public KelvinletCurve
{
public:
    TwistKelvinletCurve();

    virtual Vector3 evaluate(const Vector3& p) const;

    float force;
};

class ScaleKelvinletCurve : public KelvinletCurve
{
public:
    ScaleKelvinletCurve();

    virtual Vector3 evaluate(const Vector3& p) const;

    float force;
};

class PinchKelvinletCurve : public KelvinletCurve
{
public:
    PinchKelvinletCurve();

    virtual Vector3 evaluate(const Vector3& p) const;

    float force;
};

#endif // KELVINLET_H
