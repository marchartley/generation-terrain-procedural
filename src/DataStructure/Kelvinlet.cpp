#include "Kelvinlet.h"
#include "DataStructure/Matrix3.h"


Kelvinlet::Kelvinlet()
{

}

KelvinletPoint::KelvinletPoint()
    : Kelvinlet()
{

}

KelvinletCurve* KelvinletPoint::cloneToCurveKelvinlet() const {
    KelvinletCurve* k;
    std::string my_type = toLower(this->getShortName());
    if (my_type == "grab") k = new GrabKelvinletCurve;
    else if (my_type == "scale") k = new ScaleKelvinletCurve;
    else if (my_type == "twist") k = new TwistKelvinletCurve;
    else if (my_type == "pinch") k = new PinchKelvinletCurve;

    k->radialScale = radialScale;
    k->mu = mu;
    k->v = v;
    k->epsLimit = epsLimit;
    k->forceDir = forceDir;
    k->forceMag = forceMag;
    return k;
}

GrabKelvinlet::GrabKelvinlet()
    : KelvinletPoint()
{

}

Vector3 GrabKelvinlet::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3::origin;
    auto identity = Matrix::identity(3);
    const float epsilon = this->radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float a = 1.f / (4 * M_PI * mu);
    float b = a / (4 * (1.f - v));
    // float c = 2 / (3 * a - 2 * b);

    Matrix first = identity * float((a - b) * 1. / repsilon);
    float second = b / (repsilon3);
    Matrix last = identity * float((a * epsilon * epsilon) / (2 * repsilon3));

    Matrix rrt = rvector.toMatrix() * rvector.toMatrix().transpose();

    const auto force = first + rrt * second + last;

    auto res = Matrix::matprod(force, this->force().toMatrix());
    return Vector3::fromMatrix(res);
}

TwistKelvinlet::TwistKelvinlet()
    : KelvinletPoint()
{

}

Vector3 TwistKelvinlet::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const float epsilon = radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float epsilon2 = epsilon * epsilon;

    float a = 1.f / (4 * M_PI * mu);
    // float b = a / (4 * (1.f - v));
    float affineScalar = 1 / repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);

    auto qxr = rvector.cross(force());

    const float scaleFactor = 100.;
    return -scaleFactor * a * affineScalar * qxr;
}

ScaleKelvinlet::ScaleKelvinlet()
    : KelvinletPoint()
{

}

Vector3 ScaleKelvinlet::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const float epsilon = this->radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float epsilon2 = epsilon * epsilon;
    float affineScalar = 1/repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);
    float a = 1.f / (4 * M_PI * mu);
    float b = a / (4 * (1.f - v));
    float scale = 2 * b - a;
    const float scaleFactor = 10.f;
    Vector3 delta = rvector * (this->force() * -scaleFactor * scale * affineScalar);

    return delta;
}

PinchKelvinlet::PinchKelvinlet()
    : KelvinletPoint()
{

}

Vector3 PinchKelvinlet::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    const Vector3 r = p - pos;
    const float r2 = r.dot(r);
    const float eps = this->radialScale;
    const float re2 = r2 + eps * eps;
    const float re  = std::sqrt(re2);

    if (re > epsLimit) return Vector3();

    const float re4 = re2 * re2;

    const float a = 1.f / (4.f * float(M_PI) * mu);
    const float b = a / (4.f * (1.f - v));

    const float strength = this->forceMag;
    if (strength == 0.f) return Vector3();

    const Vector3 d = this->force() / strength;
    const float dx = d.x(), dy = d.y(), dz = d.z();

    // F = strength * (d d^T - I/3)  (symmetric, trace = 0)
    Matrix I = Matrix::identity(3);
    Matrix ddT = Matrix({
        {dx*dx, dx*dy, dx*dz},
        {dy*dx, dy*dy, dy*dz},
        {dz*dx, dz*dy, dz*dz}
    });
    Matrix F = (ddT - I * (1.f/3.f)) * strength;

    // Compute Fr and s = r^T F r
    Matrix rM = r.toMatrix();                 // 3x1
    Matrix FrM = F.product(rM);               // 3x1
    const Vector3 Fr = Vector3::fromMatrix(FrM);

    Matrix rT = rM.transpose();               // 1x3
    const float s = rT.product(F).product(rM)(0, 0);

    // Pinch Kelvinlet (regularized) displacement:
    //
    // p_eps(r)= -2a(1/re^2 + eps^2/re^4) F r
    //         + 4b[ (1/re^2) F - (s/re^4) I ] r
    //
    // which can be arranged as:
    // u = [(-2a(1/re2 + eps^2/re4) + 4b*(1/re2))] * (Fr)
    //     - 4b*(s/re4) * r
    //
    const float invRe2 = 1.f / re2;
    const float invRe4 = 1.f / re4;

    const float coeffFr =
        (-2.f * a) * (invRe2 + (eps * eps) * invRe4)
        + (4.f * b) * invRe2;

    const float coeffR = -(4.f * b) * (s * invRe4);

    return Fr * coeffFr + r * coeffR;
}


KelvinletCurve::KelvinletCurve()
    : Kelvinlet()
{

}

GrabKelvinletCurve::GrabKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 GrabKelvinletCurve::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    float closestTime = curve->estimateClosestTime(p);
    Vector3 pos = curve->getPoint(closestTime);
    Vector3 dir = curve->getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const float epsilon = this->radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float a = 1.f / (4 * M_PI * mu);
    float b = a / (4 * (1.f - v));
    // float c = 2 / (3 * a - 2 * b);

    Matrix first = identity * float((a - b) * 1. / repsilon);
    float second = b / (repsilon3);
    Matrix last = identity * float((a * epsilon * epsilon) / (2 * repsilon3));

    Matrix rrt = rvector.toMatrix() * rvector.toMatrix().transpose();

    const auto force = first + rrt * second + last;

    auto res = force * (dir * this->force()).toMatrix();
    return Vector3::fromMatrix(res);
}

TwistKelvinletCurve::TwistKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 TwistKelvinletCurve::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    float closestTime = curve->estimateClosestTime(p);
    Vector3 pos = curve->getPoint(closestTime);
    Vector3 dir = curve->getBinormal(closestTime); //curve->getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const float epsilon = radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float epsilon2 = epsilon * epsilon;

    float a = 1.f / (4 * M_PI * mu);
    // float b = a / (4 * (1.f - v));
    float affineScalar = 1 / repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);

    auto qxr = rvector.cross(dir * force());

    float scaleFactor = 100.;
    return -scaleFactor * a * affineScalar * qxr;
}

ScaleKelvinletCurve::ScaleKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 ScaleKelvinletCurve::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    float closestTime = curve->estimateClosestTime(p);
    Vector3 pos = curve->getPoint(closestTime);
//    Vector3 dir = curve->getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const float epsilon = this->radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float epsilon2 = epsilon * epsilon;
    float affineScalar = 1/repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);
    float a = 1.f / (4 * M_PI * mu);
    float b = a / (4 * (1.f - v));
    float scale = 2 * b - a;
    float scaleFactor = 10.f;
    Vector3 delta = rvector * (this->force() * -scaleFactor * scale * affineScalar);

    return delta;
}

PinchKelvinletCurve::PinchKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 PinchKelvinletCurve::evaluate(const Vector3 &p) const
{
    if (!this->valid()) return Vector3::origin;
    float closestTime = curve->estimateClosestTime(p);
    Vector3 pos = curve->getPoint(closestTime);
    Vector3 dir = curve->getDirection(closestTime) * this->force();

    Vector3 rvector =  p - pos;
    const float radius = rvector.norm();
    const float repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const float epsilon = this->radialScale;
    const float repsilon3 = repsilon * repsilon * repsilon;

    float a = 1.f / (4 * M_PI * mu);
    float b = a / (4 * (1.f - v));
    float c = 2 / (3 * a - 2 * b);

    float fx = dir.x();
    float fy = dir.y();
    float fz = dir.z();
    float mean = (fx + fy + fz) / 3.f;
    fx -= mean;
    fy -= mean;
    fz -= mean;
    Matrix F = Matrix({
                          {fx, a, b},
                          {a, fy, c},
                          {b, c, fz}
                      });

    Matrix r = rvector.toMatrix();
    Matrix first = F.product(r) * (2 * b - a) / repsilon3;
    Matrix rF = r.transpose().product(F);
    float rFr = rF.product(r)(0, 0);
    Matrix rFrI = identity * rFr;
    Matrix second = (rFrI + F * a * epsilon * epsilon).product(r) * (3.f/(2.f * repsilon3*repsilon*repsilon));

    const auto force = first - second;

    return Vector3::fromMatrix(force) * sign(this->force());
}

RelativeKelvinlet::RelativeKelvinlet(Kelvinlet* kelvinlet, const Vector3 &anchorPoint, float zoom)
    : kelvinlet(kelvinlet), anchorPoint(anchorPoint), zoom(zoom)
{}

Vector3 RelativeKelvinlet::evaluate(const Vector3& p, float angle, float scaleFactor, bool KelvinletPositionIsRelativeToAnchor) const
{
    Vector3 evaluationPoint = ((p - anchorPoint) / zoom).rotate(Vector3(0, 0, -angle)) + (KelvinletPositionIsRelativeToAnchor ? Vector3::origin : anchorPoint * zoom);
    return kelvinlet->evaluate(evaluationPoint).rotate(Vector3(0, 0, angle)) * scaleFactor;
}
