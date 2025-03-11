#include "Kelvinlet.h"
#include "DataStructure/Matrix3.h"


Kelvinlet::Kelvinlet()
{

}

KelvinletPoint::KelvinletPoint()
    : Kelvinlet()
{

}

TranslateKelvinlet::TranslateKelvinlet()
    : KelvinletPoint()
{

}

Vector3 TranslateKelvinlet::evaluate(const Vector3 &p) const
{
    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double c = 2 / (3 * a - 2 * b);

    Matrix first = identity * float((a - b) * 1. / repsilon);
    float second = b / (repsilon3);
    Matrix last = identity * float((a * epsilon * epsilon) / (2 * repsilon3));

    Matrix rrt = rvector.toMatrix() * rvector.toMatrix().transpose();

    const auto force = first + rrt * second + last;

    auto res = Matrix::matprod(force, this->force.toMatrix());
    return Vector3::fromMatrix(res);
}

TwistKelvinlet::TwistKelvinlet()
    : KelvinletPoint()
{

}

Vector3 TwistKelvinlet::evaluate(const Vector3 &p) const
{
    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const double epsilon = radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double epsilon2 = epsilon * epsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double affineScalar = 1 / repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);

    auto qxr = rvector.cross(force);

    double scaleFactor = 100.;
    return -scaleFactor * a * affineScalar * qxr;
}

ScaleKelvinlet::ScaleKelvinlet()
    : KelvinletPoint()
{

}

Vector3 ScaleKelvinlet::evaluate(const Vector3 &p) const
{
    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double epsilon2 = epsilon * epsilon;
    double affineScalar = 1/repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);
    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double scale = 2 * b - a;
    float scaleFactor = 10.f;
    Vector3 delta = rvector * (this->scale * -scaleFactor * scale * affineScalar);

    return delta;
}

PinchKelvinlet::PinchKelvinlet()
    : KelvinletPoint()
{

}

Vector3 PinchKelvinlet::evaluate(const Vector3 &p) const
{
    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double c = 2 / (3 * a - 2 * b);

    float fx = this->force.x;
    float fy = this->force.y;
    float fz = this->force.z;
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
    float rFr = rF.product(r)[0][0];
    Matrix rFrI = identity * rFr;
    Matrix second = (rFrI + F * a * epsilon * epsilon).product(r) * (3.f/(2.f * repsilon3*repsilon*repsilon));

    const auto force = first - second;

    return Vector3::fromMatrix(force);
}



KelvinletCurve::KelvinletCurve()
    : Kelvinlet()
{

}

TranslateKelvinletCurve::TranslateKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 TranslateKelvinletCurve::evaluate(const Vector3 &p) const
{
    float closestTime = curve.estimateClosestTime(p);
    Vector3 pos = curve.getPoint(closestTime);
    Vector3 dir = curve.getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double c = 2 / (3 * a - 2 * b);

    Matrix first = identity * float((a - b) * 1. / repsilon);
    float second = b / (repsilon3);
    Matrix last = identity * float((a * epsilon * epsilon) / (2 * repsilon3));

    Matrix rrt = rvector.toMatrix() * rvector.toMatrix().transpose();

    const auto force = first + rrt * second + last;

    auto res = force * (dir * this->force).toMatrix();
    return Vector3::fromMatrix(res);
}

TwistKelvinletCurve::TwistKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 TwistKelvinletCurve::evaluate(const Vector3 &p) const
{
    float closestTime = curve.estimateClosestTime(p);
    Vector3 pos = curve.getPoint(closestTime);
    Vector3 dir = curve.getBinormal(closestTime); //curve.getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const double epsilon = radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double epsilon2 = epsilon * epsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double affineScalar = 1 / repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);

    auto qxr = rvector.cross(dir * force);

    double scaleFactor = 100.;
    return -scaleFactor * a * affineScalar * qxr;
}

ScaleKelvinletCurve::ScaleKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 ScaleKelvinletCurve::evaluate(const Vector3 &p) const
{
    float closestTime = curve.estimateClosestTime(p);
    Vector3 pos = curve.getPoint(closestTime);
//    Vector3 dir = curve.getDirection(closestTime);

    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double epsilon2 = epsilon * epsilon;
    double affineScalar = 1/repsilon3 + (3 * epsilon2) / (2 * repsilon3 * repsilon * repsilon);
    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double scale = 2 * b - a;
    float scaleFactor = 10.f;
    Vector3 delta = rvector * (this->force * -scaleFactor * scale * affineScalar);

    return delta;
}

PinchKelvinletCurve::PinchKelvinletCurve()
    : KelvinletCurve()
{

}

Vector3 PinchKelvinletCurve::evaluate(const Vector3 &p) const
{
    float closestTime = curve.estimateClosestTime(p);
    Vector3 pos = curve.getPoint(closestTime);
    Vector3 dir = curve.getDirection(closestTime) * this->force;

    Vector3 rvector =  p - pos;
    const double radius = rvector.norm();
    const double repsilon = regularDistance(radius);
    if (repsilon > epsLimit) return Vector3();
    auto identity = Matrix::identity(3);
    const double epsilon = this->radialScale;
    const double repsilon3 = repsilon * repsilon * repsilon;

    double a = 1.f / (4 * M_PI * mu);
    double b = a / (4 * (1.f - v));
    double c = 2 / (3 * a - 2 * b);

    float fx = dir.x;
    float fy = dir.y;
    float fz = dir.z;
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
    float rFr = rF.product(r)[0][0];
    Matrix rFrI = identity * rFr;
    Matrix second = (rFrI + F * a * epsilon * epsilon).product(r) * (3.f/(2.f * repsilon3*repsilon*repsilon));

    const auto force = first - second;

    return Vector3::fromMatrix(force);
}
