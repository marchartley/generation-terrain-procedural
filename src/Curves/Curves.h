#ifndef CURVES_H
#define CURVES_H

#include "Polyline.h"
#include "BezierCurve.h"
#include "CatmullRomSpline.h"
#include "BSpline.h"

CatmullRomSpline toCatmullRom(const BezierCurve& curve);
CatmullRomSpline toCatmullRom(const Polyline& curve);
CatmullRomSpline toCatmullRom(const BSpline& curve);

BezierCurve toBezier(const CatmullRomSpline& curve);
BezierCurve toBezier(const Polyline& curve);
BezierCurve toBezier(const BSpline& curve);

Polyline toPolyline(const Curve& curve);
// Polyline toPolyline(const CatmullRomSpline& curve);
// Polyline toPolyline(const BezierCurve& curve);

#endif // CURVES_H
