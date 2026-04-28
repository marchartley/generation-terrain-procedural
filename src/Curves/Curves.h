#ifndef CURVES_H
#define CURVES_H

#include "Polyline.h"
#include "BezierCurve.h"
#include "CatmullRomSpline.h"

CatmullRomSpline toCatmullRom(const BezierCurve& curve);
CatmullRomSpline toCatmullRom(const Polyline& curve);

BezierCurve toBezier(const CatmullRomSpline& curve);
BezierCurve toBezier(const Polyline& curve);

Polyline toPolyline(const CatmullRomSpline& curve);
Polyline toPolyline(const BezierCurve& curve);

#endif // CURVES_H
