#ifndef CURVES_H
#define CURVES_H

#include "Polyline.h"
#include "BezierCurve.h"
#include "BSpline.h"

BSpline toCatmullRom(const BezierCurve& curve);
BSpline toCatmullRom(const Polyline& curve);

BezierCurve toBezier(const BSpline& curve);
BezierCurve toBezier(const Polyline& curve);

Polyline toPolyline(const BSpline& curve);
Polyline toPolyline(const BezierCurve& curve);

#endif // CURVES_H
