#ifndef CURVES_H
#define CURVES_H

#include "Polyline.h"
#include "BezierCurve.h"
#include "CatmullRomSpline.h"
#include "BSpline.h"

CatmullRomSpline toCatmullRom(const BezierCurve& curve);
CatmullRomSpline toCatmullRom(const Polyline& curve);
CatmullRomSpline toCatmullRom(const BSpline& curve);
CatmullRomSpline toCatmullRom(const CatmullRomSpline& curve);


BezierCurve toBezier(const CatmullRomSpline& curve);
BezierCurve toBezier(const Polyline& curve);
BezierCurve toBezier(const BSpline& curve);
BezierCurve toBezier(const BezierCurve& curve);


Polyline toPolyline(const Curve& curve, int samplesFactor = 2);
Polyline toPolyline(const Polyline& curve);

BSpline toBSpline(const CatmullRomSpline& curve);
BSpline toBSpline(const BezierCurve& curve);
BSpline toBSpline(const Polyline& curve);
BSpline toBSpline(const BSpline& curve);

std::vector<Vector3> curveIntersection(const Polyline& a,           const Polyline& b);
std::vector<Vector3> curveIntersection(const Polyline& a,           const CatmullRomSpline& b);
std::vector<Vector3> curveIntersection(const Polyline& a,           const BezierCurve& b);
std::vector<Vector3> curveIntersection(const Polyline& a,           const BSpline& b);

inline std::vector<Vector3> curveIntersection(const CatmullRomSpline& a,   const Polyline& b) { return curveIntersection(b, a); }
std::vector<Vector3> curveIntersection(const CatmullRomSpline& a,   const CatmullRomSpline& b);
std::vector<Vector3> curveIntersection(const CatmullRomSpline& a,   const BezierCurve& b);
std::vector<Vector3> curveIntersection(const CatmullRomSpline& a,   const BSpline& b);

inline std::vector<Vector3> curveIntersection(const BezierCurve& a,        const Polyline& b) { return curveIntersection(b, a); }
inline std::vector<Vector3> curveIntersection(const BezierCurve& a,        const CatmullRomSpline& b) { return curveIntersection(b, a); }
std::vector<Vector3> curveIntersection(const BezierCurve& a,        const BezierCurve& b);
std::vector<Vector3> curveIntersection(const BezierCurve& a,        const BSpline& b);

inline std::vector<Vector3> curveIntersection(const BSpline& a,            const Polyline& b) { return curveIntersection(b, a); }
inline std::vector<Vector3> curveIntersection(const BSpline& a,            const CatmullRomSpline& b) { return curveIntersection(b, a); }
inline std::vector<Vector3> curveIntersection(const BSpline& a,            const BezierCurve& b) { return curveIntersection(b, a); }
std::vector<Vector3> curveIntersection(const BSpline& a,            const BSpline& b);

#endif // CURVES_H
