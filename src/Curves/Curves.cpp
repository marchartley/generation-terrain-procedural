#include "Curves.h"

BSpline toCatmullRom(const BezierCurve &curve)
{
    return BSpline(curve.getPoints());
}

BSpline toCatmullRom(const Polyline &curve)
{
    return BSpline(curve.getPoints());
}

BezierCurve toBezier(const BSpline &curve)
{
    std::vector<Vector3> points = curve.getPoints();
    std::vector<Vector3> handles((points.size() - 1) * 2);

    for (size_t i = 0; i < points.size() - 1; i++) {
        float t0 = float(i+0) / float(points.size() - 1);
        float t1 = float(i+1) / float(points.size() - 1);

        float scale = std::pow((points[i+1]-points[i]).norm(), curve.alpha) / 3.f;

        handles[2 * i + 0] = points[i + 0] + curve.getDerivative(t0) * scale; //points[i-1], points[i], points[i+1], points[i+2], curve.alpha, 0.f) * scale;
        handles[2 * i + 1] = points[i + 1] - curve.getDerivative(t1) * scale; //points[i-1], points[i], points[i+1], points[i+2], curve.alpha, 1.f) * scale;
    }

    return BezierCurve(points, handles);
}

BezierCurve toBezier(const Polyline &curve)
{
    return BezierCurve(curve.getPoints());
}







Polyline toPolyline(const BSpline &curve)
{
    return Polyline(curve.getPoints());
}

Polyline toPolyline(const BezierCurve &curve)
{
    return Polyline(curve.getPoints());
}
