#include "Curves.h"

CatmullRomSpline toCatmullRom(const BezierCurve& curve)
{
    return CatmullRomSpline(curve.getPoints());
}

CatmullRomSpline toCatmullRom(const Polyline& curve)
{
    return CatmullRomSpline(curve.getPoints());
}

CatmullRomSpline toCatmullRom(const BSpline& curve)
{
    return toCatmullRom(toBezier(curve));
}

BezierCurve toBezier(const CatmullRomSpline& curve)
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

BezierCurve toBezier(const Polyline& curve)
{
    return BezierCurve(curve.getPoints());
}

BezierCurve toBezier(const BSpline& curve)
{
    const auto splinePoints = curve.getPoints();
    std::vector<Vector3> points(splinePoints.size() - 2);
    std::vector<Vector3> handles(points.size() * 2 - 2);
    size_t currentPointIdx = 0;
    size_t currentHandleIdx = 0;
    for (int i = 1; i < curve.numPoints() - 2; i++) {
        const auto& B0 = splinePoints[i - 1];
        const auto& B1 = splinePoints[i - 0];
        const auto& B2 = splinePoints[i + 1];
        const auto& B3 = splinePoints[i + 2];

        const auto P0 = (B0 + B1 * 4.f + B2) / 6.f;
        const auto P1 = (B1 * 2.f + B2) / 3.f;
        const auto P2 = (B1 + B2 * 2.f) / 3.f;
        const auto P3 = (B1 + B2 * 4.f + B3) / 6.f;

        if (i == 1) {
            points[currentPointIdx++] = P0;
        }
        points[currentPointIdx++] = P3;
        handles[currentHandleIdx++] = P1;
        handles[currentHandleIdx++] = P2;
    }
    return BezierCurve(points, handles);
}






Polyline toPolyline(const Curve& curve) {
    return Polyline(curve.getPath(curve.numPoints()));
}
/*Polyline toPolyline(const CatmullRomSpline& curve)
{
    return Polyline(curve.getPoints());
}

Polyline toPolyline(const BezierCurve& curve)
{
    return Polyline(curve.getPoints());
}*/

