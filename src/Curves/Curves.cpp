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
    std::vector<Vector3> handles(points.size() * 2);

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
/*
BezierCurve toBezier(const BSpline& curve)
{
    const auto splinePoints = curve.getPoints();

    const size_t nbBezierPoints = splinePoints.size() - 2;

    std::vector<Vector3> points(nbBezierPoints);
    std::vector<Vector3> handles(nbBezierPoints * 2, Vector3::origin);

    size_t currentPointIdx = 0;

    for (size_t i = 1; i < curve.numPoints() - 2; i++) {
        const auto& B0 = splinePoints[i - 1];
        const auto& B1 = splinePoints[i];
        const auto& B2 = splinePoints[i + 1];
        const auto& B3 = splinePoints[i + 2];

        const auto P0 = (B0 + B1 * 4.f + B2) / 6.f;
        const auto P1 = (B1 * 2.f + B2) / 3.f;
        const auto P2 = (B1 + B2 * 2.f) / 3.f;
        const auto P3 = (B1 + B2 * 4.f + B3) / 6.f;

        if (i == 1) {
            points[currentPointIdx++] = P0;
        }

        const size_t prevPointIdx = currentPointIdx - 1;
        const size_t nextPointIdx = currentPointIdx;

        points[currentPointIdx++] = P3;

        handles[BezierCurve::handleOut(prevPointIdx, nbBezierPoints, false)] = P1 - P0;
        handles[BezierCurve::handleIn(nextPointIdx, nbBezierPoints, false)]  = P2 - P3;
    }

    return BezierCurve(points, handles);
}
*/
BezierCurve toBezier(const BSpline& curve)
{
    BSpline tmp = curve;

    const int p = tmp.getDegree(); // 3

    // 1. Insert every internal knot until its multiplicity == degree
    for (float u : tmp.uniqueInternalKnots()) {
        int mult = tmp.knotMultiplicity(u);
        while (mult < p) {
            tmp.insertKnot(u);
            mult++;
        }
    }

    // 2. Now every knot span is one Bezier segment
    const auto& P = tmp.getPoints();
    const auto& U = tmp.getKnots();

    std::vector<Vector3> bezierPoints;
    std::vector<Vector3> handles;

    // For cubic, each Bezier segment uses 4 consecutive control points.
    // Consecutive segments share endpoints.
    for (size_t i = 0; i + 3 < P.size(); i += 3) {
        Vector3 C0 = P[i + 0];
        Vector3 C1 = P[i + 1];
        Vector3 C2 = P[i + 2];
        Vector3 C3 = P[i + 3];

        if (bezierPoints.empty()) {
            bezierPoints.push_back(C0);
            handles.push_back(Vector3::origin); // in handle of first point
            handles.push_back(C1 - C0);         // out handle
        } else {
            handles.back() = C1 - C0;           // out handle of previous point
        }

        bezierPoints.push_back(C3);

        handles.push_back(C2 - C3);             // in handle
        handles.push_back(Vector3::origin);     // out handle, filled by next segment
    }
    return BezierCurve(bezierPoints, handles);
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


BSpline toBSpline(const CatmullRomSpline& curve)
{
    return BSpline(curve.getPoints());
}

BSpline toBSpline(const BezierCurve& curve)
{
    const auto bezierPoints = curve.getPoints();

    const size_t numSegments = bezierPoints.size() - 1;

    std::vector<Vector3> splinePoints(numSegments * 4);
    std::vector<float> knots((numSegments - 1) * 3 + 2 * (3 + 1));
    for (size_t i = 0; i < numSegments; i++) {
        float t = float(i) / float(numSegments);

        const auto C0 = curve.get(i);
        const auto C1 = curve.handlePos(curve.handleOut(i));
        const auto C2 = curve.handlePos(curve.handleIn(i + 1));
        const auto C3 = curve.get(i + 1);

        if (i == 0) {
            splinePoints[0] = C0;
            knots[0] = 0;
        }
        splinePoints[3 * i + 1] = C1;
        splinePoints[3 * i + 2] = C2;
        splinePoints[3 * i + 3] = C3;
        knots[3 * i + 1] = t;
        knots[3 * i + 2] = t;
        knots[3 * i + 3] = t;
    }
    for (size_t i = knots.size() - 4; i < knots.size(); i++) {
        knots[i] = 1.f;
    }
    return BSpline(splinePoints, knots);
    /*

    std::vector<Vector3> splinePoints(numSegments + 3);

    for (size_t i = 0; i < numSegments; i++) {
        const auto& B0 = bezierPoints[i];
        const auto& B1 = curve.handlePos(curve.handleOut(i));
        const auto& B2 = curve.handlePos(curve.handleIn(i + 1));
        const auto& B3 = bezierPoints[i + 1];

        const auto P0 = B0 * 6.f - B1 * 7.f + B2 * 2.f;
        const auto P1 = B1 * 2.f - B2;
        const auto P2 = B2 * 2.f - B1;
        const auto P3 = B3 * 6.f + B1 * 2.f - B2 * 7.f;

        if (i == 0) {
            splinePoints[0] = P0;
            splinePoints[1] = P1;
            splinePoints[2] = P2;
        } else {
            splinePoints[i + 2] = P2;
        }

        splinePoints[i + 3] = P3;
    }

    return BSpline(splinePoints);
    */
}

BSpline toBSpline(const Polyline& curve)
{
    return BSpline(curve.getPoints());
}
