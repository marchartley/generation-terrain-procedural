#include "Curves.h"

#include "Utils/Collisions.h"

CatmullRomSpline toCatmullRom(const BezierCurve& curve)
{
    auto result = CatmullRomSpline(curve.getPoints());
    result.setClosed(curve.isClosed());
    return result;
}

CatmullRomSpline toCatmullRom(const Polyline& curve)
{
    auto result = CatmullRomSpline(curve.getPoints());
    result.setClosed(curve.isClosed());
    return result;
}

CatmullRomSpline toCatmullRom(const BSpline& curve)
{
    return CatmullRomSpline(curve.getPath(curve.numSegments() * 4)).setClosed(curve.isClosed());
}

BezierCurve toBezier(const CatmullRomSpline& curve)
{
    std::vector<Vector3> points = curve.getPoints();
    std::vector<Vector3> handles(points.size() * 2);

    for (size_t i = 0; i < points.size(); i++) {
        float t0 = float(i+0) / float(points.size() - (curve.isClosed() ? 0 : 1));
        float t1 = float(i+1) / float(points.size() - (curve.isClosed() ? 0 : 1));

        // float scale = std::pow((points[i+1]-points[i]).norm(), curve.alpha) / 3.f;
        float scale = 1.f / (3.f * float(points.size() - 1));


        const Vector3 hOut = curve.getDerivative(t0) * scale;
        const Vector3 hIn = -curve.getDerivative(t1) * scale;
        handles[BezierCurve::handleOut(i, points.size(), curve.isClosed())] = hOut;
        handles[BezierCurve::handleIn(i + 1, points.size(), curve.isClosed())] = hIn;
    }

    auto result = BezierCurve(points, handles);
    result.setClosed(curve.isClosed());
    return result;
}

BezierCurve toBezier(const Polyline& curve)
{
    auto result = BezierCurve(curve.getPoints());
    result.setClosed(curve.isClosed());
    return result;
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

    auto result = BezierCurve(points, handles);
    result.setClosed(curve.isClosed());
    return result;
}
*/
BezierCurve toBezier(const BSpline& curve)
{
    BSpline tmp = (curve.isClosed() ? curve.generateFakeClosedCurve() : curve);

    const int p = tmp.getDegree();
    assert(p == 3);

    // 1. Insert every internal knot until its multiplicity == degree
    auto a = tmp.getKnots()[p];
    auto b = tmp.getKnots()[tmp.getKnots().size() - p - 1];

    while (tmp.knotMultiplicity(a) < p + 1)
        tmp.insertKnot(a);

    while (tmp.knotMultiplicity(b) < p + 1)
        tmp.insertKnot(b);

    auto internal = tmp.uniqueInternalKnots(a, b);

    for (float u : internal) {
        if (u <= a || u >= b)
            continue;

        while (tmp.knotMultiplicity(u) < p)
            tmp.insertKnot(u);
    }

    std::vector<Vector3> bezierPoints;
    std::vector<Vector3> handles;
    const auto& P = tmp.getPoints();
    const auto& U = tmp.getKnots();

    for (size_t j = p; j + 1 < U.size() - p; ++j) {
        if (U[j] == U[j + 1])
            continue;

        if (U[j] < a || U[j + 1] > b)
            continue;

        Vector3 C0 = P[j - 3];
        Vector3 C1 = P[j - 2];
        Vector3 C2 = P[j - 1];
        Vector3 C3 = P[j];

        if (bezierPoints.empty()) {
            bezierPoints.push_back(C0);
            handles.push_back(Vector3::origin);
            handles.push_back(C1 - C0);
        } else {
            handles.back() = C1 - C0;
        }

        bezierPoints.push_back(C3);

        handles.push_back(C2 - C3);
        handles.push_back(Vector3::origin);
    }
    auto result = BezierCurve(bezierPoints, handles);
    result.setClosed(curve.isClosed());
    return result;
}




Polyline toPolyline(const Curve& curve, int samplesFactor) {
    auto result = Polyline(curve.getPath(curve.numPoints() * samplesFactor));
    result.setClosed(curve.isClosed());
    return result;
}
Polyline toPolyline(const Polyline& curve) {
    return curve;
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
    auto result = BSpline(curve.getPoints());
    result.setClosed(curve.isClosed());
    return result;
}

BSpline toBSpline(const BezierCurve& curve)
{
    const size_t numSegments = curve.getPoints().size() - 1;
    assert(numSegments >= 1);

    constexpr int p = 3;

    std::vector<Vector3> splinePoints(3 * numSegments + 1);
    std::vector<float> knots;
    knots.reserve(3 * (numSegments - 1) + 2 * (p + 1));

    // Start clamped knot: 0,0,0,0
    for (int i = 0; i < p + 1; ++i)
        knots.push_back(0.f);

    for (size_t i = 0; i < numSegments; ++i) {
        const Vector3 C0 = curve.get(i);
        const Vector3 C1 = curve.handlePos(curve.handleOut(i));
        const Vector3 C2 = curve.handlePos(curve.handleIn(i + 1));
        const Vector3 C3 = curve.get(i + 1);

        if (i == 0)
            splinePoints[0] = C0;

        splinePoints[3 * i + 1] = C1;
        splinePoints[3 * i + 2] = C2;
        splinePoints[3 * i + 3] = C3;

        // Internal knot with multiplicity degree = 3
        if (i > 0) {
            const float u = float(i) / float(numSegments);
            knots.push_back(u);
            knots.push_back(u);
            knots.push_back(u);
        }
    }

    // End clamped knot: 1,1,1,1
    for (int i = 0; i < p + 1; ++i)
        knots.push_back(1.f);

    // assert(knots.size() == splinePoints.size() + p + 1);

    auto result = BSpline(splinePoints, knots);
    result.setClosed(curve.isClosed());
    return result;

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
    auto result = BSpline(curve.getPoints());
    result.setClosed(curve.isClosed());
    return result;
}

CatmullRomSpline toCatmullRom(const CatmullRomSpline& curve)
{
    return curve;
}

BezierCurve toBezier(const BezierCurve& curve)
{
    return curve;
}

BSpline toBSpline(const BSpline& curve)
{
    return curve;
}

std::vector<Vector3> curveIntersection(const Polyline &a, const Polyline &b)
{
    std::vector<Vector3> intersections;
    for (size_t i = 0; i < a.numSegments(); i++) {
        for (size_t j = 0; j < b.numSegments(); j++) {
            auto collide = Collision::intersectionBetweenTwoSegments(a.get(i), a.get(i+1), b.get(j), b.get(j+1), 1e-8);
            if (collide.isValid()) intersections.push_back(collide);
        }
    }
    return intersections;
}

std::vector<Vector3> curveIntersection(const Polyline& a, const CatmullRomSpline& b)
{
    return curveIntersection(a, toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const Polyline& a, const BezierCurve& b)
{
    return curveIntersection(a, toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const Polyline& a, const BSpline& b)
{
    return curveIntersection(a, toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const CatmullRomSpline& a, const CatmullRomSpline& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const CatmullRomSpline& a, const BezierCurve& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const CatmullRomSpline& a, const BSpline& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const BezierCurve& a, const BezierCurve& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const BezierCurve& a, const BSpline& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}

std::vector<Vector3> curveIntersection(const BSpline& a, const BSpline& b)
{
    return curveIntersection(toPolyline(a, 5), toPolyline(b, 5));
}
