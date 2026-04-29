#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include "Curves/Curve.h"

class BezierCurve : public Curve
{
public:
    BezierCurve();
    BezierCurve(const std::vector<Vector3>& points);
    BezierCurve(const std::vector<Vector3>& points, const std::vector<Vector3>& handles);

    CLONE_FUNCTION(BezierCurve)
    std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    Vector3 getPoint(float x) const override;
    Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const override;
    Vector3 getDerivative(float x, bool normalize = false) const override;
    Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    float estimateClosestTime(const Vector3& pos) const override;
    Vector3 estimateClosestPos(const Vector3& pos) const override;
    float estimateSqrDistanceFrom(const Vector3& pos) const override;
    float length() const override;

    size_t getIndex(int i) { return (i + numPoints()) % numPoints(); }
    BezierCurve& setPoint(int i, const Vector3& newPos) override;

    BezierCurve& resamplePoints(int newNbPoints = -1) override;

    BezierCurve& reverseVertices() override;

    // BezierCurve simplifyByRamerDouglasPeucker(float epsilon, BezierCurve subBezierCurve = BezierCurve()) override;

    std::pair<Vector3, Vector3> AABBox() const override;

    using Curve::scale;
    BezierCurve& scale(const Vector3& factor) override;
    // BezierCurve scaled(float factor);
    // BezierCurve scaled(const Vector3& factor);

    BezierCurve& translate(const Vector3& translation) override;

    BezierCurve& removeDuplicates() override;

    size_t size() const { return numPoints(); }
    size_t numPoints() const override;

    Vector3& operator[](size_t i) override;
    const Vector3& operator[](size_t i) const override;

    std::string toString() const override;

    BezierCurve& close() override;

    BezierCurve& reset() { points.clear(); handles.clear(); return *this; }

    inline size_t handleIndex(int index) const;
    inline size_t handleIn(int pointIdx) const;
    inline size_t handleOut(int pointIdx) const;

    std::vector<Vector3> getPoints() const { return points; }
    std::vector<Vector3> getHandles() const { return handles; }

    BezierCurve& autosmooth(int pointIndex);
protected:
    std::vector<Vector3> points;
    std::vector<Vector3> handles;
};

#endif // BEZIERCURVE_H
