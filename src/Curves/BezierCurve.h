#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H

#include "Curves/Curve.h"

class BezierCurve : public Curve
{
public:
    BezierCurve();
    BezierCurve(const std::vector<Vector3>& points);
    BezierCurve(const std::vector<Vector3>& points, const std::vector<Vector3>& handles);

    virtual std::vector<Vector3> getPath(int numberOfPoints = -1) const override;
    virtual Vector3 getPoint(float x) const override;
    virtual Vector3 getPoint(float x, const Vector3& a, const Vector3& b) const override;
    virtual Vector3 getDerivative(float x, bool normalize = false) const override;
    virtual Vector3 getSecondDerivative(float x, bool normalize = false) const override;
    virtual float estimateClosestTime(const Vector3& pos, float epsilon = 1e-4, float nbChecksFactor = 2.f, float earlyExitThreshold = 1e-3) const override;
    virtual Vector3 estimateClosestPos(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const override;
    virtual float estimateSqrDistanceFrom(const Vector3& pos, bool useNativeShape = false, float epsilon = 1e-3) const override;
    virtual float length() const override;

    size_t getIndex(int i) { return (i + numPoints()) % numPoints(); }
    virtual BezierCurve& setPoint(int i, const Vector3& newPos) override;

    virtual BezierCurve& resamplePoints(int newNbPoints = -1) override;

    virtual BezierCurve& reverseVertices() override;

    // virtual BezierCurve simplifyByRamerDouglasPeucker(float epsilon, BezierCurve subBezierCurve = BezierCurve()) override;

    virtual std::pair<Vector3, Vector3> AABBox() const override;

    BezierCurve& scale(float factor);
    BezierCurve& scale(const Vector3& factor);
    // BezierCurve scaled(float factor);
    // BezierCurve scaled(const Vector3& factor);

    virtual BezierCurve& translate(const Vector3& translation) override;

    virtual BezierCurve& removeDuplicates() override;

    size_t size() const { return numPoints(); }
    virtual size_t numPoints() const override;

    virtual Vector3& operator[](size_t i) override;
    virtual const Vector3& operator[](size_t i) const override;

    virtual std::string toString() const override;

    virtual BezierCurve& close() override;

    std::vector<Vector3> getPoints() const { return points; }
    std::vector<Vector3> getHandles() const { return handles; }

protected:
    std::vector<Vector3> points;
    std::vector<Vector3> handles;
};

#endif // BEZIERCURVE_H
